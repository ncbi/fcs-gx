/*
-----------------------------------------------------------------------------
                             PUBLIC DOMAIN NOTICE
                 National Center for Biotechnology Information

  This software is a "United States Government Work" under the terms of the
  United States Copyright Act.  It was written as part of the author's official
  duties as a United States Government employees and thus cannot be copyrighted.
  This software is freely available to the public for use. The National Library
  of Medicine and the U.S. Government have not placed any restriction on its use
  or reproduction.

  Although all reasonable efforts have been taken to ensure the accuracy and
  reliability of this software, the NLM and the U.S. Government do not and
  cannot warrant the performance or results that may be obtained by using this
  software. The NLM and the U.S. Government disclaim all warranties, expressed
  or implied, including warranties of performance, merchantability or fitness
  for any particular purpose.

  Please cite NCBI in any work or product based on this material.

-----------------------------------------------------------------------------
*/
#define RANGELESS_ENABLE_TSV 1

#include "types.hpp"
#include "segment.hpp"

#include <set>
#include <cmath>
#include <smmintrin.h>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;

/////////////////////////////////////////////////////////////////////////////

// GP-33388
struct repeat_tracks_t
{
    ivls_t transposons;     // transposons from query-genome.
    ivls_t low_complexity;  // low-entropy regions.
    ivls_t ultraconserved;  // regions covered by hits from at least 5 distinct non-prok divs.
    ivls_t n_runs;          // spans of Ns.
    ivls_t combined;        // Union of the above.

    void finalize()
    {
        combined.clear();
        for (const auto p : { &transposons, &low_complexity, &ultraconserved, &n_runs }) {
            ivl_t::sort_and_merge(*p);
            combined <<= std::as_const(*p);
        }
        ivl_t::sort_and_merge(combined);

        // make non-redundant. GP-33482
        // NB: this is relected in the repeats sub-columns order: (xp,lc,co,n).
        low_complexity = ivl_t::subtract(low_complexity, ultraconserved);

        transposons    = ivl_t::subtract(transposons, low_complexity);
        transposons    = ivl_t::subtract(transposons, ultraconserved);
    }

    void clear()
    {
        for (const auto p : { &transposons, &low_complexity, &ultraconserved, &n_runs, &combined }) {
            p->clear();
        }
    }
};

/////////////////////////////////////////////////////////////////////////////

struct node_t : public ivl_t // ivl on query coordinates
{
    enum class ord_t : uint32_t{ none = uint32_t(-1) }; // 0-based ordinal into nodes_t

    tax_id_t tax_id   = {}; // sbj-tax-id
    float seg_score   = 0;  // squared l2-norm

    static node_t collapse(const node_t& a, node_t b)
    {
        auto combined_score = a.seg_score + b.seg_score;

        // If have overlap, proportionally subtract score-contribution of a or b, whichever is less.
        //
        // We weight it a little less (0.8) to give some score boost for overlaps, to promote
        // score to species that have multiple hits. E.g. NZ_BMTN01000049.1 has high identity full-coverage
        // hits to several species of prok:high GC Gram+, and lower identity full-coverage hits to primate repeats,
        // so we want the scoring model to yield higher score for the latter, even though hits
        // to proks are higher quality.
        const auto overlap = (float)ivl_t::uoverlap(a, b);
        if (overlap > 0) {
            const auto a_frac_score = a.seg_score / float(a.len);
            const auto b_frac_score = b.seg_score / float(b.len);
            combined_score -= overlap * std::min(a_frac_score, b_frac_score) * 0.8f;
        }

        const auto combined = ivl_t::collapse(a, b);
        b.pos = combined.pos;
        b.len = combined.len;
        b.seg_score = combined_score;

        return b;
    }

};
static_assert(sizeof(node_t) == 16, "");
using nodes_t = std::vector<node_t>;

static inline uint32_t operator+(node_t::ord_t ord)
{
    return (uint32_t)ord;
}

#if 0
static void print_aggregate_scores_by_tax_ids(const nodes_t& nodes)
{
    struct dat_t
    {
        len_t len = 0;
        double score = 0;
    };

    std::cerr << "\nAggregate scores by tax-id\n";
    nodes
  % fn::foldl(std::map<tax_id_t, dat_t>{}, [](auto ret, const node_t& node)
    {
        auto& dest = ret[node.tax_id];
        dest.len += node.len;
        dest.score += node.seg_score;
        return ret;
    })
  % fn::to_vector()
  % fn::sort_by L(_.second.score * -1)
  % fn::take_first(20)
  % fn::for_each([&](const auto& kv)
    {
        std::cerr << +kv.first << "\t" << kv.second.len << "\t" << kv.second.score << "\n";
    });

    std::cerr << "\n\n";
}
#endif

struct track_t
{
    taxdiv_oid_t div_oid = taxdiv_oid_t{};
    ivls_t ivls = {};
    size_t cvg_len = 0ul; // sum of interval lengths (raw coverage)
};
using tracks_t = std::vector<track_t>;


/////////////////////////////////////////////////////////////////////////////
// GP-34068, GP-34188
// For each tax-div will keep track of coverage length and intersection with top-div.
// Will use this information in the end to compute a set of primary divs (top-div + highly-overlapping)
class div_stats_t
{
public:
    struct div_stat_t
    {
        taxdiv_oid_t div_oid = {};
        std::string gx_taxdiv;

        size_t cvg_len    = 0ul;
        size_t cvg_len_nr = 0ul;           // excluding repeats
        size_t top_div_isect_len_nr = 0ul; // excluding repeats

        double get_frac() const
        {
            return (float)(top_div_isect_len_nr)/((float)cvg_len_nr + 1);
        }
    };

    using div_stats_vec_t = std::vector<div_stat_t>;

    div_stats_vec_t m_stats = {};
    uint64_t m_genome_len = 0;
    uint64_t m_aggregate_cvg_len = 0;

    div_stats_t(const tax_map_t& tax_map)
    {
        m_stats.resize(256);
        uint8_t i = 0ul;
        for (auto& x : m_stats) {
            x.div_oid = taxdiv_oid_t{ i++ };
            x.gx_taxdiv = fn::first_or_default(tax_map, L(_.second.taxdiv_oid == x.div_oid))
                                                           .second.gx_taxdiv;
        }
    }

    const div_stat_t& at(taxdiv_oid_t div_oid) const
    {
        return m_stats.at(+div_oid);
    }

    void update(const tracks_t& tracks, const ivls_t& all_repeats)
    {
        // First update the cvg_len, so that selecting top_div_track based on aggregate coverage takes it into account
        for (const auto& track : tracks) {
            m_stats[+track.div_oid].cvg_len += ivl_t::sum_lens(track.ivls);
        }

        const auto& top_div_track = // based on cvg_len accumulated so far
            *std::max_element(tracks.begin(), tracks.end(), BY(m_stats[+_.div_oid].cvg_len));

        const auto top_div_track_nr = ivl_t::subtract(top_div_track.ivls, all_repeats);

        for (const auto& track : tracks) {
            auto& dest = m_stats[+track.div_oid];
            const auto len = ivl_t::sum_lens(track.ivls);

            dest.cvg_len_nr += len - ivl_t::intersection_len(track.ivls, all_repeats);
            dest.top_div_isect_len_nr += ivl_t::intersection_len(track.ivls, top_div_track_nr);
        }
    }


    static std::string get_kingdom(const std::string& gx_taxdiv)
    {
        //VERIFY(gx_taxdiv == "NULL" || (gx_taxdiv.size() > 5 && gx_taxdiv[4] == ':'));
        return gx_taxdiv.substr(0, 4);
    };


    std::vector<taxdiv_oid_t> select_primary_divs() const
    {
        auto ret = std::vector<taxdiv_oid_t>{};

        // sort by decreasing intersection length with top-div
        auto iters = make_vec_of_iterators(m_stats)
                   % fn::sort_by(fn::by::decreasing L(_->top_div_isect_len_nr));

        const auto top_div_frac_cvg = (double)iters.front()->cvg_len / ((double)m_genome_len + 1);

        // If the top-div has high aggregate coverage, we don't need to select additional ones.
        // For low-coverage cases we need to select more divs. How to decide the min_frac() threshold?
        // The higher the aggregate coverage is, the higher the threshold should be.
        // For now setting min_frac to top_div_frac_cvg, capped at minimum of 0.3
        const auto min_frac = std::max(0.3, top_div_frac_cvg);

        size_t prev_isect_len = 0;
        for (const auto& it : iters)
            if (it->gx_taxdiv != "NULL"
                && (it == iters.front()
                    || ( it->cvg_len > 0
                         && get_kingdom(it->gx_taxdiv) == get_kingdom(iters.front()->gx_taxdiv)
                         && (it->get_frac() > min_frac)
                         && (it->top_div_isect_len_nr * 20 > iters.front()->cvg_len_nr
                             || (double)it->top_div_isect_len_nr > 0.9 * (double)prev_isect_len)
                       )
                   )
               )
        {
            ret.push_back(it->div_oid);
            prev_isect_len = it->top_div_isect_len_nr;
        }

        return ret;
    }

    void print(std::ostream& ostr, size_t top_n = 10) const
    {
        ostr << "\n#top-divs by aggregate coverage (excluding repeats)\n";
        ostr << "#\tis-primary\tcvg-len-nr\ttop-div-overlap-nr\tdiv\n";

        const auto primary_divs = select_primary_divs();

        make_vec_of_iterators(m_stats)
      % fn::unstable_sort_by(fn::by::decreasing L(_->cvg_len_nr))
      % fn::for_each([&, i = 0ul](const auto it) mutable
        {
            const auto& stat = *it;
            const bool is_primary = std::count(primary_divs.begin(), primary_divs.end(), stat.div_oid);
            if (stat.cvg_len > 0 && (i++ < top_n || is_primary)) {
                ostr << "\t" << (is_primary ? 'T' : 'F')
                     << "\t" << stat.cvg_len_nr
                     << "\t" << float(int(stat.get_frac() * 1000 + 0.5)) / 1000.0f
                     << "\t" << stat.gx_taxdiv
                     << "\n";
            }
        });
        ostr << "\n";
    }

    // e.g. {"agg-cvg": 0.056, "asserted-div": "anml:crustaceans", "inferred-primary-divs": ["anml:crustaceans", "anml:insects"]}
    std::string make_run_info_json(const std::string& asserted_div) const
    {
        const float agg_cvg = (float)this->m_aggregate_cvg_len/(float)m_genome_len;
        VERIFY(agg_cvg <= 1);
        auto ret = std::string{"{\"agg-cvg\": "} + std::to_string(agg_cvg);
        const auto is_asserted_div =
            L( _.gx_taxdiv == asserted_div || (asserted_div == "virs:viruses" && str::startswith(_.gx_taxdiv, "virs:"))); // GP-33387

        if (!str::contains(asserted_div, "unknown") && asserted_div != "") {
            if (!fn::exists_where(is_asserted_div)(m_stats) && !str::contains(asserted_div, "metagenome")) {
                std::cerr << "\nWarning: asserted div '" + asserted_div + "' is not represented in the output!\n\n";
            }
            ret += ", \"asserted-div\": \"" + asserted_div + "\"";
        }

        ret += ", \"inferred-primary-divs\": [";
        bool first = true;
        for (const taxdiv_oid_t primary_div : select_primary_divs()) {
            ret += (first ? "\"" : ", \"") + m_stats[+primary_div].gx_taxdiv + "\"";
            first = false;
        }

        ret += "]}";
        return ret;
    }
};


/////////////////////////////////////////////////////////////////////////////
static ivls_t make_primary_track(const tracks_t& tracks, const div_stats_t& div_stats)
{
    const auto primary_divs = div_stats.select_primary_divs();

    auto is_primary = [&](const track_t& track) -> bool
    {
        return primary_divs % fn::exists_where L(_ == track.div_oid);
    };

    // select top-track by coverage, rather than taking first from primary-divs above,
    // because if this query is predominantly contaminant (not in primary-divs),
    // we want to use just that; if it's in the set of primary divs then use whole set.
    const auto& top_track = *std::max_element(tracks.begin(), tracks.end(), BY(_.cvg_len));
    const bool top_track_is_primary = is_primary(top_track);

    auto ret = ivls_t{};
    for (const auto& track : tracks)
        if (track.div_oid == top_track.div_oid || (top_track_is_primary && is_primary(track)))
    {
        ret <<= track.ivls;
    }

    ivl_t::sort_and_merge(ret);
    return ret;
}

// Scan chimeric_track with sliding windw, capturing and merging high-coverage instances.
static ivls_t find_chimeric_intervals(
                       const ivls_t& chimeric_track,
                                 int window_size,
                                 int terminal_window_size,
                                 int min_pct_cvg,
                               len_t qry_len)
{
    auto ret = ivls_t{};
    size_t cvg_in_window = 0;

    auto capture_chimeric_ivl = [&](ivls_t::const_iterator first, ivls_t::const_iterator last)
    {
        VERIFY(first <= last && last < chimeric_track.cend());

        VERIFY(first->endpos() < last->pos || first == last);

        const auto ivl = ivl_t::collapse(*first, *last);
        const auto pct_cvg = 100.0 * double(cvg_in_window) / (1.0 + ivl.len);

        //std::cerr << ivl.to_string() << "; pct_cvg=" << pct_cvg << "; n=" << (last-first+1) << "\n";

        VERIFY(0 <= pct_cvg && pct_cvg < 100.0);

        // range can be less than window_size in case of terminal cases where
        // the chimeric interval is close to sequnece ends; if so, require
        // that the interval must be at least 25% of window_sivze
        const bool ok = pct_cvg < min_pct_cvg            ? false
                      : ivl.len >= window_size           ? true
                      : ivl.len < terminal_window_size   ? false
                      : ivl.endpos() < window_size       ? true
                      : ivl.pos + window_size > qry_len  ? true
                      :                                    false;
        if (ok) {
            ivl_t::push_or_merge(ret, ivl);
        }
    };

    // Scan track with window [it_tail ..it_head] (inclusive)
    auto it_tail = chimeric_track.cbegin();
    for (const auto it_head : irange{ chimeric_track.cbegin(), chimeric_track.cend() }) {
        cvg_in_window += it_head->len;

        // advance tail while keeping the window as small as possible, but above window_size
        while (it_tail < it_head && (it_head->endpos() - (it_tail+1)->pos > window_size)) {
            cvg_in_window -= it_tail->len;
            ++it_tail;
        }

        capture_chimeric_ivl(it_tail, it_head);
    }

    if (!chimeric_track.empty()) {
        while (it_tail+1 < chimeric_track.cend()) {
            cvg_in_window -= it_tail->len;
            ++it_tail;
            capture_chimeric_ivl(it_tail, chimeric_track.cend()-1);
        }
        VERIFY(cvg_in_window == (size_t)chimeric_track.rbegin()->len);
    }

    return ret;
}


static ivls_t find_div_intervals(
                 const nodes_t& all_nodes,
                  const ivls_t& all_repeats,
                    const len_t qry_len,
               const tax_map_t& tax_map,
                   div_stats_t& div_stats)
{
    // NB, when adjusting, test on GCA_009936405.1
    // which has tons of little chimeras from insects and bacteria
    // also JHUN02000516.1 and NT_187580.1 and NKLS02002111.1

    static const int s_chimeric_window          = get_env("GX_CHIMERIC_WINDOW"         , 2000);
    static const int s_chimeric_terminal_window = get_env("GX_CHIMERIC_TERMINAL_WINDOW", 500); // TODO: change to 250, GP-32590
    static const int s_chimeric_min_pct_cvg     = get_env("GX_CHIMERIC_MIN_PCT_CVG"    , 50);

    VERIFY(s_chimeric_window >= 0);
    if (s_chimeric_window == 0) {
        return ivls_t{}; // chimeric identification disabled
    }

    auto merged_all_track = ivls_t{};

    static_assert(sizeof(taxdiv_oid_t) == 1, ""); // up to 256 values
    auto tracks = std::vector<track_t>(256); // by taxdiv_oid

    // group nodes by tax-divs
    for (const auto& node : all_nodes) {

        ivl_t::push_or_merge(merged_all_track, node);

        const auto taxdiv_oid = tax_map.at(node.tax_id).taxdiv_oid;
        auto& dest = tracks.at(+taxdiv_oid);
        ivl_t::push_or_merge(dest.ivls, node);
    }

    for (const auto i : irange{ tracks.size() }) {
        auto& track = tracks[i];
        track.div_oid = taxdiv_oid_t{ (uint8_t)i };
        track.cvg_len = ivl_t::sum_lens(track.ivls);

        fn::for_each_adjacent([](const ivl_t& a, const ivl_t& b)
        {
            VERIFY(a.endpos() < b.pos);
        })(track.ivls);
    }

    div_stats.update(tracks, all_repeats);

    const ivls_t primary_track = make_primary_track(tracks, div_stats);
    const auto chimeric_track =
        ivl_t::subtract(
                ivl_t::subtract(
                    merged_all_track,
                    primary_track),
                all_repeats);

#if 0
    std::cerr << "Cvg by primary  :" << ivl_t::sum_lens(primary_track) << "\n";
    std::cerr << "Cvg by all      :" << ivl_t::sum_lens(merged_all_track) << "\n";
    std::cerr << "Cvg by chimeric :" << ivl_t::sum_lens(chimeric_track)
                                     << "; n=" << chimeric_track.size() << "\n";

    if (chimeric_track.size() < 100)
        for (const auto& ivl : chimeric_track)
    {
        std::cerr << ivl.to_string() << " ";
    }
    std::cerr << "\n";
#endif

    auto chimeric_ivls =
        find_chimeric_intervals( chimeric_track,
                               s_chimeric_window,
                               s_chimeric_terminal_window,
                               s_chimeric_min_pct_cvg,
                               qry_len);

    // Keep those where coverage by a single div is high.
    // e.g. Right end of KV700477.1~126735..168854 has hits from across many divs,
    // where div-specific coverages are low, but aggregate coverage is high -
    // we don't want to keep this interval.
    chimeric_ivls %= fn::where([&](const ivl_t& chimeric_ivl)
    {
        len_t max_cvg = 0;
        for (const auto& track : tracks) {
            const auto overlapping_ivls = ivl_t::get_overlapping_v(track.ivls, chimeric_ivl);
            const auto cvg = sum_by(overlapping_ivls, L(ivl_t::intersect(chimeric_ivl, _).len));
            max_cvg = std::max(max_cvg, cvg);
        }
        return max_cvg * 100 > chimeric_ivl.len * s_chimeric_min_pct_cvg;
    });


    if (chimeric_ivls.empty()) {
        return chimeric_ivls;
    }

    // We don't want to merge chimeric ivls that are too far away, because
    // this will complicate the downstream assessment (i.e. lower coverage
    // within the interval, and potentially coverage coming from primary-div).

    const auto len_thr = s_chimeric_window / 4;
    ivl_t::sort_and_merge(chimeric_ivls, len_thr);

    // Merge adjacent (not necessarily abutting) chimeric intervals if
    // sum_of_lengths/total_range is high (i.e. gap is *relatively* small)
    fn::for_each_adjacent([](ivl_t& a, ivl_t& b)
    {
        auto ab = ivl_t::collapse(a, b);
        if (a.len + b.len > ab.len * 0.9) {
            a = ivl_t{};
            b = ab;
        }
    })(chimeric_ivls);
    chimeric_ivls %= fn::where L(_.len > 0);

    // Extend to sequence boundaries if within len_thr
    auto& front = chimeric_ivls.front();
    auto& back  = chimeric_ivls.back();

    if (pos1_t len_l = front.pos - 1; len_l < len_thr) {
        chimeric_ivls.front().len += len_l;
        chimeric_ivls.front().pos -= len_l;
    }

    if (const pos1_t len_r = qry_len + 1 - back.endpos(); len_r < len_thr) {
        back.len += len_r;
    }

    auto primary_ivls = ivl_t::subtract(ivls_t{{ 1, qry_len}}, chimeric_ivls);

    return std::move(primary_ivls)
         % fn::append(std::move(chimeric_ivls))
         % fn::where L(_.len > 0)
         % fn::sort();
}

/////////////////////////////////////////////////////////////////////////////

template<typename T>
class counter_t : public std::map<T, size_t>
{
public:
    void add(const T& key)
    {
        ++((*this)[key]);
    }

    void remove(const T& key)
    {
        const auto it = this->find(key);
        if (it == this->end()) {
            ;
        } else if (it->second == 1) {
            this->erase(it);
        } else {
            --it->second;
        }
    }
};

/////////////////////////////////////////////////////////////////////////////

using nodes_view_t = fn::view<nodes_t::const_iterator>;

static ivls_t find_ultraconserved_intervals(
                         const size_t qry_len,
                       const nodes_t& nodes,
                     const tax_map_t& tax_map,
                               size_t num_divs = 5)
{
    const len_t stride = 20;

    // count of set bits will correspond to count of unique divs from hits
    // overlapping the window. (NB: we have 55 divs now, so can use unt64_t).
    auto track = std::vector<uint64_t>(1 + qry_len/stride);

    for (const auto& node : nodes)
        if (const auto& ti = tax_map.at(node.tax_id); (ASSERT(+ti.taxdiv_oid < 64), !ti.is_prok_or_virus()))
            for(auto pos = node.pos + stride/2; pos < node.endpos(); pos += stride)
    {
        track.at(pos / stride) |= (1ul << (uint64_t)ti.taxdiv_oid); // NB: using 1ul instead of 1 here because triggered ubsan,
    }                                                               // but shouldn't it have been promoted based on rhs type?

    auto ret = ivls_t{};
    for (auto i : irange{ 1ul, track.size() })
        if ((uint8_t)_mm_popcnt_u64(track[i]) >= num_divs)
    {
        ivl_t::push_or_merge(ret, ivl_t{pos1_t(i * stride - stride/2), stride}, stride);
    }
    return ret;
}

/////////////////////////////////////////////////////////////////////////////
// NB: 0 can be valid length (empty-fasta seq) - GP-33468
static constexpr auto k_invalid_len = len_t(-1);

static void process_ivl(
              const ivl_t div_ivl, // div-specific (chimera-free) interval on query coordinates
       const nodes_view_t nodes,   // sorted by pos
      const seq_id_str_t& qry_id,
              const len_t qry_len,
            const ivls_t& collapsed_all_track,
   const repeat_tracks_t& repeat_tracks,
            const ivls_t& xtrachr_track,
         const tax_map_t& tax_map,
            std::ostream& ostr)
{
    // for accumulating tax-specific score and coverage
    struct data_t
    {
        tax_id_t tax_id = {};
           float score = 0;
           float adj_score = 0;
        uint32_t len = 0;
    };
    auto data = std::unordered_map<tax_id_t, data_t>{};

    // for accumulating taxdiv-specific coverage
    struct taxdiv_data_t
    {
        ivl_t ivl = {};
        len_t len = 0;

        void add(const ivl_t ivl2)
        {
             this->len += ivl2.len - ivl_t::uoverlap(this->ivl, ivl2);
             this->ivl = this->ivl.len == 0 ? ivl2 : ivl_t::collapse(this->ivl, ivl2);
        }
    };
    auto taxdiv_cvg = std::unordered_map<taxdiv_oid_t, taxdiv_data_t>{};

    // accumulate data and taxdiv_cvg
    for (const auto& n : nodes)
        if (const auto ivl = ivl_t::intersect(n, div_ivl); ivl.len > 0)
    {
        {
            auto& dest  = data[n.tax_id];
            dest.tax_id = n.tax_id;
            dest.len   += ivl.len;
            dest.score += n.seg_score * float(ivl.len) / float(n.len);
        }

        const auto taxdiv_oid = tax_map.at(n.tax_id).taxdiv_oid;
        taxdiv_cvg[taxdiv_oid].add(ivl);
    }

    const auto get_overlap_len = [&](const ivls_t& ivls)
    {
        auto ret = len_t{0};
        for (const auto& ivl : ivl_t::get_overlapping_v(ivls, div_ivl)) {
            ret += ivl_t::uoverlap(ivl, div_ivl);
        }
        return ret;
    };

    ostr << qry_id;

    // query-interval suffix unless whole-query.
    // NB, to differentiate from intervals due to split-fasta upstream,
    // we use ~~ instead of ~.
    if (!(div_ivl.pos == 1 && div_ivl.len == qry_len)) {
        ostr << "~~" << div_ivl.pos << ".." << div_ivl.endpos()-1;
    }

    ostr << "\t" << div_ivl.len
         << "\t" << get_overlap_len(repeat_tracks.transposons)
         << ","  << get_overlap_len(repeat_tracks.low_complexity)
         << ","  << get_overlap_len(repeat_tracks.ultraconserved)
         << ","  << get_overlap_len(repeat_tracks.n_runs)
         << ","  << get_overlap_len(xtrachr_track)
         << "\t" << get_overlap_len(collapsed_all_track)
         << "\t|";

    static const size_t num_taxa_to_print = get_env("GX_TAXIFY_NUM_TAXA", 4ul);

    // Take values from data map;
    // Sort by decreasing score, giving boost to "synthetic" and "virus" divs.
    // Pick up to top-2 per tax-div.
    // Resize to num_taxa_to_print.
    // Format to ostr.

    // Set adj_score, giving boost to viruses and synthetic divs.
    for (auto& kv : data) {
        auto& d = kv.second;
        const auto& tax = tax_map.at(d.tax_id);

        // NB: the following check works for old and new gx-div names (GP-33245)
        const auto priority_factor = str::startswith(tax.gx_taxdiv, "synt") ? 1.5f
                                   : str::contains(tax.gx_taxdiv, "virus")  ? 1.5f
                                   :                                            1.0f;
        d.adj_score = d.score * priority_factor;
    }

    fn::cfrom(data)
  % fn::transform L(_.second)
  % fn::sort_by L(_.adj_score * -1)
  % fn::where([&, counts = std::map<taxdiv_oid_t, size_t>()](const data_t& tax) mutable
    {
        VERIFY(+tax.tax_id);
        const auto taxdiv_oid = tax_map.at(tax.tax_id).taxdiv_oid;
        return taxdiv_oid != taxdiv_oid_t::N_run && ++counts[taxdiv_oid] <= 2;
    })

#if 0
   // Experimental:
   // Take top-N by score*len.
   // Note: need to take len into account:
   // e.g. JXIK01000309.1 - the highest coverage is to bats at <2.5kb
   // but have high-scoring non-specific hits at >2.5kb
  % fn::sort_by L(_.score * (float)_.len * -1)
  % L((_.resize(num_taxa_to_print), std::move(_)))
  % fn::sort_by L(_.score * -1) // resort by descending score for output
#else
  % fn::sort_by L(_.score * -1)
  % L((_.resize(num_taxa_to_print), std::move(_)))
#endif

  % fn::for_each([&, num_printed_taxa = 0ul](const data_t tax) mutable
    {
        if (num_printed_taxa == 0) {
            ostr << "\t" << (!+tax.tax_id ? "" : tax_map.at(tax.tax_id).species);
        }

        if (!+tax.tax_id) {
            ostr << "\t\t\t\t\t\t|";
        } else {
            const auto taxdiv_oid = tax_map.at(tax.tax_id).taxdiv_oid;
            ostr << "\t" << +tax.tax_id
                 << "\t" << tax_map.at(tax.tax_id).gx_taxdiv
                 << "\t" << taxdiv_cvg.at(taxdiv_oid).len // cvg-by-div
                 << "\t" << tax.len
                 << "\t" << int(sqrt(tax.score) + 0.5)
                 << "\t|";
        }
        ++num_printed_taxa;
    });

    ostr << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

static void process_qry(
                     nodes_t& nodes,
          const seq_id_str_t& qry_id,
                  const len_t qry_len,
             repeat_tracks_t& repeat_tracks,
                      ivls_t& xtrachr_track,
             const tax_map_t& tax_map,
                 div_stats_t& div_stats,
                std::ostream& ostr)
{
    const auto execption_guard = make_exception_scope_guard([&]
    {
        std::cerr << "NB: " << GX_SOURCE_LOCATION_STR << " while processing " << qry_id << "...\n";
    });

    if (qry_id.empty()) {
        return;
    }

    VERIFY(qry_len != k_invalid_len);

    // print_aggregate_scores_by_tax_ids(nodes);

    // TODO: sorting is a bottleneck. Instead:
    // fill std::unordered_map<tax_id_t, nodes_t>;
    // sort by pos and collapse, in parallel
    // n-way merge:
    //   move into min-heap, ordered by size
    // pop two smallest, merge, push result into min-heap
    // (Note, can't use priority-queue as it does not support move-semantics;
    // top() is const, and pop() does not return the value)

    // merge by tax-id, collapsing over gaps less than 10bp
    nodes %= fn::unstable_sort_by L(std::make_pair(_.tax_id, _.pos));
    nodes %  fn::for_each_adjacent([](node_t& a, node_t& b)
    {
        if (a.tax_id == b.tax_id && ivl_t::sdist(a, b) < 10) {
            b = node_t::collapse(a, b);
            a.len = 0;
        }
    });
    nodes %= fn::where L(_.len > 0);
    nodes %= fn::unstable_sort_by L(_.pos);

    const ivls_t collapsed_all_track = [&]
    {
        auto ret = ivls_t{};
        for (const auto& node : nodes) {
            ivl_t::push_or_merge(ret, node);
        }
        return ret;
    }();

    repeat_tracks.ultraconserved = find_ultraconserved_intervals(qry_len, nodes, tax_map);
    repeat_tracks.finalize();

    ivl_t::sort_and_merge(xtrachr_track);

    // print_aggregate_scores_by_tax_ids(nodes);

    div_stats.m_genome_len += qry_len;
    div_stats.m_aggregate_cvg_len += ivl_t::sum_lens(collapsed_all_track);

    const ivls_t div_ivls = find_div_intervals(nodes, repeat_tracks.combined, qry_len, tax_map, div_stats);

    if (div_ivls.size() <= 1) { // no hits (div_ivls is empty), or no chimeras (div_ivls.size() == 1)
        process_ivl(
            ivl_t{1, qry_len},
            fn::cfrom(nodes),
            qry_id,
            qry_len,
            collapsed_all_track,
            repeat_tracks,
            xtrachr_track,
            tax_map,
            ostr);
        return;
    }

    // for each div_ivl we'll locate the span (pair of iterators into nodes)
    // of overlapping nodes, and call process_ivls on that. Note that if we
    // encounter a node tha overlaps next downstream, div_ivl, that will be
    // the beginning of the next span.
    auto next_beg = nodes.cbegin();

    for (const auto i : irange{ div_ivls.size() }) {
        const ivl_t& div_ivl = div_ivls[i];
        const ivl_t* next_div_ivl_ptr = i + 1 < div_ivls.size() ? &div_ivls[i+1] : nullptr;

        const auto b = next_beg;
        auto e = b; // advance e while overlapping div_ivl
        while (e != nodes.end() && e->pos < div_ivl.endpos()) {

            if (next_beg == b     // not yet updated next_beg in this loop
               && next_div_ivl_ptr  // have next div-ivl
               && e->endpos() > next_div_ivl_ptr->pos) // overlapping next div-ivl
            {
                next_beg = e;
            }
            ++e;
        }
        if (next_beg == b) {
            next_beg = e;
        }

        process_ivl(
            div_ivl,
            fn::from(b, e),
            qry_id,
            qry_len,
            collapsed_all_track,
            repeat_tracks,
            xtrachr_track,
            tax_map,
            ostr
        );
    }
}

/////////////////////////////////////////////////////////////////////////////
namespace GP_31225
{
    struct region_t : public ivl_t
    {
        seq_id_str_t sbj_id = {};
        len_t total_uncovered_length = 0; // length of all gaps
        float score = 0.0;

        static region_t collapse(const region_t& a, region_t b)
        {
            VERIFY(a.sbj_id == b.sbj_id);
            VERIFY(a.pos <= b.pos);

            auto combined_score = a.score + b.score;

            // If have overlap, proportionally subtract score-contribution of a or b, whichever is less.
            const auto sdist = ivl_t::sdist(a, b);
            const auto overlap = 0.0f - (float)sdist;
            if (overlap > 0) {
                const auto a_frac_score = a.score / float(a.len);
                const auto b_frac_score = b.score / float(b.len);
                combined_score -= overlap * std::min(a_frac_score, b_frac_score);
            }

            const auto combined = ivl_t::collapse(a, b);
            b.pos = combined.pos;
            b.len = combined.len;
            b.score = combined_score;
            b.total_uncovered_length = a.total_uncovered_length + std::max(sdist, 0);
            return b;
        }

        void to_stream(std::ostream& ostr) const
        {
            if (len == 0) {
                return; // empty
            }

            ostr <<         sbj_id
                 << "\t" << pos
                 << "\t" << endpos() - 1
                 << "\t" << len - total_uncovered_length
                 << "\t" << int(sqrt(score) + 0.5)
                 << "\n";
        }
    };

    // input is GX-hits TSV
    static void merge_subject_intervals(std::istream& istr, std::ostream& ostr)
    {
        namespace tsv = rangeless::tsv;

        //auto seen_sbj_ids = std::set<seq_id_str_t>{};
        auto current_rgn = region_t{};

        ConsumeMetalineHeader(istr, GX_TSV_HEADER__HITS);
        //"#q-id\tq-pos1\ts-taxid\ts-id\ts-pos1\tlen\tmatchrun_L2\n";
        for (const tsv::row_t& row : tsv::from(istr)) {
            VERIFY(row.size() == 8);

          //const auto qry_seq_id =           seq_id_str_t(row[0]);
          //const auto qry_start  =   (pos1_t) tsv::to_num(row[1]);
          //const auto sbj_tax_id = (tax_id_t) tsv::to_num(row[2]);
            const auto sbj_seq_id =           seq_id_str_t(row[3]);
            const auto sbj_start  =   (pos1_t) tsv::to_num(row[4]);
            const auto len        =    (len_t) tsv::to_num(row[5]);
            const auto l2_score   = (uint64_t) tsv::to_num(row[6]);
          //const auto num_mismch =    (len_t) tsv::to_num(row[7]);

            // ignore synthetic hits _self _N_run _repeat, _transposon, etc
            if (sbj_seq_id.at(0) == '_') {
                continue;
            }

            VERIFY(sbj_start > 0); // expecting plus-anchored on sbj
            VERIFY(sbj_start != 0);

            // expecting sorted by sbj
            if (sbj_seq_id == current_rgn.sbj_id) {
                VERIFY(sbj_start >= current_rgn.pos);
            }

            const auto ivl = ivl_t{ sbj_start, len };

            auto score = std::min( (float)l2_score * (float)l2_score,
                                   (float)ivl.len * 20.0f); // See notes in Taxify about this
            auto rgn = region_t{ ivl, sbj_seq_id, 0, score };

            if (sbj_seq_id == current_rgn.sbj_id && current_rgn.endpos() + 1000 > sbj_start) {
                current_rgn = region_t::collapse(current_rgn, std::move(rgn));
            } else {
                current_rgn.to_stream(ostr);
                current_rgn = std::move(rgn);
            }
        }

        current_rgn.to_stream(ostr);
    }
}

std::string add_db_info(std::string metaline, std::string db_path);


// Add element to the metaline
// at the beginning of the file, presumed to have enough
// extra-whitespace padding to accomodate the extra info:
// e.g.
//
// metaline  : ##[["GX taxonomy pre-analysis report",3,1], {"git-rev":"0.2.3-95-g9be365c5-dirty", "run-date":"Wed Oct  5 15:06:46 2022", "db":{"build-date":"2022-07-09", "seqs":3086610, "Gbp":708.351}}]
// json_key  : run-info
// json-value: "run-info":{"agg-cvg": 0.782690, "inferred-primary-divs": ["prok:CFB group bacteria", "prok:bacteria", "prok:g-proteobacteria", "prok:d-proteobacteria"]}
// output : ##[["GX taxonomy pre-analysis report",3,1], {"git-rev":"0.2.3-95-g9be365c5-dirty", "run-date":"Wed Oct  5 15:06:46 2022", "db":{"build-date":"2022-07-09", "seqs":3086610, "Gbp":708.351}, "run-info": {"agg-cvg": 0.782690, "inferred-primary-divs": ["prok:CFB group bacteria", "prok:bacteria", "prok:g-proteobacteria", "prok:d-proteobacteria"]}}]
static void add_metaline_content(
                      std::string metaline,
               const std::string& out_path,
               const std::string& json_key,
               const std::string& json_value)
{
    VERIFY(!str::contains(metaline, json_key));

    // verify that the out-file begins with expected meta-line
    {
        auto ifstr = open_ifstream(out_path);
        VERIFY(ifstr);
        auto str = std::string( metaline.size(), ' ');
        ifstr.read(str.data(), str.size());
        VERIFY(str == metaline);
    }

    const auto orig_metaline_size = metaline.size();

    while (!metaline.empty() && metaline.back() == ' ') {
        metaline.pop_back();
    }

    VERIFY(str::endswith(metaline, "}]"));
    metaline.resize(metaline.size() - 2);
    metaline += ", \"" + json_key + "\": " + json_value + "}]";

    VERIFY(metaline.size() <= orig_metaline_size);
    metaline += std::string(orig_metaline_size - metaline.size(), ' ');
    // NB: adding the padding to keep the final size the same as original

    auto ofstr = std::ofstream(out_path, std::ios::in | std::ios::out | std::ios::binary);
    VERIFY(ofstr);
    ofstr.write(metaline.data(), metaline.size());
}


void gx::Taxify(
             std::istream& istr,
             std::istream& taxa_istr,
        const std::string& hardmask_locs_path,  // GP-34552
        const std::string& asserted_div,        // e.g. "proks:firmicutes"
        const std::string& db_path,
        const std::string& out_path)
{
    const bool to_file = out_path != "stdout" && out_path != "/dev/stdout" && out_path != "-";
    std::ofstream ofstr = to_file ? std::ofstream{ out_path } : std::ofstream{};
    std::ostream& ostr = to_file ? ofstr : std::cout;
    VERIFY(ostr);

    if (get_env("GX_COLLAPSE_SBJ_REGIONS", false)) {
        GP_31225::merge_subject_intervals(istr, ostr);
        return;
    }

    const tax_map_t tax_map = LoadTaxa(&taxa_istr);

    auto div_stats = div_stats_t{ tax_map };

    namespace tsv = rangeless::tsv;

    auto current_qry_id    = seq_id_str_t{};
    auto qry_len           = len_t{ k_invalid_len }; // NB: not 0, GP-33468
    auto nodes             = nodes_t{};
    auto repeat_tracks     = repeat_tracks_t{};
    auto xtrachr_track     = ivls_t{}; // extrachromosomal: plastids, plasmids, mito
    auto seen_queries      = std::set<seq_id_str_t>{};

    auto process_and_reset = [&](const seq_id_str_t& next_seq_id)
    {
        process_qry(nodes, current_qry_id, qry_len, repeat_tracks, xtrachr_track, tax_map, div_stats, ostr);

        current_qry_id = next_seq_id;
        qry_len = k_invalid_len;
        nodes.clear();
        repeat_tracks.clear();
        xtrachr_track.clear();

        if (!seen_queries.insert(next_seq_id).second) {
            GX_THROW("Duplicate seq_id in in the input of `gx taxify`: " + next_seq_id);
        }
    };

    // extra-whitespace at the end of the header to add run-info in the end
    static const size_t metaline_header_reserved_size = 
        get_env("GX_METALINE_HEADER_RESERVED_SIZE", 1024ul);

    const std::string metaline =
        add_db_info(MakeMetaLine(GX_TSV_HEADER__PRE_TAXONOMY_RPT), db_path)
        + std::string(metaline_header_reserved_size, ' '); 

    // (xp,lc,co,n,xc)-len:
    //      xp: transposon-specific
    //      lc: low-complexity
    //      co: conserved across many divs
    //      n:  N-runs
    //      xc: extrachromosomal (plastid|plasmid|mito)
    ostr << metaline
         <<  "\n#seq-id\tseq-len\t(xp,lc,co,n,xc)-len\tcvg-by-all\tsep1"
             "\ttax-name-1\ttax-id-1\tdiv-1\tcvg-by-div-1\tcvg-by-tax-1\tscore-1"
                   "\tsep2\ttax-id-2\tdiv-2\tcvg-by-div-2\tcvg-by-tax-2\tscore-2"
                   "\tsep3\ttax-id-3\tdiv-3\tcvg-by-div-3\tcvg-by-tax-3\tscore-3"
                   "\tsep4\ttax-id-4\tdiv-4\tcvg-by-div-4\tcvg-by-tax-4\tscore-4\tsep5\n";

    // GP-32758
    const std::set<tax_id_t> exclude_tax_ids =
        get_env("GX_EXCLUDE_TAXA", std::string{})
      % tsv::split_on_delim(',')
      % fn::where L(_ != "")
      % fn::transform L(tax_id_t{ tsv::to_num(_) })
      % fn::to(std::set<tax_id_t>());

    // GP-34226
    // ######################################################################
    const locs_map_t exclude_locs = [&]
    {
        auto istr_opt = open_ifstream_opt(hardmask_locs_path);
        auto ret = istr_opt ? LoadLocsMap(*istr_opt) : locs_map_t{};
        for (const auto& kv : ret) {
            ivl_t::verify_regular(kv.second);
        }
        return ret;
    }();

    // true iff interval overlaps any of the exclude_locs by at least 1/3 of its length
    const auto is_on_excludelist = [&exclude_locs](const seq_id_str_t& seq_id, ivl_t ivl)
    {
        ivl.pos = ivl.pos < 0 ? flip_pos1(ivl.pos, ivl.len) : ivl.pos;
        const auto it = exclude_locs.find(seq_id);
        return it != exclude_locs.end()
            && (    
                    ivl_t::get_overlapping_v(it->second, ivl) 
                  % fn::exists_where L(ivl_t::uoverlap(_, ivl) * 3 > ivl.len)
               );
    };

    const locs_map_t extrachromosomal_locs = [&]
    {
        const auto path = str::replace_suffix(db_path, ".gxi", ".extrachromosomal.tsv");
        std::ifstream ifstr{ path }; // only present in newer dbs.
        return ifstr.good() ? LoadLocsMap(ifstr) : locs_map_t{};
    }();

    // ######################################################################


    ConsumeMetalineHeader(istr, GX_TSV_HEADER__HITS);
    //"#q-id\tq-pos1\ts-taxid\ts-id\ts-pos1\tlen\tmatchrun_L2\n";
    size_t row_num = 0;
    for (const tsv::row_t& row : tsv::from(istr)) {
        ++row_num;

        const auto execption_guard = make_exception_scope_guard([&]
        {
            std::cerr << "Exception on row: " << row_num << "\n---- fields:\n";
            for (const auto & field : row) {
                std::cerr << field << "\n";
            }
            std::cerr << "----\n";
        });

        VERIFY(row.size() == 8);

        const auto qry_seq_id =           seq_id_str_t(row[0]);
        const auto qry_start  =   (pos1_t) tsv::to_num(row[1]);
        const auto sbj_tax_id = (tax_id_t) tsv::to_num(row[2]);
        const auto sbj_seq_id =           seq_id_str_t(row[3]);
        const auto sbj_start  =   (pos1_t) tsv::to_num(row[4]);
        const auto len        =    (len_t) tsv::to_num(row[5]);
        const auto l2_score   = (uint64_t) tsv::to_num(row[6]);
      //const auto num_mismch =    (len_t) tsv::to_num(row[7]);

        VERIFY(qry_start > 0); // plus-anchored on query
        VERIFY(sbj_start != 0);

        if (qry_seq_id != current_qry_id) { // reached next query
            process_and_reset(qry_seq_id);
        }

        const auto ivl     = ivl_t{ qry_start, len };
        const auto sbj_ivl = ivl_t{ sbj_start, len };

         // self-alignment reporting query-length
        if (sbj_seq_id == "_self") {

            VERIFY(qry_start == 1 && sbj_start == 1 && qry_len == k_invalid_len);
            qry_len = len;

        } else if (sbj_seq_id == "_N_run") {

            repeat_tracks.n_runs.push_back(ivl);

        } else if (sbj_seq_id == "_transposon") {

            repeat_tracks.transposons.push_back(ivl);

        } else if (sbj_seq_id == "_low_complexity") {

            repeat_tracks.low_complexity.push_back(ivl);

        } else if (exclude_tax_ids.count(sbj_tax_id) || is_on_excludelist(sbj_seq_id, sbj_ivl)) {
            ;
        } else if (!tax_map.count(sbj_tax_id)) {

            static thread_local std::set<tax_id_t> reported;
            if (reported.insert(sbj_tax_id).second && reported.size() <= 10) {
                std::cerr << "Warning: Skipping inputs for tax-id:" << +sbj_tax_id << " (not in taxa-file).\n";
                if (reported.size() == 10) {
                    std::cerr << "Will not report further warnings of this type.\n";
                }
            }

        } else {
            VERIFY(sbj_seq_id.at(0) != '_');
            VERIFY(+sbj_tax_id);

            auto score = std::min( (float)l2_score * (float)l2_score,
                                   (float)ivl.len * 50.0f);
            // L2-norms are not additive, but their squares are, so we have to square the scores
            // as we'll be summing them.
            //
            // We cap it at len*50 for the following reason: when the identity is high enough (e.g 98%)
            // it's no longer important whether it is 98% or 99.5% for our purposes - in both cases
            // it is indicative of same or close species, so anything above 98% we'll consider "perfect".
            //
            // This corresponds to seg composed of 50-bp matchruns, len/50 of them,
            // so the sum-of-squares score is 50^2*len/50 = len*50.
            //
            // We do it to guard against preferring higher-identity hits that are really contaminant-in-reference,
            // and instead be more sensitive to high-coverage hits in correct species, even though
            // they may have lower, but high-enough identity. (GP-31561)
            //
            // NB: used to use the value of 20 instead of 50, but changing it to the latter to make
            // the high identity hits to human stand out more (See GP-32031 comment 11/19/21 08:37 AM)

            nodes.push_back(node_t{ ivl, sbj_tax_id, score });


            if (   extrachromosomal_locs.count(sbj_seq_id)
                && tax_map.at(sbj_tax_id).gx_taxdiv == asserted_div)
            {
                xtrachr_track.push_back(ivl);
            }
        }
    }

    process_and_reset(seq_id_str_t{});


    const bool verbose = get_env("GX_TAXIFY_VERBOSE", false);
    if (verbose) {
        div_stats.print(std::cerr);
    }

    if (to_file) {
        ofstr.close();
        add_metaline_content(metaline, out_path, "run-info", div_stats.make_run_info_json(asserted_div));
    } else if(verbose) {
        std::cerr << "\nNote: the output-stream of taxify is stdout - will not add run-info to the meta-line header.\n"
                  << "run-info:" << div_stats.make_run_info_json(asserted_div) << "\n\n";
    }
}
