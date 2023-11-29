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
#define RANGELESS_FN_ENABLE_PARALLEL 1
#define RANGELESS_ENABLE_TSV 1

#include "align.hpp"
#include "types.hpp"
#include "serial_util.hpp"
#include "index.hpp"
#include "seq_info.hpp"
#include "segment.hpp"
#include "tmasker.hpp"

#include "subbyte_array.hpp"
#include "small_index.hpp"
#include "ext/thread_pool.hpp"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <set>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;

namespace tsv = rangeless::tsv;

// return geometric mean (ignoring values <= 0)
template<typename Xs, typename F>
double gmean_by(const Xs& xs, const F& fn)
{
    size_t n = 0;
    double sum_logs = 0;
    for (const auto& x : xs) {
        auto y = fn(x);
        if (y > 0) {
            ++n;
            sum_logs += log(1 + (double)y);
        }
    }
    return n == 0 ? 0 : exp(sum_logs / (double)n) - 1;
}


/////////////////////////////////////////////////////////////////////////////
static ivls_t get_low_complexity_mask_impl(
                    const iupacna_seq_t& seq,
                            const size_t window_size,
                            const double min_entropy, // for hexamers
                                  ivls_t ret = {})
{
    if (min_entropy <= 0 || seq.empty()) {
        return ret;
    }

    // hexamer -> count-in-window
    static thread_local std::vector<int32_t> counts; // signed, because will test whether >0
    counts.clear();
    counts.resize(4096); // 4^6 12-bit hexamers

    // circular buffer of last window_size hexamers
    static thread_local std::vector<uint32_t> hexamers;
    hexamers.clear();
    hexamers.resize(window_size);

    double entropy = 0; // Shannon entropy of distribution of hexamers in the window,
                        // -1 * sum p[i]*log2(p[i])

    const auto p_log2p = [&](int count) -> double
    {
        const auto p = (double)count / (double)window_size;
        return count == 0 ? 0 : p * log2(p);
    };

    const auto whole_ivl = ivl_t{ 1, (len_t)seq.size() };
    uint32_t head_hxm = 0;

    for (const auto i : irange{ std::max(seq.size(), window_size) }) {

        // To deal with queries shorter than window_size, we'll
        // consider them as concatenation of copies of the query
        // amounting to window_size.

        const auto na2 = (uint8_t)as_na2(na_t(seq[i % seq.size()]), i);

        head_hxm = ((head_hxm << 2) | na2) & Ob1x(12);

        const auto tail_hxm = hexamers[i % window_size];

        // Update the entropy-delta for the two hexamers, the one leaving the window
        // (whose count is decremented) and the one incoming (whose count is incremented).
        if (auto& count = counts[tail_hxm]; i >= window_size && count > 0) {
            entropy += p_log2p(count) - p_log2p(count - 1);
            --count;
        }

        if (auto& count = counts[head_hxm]; true) {
            entropy += p_log2p(count) - p_log2p(count + 1);
            ++count;
        }

        hexamers[i % window_size] = head_hxm;

        if (i >= window_size-1 && entropy < min_entropy) {
            auto ivl = ivl_t{ pos1_t(i + 2 - window_size), len_t(window_size) };
                                       //^ +2 to make 1-based and convert to endpos
            ivl_t::push_or_merge(ret, ivl_t::intersect(ivl, whole_ivl));
        }
    }

    // Erode intevals by few bp to better localize to low-complexity-only,
    // unless they extend all the way to sequence bounds.
    static const len_t erode_len = std::min((len_t)window_size / 3, get_env("GX_LOW_COMPLEXITY_ERODE_LEN", 30));
    for (auto& ivl : ret) {
        if (ivl.pos > 1 && ivl.len > erode_len) {
            ivl.pos += erode_len;
            ivl.len -= erode_len;
        }
        if (ivl.endpos()-1 < (int64_t)seq.size() && ivl.len > erode_len) {
            ivl.len -= erode_len;
        }
    }

    return ret;
}

/////////////////////////////////////////////////////////////////////////////
// GP-32021
// Most simple repeats have short periods e.g. {mono,di,tri,hexa}-mers, but some have long ones,
// e.g. example in the ticket has period of 34bp: TGTAGAGAAGGTCTGACACTGGAGGAGGACTATA
// These occur only twice per 100bp window, and so are not detected as repeats.
//
// To deal with this, do another pass with wider window and higher stringency.
// We can't just always use the wider window because we'll miss regions under 500bp,
// and so have to do two passes.
//
// This affects sensitivity, so instead of merging the results of the narrow and wide window runs,
// we only use get_low_complexity_mask_narrow results for masking, and augment the wide-window pass
// results for reporting purposes.

namespace gx
{

ivls_t get_low_complexity_mask_narrow(const iupacna_seq_t& seq); // NB: also used in make_db.cpp
ivls_t get_low_complexity_mask_narrow(const iupacna_seq_t& seq)
{
    static const size_t window_size = get_env("GX_LOW_COMPLEXITY_WINDOW_NARROW",        50ul);
    static const double min_entropy = get_env("GX_LOW_COMPLEXITY_MIN_ENTROPY_NARROW",   4.5); // higher=more masking
    return get_low_complexity_mask_impl(seq, window_size, min_entropy);
}

static ivls_t get_low_complexity_mask_wide(const iupacna_seq_t& seq)
{
    static const size_t window_size = get_env("GX_LOW_COMPLEXITY_WINDOW_WIDE",          500ul);
    static const double min_entropy = get_env("GX_LOW_COMPLEXITY_MIN_ENTROPY_WIDE",     6.5); // higher=more masking
    return get_low_complexity_mask_impl(seq, window_size, min_entropy);
}

}

/////////////////////////////////////////////////////////////////////////////

struct scores_t
{
    uint64_t matchruns_l2 = 0; // l2 norm of matchruns
    uint32_t num_mismatches = 0;
    uint32_t longest_matchrun = 0;

    scores_t() = default;

    scores_t(   const fasta_seq_t& qry_seq,
                  const sbj_seq_t& sbj_seq,
                  const segment_t& seg)
    {
        if (abs(seg.s) >= k_repeat_marker) {
            return; // special value - not a real seg
        }

        VERIFY(sbj_seq.size() > 0 || seg.len == 0);

        auto update = [&, last_mismatch_pos = int32_t(-1)](int32_t mm_pos) mutable
        {
            auto matchrun_len = uint32_t(mm_pos - last_mismatch_pos - 1); // 0 when adjacent
            matchruns_l2 += matchrun_len * matchrun_len;
            last_mismatch_pos = mm_pos;
            longest_matchrun = std::max(longest_matchrun, matchrun_len);
        };

        for (const int32_t i : irange{ seg.len }) // for each mismatch
            if (qry_seq.at1(seg.q + i) != sbj_seq.at1(seg.s + i))
        {
            num_mismatches++;
            update(i);
        }
        update(seg.len); // final matchurn

        matchruns_l2 = uint64_t(sqrt((double)matchruns_l2) + 0.5);
   }
};

/////////////////////////////////////////////////////////////////////////////

struct prelim_align_ret_t
{
    segments_t segs;
        ivls_t low_complexity;
        ivls_t transposons;
};

static prelim_align_ret_t prelim_align(const fasta_seq_t& qry, const CIndex& index, const CTmasker& tmasker)
{
    auto elapsed = timer{};
    const size_t qry_len = qry.seq.size();

    /////////////////////////////////////////////////////////////////////////
    // We need to do multiple passes over hitses, so we precompute these.
    using nodes_view_t = CIndex::nodes_view_t;
    static thread_local std::vector<nodes_view_t> hitses{};
    static thread_local std::vector<bool> is_flippeds{}; // whether i'th query-hmer is flipped
    {
        hitses.clear();
        hitses.resize(qry_len);

        is_flippeds.clear();
        is_flippeds.resize(qry_len);

        for (auto kmer = kmer_ci_t(qry.seq, CIndex::k_word_tlen); kmer; ++kmer) {
            const auto hmer     = CIndex::hmer38_t{ kmer.buf };
            hitses[kmer.i]      = index.at(hmer);
            is_flippeds[kmer.i] = hmer.is_flipped;
        }
    }

    static thread_local std::ostringstream ostr{};
    ostr.str("");
    ostr << "Chunk length:" << qry.seq.size() << "\n";
    ostr << "Initialized hitses in " << float(elapsed) << "s.\n";

    /////////////////////////////////////////////////////////////////////////
    // High occurrence intervals, having more than 10k hits in 20-bp window.
    // TODO: pass them downstream to the taxify stage?
    ivls_t high_occ_ivls = [&]
    {
        auto ret = ivls_t{};

        static const auto win_size           = get_env("GX_HI_OCC_WINDOW", 20ul);
        static const auto max_hits_in_window = get_env("GX_HI_OCC_HITS", 10000ul);
        VERIFY(win_size <= CIndex::k_word_tlen);

        size_t count_in_window = 0;
        for (const auto i : irange{ CIndex::k_word_tlen, hitses.size()} ) {
            count_in_window += hitses[i].size();
            count_in_window -= i < win_size ? 0 : hitses[i - win_size].size();

            if (count_in_window > max_hits_in_window) {
                const auto len = len_t(CIndex::k_word_tlen + win_size);
                const auto pos = std::max(1, pos1_t(i + 1) - len);
                ivl_t::push_or_merge(ret, ivl_t{ pos, len });
            }
        }

        return ret;
    }();

    ivls_t low_complexity_ivls = [&]
    {
        ivls_t ret = get_low_complexity_mask_narrow(qry.seq);
        ivl_t::sort_and_merge(ret, CIndex::k_word_tlen);
        ret %= fn::where L(_.len > CIndex::k_word_tlen); // drop singletons
        return ret;
    }();

    auto transposon_ivls = tmasker.find_repeats(qry.seq);


    static const bool enable_repeat_filtering = get_env("GX_ENABLE_REPEAT_FILTERING", true);
    if (enable_repeat_filtering)
        for (const auto ivls_ptr : { &transposon_ivls, &low_complexity_ivls, &high_occ_ivls })
            for (const auto& ivl : *ivls_ptr)
                for (const auto pos1 : irange{ ivl.pos + CIndex::k_word_tlen - 1,
                                               std::min(ivl.endpos(), pos1_t(hitses.size() + 1)) })
    {
        hitses.at(pos1 - 1) = nodes_view_t{};
    }

    /////////////////////////////////////////////////////////////////////////

    // NB: experimental; disabling by default for now.
    static const bool enable_augment_singletons  = get_env("GX_ENABLE_AUGMENT_SINGLETONS", false);
    static const bool enable_singleton_filtering = get_env("GX_ENABLE_SINGLETON_FILTERING", true);


    // Count of hits per subject seq-id-hash; used in enable-augment-singletons stage
    // std::unordered_map affects performance, so using just plain vector and disregarding collisions.
    static thread_local std::vector<uint32_t> h_oid2count{};
    h_oid2count.clear();
    h_oid2count.resize(4096);

    // Initialize segments from hits.
    //
    // NB: used to do this in multiple passes, for each stride-phase separately.
    // Doing it as a single pass has better sensitivity e.g. euglenoids in GCA_024586455.1
    segments_t segs{};
    {
        segs.clear();
        segs.reserve(qry_len);

        // Cap the size of the filter at 10x the query size as otherwise it's too much overhead for short queries.
        static thread_local CSingletonFilter filter(nullptr, 0UL);
        const size_t filter_size = std::min(100000UL, qry.seq.size() * 10);
        filter.reset(&segs, filter_size); // will push into filter, which will fill segs.

        // geometric mean of count of hits per position, excluding repeat-specific (size() == 0 after above)
        const double gmean_occ = gmean_by(hitses, L(_.size()));

        // to further filter hits that slip by repeat-filtering
        const auto max_occ_thr = (size_t)std::max(20.0, gmean_occ * 4);

        for (size_t i = 0; i < qry_len; ++i) {
            const auto q_pos = as_ivl(i, CIndex::k_word_tlen, is_flippeds[i]).pos;
            const nodes_view_t hits = hitses[i];

            for (const auto& h : hits)
                if (h.pos != k_repeat_marker && hits.size() < max_occ_thr)
            {
                auto seg = segment_t{ q_pos, h.pos, h.seq_oid, CIndex::k_word_tlen };
                seg.make_s_fwd();

                if (enable_augment_singletons) {
                    ++h_oid2count[uint64_hash(+seg.s_oid) % h_oid2count.size()];
                }

                if (enable_singleton_filtering) {
                    filter.push(seg);
                } else {
                    segs.push_back(seg);
                }
            }
        }
        filter.finalize();
    }


    /////////////////////////////////////////////////////////////////////

    // Need to do another separate round of coalescing because the filter
    // can't do it 100% due to segs being expunged from the hotlist prematurely.
    Coalesce(segs);
    DropShadowedOnSbj(segs);

    // Another round of singleton dropping singleton for the same reason
    // as Coalesce above. We do DropShadowedOnSbj first so if there's a
    // stack of repeats on subject between neighbors, we remove it so that
    // there's better chance of neighborship being detected.
    if (enable_singleton_filtering && segs.size() > 10) {
        DropSingletons(segs, CIndex::k_word_tlen);
    }

    // Do another pass over the initial hits.
    // For top-N subject seq-ids, by count of hits, augment the hits
    // that were dropped out by the filtering steps.
    // This picks up more relevant alignments for distant taxa, e.g.
    // euglenoids hits for JAMJQZ010001831.1; also significantly improves
    // coverage by insects for cockroach, and coverage by other fungal divs
    // for GCA_001574975.1 (monoblepharidomycetes).
    if (enable_augment_singletons && !segs.empty()) {
        const auto min_count = fn::first_or_default(fn::cfrom(h_oid2count) % fn::take_top_n(20));

        for (size_t i = 0; i < qry_len; ++i) {
            const auto q_pos = as_ivl(i, CIndex::k_word_tlen, is_flippeds[i]).pos;
            for (const auto& h : hitses[i])
                if (h_oid2count[uint64_hash(+h.seq_oid) % h_oid2count.size()] > min_count)
                // NB: >, not >=, as there may be arbitrarily many
            {
                segs.push_back(segment_t{ q_pos, h.pos, h.seq_oid, CIndex::k_word_tlen });
            }
        }

        Coalesce(segs);
    }


#if 0
    // For now this is disabled because it may leave a short query entirely unaligned, e.g.
    // NW_016681909.1 vs. NW_017548020.1
    //
    // Need to develop a better algorithm for filtering the repeats, that will
    // retain only the "best" repeat-specific segs.
    // e.g.
    /*
        1) For every point at k-bp stride (e.g. k=50, k < word-len) compute the count of overlapping
        segs and the longest overlapping seg.
        2) Identify repeat-specific segs (whose every overlapping k-point is covered with high depth)
        3) Accept segs that are not repeat-specific, OR is-longest for some k-point.
        4) Accept segs that are repeat-specific AND are syntenic neighbors to an already accepted seg.
    */

    FilterOutRepeatsOnSbj(segs, sbj_infos); // Does this need to be per-taxon?
    ostr << "After filtering repeats on sbj:" << segs.size() << "\n";
#endif

    /////////////////////////////////////////////////////////////////////////

    static const bool do_print = get_env("GX_VERBOSE", false);
    if (do_print) {
        static std::mutex mut;
        std::lock_guard<std::mutex> g{ mut };
        std::cerr << ostr.str(); // Note: this will mess with std::ostream if multithreading.
    }

    for (auto& seg : segs) {
        seg.make_q_fwd();
    }

    // Augment wide-window-low-complexity and high-occurrence intervals for reporting.
    low_complexity_ivls <<= get_low_complexity_mask_wide(qry.seq);
    low_complexity_ivls <<= std::as_const(high_occ_ivls);
    ivl_t::sort_and_merge(low_complexity_ivls);

    // Make transposons-track pure (subtract low-complexity)
    transposon_ivls = ivl_t::subtract(transposon_ivls, low_complexity_ivls);

    return { std::move(segs), std::move(low_complexity_ivls), std::move(transposon_ivls) };
}



//////////////////////////////////////////////////////////////////////////////////

static std::string format_output(
                     const fasta_seq_t& inp_chunk,
                     const seq_infos_t& sbj_infos,
    std::function<sbj_seq_t(seq_oid_t)> get_sbj_seq,
              const prelim_align_ret_t& prelim_align_result)
{
    static thread_local std::stringstream ostr;
    ostr.clear();
    ostr.str("");

    // When filtering-out self-alignments (query and subject ids are the same),
    // disregard ~ -suffix, e.g. NW_023312738.1~268840..270839
    const std::string no_suffix_seq_id = tsv::split_on_delim{'~'}(inp_chunk.seq_id)[0];

    static const bool s_report_self_alignments = get_env("GX_REPORT_SELF_ALIGNMENTS", false);
    static const float s_min_pct_idty = get_env("GX_MIN_PCT_IDENTITY", 40.0f);

    for (const auto& seg : prelim_align_result.segs) {

        const auto& sbj_info = sbj_infos.at(seg.s_oid);

        if (!s_report_self_alignments && no_suffix_seq_id == sbj_info.get_seq_id()) {
            continue;
        }

        const auto sbj_seq   = get_sbj_seq(seg.s_oid);
        const auto scores    = scores_t(inp_chunk, sbj_seq, seg);

        if ((float)scores.num_mismatches > (float)seg.len * (1.0 - (s_min_pct_idty / 100.0f))) {
             continue;
        }

        ostr         << inp_chunk.seq_id << "\t" << seg.q
             << "\t" << +sbj_info.tax_id
             << "\t" << sbj_info.get_seq_id() << "\t" << seg.s
             << "\t" << seg.len
             << "\t" << scores.matchruns_l2
             << "\t" << scores.num_mismatches
             << "\n";
    }

    // report N-runs
    {
        auto n_run = ivl_t{}; // in chunk-local coords
        const auto report_n_run = [&]
        {
            if (n_run.len >= 10) {
                ostr         << inp_chunk.seq_id << "\t" << size_t(n_run.pos) + inp_chunk.offset
                     << "\t" << 0  // tax-id
                     << "\t" << "_N_run" << "\t" << 1
                     << "\t" << n_run.len
                     << "\t" << 0  // matchruns_l2
                     << "\t" << 0  // mismatches
                     << "\n";
            }
        };

        for (const auto i : irange{ inp_chunk.seq.size() }) {
            const auto c = inp_chunk.seq[i];
            if (c != 'N' && c != 'n') {
                continue;
            } else if (n_run.endpos() == int32_t(i+1)) {
                n_run.len++; // extend n-run
                continue;
            }
            report_n_run();
            n_run.pos = int32_t(i+1);
            n_run.len = 1;
        }
        report_n_run();
    }

    // report repeats
    for (const auto ivls_p : { &prelim_align_result.low_complexity, &prelim_align_result.transposons })
        for (const ivl_t& ivl : *ivls_p)
    {
        ostr         << inp_chunk.seq_id << "\t" << size_t(ivl.pos) + inp_chunk.offset
             << "\t" << 0  // tax-id
             << "\t" << (ivls_p == &prelim_align_result.low_complexity ? "_low_complexity" : "_transposon")
             << "\t" << 1  // Theres's no subject-sequence, so the coordinate is a sentinel.
             << "\t" << ivl.len
             << "\t" << 0  // matchruns_l2
             << "\t" << 0  // mismatches
             << "\n";
    }

    return ostr.str();
}

/////////////////////////////////////////////////////////////////////////////

static void align_round2(
        const tax_map_t&                    taxa,
        const seq_infos_t&                  sbj_infos,
        const fasta_seq_t&                  inp_chunk,
        std::function<sbj_seq_t(seq_oid_t)> get_sbj_seq,
        segments_t&                         segs1)
{
    static const auto min_qry_len_for_round2 = get_env("GX_ROUND2_MIN_QRY_LEN", 300ul);
    const bool qry_too_short_for_round2 = inp_chunk.seq.size() < min_qry_len_for_round2;

    // 0 - disable; 1 - for select taxa; 2 - always
    static const auto enable_round2 = get_env("GX_ENABLE_ROUND2", 1U);
    if (!enable_round2 || segs1.empty() || qry_too_short_for_round2) {
        return;
    }

#if 0
    // if query is under 1kb and high-coverage - skip it too
    if (inp_chunk.seq.size() < min_qry_len_for_round2 * 2) {
        segs1 %= fn::unstable_sort_by L(_.q);
        const auto raw_cvg_len = get_raw_cvg_len(segs1);

        if(raw_cvg_len * 10 > (len_t)inp_chunk.seq.size() * 9) {
            return;
        }
    }
#endif

    const auto round2_taxa = SelectTaxaForRound2(segs1, sbj_infos, taxa);

    VERIFY(std::is_sorted(round2_taxa.begin(), round2_taxa.end()));

#if 0
    std::cerr << "Round-2 taxa: ";
    for (const auto t : round2_taxa) {
        std::cerr << " " << t;
    }
    std::cerr << "\n";
#endif

    static thread_local CSmallIndex query_index{};
    query_index = MakeQueryIndex(inp_chunk, std::move(query_index));

    // Orient and sort on subject
    for (auto& seg : segs1) {
        seg.make_s_fwd();
    }
    segs1 %= fn::unstable_sort_by(get_sbj);

    static thread_local segments_t segs2{};
    segs2.clear();
    for_each_group_by(segs1, L(_.s_oid), [&](segs_view_t sv)
    {
        const auto sbj_oid = sv.begin()->s_oid;
        const auto sbj_seq = get_sbj_seq(sbj_oid);
        const auto& sbj_inf = sbj_infos.at(sbj_oid);

        const bool taxon_in_scope_for_round2 =
            std::binary_search(round2_taxa.begin(), round2_taxa.end(), sbj_inf.tax_id);

        if (enable_round2 == 1 && !taxon_in_scope_for_round2) {
            return;
        }

        static thread_local segments_t tmp_segs{};
        tmp_segs.clear();
        tmp_segs = SeedRound2(query_index, inp_chunk, sbj_seq, sbj_oid, sv, std::move(tmp_segs));
        segs2 <<= std::as_const(tmp_segs);
    });

    if (!segs2.empty()) {
        segs1 <<= std::as_const(segs2);
        Coalesce(segs1, 5);
    }
}

/////////////////////////////////////////////////////////////////////////////
struct result_for_chunk_t
{
    fasta_seq_t qry_chunk;
     segments_t segs;
    std::string formatted_lines;
         double elapsed;
};
using results_for_qry_t = std::vector<result_for_chunk_t>;

static auto align_chunk( const tax_map_t& taxa,
                       const seq_infos_t& sbj_infos,
                            const CIndex& index,
                          const CTmasker& tmasker,
                        const char* const mmapped_seq_db_data,
                              fasta_seq_t inp_chunk) -> result_for_chunk_t
{
    const auto execption_guard = make_exception_scope_guard([&]
    {
        std::cerr << "NB: " << GX_SOURCE_LOCATION_STR << ": while processing "
                  << inp_chunk.defline << " @" << inp_chunk.offset << "...\n";
    });

    const auto get_sbj_seq = [&](seq_oid_t sbj_oid)
    {
        const seq_info_t& si = sbj_infos.at(sbj_oid);
        VERIFY(si.length > 0);
        return sbj_seq_t{ si, mmapped_seq_db_data };
    };

    const auto elapsed = timer{};

    auto prelim_align_result = prelim_align(inp_chunk, index, tmasker);
    auto& segs = prelim_align_result.segs;

    // adjust query position from query-chunk coordinates into query-coordinates
    for (auto& seg : segs) {
        seg.make_q_fwd();
        seg.q += (int32_t)inp_chunk.offset;
    }

    // GP-33150. Related: GP-32758
    // Filter out hits to unwanted taxa from dowstream processing.
    {
        static const std::set<tax_id_t> exclude_tax_ids =
            get_env("GX_ALIGN_EXCLUDE_TAXA", std::string{}) // comma-delimited list of ints
          % tsv::split_on_delim(',')
          % fn::where L(_ != "")
          % fn::transform L(tax_id_t{ tsv::to_num(_) })
          % fn::to(std::set<tax_id_t>());

        if(!exclude_tax_ids.empty()) {
            segs %= fn::where([&](const segment_t& seg)
            {
                const auto tax_id = sbj_infos.at(seg.s_oid).tax_id;
                return ! exclude_tax_ids.count(tax_id);
            });
        }
    }

    auto extend_and_filter = [&](bool gapped_extend)
    {
        // orient on query and sort by subject, as required by UngappedSegsInPlace
        for (auto& seg : segs) {
            seg.make_q_fwd();
        }
        segs %= fn::unstable_sort_by(get_sbj);

        for_each_group_by(segs, L(_.s_oid), [&](segs_view_t sv)
        {
            const auto sbj_oid = sv.begin()->s_oid;
            const auto sbj_seq = get_sbj_seq(sbj_oid);
            UngappedExtendSegsInPlace(inp_chunk, sbj_seq, sv);
        });

        segs %= fn::unstable_sort_by(get_sbj); // extensions may have invalidated the order

        if (gapped_extend) {
            GappedExtend(inp_chunk, get_sbj_seq, segs);
            segs %= fn::unstable_sort_by(get_sbj); // restore order including auxiliary segs from GappedExtend
        }

        VERIFY(std::is_sorted(segs.begin(), segs.end(), by_sbj));

        // Filter to best-chains. NB this invalidates the sort order
        for_each_group_by( segs,
                           L(sbj_infos.at(_.s_oid).tax_id), // group by taxon
                           L(FilterToBestChain(_, 1000ul/*max-backtrack*/)));
        segs %= fn::where L(_.flags != 0); // keep segs marked as belonging to best-chain

        for (auto& s : segs) {
            s.flags = 0;
            VERIFY(s.q > 0);
        }
    };

    extend_and_filter(false); // extend-and-filter to get better signal for round-2 taxa selection
    segs %= fn::where L(_.len > CIndex::k_word_tlen);
    align_round2(taxa, sbj_infos, inp_chunk, get_sbj_seq, segs);
    extend_and_filter(true);

    // Drop crud that made it through the singleton filter but couldn't be extended.
    segs %= fn::unstable_sort_by(get_sbj);
    for_each_group_by( segs,
                       L(sbj_infos.at(_.s_oid).tax_id), // group by taxon
                       [&](auto group)
    {
        const auto total_len = sum_by(group, L(_.len));

        // length of extensions (subtract the seed length; ignore short)
        const auto total_ext_len = sum_by(group, L(std::max(0, _.len - CIndex::k_word_tlen)));

        if (total_ext_len < 50 && (size_t)total_len * 10 < inp_chunk.seq.size())
            for(auto& seg : group)
        {
            seg.flags = 1;
        }
    });
    segs %= fn::where L(_.flags != 1);

    //////////////////////////////////////////////////////////////////////////////////

    std::string formatted_output = format_output(inp_chunk, sbj_infos, get_sbj_seq, prelim_align_result);
    return { std::move(inp_chunk), std::move(segs), std::move(formatted_output), elapsed };
}

/////////////////////////////////////////////////////////////////////////////

std::string add_db_info(std::string metaline, std::string db_path);

std::string add_db_info(std::string metaline, std::string db_path)
{
    // NB: for not using open_ifstr_opt, as for now existence of .meta.jsonl is optional
    errno = 0;
    std::ifstream ifstr{ str::replace_suffix(db_path, ".gxi", ".meta.jsonl") };
    if (!ifstr) {
        errno = 0;
        return metaline;
    }

    std::string line;
    std::getline(ifstr, line);
    VERIFY(!line.empty() && line.front() == '{' && line.back() == '}');

    VERIFY(!metaline.empty() && str::endswith(metaline, "}]"));
    metaline.insert(metaline.size() - 2, ", \"db\":" + line);
    return metaline;
}


void gx::ProcessQueries(  const std::string& db_path,
                          const std::string& taxa_path,          // can be empty
                          const std::string& repeats_fasta_path, // can be empty
                               std::istream& fasta_istr,
                               std::ostream& ostr)
{
    const std::string gxs_file_path = str::replace_suffix(db_path, ".gxi", ".gxs");
    const std::string_view mmapped_seq_db = ser::mmap(gxs_file_path);
    const std::string_view mmapped_db = ser::mmap(db_path);

    static const auto prefetch_mode = get_env("GX_PREFETCH", 1U); // 0=no; 1=autodetect; 2=yes
    if(prefetch_mode) {
        const bool force = prefetch_mode == 2;
        ser::prefetch_mmapped_pages(gxs_file_path, mmapped_seq_db, force);
        ser::prefetch_mmapped_pages(db_path,       mmapped_db,     force);
    }

    ser::memistream istr(mmapped_db.data(), mmapped_db.size());

    auto t = timer{};
    const auto sbj_infos = seq_infos_t(ser::from_stream(istr));
    const auto index     = CIndex(istr);
    const auto tmasker   = CTmasker(repeats_fasta_path);

    /////////////////////////////////////////////////////////////////////////

    const tax_map_t taxa = LoadTaxa(open_ifstream_opt(
                !taxa_path.empty() ? taxa_path
                                   : str::replace_suffix(db_path, ".gxi", ".taxa.tsv")).get());

    // Verify the size of mmapped_db (.gxi) and mmapped_seq_db (.gxs), e.g. that files are not truncated.
    {
        if ((size_t)istr.tellg() != mmapped_db.size()) {
            GX_THROW("Unexpected error while reading " + db_path
                   + "; file size is " + std::to_string(mmapped_db.size()) + " bytes.");
        }

        // Find seq-info with highest offset_in_seq_db (last in mmapped_seq_db)
        // NB: not necessarily last in sbj_infos.
        VERIFY(!sbj_infos.empty());
        const auto& last_sbj_info = *std::max_element(sbj_infos.begin(), sbj_infos.end(), BY(_.offset_in_seq_db));

        const auto last_sbj_seq  = sbj_seq_t{ last_sbj_info, mmapped_seq_db.data() };
        const auto expected_size = last_sbj_info.offset_in_seq_db + last_sbj_seq.data().size_bytes();

        if (expected_size != mmapped_seq_db.size()) {
            GX_THROW("Expected " + gxs_file_path + " of " + std::to_string(expected_size)
                              + " bytes. Actual size is " + std::to_string(mmapped_seq_db.size()) + ".");
        }
    }

    const auto align_chunk_wr = [&](fasta_seq_t chunk) -> result_for_chunk_t
    {
        return align_chunk(taxa, sbj_infos, index, tmasker, mmapped_seq_db.data(), std::move(chunk));
    };

    /////////////////////////////////////////////////////////////////////////

    size_t query_count = 0;
    size_t total_qry_bases = 0;
    auto output_result =
        [             &ostr
        ,      &query_count // for updating
        ,  &total_qry_bases
        ,           prev_id = seq_id_str_t()
        ,   partial_qry_len = 0UL
        , partial_qry_N_len = 0UL ](result_for_chunk_t res) mutable
    {
        total_qry_bases += res.qry_chunk.seq.size();

        const auto execption_guard = make_exception_scope_guard([&]
        {
            std::cerr << "NB: " << GX_SOURCE_LOCATION_STR << " while processing "
                      << res.qry_chunk.defline << " @" << res.qry_chunk.offset << "...\n";
        });

        // std::cerr << res.qry_chunk.seq_id << "@" << res.qry_chunk.offset << "\t" << res.elapsed << std::endl;

        // reached next query - print terminator row for the previous query reporting the whole-query self-alignment
        // for query-length reporting. (we only know length when we reach the last query chunk).
        if (res.qry_chunk.seq_id != prev_id && !prev_id.empty()) {
            ostr         << prev_id << "\t" << 1 // q-pos
                 << "\t" << 0
                 << "\t" << "_self" << "\t" << 1 // s-pos
                 << "\t" << partial_qry_len   // length
                 << "\t" << partial_qry_len   // L2
                 << "\t" << partial_qry_N_len // report conut of Ns in mismatches column
                 << "\n";
            ++query_count;
        }

        // reset query-len accumulators and prev_id upon reaching the new query
        if (res.qry_chunk.seq_id != prev_id) {
            partial_qry_len = 0;
            partial_qry_N_len = 0;
            prev_id = res.qry_chunk.seq_id;
        }

        // update qry_N_len and qry_len
        const auto overlap_size = int32_t(partial_qry_len - res.qry_chunk.offset);

        VERIFY(overlap_size >= 0);
        VERIFY(overlap_size <= (int32_t)res.qry_chunk.seq.size());
        VERIFY(overlap_size <= (int32_t)partial_qry_len);

        partial_qry_N_len += (size_t)std::count( res.qry_chunk.seq.begin() + overlap_size,
                                                 res.qry_chunk.seq.end(), 'N');
        partial_qry_len = res.qry_chunk.offset + res.qry_chunk.seq.size();

        ostr << res.formatted_lines << std::flush;
    };

    /////////////////////////////////////////////////////////////////////////
    // Cap at 48 cores because of diminishing returns. If three's more CPUs than that
    // it is better to split the inputs and launch multiple processes.
    static const auto num_cores = get_env("GX_NUM_CORES", std::min(48U, std::thread::hardware_concurrency() - 1));

    ostr << add_db_info(MakeMetaLine(GX_TSV_HEADER__HITS), db_path)
         << "\n#q-id\tq-pos1\ts-taxid\ts-id\ts-pos1\tlen\tmatchruns-L2\tmismatches" << std::endl;

    static const auto k_chunk_stride  = get_env("GX_QUERY_CHUNK_STRIDE", 100000UL);
    static const auto k_chunk_overlap = get_env("GX_QUERY_CHUNK_OVERLAP", 100UL);

    auto fasta_reader = MakeFastaReader(fasta_istr, k_chunk_stride, k_chunk_overlap);

    t = timer{};

#if 1   //using thread-pool
    const auto capacity = num_cores*4;
    gx::thread_pool executor{ num_cores, capacity };

#else   // uing std::async
    const auto capacity = num_cores;
    const auto executor = [](auto job)
    {
        return std::async(std::launch::async, std::move(job));
    };
#endif

    auto num_jobs = std::atomic_size_t{};

    if (num_cores == 0) {
        // In num_cores==0 will do the entire processing using the simple loop, in-this-thread,
        // without setting up the pipeline to get simpl(er) stack traces when debugging.
        for (fasta_seq_t qry_chunk : fasta_reader) {
            output_result(align_chunk_wr(std::move(qry_chunk)));
        }
    } else {
        size_t partial_sum = 0;

        std::move(fasta_reader)
#if 0   // one async job per input fasta-chunk
      % fn::transform_in_parallel(align_chunk_wr, std::ref(executor)).queue_capacity(capacity)
#else
        // For genomes with many small scaffolds, the above is not very efficitent,
        // because if sequences are of highly variable size, then the thread pool's
        // queue may saturate and block while waiting on a large job. Addditionally,
        // the parallelization-related overhead becomes significant.
        //
        // To deal with that, we will group the input sequences into batches
        // of aggregate length at most k_chunk_stride, and process each batch
        // by a single async-job.
        //
        // This yields 30% increase in throughput for GCA_900067735.1
      % fn::group_adjacent_by L((partial_sum += _.seq.size()) / k_chunk_stride)

      % fn::transform_in_parallel([&](std::vector<fasta_seq_t> batch)
        {
            ++num_jobs;
            return std::move(batch)
                 % fn::transform(align_chunk_wr)
                 % fn::to_vector();

        }, std::ref(executor)).queue_capacity(capacity)

      % fn::concat() // ungroup batches of results
#endif
      % fn::for_each(std::ref(output_result)); // NB: by ref, because stateful and mutable, and we want
    }                                          // the terminator-call below to refer to same instance.

    output_result(result_for_chunk_t{});       // Report "_self" row for last query.

    std::cerr << "Processed " << query_count << " queries, "
              << float(total_qry_bases)/1e6f << "Mbp in " << float(t)
              << "s. (" << float(total_qry_bases)/1e6f/float(t) << "Mbp/s);"
              << " num-jobs:" << num_jobs
              << std::endl;
}

