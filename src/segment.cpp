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
#include "segment.hpp"
#include "types.hpp"

using namespace gx;

namespace fn = rangeless::fn;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;


const ivl_t ivl_t::s_whole_ivl = ivl_t{ 1, 2000000000 };

auto gx::LoadLocsMap(std::istream& istr) -> locs_map_t
{
    auto ret = locs_map_t{};

    ConsumeMetalineHeader(istr, GX_TSV_HEADER__LOCS);
    namespace tsv = rangeless::tsv;
    for (const tsv::row_t& row : tsv::from(istr)) {
        VERIFY(row.size() >= 3);
        const auto id = seq_id_str_t{ row[0] };

        if (row[1] == "." && row[2] == ".") {
            ret[id].push_back(ivl_t::s_whole_ivl);
        } else {
            int32_t from1 = tsv::to_num(row[1]);
            int32_t to1   = tsv::to_num(row[2]);
            VERIFY(from1 <= to1);
            VERIFY(from1 != 0); // For now require valid from..to intervals
            ret[id].push_back(ivl_t{ pos1_t(from1), to1 + 1 - from1 });
        }
    }

    for (auto& kv : ret) {
        ivl_t::sort_and_merge(kv.second);
    }

    return ret;
}


void gx::Coalesce(segments_t& segs, len_t gap_thr)
{
    for (auto& s : segs) {
        VERIFY(s.valid());
        s.make_q_fwd(); // shouldn't matter whether q-fwd or s-fwd,
                        // as long as all are in same orientation.
    }

    // NB: diag as int64_t to prevent signed-int-overflow. GP-34495
    segs %= fn::unstable_sort_by L(std::make_tuple(_.s_oid, int64_t(_.s) - _.q, _.q)); // by sbj-id, diag, qry

    for (const auto i : irange{ 1, segs.size() }) {
        auto& l = segs[i-1];
        auto& r = segs[i];

        if (segment_t::are_coalescible(l, r) && l.q_end() + gap_thr >= r.q) {
            l.len = r.q_end() - l.q;
            r = l;     // assign it to the right's so can keep coalescing
            l.len = 0; // mark as garbage (will drop these below).
        }
    }

    segs %= fn::where L(_.len > 0);
}

/////////////////////////////////////////////////////////////////////////////

void gx::DropShadowedOnSbj(segments_t& segs)
{
    for (auto& s : segs) {
        VERIFY(s.valid());
        s.make_s_fwd();
    }

    segs %= fn::unstable_sort_by(get_sbj);

    for (const auto i : irange{ 1, segs.size() }) {
        auto& l = segs[i-1];
        auto& r = segs[i];

        if (l.s_oid == r.s_oid
            &&  l.s_end()   > r.s_end()
            &&  l.len       > r.len * 2)
        {
            r = l;
            l.len = 0; // mark as garbage (will drop these below).
        }
    }

    // drop hits from repeats on query (stacked on subject after sorting)
    for_each_group_by(segs, get_sbj, [&](auto subsegs)
    {
        const auto n = subsegs.end() - subsegs.begin();
        if ((n >= 32 && subsegs.begin()->len <= 110)
            || (n >= 64 && subsegs.begin()->len <= 150))
        {
            for (auto& seg : subsegs) {
                seg.len = 0;
            }
        }
    });

    segs %= fn::where L(_.len > 0);
}

/////////////////////////////////////////////////////////////////////////////

void gx::DropSingletons(
               segments_t& segs,
                     len_t word_len,
                     len_t diag_dist_thr,
                     len_t dist_thr)
{
    for (auto& seg : segs) {
        VERIFY(seg.valid());
        seg.flags = 0;
        seg.make_q_fwd(); // this spreads hits on both strands on subject,
                          // so neighbrs, are closer together.
    }

    if (!std::is_sorted(segs.begin(), segs.end(), by_sbj)) {
        segs %= fn::unstable_sort_by(get_sbj);
    }

    for (const auto r : irange{ 1, segs.size() }) {
        for (size_t l = r - 1; l + 16 > r && l < segs.size(); l--) {
            VERIFY(segs[l] != segs[r]);

            if (abs(segs[l].diag() - segs[r].diag()) <= diag_dist_thr
                && segs[l].s_end() + dist_thr >= segs[r].s
                && segs[l].s_oid == segs[r].s_oid)
            {
                segs[l].flags = 1;
                segs[r].flags = 1;
                break;
            }
        }
    }

    // keep those that had neighbors, or those that were coalesced prior
    segs %= fn::where L(_.flags != 0 || _.len > word_len);

    for (auto& s : segs) {
        s.flags = 0;
    }
}


void gx::FilterToHotWindows(segments_t& segs, len_t diag_w, len_t antidiag_w)
{
    VERIFY(diag_w > 0);
    VERIFY(antidiag_w > 0);

    VERIFY(false); // disabling for now

    static thread_local std::vector<int32_t> counts;
    counts.clear();
    counts.resize(65536);

    auto get_hash = [&](seq_oid_t seq_oid, int32_t diag, int32_t antidiag)
    {
        const auto h1 = diag / diag_w;
        const auto h2 = antidiag / antidiag_w;

        return uint64_hash(uint64_t(seq_oid) ^ uint64_hash((uint64_t(h1) << 32) | uint64_t(h2))) % counts.size();
    };

    for (const auto& s : segs) {
        counts[get_hash(s.s_oid, s.diag()           , s.antidiag())]                += s.len;
        counts[get_hash(s.s_oid, s.diag() + diag_w/2, s.antidiag())]                += s.len;
        counts[get_hash(s.s_oid, s.diag()           , s.antidiag() + antidiag_w/2)] += s.len;
        counts[get_hash(s.s_oid, s.diag() + diag_w/2, s.antidiag() + antidiag_w/2)] += s.len;
    }

    //const auto count_thr = max_count / 32 + 1;

    const auto count_thr = 4000;


    // Instead, determine count_thr based on top-n bins.
    // Each bin spans antidiag_w of query-sequence

    /*
    // todo: get harmonic mean over non-zero counts to get the baseline
    const auto count_thr = [&]
    {
        double acc = 0;
        for (const auto& x : counts) {
            acc += 1/(x+1);
        }
        const double h_mean = (counts.size()/acc)-1;

        return size_t(h_mean*4);
    }();
    */

    //std::cerr << "Max-count;" << max_count << "; thr:" << count_thr << "\n";

    segs %= fn::where([&](const segment_t& s)
    {
        return
            counts[get_hash(s.s_oid, s.diag()           , s.antidiag())]                >= count_thr
         || counts[get_hash(s.s_oid, s.diag() + diag_w/2, s.antidiag())]                >= count_thr
         || counts[get_hash(s.s_oid, s.diag()           , s.antidiag() + antidiag_w/2)] >= count_thr
         || counts[get_hash(s.s_oid, s.diag() + diag_w/2, s.antidiag() + antidiag_w/2)] >= count_thr;
    });
}




// todo: take tax-infos and process per-taxon.
//
// Using a rolling window of size w_elems on query: if the distance between first and last hit is
// small (less then w_span_thr), it means we have a stacking of hits on query (repeats), so we'll
// filter those out.
void gx::FilterOutRepeatsOnSbj(segments_t& all_segs, const seq_infos_t& sbj_infos, uint32_t w_elems_unsigned, len_t w_span_thr)
{
    const auto w_elems = (int32_t)w_elems_unsigned;
    VERIFY(w_elems > 0);

    for (auto& seg : all_segs) {
        VERIFY(seg.valid());
        seg.make_q_fwd();
        seg.flags = 0;
    }

    if (!std::is_sorted(all_segs.begin(), all_segs.end(), by_sbj)) {
        all_segs %= fn::unstable_sort_by(get_sbj);
    }

    // After sotring by subject, all segs are collated by tax-ids.
    // Processing per-taxon:

    for_each_group_by(all_segs, L(sbj_infos.at(_.s_oid).tax_id), [&](segs_view_t segs)
    {
        static thread_local std::vector<segs_view_t> windows; // windows containing repeats
        windows.clear();

        if (segs.end() < segs.begin() + w_elems) {
            return;
        }

        std::sort(segs.begin(), segs.end(), BY(_.q));

        for (const auto it : irange{ segs.begin() + w_elems, segs.end() }) {

            const auto& l = *(it - w_elems);
            const auto& r = *it;

            if (r.q - l.q > int64_t(w_span_thr)) {
                continue;
            }

            auto w = segs_view_t{ it - w_elems, it + 1 };

            // extend the existing window if possible
            if (!windows.empty() && windows.back().end() >= w.begin()) {
                windows.back() = segs_view_t{ windows.back().begin(), w.end() };
            } else {
                windows.push_back(w);
            }
        }

        for (const auto& w : windows) {
            VERIFY(w.begin() < w.end());

            // compute average seg-length in this repeat-window
            int64_t tot_len = 0;
            for (const auto& seg : w) {
                tot_len += seg.len;
            }
            const auto avg_len = double(tot_len) / double(w.end() - w.begin());

            for (auto& seg : w)
                if (seg.len < avg_len * 4)
            {
                seg.flags = 1; // mark as junk
            }
        }
    });


    all_segs %= fn::where L(_.flags == 0);


    // Per-taxon.

    // Make-plus on query; sort by query.
    //
    // Using 100-element rolling-window.
    // If distance between last and first in window on query is small,
    // treat the segs in window as repeat-specific.

}


/////////////////////////////////////////////////////////////////////////////
void gx::FilterToBestChain  (const segs_view_t segs,
                                  const size_t max_backtrack_segs,
                                   const len_t max_backtrack_dist,
                                   const float rearrangement_penalty)
{
    const auto num_segs = size_t(segs.end() - segs.begin());

    if (num_segs == 0) {
        return;
    }

    for (auto& seg : segs) {
        VERIFY(seg.flags == 0);
        VERIFY(seg.valid());
        seg.make_q_fwd();
    }

    std::sort(segs.begin(), segs.end(), BY(_.q));

    ///////////////////////////////////

    auto seg_at = [&](size_t i) -> segment_t&
    {
        return *(segs.begin() + (long)i);
    };

    const auto max_antidiag_dist = max_backtrack_dist*2;


    auto are_consistent = [&](const segment_t& a, const segment_t& b)
    {
        return segment_t::are_consistent(a, b)
            && a.antidiag() + max_antidiag_dist > b.antidiag();
    };


    // [0..1]; approaches 1 if a and b are close on diag.
    // precondition: are_consistent(a,b)
    auto relative_diag_score = [](const segment_t& a, const segment_t& b) -> float
    {
        auto delta_diag = float(std::abs(a.diag() - b.diag()));
        auto delta_antidiag = float(b.antidiag() + b.len*2 - a.antidiag()) / 2;
        return std::max(0.0f, 1.0f - (delta_diag / (delta_antidiag + 1000.0f)));
    };

    // precondition: are_consistent(a,b)
    auto nonoverlapping_b_len = [](const segment_t& a, const segment_t& b)
    {
        return std::min( b.len,
                         std::min( b.q_end() - a.q_end(),
                                   b.s_end() - a.s_end()));
    };

    auto delta_score = [&](const segment_t& a, const segment_t& b) -> float
    {
        return (float)nonoverlapping_b_len(a,b) * (1 + relative_diag_score(a,b));
    };

    ///////////////////////////////////
    // Chain segs

    struct node_t
    {
        float score = 0;
        uint32_t prev = uint32_t(-1); // -1 means none

        void set_if_better_score(float sc, uint32_t pr)
        {
            if (sc > this->score) {
                this->score = sc;
                this->prev = pr;
            }
        }
    };

    static thread_local std::vector<node_t> nodes;
    nodes.clear();
    nodes.resize(num_segs);

    uint32_t i_maxnode = 0;


    for (const auto ib : irange{ 1u, (uint32_t)nodes.size() }) {

        const segment_t& b = seg_at(ib);

        node_t& node_b = nodes[ib];

        // Chain to max-node, paying rearrangement-penalty.
        if (ib == i_maxnode) {
            node_b.score = float(b.len);
        } else {
            // TODO: do not pay rearrangement penaly segs are consistent,
            // e.g. if we have a high stack of alignments on query,
            // greater than max_backtrack_segs, and we can't reach the
            // best previous-node by backtracking, we will always try it anyway.


            const segment_t& maxnode_seg = seg_at(i_maxnode);

            const bool is_q_consistent = maxnode_seg.q < b.q && maxnode_seg.q_end() < b.q_end();
            const auto nonoverlapping_b_len_on_q = std::max(0, std::min(b.len, b.q_end() - maxnode_seg.q_end()));

            // suppose on q we have:
            //         ------------------------------------------------         maxnode-seg
            // ...
            //            ------- --------     ---------- ---------- --------   cruft segs elsewhere on s
            //            ^segment b
            // If we chain b to the maxnode paying the rearrangement penalty,
            // the nodes following b will pay for the penalty and we will end up with
            // the maxnode and the cruft chained to it. Rather, if the maxnode seg
            // covers b, we don't want to chain b, and prevent any following segs
            // to be chained to b; so if not is-q-consistent
            node_b.score = nodes[i_maxnode].score
                         + rearrangement_penalty  // for chaining directly to maxnode
                         + (float)nonoverlapping_b_len_on_q
                         + (-1e9f * !is_q_consistent);
            node_b.prev = i_maxnode;
        }

        // Try to chain to an upstream neighbor.
        for (uint32_t ia = ib - 1; ia < nodes.size() && ia + max_backtrack_segs > ib; --ia) {
            const segment_t& a = seg_at(ia);

            VERIFY(a.q <= b.q);

            if (a.q_end() + int32_t(max_backtrack_dist) < b.q) {
                break; //backtracked far-enough
            }

            if (!are_consistent(a, b)) {
                continue;
            }

            const node_t& node_a = nodes[ia];

            node_b.set_if_better_score(node_a.score + delta_score(a, b), ia);
        }

        if (node_b.score > nodes[i_maxnode].score) {
            i_maxnode = ib;
        }
    }

    ///////////////////////////////////
    // Mark the segs belonging to the filetered chain
    for (auto i = i_maxnode; i != uint32_t(-1); i = nodes[i].prev) {
        seg_at(i).flags = 1;
    }
}


void ivl_t::verify_regular(const ivls_t& ivls)
{
    const ivl_t* prev = nullptr;
    for (const auto& ivl: ivls) {
        VERIFY(ivl.pos > 0);
        VERIFY(!prev || prev->endpos() <= ivl.pos);
        prev = &ivl;
    }
}


auto ivl_t::get_overlapping_v(const ivls_t& ivls, const ivl_t ivl) -> fn::view<ivls_t::const_iterator>
{
    const auto e = std::lower_bound(ivls.begin(), ivls.end(), ivl_t{ivl.endpos(), len_t(0)});
    auto b = e;
    while (ivls.begin() < b && ivl.pos < std::prev(b)->endpos()) {
        --b;
    }
    return fn::from(b, e);
}


void ivl_t::push_or_merge(ivls_t& ivls, ivl_t ivl, len_t max_dist)
{
    if (ivl.len == 0) {
        return;
    }

    const auto prev = !ivls.empty() ? &ivls.back() : nullptr;

    VERIFY(!prev || prev->pos <= ivl.pos);

    if (prev && ivl_t::sdist(*prev, ivl) <= max_dist) { // overlapping, or within max_dist
        *prev = ivl_t::collapse(*prev, ivl);
    } else {
        ivls.push_back(ivl);
    }
}


void ivl_t::sort_and_merge(ivls_t& ivls, int32_t dist_thr)
{
    ivls %= fn::unstable_sort_by L(_.pos);

    ivl_t* prev = nullptr;
    for (auto& ivl : ivls) {
        if (!prev) {
            prev = &ivl;
            continue;
        }

        VERIFY(prev->pos <= ivl.pos);
        VERIFY((prev->pos < 0) == (ivl.pos < 0));

        if (prev->len > 0 && ivl.len > 0 && ivl_t::sdist(*prev, ivl) <= dist_thr) {
            ivl = ivl_t::collapse(*prev, ivl);
            VERIFY(ivl.valid());
            *prev = ivl_t{};
        }
        prev = &ivl;
    }
    ivls %= fn::where L(_ != ivl_t{});
}

ivls_t ivl_t::invert(ivls_t ivls)
{
    if (ivls.empty()) {
        return ivls;
    }

    ivl_t::verify_regular(ivls);

    const auto span_len = ivls.back().endpos() - ivls.front().pos;
    const auto orig_len = ivl_t::sum_lens(ivls);

    ivls.push_back(ivl_t{ ivls.back().endpos(), 0});
    for (size_t i = ivls.size() - 2; i > 0; --i) {
        ivls[i].len = ivls[i].pos - ivls[i-1].endpos();
        ivls[i].pos = ivls[i-1].endpos();
    }
    ivls.front().len = 0;

    // Drop the zero-length placeholders at front at back if the are not necessary
    if (ivls.size() >= 2 && ivls[ivls.size() - 2].endpos() == ivls.back().endpos()) {
        ivls.pop_back();
    }

    if (ivls.size() >= 2 && ivls[1].pos == ivls.front().pos) {
        ivls.erase(ivls.begin());
    }

    const auto inv_len = ivl_t::sum_lens(ivls);
    VERIFY(orig_len + inv_len == (size_t)span_len);

    return ivls;
}


ivls_t ivl_t::invert(ivls_t ivls, ivl_t whole_ivl)
{
    VERIFY(whole_ivl.pos >= 0);

    if (ivls.empty()) {
        ivls.push_back(whole_ivl);
        return ivls;
    }

    ivl_t::sort_and_merge(ivls);
    ivls = ivl_t::invert(std::move(ivls));

    // front and back are 0-length boundaries; extend them all the way.
    {
        VERIFY(ivls.front().len == 0 && ivls.back().len == 0);
        VERIFY(ivls.front().pos >= 1); // verify that in plus-orientation
        VERIFY(ivls.back().pos >= 1);

        ivl_t& first = ivls.front();
        if (whole_ivl.pos < first.pos) {
            const auto d = first.pos - whole_ivl.pos;
            first.pos -= d;
            first.len += d;
        }

        if (whole_ivl.endpos() > ivls.back().endpos()) {
            ivls.back().len = whole_ivl.endpos() - ivls.back().pos;
        }
    }

    // constrain to inp-chunk boundaries
    for (auto& ivl : ivls) {
        ivl = ivl_t::intersect(ivl, whole_ivl);
    }

    ivls %= fn::where L(_.len > 0);

    return ivls;
}


ivls_t ivl_t::intersect(const ivls_t& lhs, const ivls_t& rhs)
{
    // TODO: iterate ivls having fewer elements.
    // Advance to starting element in the other set by binary search.

    ivl_t::verify_regular(lhs);
    ivl_t::verify_regular(rhs);

    auto ret = ivls_t{};
    ret.reserve(std::min(lhs.size(), rhs.size()));

    auto rhs_it = rhs.begin();

    for (const ivl_t& l_ivl : lhs) {
        while (rhs_it != rhs.end() && rhs_it->endpos() <= l_ivl.pos) {
            ++rhs_it;
        }

        if (rhs_it == rhs.end()) {
            break;
        }

        for (auto it = rhs_it; it != rhs.end() && it->pos < l_ivl.endpos(); ++it) {
            auto i_ivl = ivl_t::intersect(l_ivl, *it);
            if (i_ivl.len > 0) {
                ret.push_back(i_ivl);
            }
        }
    }

    ivl_t::verify_regular(ret);

    return ret;
}

len_t ivl_t::intersection_len(const ivls_t& lhs, const ivls_t& rhs)
{
    auto rhs_it = rhs.begin();
    auto ret = len_t{0};

    for (const ivl_t& l_ivl : lhs) {
        while (rhs_it != rhs.end() && rhs_it->endpos() <= l_ivl.pos) {
            ++rhs_it;
        }

        if (rhs_it == rhs.end()) {
            break;
        }

        for (auto it = rhs_it; it != rhs.end() && it->pos < l_ivl.endpos(); ++it) {
            ret += ivl_t::intersect(l_ivl, *it).len;
        }
    }
    return ret;
}



ivls_t ivl_t::subtract(const ivls_t& lhs, const ivls_t& rhs)
{
    if (lhs.empty() || rhs.empty()) {
        return lhs;
    }

    ivl_t::verify_regular(lhs);
    ivl_t::verify_regular(rhs);

    const auto lhs_whole_ivl = ivl_t{ lhs.front().pos, lhs.back().endpos() - lhs.front().pos };

    static thread_local ivls_t irhs;
    irhs.clear();
    irhs.insert(irhs.end(), rhs.begin(), rhs.end());
    irhs = ivl_t::invert(std::move(irhs), lhs_whole_ivl);

    return ivl_t::intersect(lhs, irhs);
}


static const bool test1 = []
{
    /*
     * ----   ------ ------  ------   # lhs
     *   ----     ---- -              # rhs
     * ^^    -----    - ^^^^^^^^^^^   # irhs (inverted); ^s are extensions to lhs_whole_ivl
     * --     ----    - ---  ------   # lhs intersected with irhs
    */
    const ivls_t lhs{ {8,4},  {15,6}, {22,6}, {30,6} };
    const ivls_t rhs{ {10,4}, {19,4}, {24,1} };
    const ivls_t e_result{ {8,2},  {15,4}, {23,1}, {25,3}, {30,6} };
    const auto a_result = ivl_t::subtract(lhs, rhs);
    VERIFY(a_result == e_result);

    // test get_overlapping_v
    const auto overlapping = ivl_t::get_overlapping_v(e_result, ivl_t{16, 10});
    VERIFY( overlapping.begin() == e_result.begin() + 1);
    VERIFY( overlapping.end() == e_result.begin() + 4);

    // test invert
    const auto irhs = ivl_t::invert(rhs);
    VERIFY((irhs == ivls_t{ {10,0}, {14,5}, {23,1}, {25,0} }));
    VERIFY(ivl_t::invert(irhs) == rhs);

    return true;
}();
