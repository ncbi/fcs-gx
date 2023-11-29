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
#pragma once

#include "util.hpp"
#include "seq_info.hpp"

namespace gx
{

struct ivl_t;
using ivls_t = std::vector<ivl_t>;

struct ivl_t
{
    pos1_t pos = 0;
    len_t len = 0;

    ivl_t(pos1_t pos_ , len_t len_)
         :   pos(pos_),  len(len_)
    {}

    ivl_t() = default;

    // NB: used to be named end(), but that name is commonly used with c++iterators
    pos1_t endpos() const
    {
        return pos1_t(pos + len);
    }

    bool operator==(const ivl_t& other) const
    {
        return pos == other.pos && len == other.len;
    }

    bool operator!=(const ivl_t& other) const
    {
        return !(this->operator==(other));
    }

    bool operator<(const ivl_t& other) const
    {
        return pos != other.pos ? pos < other.pos
                                : len < other.len;
    }

    bool valid() const
    {
        return (pos > 0 || (pos == 0 && len == 0) || (pos < 0 && pos + len <= 0));
    }

    void validate() const
    {
        if (!valid()) {
            GX_THROW("Invalid interval:" + this->to_string());
        }
    }

    std::string to_string() const
    {
        return    "[" + std::to_string(pos)
               + ".." + std::to_string(endpos()-1)
               + "]len:" + std::to_string(len);
    }

    /////////////////////////////////////////////////////////////////////////

    // signed distance negative means overlap
    static int32_t sdist(const ivl_t& a, const ivl_t& b)
    {
        VERIFY((a.pos < 0) == (b.pos < 0));
        const auto span = std::max(a.endpos(), b.endpos()) - std::min(a.pos, b.pos);
        return span - (int32_t)a.len - (int32_t)b.len;
    }

    // unsigned overlap (0 if no overlap)
    static len_t uoverlap(const ivl_t& a, const ivl_t& b)
    {
        auto d = ivl_t::sdist(a, b);
        return d >= 0 ? 0 : 0 - d;
    }

    static ivl_t collapse(const ivl_t& a, const ivl_t& b)
    {
        VERIFY(a.valid() && b.valid() && (a.pos < 0) == (b.pos < 0));
        const auto minpos = std::min(a.pos, b.pos);
        const auto maxend = std::max(a.endpos(), b.endpos());
        return ivl_t(minpos, maxend - minpos);
    }

    static ivl_t intersect(const ivl_t& a, const ivl_t& b)
    {
        VERIFY(a.valid() && b.valid() && (a.pos < 0) == (b.pos < 0));
        const auto maxpos = std::max(a.pos, b.pos);
        const auto minend = std::min(a.endpos(), b.endpos());
        return maxpos < minend ? ivl_t(maxpos, minend - maxpos)
                               : ivl_t();
    }

    /////////////////////////////////////////////////////////////////////////

    // Return [begin, end) into ivls for ivls overlapping this ivl.
    // Precondition: ivls' starts and ends are ordered.
    static auto get_overlapping_v(const ivls_t& ivls, const ivl_t ivl) -> fn::view<ivls_t::const_iterator>;

    // Collapse adjacent intervals if the distance between them is less than or equal to dist_thr.
    static void sort_and_merge(ivls_t& ivls, int32_t dist_thr = 0);

    // Collapse into last interval if overlapping within max_dist; otherwise push
    static void push_or_merge(ivls_t& ivls, ivl_t ivl, int max_dist = 0);

    // throw unless that all intervals are in plus-orientation, sorted, and non-overlapping
    static void verify_regular(const ivls_t& ivls);

    template<typename Ivls> // so it works for nodes_t
    static size_t sum_lens(const Ivls& ivls)
    {
        size_t ret = 0ul;
        for (const auto& ivl : ivls) {
            ret += (size_t)ivl.len;
        }
        return ret;
    }

    // precondition: ivls are sorted and non-overlapping.
    // returns ivls over gaps;
    // first and last are 0-length, capturing original boundaries
    // (except if they are abutting to inverted intervals and not necessary,
    // i.e. invert(invert(segs)) == segs.
    //
    // e.g. input: [ [10,20), [30,40) ]
    //     output: [ [10,10), [20,30), [40,40) ]
    static ivls_t invert(ivls_t ivls);

    // do the above, extend the boundaries upstream and downstream to whole_ivl's boundaries.
    // Constrain to whole_ivl, drop empty ivls.
    static ivls_t invert(ivls_t ivls, ivl_t whole_ivl);

    // precondition for intersect and subtract: everything is sorted and non-overlapping
    static ivls_t intersect(const ivls_t& ivls, const ivls_t& other);
    static ivls_t  subtract(const ivls_t& ivls, const ivls_t& other);

    static len_t intersection_len(const ivls_t& ivls, const ivls_t& other);

    static const ivl_t s_whole_ivl;
};

static inline ivl_t as_ivl(size_t stop_pos0, size_t len, bool is_flipped)
{
    auto start_pos = 
        is_flipped ? (int32_t(stop_pos0) + 1) * -1 // stop-pos becomes start-pos in opposite orientation
                   : (int32_t(stop_pos0) + 1) + 1 - (len_t)len; 
                                      // + 1 to make 1-based
                                      //      + 1 to convert stop-pos to exclusive end-pos
                                      //            subtract len to convert to start-pos
    return ivl_t{ start_pos, (len_t)len };
}


using locs_map_t = std::map<seq_id_str_t, ivls_t>;
auto LoadLocsMap(std::istream& istr_ptr) -> locs_map_t;

/////////////////////////////////////////////////////////////////////


struct segment_t
{
       pos1_t q = 0;
       pos1_t s = 0;
        len_t len = 0;
    seq_oid_t s_oid : 28;
      uint8_t flags : 4;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
     segment_t( pos1_t q_,
                pos1_t s_,
             seq_oid_t s_oid_,
                 len_t len_)
         : q(q_), s(s_), len(len_), s_oid(s_oid_), flags(0)
     {}
#pragma GCC diagnostic pop

     segment_t() = default;

     pos1_t q_end() const { return q + len; }
     pos1_t s_end() const { return s + len; }

    int32_t diag() const { return q - s; }
    int32_t antidiag() const { return q + s; }

    bool valid() const
    {
        // if positive - ok
        // if zero, must be zero-length
        // if negative, must not exceed the 0-bound

        return (q > 0 || (q == 0 && len == 0) || (q < 0 && q + len <= 0))
            && (s > 0 || (s == 0 && len == 0) || (s < 0 && s + len <= 0));
    }

    void validate() const
    {
        if (!valid()) {
            GX_THROW("Invalid segment:" + this->to_string());
        }
    }


    bool operator==(const segment_t& other) const
    {
        return     q == other.q
            &&     s == other.s
            && s_oid == other.s_oid
            &&   len == other.len;
    }

    bool operator!=(const segment_t& other) const
    {
        return !(*this == other);
    }


    std::string to_string() const
    {
        return   "{ q:" + std::to_string(q)
             +   ", s:" + std::to_string((uint32_t)+s_oid) // sign-promo warning without the cast (??)
             +      "@" + std::to_string(s)
             + ", len:" + std::to_string(len)
             +   ", d:" + std::to_string(diag())
             +  ", fl:" + std::to_string(+flags)
             + " }";
    }

    /// NB: (un)diagonalizing may potentially result in signed integer overflow.
    /// May consider removing the two methods below, and have diag() antidiag()
    /// return int64_t

    segment_t& flip_if(bool do_flip)
    {
        // Note: minus/minus configuration is also valid, e.g.
        // we want antisense configuration, and the genomic
        // strand follows.
        if (do_flip) {
            q = flip_pos1(q, len);
            s = flip_pos1(s, len);
        }
        return *this;
    }

    segment_t& flip()
    {
        return flip_if(true);
    }

    segment_t& make_q_fwd()
    {
        if (q < 0) {
            flip();
        }
        return *this;
    }

    segment_t& make_s_fwd()
    {
        if (s < 0) {
            flip();
        }
        return *this;
    }

    segment_t& coalesce(const segment_t& other)
    {
        VERIFY(are_coalescible(*this, other));
        VERIFY(flags == other.flags); // so there's no shenanigans

        const auto max_q = std::max(q_end(), other.q_end());
        q = std::min(q, other.q);
        s = std::min(s, other.s);
        len = max_q - q;

        VERIFY(valid() && diag() == other.diag());
        return *this;
    }

    static bool are_coalescible(const segment_t& a, const segment_t& b)
    {
        return a.diag()  == b.diag()
            && a.s_oid   == b.s_oid
            && (a.q < 0) == (b.q < 0)
            && (a.s < 0) == (b.s < 0);
    }

    static bool are_q_consistent(const segment_t& a, const segment_t& b)
    {
        return  a.q        <  b.q
            &&  a.q_end()  <  b.q_end()
            &&  (a.q < 0) == (b.q < 0);
    };

    static bool are_s_consistent(const segment_t& a, const segment_t& b)
    {
        return  a.s_oid  ==  b.s_oid
            &&  a.s       <  b.s
            &&  a.s_end() <  b.s_end()
            && (a.s < 0) == (b.s < 0);
    };

    static bool are_consistent(const segment_t& a, const segment_t& b)
    {
        return  are_q_consistent(a, b) && are_s_consistent(a,b);
    };

    // negative means overlap
    static int32_t get_q_dist(const segment_t& a, const segment_t& b)
    {
        VERIFY((a.q < 0) == (b.q < 0));
        const auto span = std::max(a.q_end(), b.q_end()) - std::min(a.q, b.q);
        return span - a.len - b.len;
    };

    // negative means overlap
    static int32_t get_s_dist(const segment_t& a, const segment_t& b)
    {
        VERIFY((a.s < 0) == (b.s < 0));
        const auto span = std::max(a.s_end(), b.s_end()) - std::min(a.s, b.s);
        return span - a.len - b.len;
    };


    static int32_t get_min_dist(const segment_t& a, const segment_t& b)
    {
        return std::min( get_q_dist(a, b),
                         get_s_dist(a, b));
    }
};
static_assert(sizeof(segment_t) == 16, "");

using segments_t = std::vector<segment_t>;
using segs_view_t = fn::view<segments_t::iterator>;

static inline const auto get_sbj = [](const segment_t& seg)
{
    return std::make_tuple(seg.s_oid, seg.s, seg.len);
};
static inline const auto by_sbj = fn::by::make_comp(get_sbj);

void Coalesce(segments_t& segs, len_t gap_thr = 100); // TODO: coalesce less aggressively?
void DropShadowedOnSbj(segments_t& segs);

void DropSingletons(
               segments_t& segs,
                     len_t word_len,
                     len_t diag_dist = 5000,
                     len_t antidiag_dist = 50000);


void FilterToHotWindows(segments_t& segs, len_t diag_w = 100000, len_t antidiag_w = 2000000);


void FilterOutRepeatsOnSbj(segments_t& segs, const seq_infos_t& sbj_infos, uint32_t num_elems = 100, len_t span_thr = 1000);

void FilterToBestChain(const segs_view_t segs,
                            const size_t max_backtrack_segs = 100,
                             const len_t max_backtrak_dist = 500000,
                             const float rearrangement_penalty = -500);

/*
Coalesce overlapping segments and drop singletons on-the-fly.
Keep hot-list (m_segs) addressable by seq-id-ord and fat-diag).

On-the-fly filtering prevents most of spurious seeds from materializing
beyond temporarily existing in the hotlist, and so cut downs on the
memory usage and speeds-up the downstream filtering stages.

Keep the incoming segs in the fixed-size hotlist.
If can coalesce into the existing seg, do that and mark it as keeper.
Else if neighbors with an existing seg, mark both as keepers.

Then replace the existing seg with the incoming one.

If the expunged seg was marked as a keeper, or is recent (not too far from
the incoming seg on q), then push it into dest-segs.
*/
class CSingletonFilter
{
public:
    CSingletonFilter(segments_t* dest, size_t size)
        : m_dest_segs_ptr(dest)
        , m_segs(size)
    {}

    void reset(segments_t* dest, size_t size)
    {
        m_dest_segs_ptr = dest;

        m_segs.clear();
        m_segs.resize(size);
    }

    void push(segment_t seg)
    {
        VERIFY(!seg.flags);

        seg.make_q_fwd(); // or s_fwd; just need to be consistent.

        // if seeds are witihn these distances, treat them as neighbors.
        // (These are the same as default parms in DropSingletons
        static const int k_diag_window = 5000;
        static const int k_antidiag_window = 50000;

        const auto are_neighbors = [&](const segment_t& l, const segment_t& r)
        {
            VERIFY(l.q <= r.q);

            return l.s_oid == r.s_oid
                && abs(l.diag() - r.diag()) <= k_diag_window
                && l.antidiag() + (l.len*2) + k_antidiag_window >= r.antidiag()
                && l.q != r.q;
            // Note: lengths double in diagonalized representation.
        };

        // prev-seg is the seg in the hotlist that is within k_diag_window on diag.
        // To get the "fat-diag" we divide by k_diag_window and round up or down.
        const auto get_prev_seg = [this](segment_t seg_, bool round_down) -> segment_t&
        {
            const auto k = round_down ? 0 : (k_diag_window/2);
            const auto fat_diag = (seg_.diag() + k) / k_diag_window;
            const auto i = uint64_hash((uint64_t(fat_diag) << 32)
                                       | (uint64_t(seg_.s_oid))) % m_segs.size();
            return m_segs[i];
        };

        segment_t& prev = get_prev_seg(seg, true);

        if (segment_t::are_coalescible(prev, seg) && prev.q_end() + 10 >= seg.q) {
            prev.len = seg.q_end() - prev.q;
            prev.flags = 1;
            return;
        } else if (are_neighbors(prev, seg)) {
            prev.flags = 1;
            seg.flags = 1;
        } else if (segment_t& other_prev = get_prev_seg(seg, false); are_neighbors(other_prev, seg)) {
            other_prev.flags = 1;
            seg.flags = 1;
        }

        VERIFY(prev.q <= seg.q);
        if (prev.flags == 1 || prev.q + k_antidiag_window/2 >= seg.q) {
            prev.flags = 0;
            m_dest_segs_ptr->push_back(prev);
        }
        prev = seg;
    }

    void finalize()
    {
        for (auto& seg : m_segs)
            if (seg.flags == 1)
        {
            seg.flags = 0;
            m_dest_segs_ptr->push_back(seg);
        }

        this->reset(m_dest_segs_ptr, m_segs.size());
    }

private:
    segments_t* m_dest_segs_ptr;
    segments_t m_segs; //fixed-size hot-list (hash-table)
};

}
