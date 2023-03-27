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
#include <iostream>

#define RANGELESS_FN_ENABLE_PARALLEL 1
#define RANGELESS_ENABLE_TSV 1

#include "align.hpp"
#include "types.hpp"
#include "util.hpp"
#include "seq_info.hpp"
#include "segment.hpp"
#include "subbyte_array.hpp"

#include <cmath>
#include <memory>
#include <algorithm>
#include <smmintrin.h>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;

///////////////////////////////////////////////////////////////

// domain is [0..size)
// Return: best endpos
template<typename F>
size_t ungapped_extend_impl(size_t size, F get_score_at, double max_dropoff)
{
    VERIFY(max_dropoff > 0);

    using score_t = decltype(get_score_at(0));
    score_t score(0), max_score(0);

    size_t best_endpos(0); // [0..size)

    int64_t run_size = 0; // run of positive scores
    for (size_t i = 0; i < size && score + max_dropoff >= max_score; i++) {
        const auto sc = get_score_at(int32_t(i));

        run_size = (run_size * (sc > 0)) + 1;
        score += sc * (score_t)run_size;

        // Extending with strict params results in underextended segs,
        // whereas extending with relaxed params results in false
        // extensions into low-identity flanks.
        //
        // The quadratic scoring extension works better than
        // traditional linear scheme where a matchrun receives
        // a score proportional to it's length.
        //
        // Instead of giving a constant positive match score for
        // each match, we give the match score multiplied by matchrun
        // length so far, such a matchrun of length N contributes
        // (1+2+..+N)*match_score overall.
        //
        // This works because a single match does not carry a
        // strong information about whether we're in aligned segment
        // or in unaligned sequence, whereas having a longer k-mer
        // matchrun has stronger signal that we're on a right track.

        // note: not score>=max_score, to avoid extending into N-runs
        if (score > max_score) {
            max_score = std::max(max_score, score);
            best_endpos = i + 1;
        }
    }

    return best_endpos;
}

///////////////////////////////////////////////////////////////

#if 0
class CSmallIndex
{
};

struct query_data_t
{
    fasta_seq_t seq;
    CSmallIndex index;
};
#endif


using ids_pair_t = std::pair<seq_id_str_t, seq_id_str_t>;

void gx::UngappedExtendSegsInPlace(
                  const fasta_seq_t& qry_seq,
                    const sbj_seq_t& sbj_seq,
                         segs_view_t segs)
{
    // TODO: an input-seg can be suprious.
    // If it's (56 + stride*N) in length, check it's identity and drop if too low.

    if (segs.empty()) {
        return;
    }

    // VERIFY that segs are oriented on query and sorted by subject, and all on same subject
    if (segs.end() > segs.begin() + 1) {
        const segment_t* prev = nullptr;
        for (const auto& seg : segs) {
            VERIFY(seg.q > 0);
            VERIFY(!prev || prev->s <= seg.s);
            VERIFY(!prev || prev->s_oid == seg.s_oid);
            prev = &seg;
        }
    }

    const len_t k_ext_cap = 2000;

    // TODO: To establish the extension params we need to assess the identity near the segment's ends.
    // But for now: match: +2; mismatch: -3; dropoff: -30.

    const pos1_t qry_end_on_plus  = pos1_t(qry_seq.offset + qry_seq.seq.size() + 1); // +1 to make 1-based
    const pos1_t qry_end_on_minus = pos1_t(qry_seq.offset + 1) * -1 + 1; // make 1-based, flip, make non-inclusive-end
    const pos1_t sbj_end_on_plus  = pos1_t(sbj_seq.size() + 1);
    const pos1_t sbj_end_on_minus = 0;

    const auto extend_seg = [&](segment_t& seg, const segment_t* next_seg, float mismatch_score, float x_dropoff_score)
    {
        VERIFY(seg.valid());

        const len_t max_ext_len = std::min(k_ext_cap, [&]
        {
            // can't extend past the end of sequence boundary
            pos1_t q_max_endpos = seg.q < 0 ? qry_end_on_minus : qry_end_on_plus;
            pos1_t s_max_endpos = seg.s < 0 ? sbj_end_on_minus : sbj_end_on_plus;

            // can't extend past the downstream segment, if in consistent configuration and on the same diag
            if (next_seg && segment_t::are_consistent(seg, *next_seg) && seg.diag() == next_seg->diag()) {
                q_max_endpos = std::min(q_max_endpos, next_seg->q);
                s_max_endpos = std::min(s_max_endpos, next_seg->s);
            }

            return q_max_endpos <= seg.q_end()
                || s_max_endpos <= seg.s_end() ? 0 : std::min(q_max_endpos - seg.q_end(),
                                                              s_max_endpos - seg.s_end());
        }());


        auto get_score_at = [&](int i) // i is position downstream of seg relative to seg's end.
        {
            const na_t q = qry_seq.at1(seg.q_end() + i);
            const na_t s = sbj_seq.at1(seg.s_end() + i);
            return q == s && q != 'N' ? 1.0f : mismatch_score;
        };

        VERIFY(max_ext_len >= 0);
        const auto best_ext = ungapped_extend_impl((size_t)max_ext_len, get_score_at, x_dropoff_score);
        seg.len += (int32_t)best_ext;

        VERIFY(seg.valid());

#if 0
        std::cout << seg.q << "\t" << seg.s << "\t" << seg.len << "\t" << best_ext << "\t" << max_ext_len<<  "\n";

        for (auto which : { -1, 0, 1 }) {
            const auto start = which == -1 ? seg.q   : seg.s;
            const auto& seq  = which == -1 ? qry_seq : sbj_seq;

            for (const auto i : irange{ seg.len }) {
                if (which == 0) {
                    auto is_match = get_score_at(i - seg.len) > 0;
                    std::cout << (is_match ? '|' : ' ');
                } else {
                    std::cout << seq.at1(start + i);
                }
            }

            std::cout << "\n";
        }
        std::cout << "\n\n";
#endif
    };

    const auto flip_segs = [&]
    {
        std::reverse(segs.begin(), segs.end());
        for (auto& seg : segs) {
            seg.flip();
        }
    };

    // to extend segs upstream, we flip them, extend downstream, and flip-back

    // With the current implementation it's possible that given adjacent segs
    // a and b we overextend a, while it would be better to extend b upstream, but
    // we stop extending while we reach the overextended boundary of a.
    //
    // To deal with this, we do several rounds of extensions with different
    // mismatch score - first extend conservatively, targeting the high-identity sequence,
    // and then follow-up with more relaxed scoring.
    //
    // NB: this is about 15% slower than single-round extensions.
    //

    for (float mismatch_score : { -9.0f, -4.0f, -2.0f }) // -9 targets 90% identity; -4 targets 80% identity -2 targets 60% identity
        for (size_t reverse_orientation : { false, true })
    {
        if (reverse_orientation) {
            flip_segs();
        }

        for (segment_t& seg : segs) {
            const segment_t* next_seg = &seg == &*(segs.end() - 1) ? nullptr : (&seg + 1);

            static const float x_dropoff_score = 20;
            extend_seg(seg, next_seg, mismatch_score, x_dropoff_score);
        }

        if (reverse_orientation) {
            flip_segs();
        }
    }
}


///////////////////////////////////////////////////////////////////////////


// a and b are 64-bit 2-bit-coding 32-mers (tail of the seq in LSBs)
//
// Return the best-shift distance between these, (maximizing the number of matches
// and penalizing for shift-amount) in a sub-k-mer of specified length located
// in MSBs of the overlap.
//
// E.g. consider two 32-mers below, where we want to locate best-shit for a 24-mer
//
// a = TTTTTTTTACGTAACCGGTTAAACCCGGGTTT
//          |    |   |     | |          // mismatches
// b =      ATTACATAAACGGTTGACCCCGGGTTTTTTTT
//          ^^^^^^^^^^^^^^^^^^^^^^^^^^^ // overlap
//          ^^^^^^^^^^^^^^^^^^^^^^^^    // best-matching 24-mer in MSBs of the overlap
// num-matches: 24-5= 19; best-shift of B relative to A = +5
//
// or the other way arround:
// a =      ATTACATAAACGGTTGACCCCGGGTTTTTTTT
//          |    |   |     | |
// b = TTTTTTTTACGTAACCGGTTAAACCCGGGTTT
//          ^^^^^^^^^^^^^^^^^^^^^^^^  // best matching 24-mer
// num-matches: 24-5= 19; best-shift of B relative to A = -5
//
// a = 0b1111111111111111000110110000010110101111000000010101101010111111
// b = 0b0011110001001100000001101011111000010101011010101111111111111111

struct best_shift_t
{
    int8_t shift = 0;
    int8_t num_matches = 0;

    best_shift_t() = default;

    best_shift_t(int8_t shift_, int8_t num_matches_)
        : shift(shift_)
        , num_matches(num_matches_)
    {}

    best_shift_t(uint64_t a, uint64_t b, uint8_t len)
    {
        auto ab = get_best_shift(a, b, len);
        auto ba = get_best_shift(b, a, len);
        ba.shift = int8_t(0 - ba.shift);
        *this = std::max(ab, ba);
    }

    bool operator<(const best_shift_t& other) const
    {
        return (this->num_matches - abs(this->shift))
             < (other.num_matches - abs(other.shift));
    }

private:
    static int8_t get_num_matches_in_msbs(uint64_t a, uint64_t b, uint8_t len)
    {
        a ^= ~b;                   // bitwise-eq
        a &= (a << 1);             // upper bits in 2-bit frames are 1 iff frame is 0b11
        a &= 0xAAAAAAAAAAAAAAAAUL; // zero-out lower bits in 2-bit frames
        a >>= 64 - (2 * len);      // shift-out irrelevant frames (2*len LSBs)
        return (int8_t)_mm_popcnt_u64(a);  // return count of frames that had 0b11 in them
    }

    static best_shift_t get_best_shift(uint64_t a, uint64_t b, uint8_t len)
    {
        const auto max_shift = int8_t(32 - len);
        auto ret = best_shift_t{};
        for (const int8_t i : irange{ max_shift }) {
            auto shift = best_shift_t{ i, get_num_matches_in_msbs(a, b, len) };
            ret = std::max(ret, shift);
            a <<= 2;
        }
        return ret;
    }
};


static const bool verify_best_shift = []
{
    const auto a = 0b1111111111111111000110110000010110101111000000010101101010111111;
    const auto b = 0b0011110001001100000001101011111000010101011010101111111111111111;

    const auto ab = best_shift_t(a, b, 24);
    const auto ba = best_shift_t(b, a, 24);

    VERIFY(ab.shift ==  5 && ab.num_matches == 19);
    VERIFY(ba.shift == -5 && ba.num_matches == 19);
    return true;
}();


segment_t gx::FindDownstreamContinuation(
                         const fasta_seq_t& qry_seq,
                           const sbj_seq_t& sbj_seq,
                            const segment_t inp_seg,
                                    uint8_t kmer_len)
{
    VERIFY(kmer_len <= 30);

    // TODO: verify that things work as expected when offset is not 0
    // todo: what if best_shift_t result is not optimal? Try top-2 instead?

    VERIFY(inp_seg.valid());

    // bail if don't have at least 32bp downstream to fill q_32mer and s_32mer
    {
        auto seg = inp_seg;
        seg.len += 32;

        seg.make_s_fwd();
        const bool s_ok = seg.s_end() <= 1 + int64_t(sbj_seq.size());

        seg.make_q_fwd();
        const bool q_ok = seg.q       >= 1 + int64_t(qry_seq.offset)
                       && seg.q_end() <= 1 + int64_t(qry_seq.offset+qry_seq.seq.size());

        if (!seg.valid() || !q_ok || !s_ok) {
            return segment_t{};
        }
    }

    const auto at1_2bit = [&](const auto& seq, pos1_t i)
    {
        const auto nt = seq.at1(i);
        return nt == 'A' ? 0
             : nt == 'C' ? 1
             : nt == 'G' ? 2
             : nt == 'T' ? 3
             :             uint64_hash(uint64_t(i)) % 4;
    };

    const best_shift_t best_shift = [&]
    {
        uint64_t q_32mer = 0;            // downstream of seg on q (all-As)
        uint64_t s_32mer = ~uint64_t(0); // downstream of seg on s (all-Ts)

        for (int i : irange{ 32 }) {
            q_32mer = (q_32mer << 2) | at1_2bit(qry_seq, inp_seg.q_end() + i);
            s_32mer = (s_32mer << 2) | at1_2bit(sbj_seq, inp_seg.s_end() + i);
        }

        return best_shift_t{ q_32mer, s_32mer, kmer_len };
    }();

    auto seg = segment_t{ inp_seg.q_end(), inp_seg.s_end(), inp_seg.s_oid, kmer_len };
    (best_shift.shift > 0 ? seg.q : seg.s) += abs(best_shift.shift);

    VERIFY(seg.valid());

    // truncate the to 1bp in the middle before ungapped-extension,
    // because it may be overextended already (i.e. less than kmer_len)
    const auto trunc_len = kmer_len/2;
    seg.q += trunc_len;
    seg.s += trunc_len;
    seg.len = 1;

    static thread_local segments_t segs; // for calling UngappedExtend
    segs.clear();
    segs.push_back(seg);
    UngappedExtendSegsInPlace(qry_seq, sbj_seq, fn::from(segs));
    VERIFY(segs.size() == 1);
    return segs.front();
}

/////////////////////////////////////////////////////////////////////////////

void gx::GappedExtend(    const fasta_seq_t& qry_seq,
   const std::function<sbj_seq_t(seq_oid_t)> get_sbj_seq,
                                 segments_t& inp_segs,
                                 const len_t min_len)
{
    VERIFY(min_len >= 0);

    VERIFY(std::is_sorted(inp_segs.begin(), inp_segs.end(), by_sbj));

    static thread_local segments_t out_segs;
    out_segs.clear();

    if (inp_segs.empty()) {
        return;
    }

    // Repeatedly FindDownstreamContinuation until can't find large-enough seg,
    // or ran out of neighborhood.
    const auto extend_downstream = [&]( segment_t from_seg,
                                      const len_t nbr_size) -> const segments_t&
    {
        static thread_local segments_t segs;
        segs.clear();

        const auto max_endpos = pos1_t(from_seg.s_end() + nbr_size);
        const sbj_seq_t sbj_seq = get_sbj_seq(from_seg.s_oid);

        while (from_seg.s_end() + min_len < max_endpos) {
            auto next_seg = FindDownstreamContinuation(qry_seq, sbj_seq, from_seg);

            if (next_seg.len < 10) {
                break;
            }

            VERIFY(next_seg.s_oid == from_seg.s_oid);
            segs.push_back(next_seg);
            from_seg = next_seg;
        }
        return segs;
    };

    /////////////////////////////////////////////////////////////////////////
    static const len_t max_neighborhood_size = 10000;

    // Find upstream extensions from the first seg, and downstream of the last
    out_segs <<= extend_downstream(segment_t{ inp_segs.front() }.flip(), max_neighborhood_size);
    out_segs <<= extend_downstream(inp_segs.back(), max_neighborhood_size);

    for (const auto i : irange{ 1, inp_segs.size() }) {
        const segment_t& a = inp_segs[i-1];
        const segment_t& b = inp_segs[i];

        // if a and b are in consistent configuration,
        // extend from a up-to b and vice-versa, otherwise up to max_neighborhood_size
        const len_t nbr_size =
            segment_t::are_consistent(a,b) ? 10 + std::max(0, segment_t::get_min_dist(a, b))
                                           : max_neighborhood_size;

        // todo: extend the search neighborhood by a few bases upstream?
        out_segs <<= extend_downstream(a                    , nbr_size);
        out_segs <<= extend_downstream(segment_t{ b }.flip(), nbr_size);
    }

    for (auto& seg : out_segs) {
        seg.make_q_fwd();
    }

    inp_segs <<= std::as_const(out_segs);

    Coalesce(inp_segs, 5);
}
