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
#include "align.hpp"
#include "types.hpp"
#include "seq_info.hpp"
#include "segment.hpp"
#include "small_index.hpp"

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;

// returns reference to static thread-local instance of CSmallIndex initialized from qry.
CSmallIndex gx::MakeQueryIndex(const fasta_seq_t& qry, CSmallIndex index)
{
    index.clear();

    const auto bits = bool_na_view_t{ qry.seq };
    auto buf = skip_mod3_buf_t::init_from(bits, CSmallIndex::k_word_tlen);

    for (const auto i : irange{ CSmallIndex::k_word_tlen, bits.size() }) {
        buf <<= bits[i];
        const auto minword = CSmallIndex::minword20_t{ buf };
        const auto pos1 = CSmallIndex::make_word_start_pos1(i + qry.offset, minword.is_flipped);
        index.insert(minword, pos1);
    }

    if (!qry.seq.empty()) {
        index.finalize(); // aovid expensive finalize() call if qry was empty
    }

    return index;
}

/////////////////////////////////////////////////////////////////////////////
// Return reference to static thread-local instance of segments_t
std::vector<tax_id_t> gx::SelectTaxaForRound2(
                                const segments_t& segs,
                               const seq_infos_t& sbj_infos,
                                 const tax_map_t& tax_map)
{
    // will aggregate segments statistics per-taxon
    struct tax_datum_t
    {
        tax_id_t tax_id       = {};
        int32_t  n            = 0;
        int64_t  sum_lens     = 0;
        int32_t  max_len      = 0;
        double   sum_lens_sq  = 0;
        bool     is_candidate = false;
    };

    std::vector<tax_datum_t> taxa = [&]
    {
        auto m = std::map<tax_id_t, tax_datum_t>{};

        for (const auto& seg : segs) {
            const auto& sbj_info = sbj_infos.at(seg.s_oid);

            auto& dest = m[sbj_info.tax_id];

            dest.tax_id       = sbj_info.tax_id;
            dest.max_len      = std::max(dest.max_len, int32_t(seg.len));
            dest.n           += 1;
            dest.sum_lens    += seg.len;
            dest.sum_lens_sq += seg.len * seg.len;
        }

        return std::move(m) % fn::transform L(_.second) % fn::to_vector();
    }();

#if 0
    {{
         auto out = std::ofstream{"tmp.tsv"};
         for (const auto& t : taxa) {
            out << +t.tax_id << "\t" << t.max_len << "\t" << t.n << "\t" << t.sum_lens << "\t" << t.sum_lens_sq << "\n";
         }
    }}
#endif

    auto mark_candidates_by = [&](auto key_fn)
    {
        taxa %= fn::unstable_sort_by(fn::by::decreasing(key_fn)); // best ones at front
        size_t n = 0; // count of marked so far.

        // Will mark up-to 3. The 3rd-best must be from a different tax-group.
        // because arbitrarily many top-taxa may be from close cross-species,
        // and for contamination-assessment purposes we also need next-best tax-group.
        for (auto& t : taxa) {
            if (n == 3) {
                break;
            } else if (  n == 0
                     ||  n == 1
                     || (n == 2 && tax_map.at(t.tax_id).taxdiv_oid
                                != tax_map.at(taxa.front().tax_id).taxdiv_oid))
            {
                t.is_candidate = true;
                n++;
            }
        }
    };

    mark_candidates_by L(_.n);           // L_0 norm
    mark_candidates_by L(_.sum_lens);    // L_1 norm
    mark_candidates_by L(_.sum_lens_sq); // L_2 norm (don't need to square-root)
    mark_candidates_by L(_.max_len);     // L_infinity norm

    mark_candidates_by L(double(_.sum_lens)    / double(_.n)); // mean-len
    mark_candidates_by L(double(_.sum_lens_sq) / double(_.n)); // RMS-len (don't need to square-root)

    auto ret = std::move(taxa)
         % fn::where L(_.is_candidate)
         % fn::transform L(_.tax_id)
         % fn::to_vector()
         % fn::unstable_sort(); // will binary-search into it

#if 0
    std::cerr << "Round-2 taxa:";
    for (const auto& x : ret) {
        std::cerr << " " << +x;
    }
    std::cerr << "\n";
#endif

    return ret;
}

///////////////////////////////////////////////////////////////////////////////

segments_t gx::SeedRound2(
        const CSmallIndex& index,
        const fasta_seq_t& qry,
          const sbj_seq_t& sbj_seq,
           const seq_oid_t sbj_oid,
        const segs_view_t& segs1,
                segments_t segs)
{
    segs.clear();

    static const auto nbr_size_abs_cap = get_env("GX_ROUND2_NBR_SIZE_ABS_CAP", 100000);
    static const auto nbr_size_rel_cap = get_env("GX_ROUND2_NBR_SIZE_REL_CAP", 2.0);
    const auto neighborhood_size = std::min(nbr_size_abs_cap, len_t(nbr_size_rel_cap * (double)qry.seq.size()));
    // neighborhood_size is the extent of the search neighborhood on subject sequneces around
    // preliminary segs. Some genomes are very fragmented and are composed of large number
    // of short sequences, so we can't spend search 100kb around every short query -
    // this becomes performance bottleneck, hence the need for relative cap.
    //
    // There's a use-case of aligning transcripts, where we may need to
    // search far enough to find neighboring exons missed in preliminary search,
    // hence making rel_cap configurable.


    // calculate subject-intervals over which we will search for alignments.
    const ivls_t& sbj_ivls = [&]
    {
        static thread_local ivls_t ivls;
        ivls.clear();

        const auto whole_s_ivl = ivl_t(pos1_t(1), (len_t)sbj_seq.size());

        if (segs1.empty()) {
            if (sbj_seq.size() <= 5000) {
                // will try to align small-query even if had no seeds
                ivls.push_back(ivl_t(pos1_t(1), (len_t)sbj_seq.size()));
            }
            return ivls;
        }

        VERIFY(std::is_sorted(segs1.begin(), segs1.end(), by_sbj));

        auto prev_seg = segment_t();
        for (const auto& seg : segs1) {
            VERIFY(seg.s > 0); // s-fwd

            if (!prev_seg.q) {
                // First interval to the left of the first seg
                ivls.push_back(ivl_t(seg.s - neighborhood_size, neighborhood_size));

                prev_seg = seg;
                continue;
            }

            auto dist = segment_t::get_s_dist(prev_seg, seg);
#if 0
            if (prev_seg.q != 0 && (prev_seg.q < 0) == (seg.q < 0)) {
                dist = std::min(dist, segment_t::get_q_dist(prev_seg, seg));
            }
#endif

            //std::cerr << prev_seg.to_string() << " -- " << seg.to_string() << "; dist=" << dist << "\t\t";
            if (dist < 30) {
                ; // too close to seed inbetween
                //std::cerr << "Skipping \n";
            } else if (dist < 2 * int64_t(neighborhood_size)) {
                // will process one interval between the segs.
                ivls.push_back(ivl_t(prev_seg.s_end(), dist));
                //std::cerr << ivls.back().to_string() << "\n";
            } else {
                // will process two intervals: after prev and before current
                ivls.push_back(ivl_t(prev_seg.s_end(), neighborhood_size));
                ivls.push_back(ivl_t(seg.s - neighborhood_size, neighborhood_size));
                //std::cerr << ivls[ivls.size() - 2].to_string() << ", " << ivls.back().to_string() << "\n";
            }

            prev_seg = seg;
        }

        // Last interval to the right of the last seg
        if (prev_seg.valid()) {
            ivls.push_back(ivl_t(prev_seg.s_end(), neighborhood_size));
        }

        // NB: at this point intervals may be overextended past the subj-seq boundaries.

        static const int k_ext = 50;
        for (auto& ivl : ivls) {
            // extend both ways.
            ivl.pos -= k_ext;
            ivl.len += k_ext*2;

            VERIFY(ivl.endpos() >= 0);

            // if extended to the left past the origin, truncate
            if (ivl.pos <= 0) {
                ivl.len -= abs(ivl.pos);
                ivl.pos = 1;
            }

            VERIFY(ivl.endpos() >= 0);

            // truncate to subject length if extended too far to the right
            ivl = ivl_t::intersect(whole_s_ivl, ivl);
        }

        ivls %= fn::where L(_.len > 0);
        ivl_t::sort_and_merge(ivls, k_ext);

        return ivls;
    }();

    /////////////////////////////////////////////////////////////////////////////////////

    static thread_local CSingletonFilter filter(nullptr, 0);
    size_t filter_size = std::min(10000UL, qry.seq.size()*10);
    filter.reset(&segs, filter_size);

    // view of subj seq in interval's coordinate system.
    struct bits_view_t
    {
        const sbj_seq_t& sbj_seq;
        const ivl_t ivl;

        bool_na_t operator[](const size_t i) const
        {
            return sbj_seq.at1_as_bool_na_t(pos1_t(ivl.pos + (int32_t)i));
        }

        size_t size() const
        {
            return size_t(ivl.len);
        }
    };

    size_t num_pushed = 0;

    for (const auto& ivl : sbj_ivls) {
        const bits_view_t bits{ sbj_seq, ivl };

        auto buf = skip_mod3_buf_t::init_from(bits, CSmallIndex::k_word_tlen);
        for (const auto i : irange{ CSmallIndex::k_word_tlen, bits.size() }) {
            buf <<= bits[i];

            const auto minword = CSmallIndex::minword20_t{ buf };
            const auto hits    = index.at(minword);
            const auto sbj_pos = CSmallIndex::make_word_start_pos1(size_t(ivl.pos-1) + i, minword.is_flipped); // pos-1 to make 0-based

            if (hits.size() <= 10)
                for (const auto& h : hits)
            {
                auto seg = segment_t{ h.pos, sbj_pos, sbj_oid, CSmallIndex::k_word_tlen };

                std::swap(seg.q, seg.s); // filter expects increasing by query, will swap-back in the end

                seg.make_q_fwd();

                ++num_pushed;
                filter.push(seg); // This may evict an older seg from the filter and push_back it into segs;
                                  // we examine it below to decide whether to keep it.

                while (!segs.empty() && segs.back().len < 50) { // 99% of hits is noise- need to filter by length
                    segs.pop_back();
                }
            }

            // TODO: else: search around diags in the lower-left and upper-right corners of the search-space.
            // (but we don't know what the search-space because we collapsed the intervals)

        }
    }

    filter.finalize();

    Coalesce(segs, 5);

    for (auto& seg : segs) {
        std::swap(seg.q, seg.s); // undo the swap we did before
        seg.make_q_fwd();
        VERIFY(seg.valid());
    }

    segs %= fn::where L(_.len >= 50);

    if (1) {
        FilterToBestChain(segs_view_t{ segs.begin(), segs.end() }, 100);
        segs %= fn::where L(_.flags != 0); // keep segs marked as belonging to best-chain
        for (auto& s : segs) {
            s.flags = 0;
        }
    }

    (void)index;
    (void)qry;
    (void)sbj_seq;
    (void)segs1;

#if 0
    std::cerr << "Round-1 segs: " << (segs1.end() - segs1.begin()) << "; len=" << sum_by(segs1, L(_.len))
              << "; round-2 segs: " << segs.size() << "; len=" << sum_by(segs, L(_.len))
              << "; ivls:" << sbj_ivls.size() << "; len=" << sum_by(sbj_ivls, L(_.len))
              << "; num_pushed:" << num_pushed << "\n";
#endif

    return segs;
}

/////////////////////////////////////////////////////////////////////////////
