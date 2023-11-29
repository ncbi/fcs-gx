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

class CIndex
{
public:
    struct hmer38_t : public hmer_t<uint64_t, 38>
    {
        hmer38_t(uint64_t buf)
            : hmer_t<uint64_t, 38>(buf)
        {
#if 0
            // Experimental
            //
            // Find min hash-neighbor within 1-bit edit-distance.
            // If w is already hash-minimal in the neighborhood, we're done.
            // Otherwise change w to the neighbor and try again.
            //
            // I.e. transitively reach local hash-minimal neighbor.

            static const auto s_enable_hmer_grouping = get_env("GX_ENABLE_HMER_GROUPING", false);

            if (s_enable_hmer_grouping)
                while (true)
            {
                auto best = std::make_pair(w, uint64_hash(w));

                for (size_t i = 0; i < 38; i++) {
                    const auto modified_w = w ^ (1ULL << i);
                    const auto h = uint64_hash(modified_w);

                    if (h < best.second) {
                        best = std::make_pair(modified_w, h);
                    }
                }

                if (best.first == w) {
                    break;
                } else {
                    w = best.first;
                }
            }
#endif
            // Rotate the bits to put the middle bits in LSBs, such that the subkey
            // is determined by bits from the middle of the word rather than the end,
            // so that if we ignore subkey when we need sensitivity boost, we don't
            // nead to think about truncating the word template-length accordingly.
            w = ((w & 0b00000000000000000000000111111111111111) << 23)  //(15+8)
              | ((w & 0b11111111111111111111111000000000000000) >> 15); //(38-(15+8))
            //                         ^^^^^^^^ // subkey will be constructed from these positions
        }
    };

    static constexpr int k_word_tlen = 38 / 2 * 3 - 1; // template: 0b11011011011011011011011011011011011011011011011011011011
                                    // ^^         - bitwidth of hmer38_t above
                                    // ^^^^^^     - number of codons
                                    // ^^^^^^^^^^ - number of bases
                                    //            ^^^ exclude last codon's third position to make the template symmetric (reverse-complementable)

    static const uint8_t k_stride = 11;
    static_assert(k_stride % 3 != 0, ""); // so that we sample coding words in all three phases

    struct node_t
    {
          seq_oid_t seq_oid; // subject-seq-id of indexed h-mer
             pos1_t pos;     // position of indexed h-mer
            uint8_t subkey8; // lower bits of the hmer38_t::w - will use std::equal_range to locate within sub-bucket.

        bool operator<(const node_t& other) const
        {
            return subkey8 != other.subkey8 ? subkey8 < other.subkey8
                 : seq_oid != other.seq_oid ? seq_oid < other.seq_oid
                 :                                pos < other.pos;
        }

        bool operator==(const node_t& other) const
        {
            return       subkey8 == other.subkey8
                &&       seq_oid == other.seq_oid
                &&           pos == other.pos;
        }
    } __attribute__((packed));
    static_assert(sizeof(node_t) == 9, "");


    using nodes_t      = std::vector<node_t>;
    using nodes_view_t = fn::view<const node_t*>;

public:
    ///////////////////////////////////////////////////////////////////

    CIndex()
    {
        const auto scope_guard = make_exception_scope_guard([&]
        {
            std::cerr << "NB: Not enough RAM for creating gx-db. Require 24GiB + 1 byte/bp.\n";
        });

        m_buckets.resize(k_num_buckets);
    }

    CIndex(std::istream& istr)
    {
        from_stream(istr);
    }

    CIndex(CIndex&&) = default;
    CIndex& operator=(CIndex&&) = default;

    void to_stream(std::ostream& ostr) const;
    void from_stream(std::istream& istr);

    /////////////////////////////////////////////////////////////////////////
    // the pos1 is expected to be negative iff hmer is flipped
    void insert(const hmer38_t hmer, seq_oid_t seq_oid, pos1_t pos);

    void finalize(const seq_infos_t&);

    nodes_view_t at(hmer38_t) const;
    /////////////////////////////////////////////////////////////////////////

private:
    static constexpr size_t k_num_buckets = 1 << 30;

    std::vector<nodes_t> m_buckets = {};  // empty in mapped-mode.
                                          // pointers below are nullptr in owning mode.

    // The subcounts arrays store the cumulative count of nodes
    // in all buckets [0..i].
    //
    // To save space, instead of using a single 2^30 * 8byte array (8GiB),
    // we will use two arrays m_subcounts64 and m_subcounts32,
    // of sizes 2^30/65536 and 2^30 respectively
    // such that subcount[i] == m_subcounts64[i/65536] + m_subcount32[i], and
    // m_subcount32 stores subcounts relative to the corresponding m_subcounts64

    static constexpr size_t k_subcounts64_stride = 65536;

    const uint64_t* m_subcounts64 = nullptr;
    const uint32_t* m_subcounts32 = nullptr;
    const node_t*   m_nodes_ptr   = nullptr;
    bool            m_finalized   = false;

    uint64_t x_subcount(size_t i) const
    {
        return m_subcounts32[i] + m_subcounts64[i / k_subcounts64_stride];
    }
};

void TestDb(const fasta_seq_t& s);


} // namespace gx
