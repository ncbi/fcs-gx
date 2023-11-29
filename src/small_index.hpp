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


// Small-key index for sensitive lookups.
// Will be used to index query-chunks on-the-fly, and then use
// subject-sequence in the neighborhoods of the preliminary low-sensitivity hits
// for querying into small-index to find additional hits.

class CSmallIndex
{
public:

#if 1
    using hmer20_t = hmer_t<uint32_t, 20>;
    static const int k_word_tlen = 29; // template: 0b11011011011011011011011011011
#else
    // GP-35566 - stalling this change pending publication
    //
    // NB: Insert extra 6 don't-care positions (two codons in two places),
    // making template length 29+6=35; to allow arbitrarily-spaced mismatches
    // (not just separated by a distance of multiple-of-3).
    static const int k_word_tlen = 35; 
    struct hmer20_t : public hmer_t<uint32_t, 20>
    {
        hmer20_t(uint64_t buf) 
            : hmer_t<uint32_t, 20>(
                    (buf & 0b00000000000000000000000000011011011)
                 |  (buf & 0b00000000000011011011011000000000000) >> 3
                 |  (buf & 0b11011011000000000000000000000000000) >> 6
              )         // 0b110110110---110110110110---11011011 - template
        {}   
    };
#endif

    struct node_t
    {
         pos1_t pos;
        uint8_t subkey8;

        bool operator<(const node_t& other) const
        {
            return subkey8 != other.subkey8 ? subkey8 < other.subkey8
                                            :     pos < other.pos;
        }

        bool operator==(const node_t& other) const
        {
            return     pos == other.pos
                && subkey8 == other.subkey8;
        }

    } __attribute__((packed));
    static_assert(sizeof(node_t) == 5, "");

    /////////////////////////////////////////////////////////////////////////

private:
    using nodes_t      = std::vector<node_t>;
    using nodes_view_t = fn::view<nodes_t::const_iterator>;
    using offsets_t    = std::array<uint32_t, 256>;

    struct bucket_t
    {
        offsets_t offsets = {}; // at position k contains count of nodes having subkey8 <= k.
        nodes_t nodes = {};
    };


public:
    ///////////////////////////////////////////////////////////////////

    CSmallIndex()
      : m_buckets(1 << 12)
    {}

    CSmallIndex(CSmallIndex&&) = default;
    CSmallIndex& operator=(CSmallIndex&&) = default;

    /////////////////////////////////////////////////////////////////////////
    // the pos1 is expected to be negative iff hmer is flipped
    void insert(const hmer20_t hmer, pos1_t pos);

    void finalize();

    void clear()
    {
        for (auto& b : m_buckets) {
            b.nodes.clear();
            b.offsets = offsets_t{};
        }
        m_finalized = false;
    }

    bool finalized() const
    {
        return m_finalized;
    }

    nodes_view_t at(hmer20_t) const;

    /////////////////////////////////////////////////////////////////////////
    std::vector<bucket_t> m_buckets = {}; // 4096 in size; upper 12 bits of key address the bucket
                                          // lower 8 bits address into buckets' offsets array
                     bool m_finalized = false;
};

} // namespace gx
