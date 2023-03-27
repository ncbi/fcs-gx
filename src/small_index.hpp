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
    struct minword20_t
    {
        minword20_t() : w(0), is_flipped(0)
        {}

        minword20_t(const skip_mod3_buf_t& b)
        {
            uint64_t buf = b.get();
            buf &= Ob1x(20);

            // Select orientation.
            auto flipped_buf = revcomp_bits(buf, 20); // number of bits in the mask
            this->is_flipped = flipped_buf < buf;
            this->w = uint32_t(is_flipped ? flipped_buf : buf);
        }

        uint32_t w;
        bool is_flipped;
    };

    static const int k_word_tlen = 29; // template: 0b11011011011011011011011011011

    // convert inclusive 0-based inclusive positional stop-pos of a word
    static inline pos1_t make_word_start_pos1(size_t stop_pos, bool is_flipped)
    {
        return is_flipped ? (int32_t(stop_pos) + 1) * -1 // stop-pos becomes start-pos in opposite orientation
                          : (int32_t(stop_pos) + 1) + 1 - k_word_tlen; // -1 to convert to endpos, and convert to startpos
    }

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
    using nodes_t = std::vector<node_t>;

    struct nodes_view_t
    {
        using iterator = nodes_t::const_iterator;

        iterator b = {};
        iterator e = {};

        const iterator begin() const
        {
            return b;
        }

        const iterator end() const
        {
            return e;
        }

        size_t size() const
        {
            return size_t(e - b);
        }
    };

    using offsets_t = std::array<uint32_t, 256>;

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
    // the pos1 is expected to be negative iff minword is flipped
    void insert(const minword20_t minword, pos1_t pos);

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

    nodes_view_t at(minword20_t) const;

    /////////////////////////////////////////////////////////////////////////
    std::vector<bucket_t> m_buckets = {}; // 4096 in size; upper 12 bits of key address the bucket
                                          // lower 8 bits address into buckets' offsets array
                     bool m_finalized = false;
};

} // namespace gx
