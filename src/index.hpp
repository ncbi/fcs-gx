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


// todo: if we maintain counts on-line while indexing, we can
// insert a value into the bucket in O(1) with up-to 16 swaps.
// Can use the upper bits for extra 4-bits of key.
// I.e. use 30-bits for addressing into the main table, another
// four bits for addressing into the subbucket, and carry four
// bits in MSBs of value (will need to use 32-bit counters).

// Index 60-bit values using 38-bit key.
//
// For our purposes the 38-bit key is a word, or rather a hash of a fingerprint
// from a window in a dna sequence, value is a seq-id (28-bit id-ordinal + 32-bit position)
class CIndex
{
public:
    // 38-bit oriented word, selected as min(w, rev_comp(w))h
    struct minword38_t
    {
        minword38_t() : w(0), is_flipped(0)
        {}

        minword38_t(const skip_mod3_buf_t& b)
        {
            uint64_t buf = b.get();
            buf &= Ob1x(38);

            // Select orientation.
            auto flipped_buf = revcomp_bits(buf, 38); // number of bits in the mask
            this->is_flipped = flipped_buf < buf;
            this->w = is_flipped ? flipped_buf : buf;

            // Rotate the bits to put the middle bits in LSBs, such that the subkey
            // is determined by bits from the middle of the word rather than the end,
            // so that if we ignore subkey when we need sensitivity boost, we don't
            // nead to think about truncating the word template-length accordingly.
            w = ((w & 0b00000000000000000000000111111111111111) << 23)  //(15+8)
              | ((w & 0b11111111111111111111111000000000000000) >> 15); //(38-(15+8))
            //                         ^^^^^^^^ // subkey will be constructed from these positions
        }

        uint64_t w;
        bool is_flipped;
    };

    static constexpr int k_word_tlen = 38 / 2 * 3 - 1; // template: 0b11011011011011011011011011011011011011011011011011011011
                                    // ^^         - bitwidth of minword38_t above
                                    // ^^^^^^     - number of codons
                                    // ^^^^^^^^^^ - number of bases
                                    //            ^^^ exclude last codon's third position to make the template symmetric (reverse-complementable)

    // convert inclusive 0-based positional stop-pos of a word
    static inline pos1_t make_word_start_pos1(size_t stop_pos, bool is_flipped)
    {
        return is_flipped ? (int32_t(stop_pos) + 1) * -1
                          : (int32_t(stop_pos) + 1) + 1 - k_word_tlen;
    }

    static const uint8_t k_stride = 10;
    static_assert(k_stride % 3 != 0, ""); // so that we sample coding words in all three phases

    struct node_t
    {
          seq_oid_t seq_oid;
             pos1_t pos;
            uint8_t subkey8; // will use std::equal_range to locate within sub-bucket.

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
    using nodes_t = std::vector<node_t>;


    static_assert(sizeof(node_t) == 9, "");

    /////////////////////////////////////////////////////////////////////////

    struct nodes_view_t
    {
        const node_t* b = nullptr;
        const node_t* e = nullptr;

        const node_t* begin() const
        {
            return b;
        }

        const node_t* end() const
        {
            return e;
        }

        size_t size() const
        {
            return size_t(e - b);
        }
    };

public:
    ///////////////////////////////////////////////////////////////////

    CIndex()
    {
        const auto scope_guard = make_exception_scope_guard([&]
        {
            std::cerr << "NB: Not enough RAM. Require 24GiB + 1 byte/bp.\n";
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
    // the pos1 is expected to be negative iff minword is flipped
    void insert(const minword38_t minword, seq_oid_t seq_oid, pos1_t pos);

    void finalize(const seq_infos_t&);

    nodes_view_t at(minword38_t) const;
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
