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
#include "segment.hpp"

namespace gx
{

/*
 * Note: we are really interested in transposons, but not tandem repeats.
 *
 * Normally we would count the numer of times a word switches orientation
 * in a scan of a genome, rather than just number of instances of a word,
 * but this necessitates some complexity when doing it in parallel deterministically.
 *
 * (e.g. each input sequence chunk is assigned an ordinal; several threads
 * process chunks having (inpt_chunk_ord % thread_ord == 0), initializing
 * their own instance of tmasker, and then finally thread-specific instances
 * of tmasker are merged into the master instance.
 *
 * For now we'll just consider raw counts of a word
 */

class CTmasker
{
public:
    struct hmer30_t // not simply uisng hmer_t<uint32_t, 30> because don't need is_flipped field
    {
        uint32_t w = 0;

        hmer30_t(uint64_t buf)
            : w{ hmer_t<uint32_t, 30>(buf).w }
        {}
    };

    static const int k_word_tlen = 44; // template: 0b11011011011011011011011011011011011011011011

public:
    ///////////////////////////////////////////////////////////////////

    CTmasker() = default;

    CTmasker(CTmasker&&) = default;
    CTmasker& operator=(CTmasker&&) = default;

    CTmasker(const std::string& fasta_path); // can be empty
    CTmasker(std::istream& fasta_istr);

#if 1
    using counts_t = std::vector<uint8_t>;
#else
    struct count_t // if we wanted to collect transposon-specific counts
    {
        uint8_t val : 7;
        bool flipped : 1;
        // we are only interested in transposons, so to differentiate
        // from tandem repeat arrays we will only count number of times
        // a word flips orientation.
    };
    static_assert(sizeof(count_t) == 1);
    using counts_t = std::vector<count_t>;
#endif


    void insert(const hmer30_t hmer);


    uint8_t at(const hmer30_t hmer) const
    {
        return m_counts[hmer.w];
    }

    const counts_t& counts() const
    {
        return m_counts;
    }

    float baseline() const
    {
        return m_baseline;
    }

    bool collected_stats() const
    {
        return m_baseline > 1.0f;
    }

    gx::ivls_t find_repeats(const iupacna_seq_t& inp_chunk) const;
    void process_fasta(std::istream& fasta_istr, std::ostream& ostr) const;

private:
    static constexpr uint8_t k_count_cap = 250;

    float m_baseline = 0.0f; // expected number of occurrences of a non-frequent word
    counts_t m_counts = {}; // empty or size = 2^30
};


} // namespace gx
