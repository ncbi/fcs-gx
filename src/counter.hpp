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

#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <cstddef>
#include <mutex>
#include <cmath>

namespace gx
{

class CCounter
{
public:
    CCounter(size_t bitwidth, uint8_t low_cap = 255)
    {
        clear(bitwidth, low_cap);
    }

    CCounter(CCounter&&) = default;
    CCounter& operator=(CCounter&&) = default;

    void insert(uint64_t key);
    size_t at(uint64_t key) const;

    struct stats_t
    {
        size_t n;
        size_t sum;
        size_t sum_of_squares;
        double sum_of_reciprocals;

        double mean()          const { return n == 0 ? 0 : double(sum)/double(n); }
        double rms()           const { return n == 0 ? 0 : std::sqrt(double(sum_of_squares)/double(n)); }
        double harmonic_mean() const { return n == 0 ? 0 : double(n) / sum_of_reciprocals; }
    };
    stats_t get_stats() const;

    void clear(size_t bitwidth, uint8_t low_cap)
    {
        m_low_cap = low_cap;
        m_high_counts = high_counts_t();
        m_low_counts = low_counts_t(1ull << bitwidth);
    }

public:
    using low_counts_t = std::vector<uint8_t>; // of size 2^bitwidth
    using high_counts_t = std::unordered_map<uint64_t, uint32_t>;

    float compute_baseline() const;

    // when a low_count for some key is about to reach 256,
    // will set it to 0 and increment the correspnding high-count.
    uint8_t m_low_cap;
    low_counts_t m_low_counts;
    high_counts_t m_high_counts; 
};


} // namespace gx
