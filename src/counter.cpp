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
#include "counter.hpp"
#include <mutex>
#include <array>
#include <algorithm>

using namespace gx;

void gx::CCounter::insert(uint64_t key)
{
    using lock_t = std::mutex;
    static std::array<lock_t, 4096> s_mutexes; // to minimize contention; size determined empirically.
    std::unique_lock<lock_t> lock{ s_mutexes[key % s_mutexes.size()] };

    if (m_low_counts.at(key) < m_low_cap) {
        ++m_low_counts.at(key);
        return;
    }

    m_low_counts.at(key) = 0;

    lock.unlock();

    static std::mutex s_mut;
    const std::lock_guard<std::mutex> lg{ s_mut };
    ++m_high_counts[key];
}

size_t gx::CCounter::at(uint64_t key) const
{
    const auto it = m_high_counts.find(key);
    return m_low_counts.at(key) + (it == m_high_counts.end() ? 0ul : it->second * m_low_cap);
}

gx::CCounter::stats_t gx::CCounter::get_stats() const
{
    stats_t ret{ 0, 0, 0, 0.0 };

    for (const auto x : m_low_counts) {
        ret.n += x > 0;
        ret.sum += x;
        ret.sum_of_squares += x * x;

        ret.sum_of_reciprocals += x == 0 ? 0 : 1.0 / x;
    }

    for (const auto& kv : m_high_counts) {
        const auto lo = m_low_counts.at(kv.first);
        const auto hi = kv.second * m_low_cap;

        ret.n += lo == 0; // as not to double-count
        ret.sum += hi;

        ret.sum_of_squares -= lo * lo; // undo the contribution from the first loop
        ret.sum_of_squares += (hi + lo) * (hi + lo);

        ret.sum_of_reciprocals -= lo == 0 ? 0 : 1.0 / lo;
        ret.sum_of_reciprocals += 1.0 / double(hi + lo);
    }

    return ret;
}
