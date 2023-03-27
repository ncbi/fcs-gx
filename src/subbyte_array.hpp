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
#include <limits>
#include <vector>
#include "util.hpp"

template<int N>
class subbyte_array
{
public:
    typedef uint8_t value_type;

    static_assert(N == 1 || N == 2 || N == 4, "");

    subbyte_array(size_t size = 0)
    {
        m_data = m_storage.data();
        resize(size); // this will repoint m_data as appropriate
    }

    // non-owning read-only-view mode.
    subbyte_array(const uint8_t* data, size_t size)
    {
        VERIFY(data != nullptr || size == 0);

        m_size = size;
        m_data = data;
    }

    size_t size() const
    {
        return m_size;
    }

    size_t size_bytes() const
    {
        return s_size_in_bytes(m_size);
    }

    static size_t s_size_in_bytes(size_t size)
    {
        return (N * size + 7)/8;
    }

    const uint8_t* data() const
    {
        return m_data;
    }

    void clear()
    {
        *this = subbyte_array();
    }

    void resize(size_t size)
    {
        VERIFY(m_data == m_storage.data()); // must be in owning mode

        m_size = size;
        m_storage.resize(s_size_in_bytes(size));

        m_data = m_storage.data();
    }

    inline uint8_t get(size_t pos) const
    {
        const size_t i = pos >> k_lowbits_shift;

        ASSERT(i < size_bytes());

        const uint8_t w = m_data[i];
        pos = (pos & k_lowbits_mask) * N;
        return uint8_t((w >> pos) & k_value_mask);
    }

    uint8_t operator[](size_t pos) const
    {
        return get(pos);
    }

    void set(size_t pos, uint8_t val)
    {
        VERIFY(m_data == m_storage.data()); // must be in owning mode

        const size_t i = pos >> k_lowbits_shift;

        VERIFY(i < m_storage.size());
        VERIFY(val <= k_value_mask);

        pos = (pos & k_lowbits_mask) * N;
        uint8_t& w = m_storage[i];
        w &= uint8_t(~(k_value_mask << pos)); // clear target bits
        w |= uint8_t(val << pos);             // fill target bits with value
    }

protected:
    // number of bits to address the value within containing uint8_t
    static constexpr uint8_t k_lowbits_shift =
        N == 1  ? 3
      : N == 2  ? 2
      : N == 4  ? 1
      :           0 ;

    static constexpr uint8_t k_lowbits_mask = ~(~0U << k_lowbits_shift);
    static constexpr uint8_t k_value_mask   = ~(~0U << N);

    /////////////////////////////////////////////////////////////////////

                  size_t m_size = 0UL;
    std::vector<uint8_t> m_storage = {};
          const uint8_t* m_data = nullptr; // into m_storage.data(), or external
};
