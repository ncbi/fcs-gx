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
#include <iostream>
#include <sstream>
#include <iterator>
#include <type_traits>

#ifdef __DEPRECATED
#   undef __DEPRECATED
#   include <strstream>  // for memistream
#   define __DEPRECATED
#else
#   include <strstream>
#endif

#include "util.hpp"

class ser
{
public:
    static std::string_view mmap(const std::string& path);

    // Access a random sample first and count page-faults.
    // If had page-faults, prefault all pages in file by accessing every page.
    static void prefetch_mmapped_pages(const std::string& filename, std::string_view data, bool force);

    // Verify that last 8 bytes in the file equal to magic_constant_eof; throw otherwise
    static void validate_eof(const std::string& filename, uint64_t magic_constant_eof);

    // Performance hack: without it std::getline (inside tsv_reader) is very slow
    // https://github.com/tmikov/linebench/blob/master/pbench.cxx/https://github.com/tmikov/linebench/blob/master/pbench.cxx
    static void prepare_standard_streams()
    {
        std::ios_base::sync_with_stdio(false);
        std::cin.sync_with_stdio(false);
        std::cout.sync_with_stdio(false);

        std::cin.imbue(std::locale::classic());
        std::cout.imbue(std::locale::classic());

        std::cin.tie(0);

        std::cin.exceptions(std::ios_base::goodbit);
        std::cout.exceptions(std::ios_base::goodbit);
    }


    ///////////////////////////////////////////////////////////////////////

    // Will align data blocks to cache lines.
    // See: GetPopcountInWindow
    static constexpr size_t k_align = 64;

    // basically istrstream but with exposed
    // ::begin() such that we can convert a
    // streampos to an address. This is used
    // for exposing a memory-mapped-file as
    // istream.
    struct memistream : virtual std::istrstream
    {
        memistream(const char* begin,
                   const size_t size)
            : istrstream(begin, std::streamsize(size))
            , m_begin(begin)
            , m_size(size)
        {
            VERIFY(reinterpret_cast<uintptr_t>(begin) % k_align == 0);
            // Verify the memory alignment of the whole memory block.
            // (memory mapped files are memory-paged aligned (4Kb).
            // The data in the stream will have proper padding
            // inserted as necessary (see read() and write() functions),
            // such that the memory-mapped POD-arrays from memistream
            // will have proper alignment in memory.
        }

        const char* begin() const
        {
            return m_begin;
        }

        size_t size() const
        {
            return m_size;
        }
    private:
        const char* m_begin;
        const size_t m_size;
    };

    ///////////////////////////////////////////////////////////////////////
    //
    //
    // Serialization helperes for POD types and vectors thereof


// https://www.google.com/#q=gcc+is_trivially_copyable
#if __GNUG__
#   define IS_TRIVIALLY_COPYABLE(T) __has_trivial_copy(T)
#else
#   define IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#endif

#define REQUIRE_TRIVIALLY_COPYABLE(T) \
    class = typename std::enable_if<IS_TRIVIALLY_COPYABLE(T)>::type

    template<typename T, REQUIRE_TRIVIALLY_COPYABLE(T)>
    static void to_stream(std::ostream& ostr, const T* ptr, size_t n = 1)
    {
        VERIFY(ostr.write(reinterpret_cast<const char*>(ptr), std::streamsize(sizeof(T)*n)));
        VERIFY(!ostr.fail());
    }

    template<typename T, REQUIRE_TRIVIALLY_COPYABLE(T)>
    static void to_stream(std::ostream& ostr, const std::vector<T>& v)
    {
        uint64_t size = v.size();
        VERIFY(size == v.size());
        to_stream(ostr, &size);
        to_stream(ostr, v.data(), v.size());
    }

    // helper for reading trivially-copyable POD types and vectors thereof,
    // serialized with to_stream above.
    //
    // const auto d = (double)from_stream(istr);
    // const auto v = (std::vector<my_pod_t>)from_stream(istr);
    struct from_stream
    {
        std::istream& m_istr;

        template<typename IStream>
        from_stream(IStream&& istr) : m_istr(istr) {}

        from_stream& operator=(const from_stream& other) = delete;

        // methods are rvalue-specific, because istr passed to constructor could be a temporary

        template<typename T, REQUIRE_TRIVIALLY_COPYABLE(T)>
        void to(T* ptr, size_t n = 1) &&
        {
            VERIFY(m_istr.read(reinterpret_cast<char*>(ptr), std::streamsize(sizeof(T)*n)));
            VERIFY(!m_istr.fail());
            VERIFY(size_t(m_istr.gcount()) == sizeof(T)*n);
        }

        template<typename T, REQUIRE_TRIVIALLY_COPYABLE(T)>
        operator T() &&
        {
            T t;
            std::move(*this).to(&t);
            return t;
        }

        template<typename T, REQUIRE_TRIVIALLY_COPYABLE(T)>
        operator std::vector<T>() &&
        {
            const uint64_t size = from_stream(m_istr);
            uninitialized_data_vector<T> v(size);
            std::move(*this).to(v.data(), size);
            return std::move(v);  // clang suggests a move here
        }
    };


    // throws if source pointer does not have proper alignment
    template<typename T>
    static inline T* safe_reinterpret_cast(void* p)
    {
        VERIFY(reinterpret_cast<uintptr_t>(p) % alignof(T) == 0);
        return reinterpret_cast<T*>(p);
    }

    template<typename T, REQUIRE_TRIVIALLY_COPYABLE(T)>
    static const T* from_memistream(memistream& mistr, size_t size)
    {
        auto ptr = mistr.begin() + mistr.tellg();
        VERIFY(size_t(ptr) % alignof(T) == 0);

        // simulate read(...) skipping over the memory-mapped bytes
        const size_t pos = mistr.tellg() + (int64_t)(size * sizeof(T));
        VERIFY(pos <= mistr.size());
        mistr.seekg(pos);

        return reinterpret_cast<const T*>(ptr);
    }


    ///////////////////////////////////////////////////////////////////////
    // add padding when necessary for memory-alignment in memory-mapped files

    static void skip_padding(std::istream& istr)
    {
        VERIFY(istr);

        // skip the padding
        const size_t padding_size = k_align - (size_t(istr.tellg()) % k_align);

        istr.seekg(istr.tellg() + (std::streamoff)padding_size);
    }

    ///////////////////////////////////////////////////////////////////////

    static void write_padding(std::ostream& ostr)
    {
        const auto padding_size = k_align - (size_t(ostr.tellp()) % k_align);

        static const std::array<char, k_align> zeroes{{0}};

        ostr.write(zeroes.data(), (std::streamsize)padding_size);

        VERIFY(!ostr.fail());
        VERIFY(size_t(ostr.tellp()) % k_align == 0);
    }

private:
    template<typename T>
    class uninitialized_data_vector : public std::vector<T>
    {
    public:
        uninitialized_data_vector(size_t n)
        {
#ifdef __GNUC__
            this->reserve(n);

            // resize without value-initializing the data
            std::vector<T>::_M_impl._M_finish =
                std::vector<T>::_M_impl._M_end_of_storage;
#else
            this->resize(n);
#endif
        }
    };
};
