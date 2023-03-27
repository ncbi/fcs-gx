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
#include <string>
#include <cctype>
#include <map>
#include "util.hpp"

namespace gx
{

// single IUPAC nucleotide 'A', 'C', 'G', 'T', 'N', etc
enum class na_t : char {};

static inline bool operator==(const na_t a, char b)
{
    return (char)a == b;
}

static inline bool operator!=(const na_t a, char b)
{
    return (char)a != b;
}

enum class na2_t : uint8_t { A, C, G, T };
static inline na2_t as_na2(na_t na, uint64_t rand)
{
    return na2_t(
           na == 'A' ? 0u
         : na == 'C' ? 1u
         : na == 'G' ? 2u
         : na == 'T' ? 3u
         :             uint64_hash(rand) % 3u);
}


///////////////////////////////////////////////////////////////////////
// iupacna seq
class iupacna_seq_t : public std::string
{
public:
    iupacna_seq_t() = default;
    explicit iupacna_seq_t(std::string s) : std::string{ std::move(s) }
    {}

    na_t at(size_t i) const
    {
        return (na_t)std::string::at(i);
    }
};

///////////////////////////////////////////////////////////////////////
// seq-id string
class seq_id_str_t : public std::string
{
public:
    seq_id_str_t() = default;
    explicit seq_id_str_t(std::string s) : std::string{ std::move(s) }
    {}
};

///////////////////////////////////////////////////////////////////////


#if defined(__GNUC__) && (__GNUC__ < 10)
using seq_oid_t = uint32_t;
#else
// We would prefer to use the following instead:

enum class seq_oid_t : uint32_t{};
static inline uint32_t operator+(seq_oid_t seq_oid)
{
    return (uint32_t)seq_oid;
}

// It compiles with unremovable false-positive warning in segment.hpp
// 'gx::segment_t::s_oid' is too small to hold all values of 'enum class gx::seq_oid_t'
// due to GCC bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51242 (fixed in GCC-10)

// So in we will use the weakly typed seq_oid_t = uint32_t,
// but occasionally we want to complie in enum-class mode
// to ensure type-correctness.
//
// Alternatively we could use uint32_t only within segment_t, but then
// it results in explicit conversions mess.
#endif

enum class tax_id_t : uint32_t{};
static inline uint32_t operator+(tax_id_t tax_id)
{
    return (uint32_t)tax_id;
}

///////////////////////////////////////////////////////////////////////

// NB: these are assigned at runtime when loading tax_map_t based on
// ordinal position in ordering of taxdiv-names. I.e. these should not
// persist in files. At the time of this writing, we have 76 tax-divs
// in the big-index (26 are prok-specific).
enum class taxdiv_oid_t : uint8_t{ empty = 0, N_run = 255 };
static inline uint8_t operator+(taxdiv_oid_t taxdiv_oid)
{
    return (uint8_t)taxdiv_oid;
}

struct taxon_t
{
    std::string species;
    std::string common_name;
    std::string blast_taxdiv;   // e.g. moths
    std::string gx_taxdiv;      // e.g. anml:insects
    taxdiv_oid_t taxdiv_oid;    // NB: 0 iff tax_id=0, 255 for N_runs

    bool is_prok_or_virus() const
    {
        // NB: the following check works for old and new gx-div names (GP-33245)
        return str::startswith(gx_taxdiv, "prok")  // "prok|" or "prok:"
            || str::startswith(gx_taxdiv, "arch")  // archaea, crenarchaeotes, euryarchaeotes
            || str::contains(gx_taxdiv, "virus");
    }
};
using tax_map_t = std::unordered_map<tax_id_t, taxon_t>;  // NB: plain std::map is perf-bottleneck in taxify.cpp

// if nullptr, return tax_map_t containing (tax_id_t{}, taxon_t { "NULL", "NULL", "NULL", "NULL", taxdiv_oid_t{}})
tax_map_t LoadTaxa(std::istream* istr);

///////////////////////////////////////////////////////////////////////


// 1-based signed coordinate.
// +n and -n address the same base-pair on plus and minus strand respectively.
// 0 is invalid pos.
// NB: this is a weak-typedef because we need to do arithmetic ops on it.
using pos1_t = int32_t;

// Length of an sequence or an interval. Expected to be non-negative.
// (Not using uint32_t to avoid mixed-sign arithmetic with pos1_t)
using len_t = int32_t;

static inline pos1_t flip_pos1(pos1_t start, len_t len)
{
    return 1 - start - len; // == (start+len-1)*-1, i.e. change start to stop and flip.
}

// convert to absolute 0-based array-pos
static inline size_t as_pos0(pos1_t pos)
{
    VERIFY(pos != 0);
    return (size_t)abs(pos) - 1;
}

/////////////////////////////////////////////////////////////////////////////
// Existence of an indexed word with this position shall indicate that
// there exists a set of frequent words that were not indexed due to
// being determined as belonging to a repeat-sequence.
static const pos1_t k_repeat_marker = 1000000000;


/////////////////////////////////////////////////////////////////////////////

struct fasta_seq_t
{
     seq_id_str_t seq_id;
        seq_oid_t seq_oid;
      std::string defline;
           size_t offset; // seq starting at this offset on this seq-id
    iupacna_seq_t seq;

    static inline na_t up(na_t a)
    {
        return (na_t)std::toupper((char)a);
    }

    na_t at1(pos1_t pos) const
    {
        auto na = up(seq.at(as_pos0(pos) - offset)); // into this->seq
        return pos > 0 ? na : s_complement(na);
    }

    // in local chunk coordinates, ignoring offset
    na_t at1_local(pos1_t pos) const
    {
        auto na = up(seq.at(as_pos0(pos)));
        return pos > 0 ? na : s_complement(na);
    }

    static na_t s_complement(na_t);
};

/////////////////////////////////////////////////////////////////////////////
// 1-bit reduced-alphabet nucleotide
enum class bool_na_t : bool { AG=false, CT=true };

struct skip_mod3_buf_t
{
    uint64_t bufs[3] = { 0, 0, 0 };
    uint8_t phase = 0;

    // After pushing the bit into LSB, get() will yield 64-bit word with
    // every 3rd bit dropped, starting at 3rd LSB
    //
    // ....x-xx-xx-xx
    //      ^  ^  ^  dropped bits
    //              ^ LSB just-pushed
    //
    // // E.g. after pushing bits...
    // ...1011001110001111000011111 // inputs
    // ...1-11-01-10-01-11-00-11-11 // every 3rd bit is skipped
    //         ...11101100111001111 // result
    void operator<<=(bool_na_t na)
    {
        bufs[0] = (bufs[0] << 1) | (bool)na;
        bufs[1] = (bufs[1] << 1) | (bool)na;
        bufs[2] = (bufs[2] << 1) | (bool)na;

        bufs[phase] >>= 1;
        phase = uint8_t((phase + 1) % 3);

        // after moving to the next phase, the corresponding buf's lowest two bits
        // contain the last two bits that were pushed.
        //
        // static thread_local bool prev_bit = false;
        // VERIFY((bufs[phase] & 3) == ((uint64_t(prev_bit) << 1) | uint64_t(bit)));
        // prev_bit = bit;
    }

    template<typename Seq>
    static skip_mod3_buf_t init_from(const Seq& seq, size_t len, size_t offset = 0)
    {
        static_assert(std::is_same<bool_na_t, decltype(seq[0])>::value, "");

        skip_mod3_buf_t ret{};
        len = std::min(offset + len, seq.size());
        for (const auto i : irange{ offset, len }) {
            ret <<= seq[i];
        }
        return ret;
    }

    uint64_t get() const
    {
        return bufs[phase];
    };
};

/////////////////////////////////////////////////////////////////////////////
struct bool_na_view_t
{
public:
    bool_na_view_t(const iupacna_seq_t& seq)
        : m_seq{ seq }
        , m_hash_salt{ seq.empty() ? 0 : uint64_hash(&seq[0], seq.size()) }
    {}

    // We could just make the hash-salt based on seq size alone, but that will be a problem
    // if we have two different seqs of the same size, and they begin with a run
    // of Ns - we don't want the pseudorandomization of those Ns to be the same.

    bool_na_t operator[](const size_t i) const
    {
        const na_t c = m_seq.at(i);
        return bool_na_t(
               c == 'A' || c == 'G' || c == 'a' || c == 'g' ? false
             : c == 'C' || c == 'T' || c == 'c' || c == 't' ? true
             : bool(1 & uint64_hash(i ^ m_hash_salt))); // pseudorandomly
    }

    size_t size() const
    {
        return m_seq.size();
    }

private:
    const iupacna_seq_t& m_seq;
    uint64_t m_hash_salt;
};

} // namespace gx
