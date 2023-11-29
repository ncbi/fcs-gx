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

// A "cursor" that scans over 1-bit ({CT, AG} alphabet) sequence and maintains the 
// recently seen bits in a 64-bit buf (current bit in lsb).
// for (auto kmer = kmer_ci_t(qry.seq, CIndex::k_word_tlen); kmer; ++kmer) { ... }
template<typename Seq>
struct kmer_ci_t
{
    const Seq& seq;
    size_t i = size_t(-1); // 0-based position into seq, corresponding to LSB in buf
    uint64_t buf = 0;
    uint64_t m_hash_salt = 0;

    kmer_ci_t(Seq&& seq_, size_t prefill_len, size_t offset = 0) = delete; // prevent binding to temporaries

    kmer_ci_t(const Seq& seq_, size_t prefill_len, size_t offset = 0)
        : seq{ seq_ }
        , i{ offset - 1 }
        , m_hash_salt{ uint64_hash(seq.size()) }
    {
        if constexpr (std::is_same_v<Seq, iupacna_seq_t>) {
            // Hash over the entire seq, as otherwise if we have a sequence
            // in db of the same length as this, and they both begin with an N-run,
            // we don't want them pseudorandomized to the same bits, but don't 
            // do that for non-zero offset cases (dealing with a sub-interval, 
            // and don't want to hash over the whole seq every time).
            m_hash_salt ^= offset != 0 ? uint64_hash(offset) : uint64_hash(seq.data(), seq.size());
        }

        do {
            ++(*this);
        } while (i < offset + prefill_len);
    }

    explicit operator bool() const
    {
        return i < seq.size();
    }

    void operator++()
    {
        ++i;
        buf <<= 1;

        if (i < seq.size()) {
            const auto na = seq[i];
            buf |= na == 'C' || na == 'c' || na == 'T' || na == 't' ? 1
                 : na == 'A' || na == 'a' || na == 'G' || na == 'g' ? 0
                 : bool(1 & uint64_hash(i ^ m_hash_salt)); // pseudorandomly
        }
    }
};



// h-mer: orientation-invariant hmer with coding-template (11-11-...-11) applied
template<typename UintT, size_t BitWidth>
struct hmer_t 
{
    static_assert(std::is_same_v<UintT, uint64_t> || std::is_same_v<UintT, uint32_t>, "");
    static_assert(BitWidth <= sizeof(UintT) * 8);
    static_assert(BitWidth <= 42); // after dropping every third bit from a 64-bit kmer.

    UintT w = 0; 
    bool is_flipped = false;

    hmer_t(uint64_t kmer) // from kmer_ci_t::buf
    {
        kmer = drop_every_3rd_bit(kmer) & Ob1x(BitWidth);
        const auto flipped_kmer = revcomp_bits(kmer, BitWidth);
        this->is_flipped = flipped_kmer < kmer;
        this->w = UintT(is_flipped ? flipped_kmer : kmer);
    }
};

} // namespace gx
