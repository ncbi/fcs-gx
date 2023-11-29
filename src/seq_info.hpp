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
#include "types.hpp"
#include "subbyte_array.hpp"
#include <map>

namespace gx
{

// describes a sequence in subjects-db
struct seq_info_t
{
               seq_oid_t seq_oid = seq_oid_t{}; // internal ordinal
                tax_id_t tax_id  = tax_id_t{};
                uint32_t length  = 0;
    std::array<char, 40> seq_id  = {};             // as array to keep seq_info_t as POD for serialization
                uint64_t offset_in_seq_db = 0;     // offset, in bytes, where the seq begins in seq-db (.gsi)

    seq_info_t() = default;

    seq_info_t(   seq_oid_t aseq_oid,
                   tax_id_t atax_id,
                   uint32_t alen,
         const std::string& aseq_id,
                   uint64_t aoffset_in_seq_db)

        : seq_oid(aseq_oid)
        , tax_id(atax_id)
        , length(alen)
        , seq_id() // set below
        , offset_in_seq_db(aoffset_in_seq_db)
    {
        set_seq_id(aseq_id);
    }

    void set_seq_id(const std::string& id)
    {
        VERIFY(id.size()+1 <= seq_id.size());
        strncpy(seq_id.data(), id.data(), id.size());
    }

    const char* get_seq_id() const
    {
        return seq_id.data();
    }

    bool operator<(const seq_info_t& other) const
    {
        return seq_oid < other.seq_oid;
    }

    bool operator==(const seq_info_t& other) const
    {
        return seq_oid == other.seq_oid;
    }
};
static_assert(sizeof(seq_info_t) == 64);

#if 0
using seq_infos_t = std::vector<seq_info_t>;
#else

// since we can't add an overload at(seq_oid_t) to a std::vector, will use a simple wrapper class

class seq_infos_t : public std::vector<seq_info_t>
{
public:
    using base_type = std::vector<seq_info_t>;
    seq_infos_t() = default;

    seq_infos_t(base_type b)
      : base_type(std::move(b))
    {}

    const seq_info_t& at(seq_oid_t i) const
    {
        return base_type::at((size_t)i);
    }
};
#endif


static inline std::map<seq_id_str_t, seq_oid_t> make_id2oid_map(const seq_infos_t& seq_infos)
{
    std::map<seq_id_str_t, seq_oid_t> id2oid;
    for (const auto& si : seq_infos) {
        id2oid[seq_id_str_t{ si.get_seq_id() }] = si.seq_oid;
    }
    return id2oid;
}

//////////////////////////////////////////////////////////////////////////////////

struct sbj_seq_t
{
    // view mode
    sbj_seq_t(const seq_info_t& si, const void* seq_db_data)
        : m_data{ reinterpret_cast<const uint8_t*>(seq_db_data) + si.offset_in_seq_db, si.length }
    {}

    // owning mode
    sbj_seq_t (const std::string& seq/*iupac*/)
        : m_data(seq.size())
    {
        const uint64_t hash_salt = seq.empty() ? 0 : uint64_hash(&seq[0], seq.size());
        for (const auto i : irange{ seq.size() }) {
            const auto na = seq[i];
            const uint8_t na2 = na == 'A' || na == 'a' ? 0
                              : na == 'C' || na == 'c' ? 1
                              : na == 'G' || na == 'g' ? 2
                              : na == 'T' || na == 't' ? 3
                              : 3 & uint64_hash(i ^ hash_salt); // pseudorandomly
            m_data.set(i, na2);
        }
    }

    sbj_seq_t() = default;

    // make move only - don't expect copy-assignment to be needed for now

    sbj_seq_t(const sbj_seq_t&) = delete;
    sbj_seq_t& operator=(const sbj_seq_t&) = delete;

    sbj_seq_t(sbj_seq_t&&) = default;
    sbj_seq_t& operator=(sbj_seq_t&&) = default;

    //////////////////////////////////////////////

    size_t size() const
    {
        return m_data.size();
    }

    void resize(size_t size)
    {
        m_data.resize(size);
    }

    na_t at1(pos1_t pos) const
    {
        static const char* bases_p = "ACGT";
        static const char* bases_m = "TGCA";
        VERIFY(pos != 0);
        const auto na2bit = m_data.get(as_pos0(pos));
        return na_t(pos > 0 ? bases_p[na2bit] : bases_m[na2bit]);
    }

    using subbyte_array_t = subbyte_array<2>;

    const subbyte_array_t& data() const
    {
        return m_data;
    }

private:
    subbyte_array<2> m_data;
};


} // namespace gx
