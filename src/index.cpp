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
#define RANGELESS_FN_ENABLE_PARALLEL 1
#include "index.hpp"
#include "types.hpp"
#include "seq_info.hpp"
#include "serial_util.hpp"
#include <mutex>
#include <math.h>
#include <algorithm>

using namespace gx;

using rangeless::fn::operators::operator%=;
using rangeless::fn::operators::operator%;


// 38-bit key
struct key38_t
{
    uint32_t key30; // index into bucket (m_buckets)
    uint8_t  key8;  // attached to the value; hits are looked-up within sub-bucket with std::equal_range

    key38_t(uint64_t a = 0)
      : key30(a & Ob1x(30))
      , key8((a >> 30) & Ob1x(8))
    {
        VERIFY(a >> 38 == 0);
    }
};


void gx::CIndex::insert(const minword38_t minword, seq_oid_t seq_oid, pos1_t pos)
{
    ASSERT(!m_finalized);

    VERIFY(minword.is_flipped == (pos < 0));

    const key38_t key{ minword.w };

    auto& nodes = m_buckets[key.key30];

    // instead of locking the entire m_buckets, we only lock the subset of buckets
    // that is synchronized by a corresponding mutex.
    using lock_t = std::mutex; //atomic_mutex;
    static std::array<lock_t, 4096> s_mutexes;
    const std::lock_guard<lock_t> lg{ s_mutexes[key.key30 % s_mutexes.size()] };

    nodes.push_back(node_t{ seq_oid, pos, key.key8 });
    m_finalized = false;
}


/////////////////////////////////////////////////////////////////////////////

void gx::CIndex::finalize(const seq_infos_t& seq_infos)
{
    ASSERT(!m_finalized);

    std::atomic_size_t nodes_total_pre_filter = {};
    std::atomic_size_t nodes_total_post_filter = {};

    static const auto num_threads = get_env("GX_NUM_CORES", std::min(32U, std::thread::hardware_concurrency() - 1));

    for_each_in_parallel(m_buckets, num_threads, [&](nodes_t& bucket_nodes)
    {
        if (bucket_nodes.empty()) {
            return;
        }

        bucket_nodes %= fn::unstable_sort();
        bucket_nodes %= fn::unique_adjacent(); // in case a seq was indexed twice

        // Keep few per-taxon.
        //
        // bucket nodes are sorted by:   (key8, value);
        // expanding value (packed_pos): (key8, seq_oid, seq_pos).
        //
        // seq-oids are collated by taxa,
        // threfore nodes are collated by (key8, tax-id).

        for_each_group_by( bucket_nodes,
                           L(std::make_pair(_.subkey8, seq_infos.at( _.seq_oid ).tax_id)),
                           [&](const auto subnodes_v) // having the same subkey8 and tax-id
        {
            // Keep first few nodes for the key; mark rest for deletion by setting pos=0.
            for(const auto it : irange{ subnodes_v.begin() + 4, subnodes_v.end() }) {
                it->pos = 0;
            }
        });

        nodes_total_pre_filter += bucket_nodes.size();
        bucket_nodes %= fn::where L(_.pos != 0);
        nodes_total_post_filter += bucket_nodes.size();
    });

    m_finalized = true;

    std::cerr << "Nodes pre-filter:" << nodes_total_pre_filter
              << "; nodes post-filter:" << nodes_total_post_filter
              << ".\n";
}

/////////////////////////////////////////////////////////////////////////

#if 0
// Disabling. Does not seem to be any better than std::equal_range when multithreaded.

template<typename Iterator, typename T, typename Compare>
auto my_equal_range( Iterator first,
                     Iterator last,
                     const T& value,
                      Compare comp) -> std::pair<Iterator, Iterator>
{
    auto len = std::distance(first, last);

    while (len > 0) {
        auto half = len >> 1;
        Iterator middle = first;
        std::advance(middle, half);

        //https://www.daemon-systems.org/man/__builtin_prefetch.3.html
        // prefetch &*middle for the next round here.
        __builtin_prefetch(&*(middle+1 + ((len - half - 1) >> 1)), 0, 1); //https://www.daemon-systems.org/man/__builtin_prefetch.3.html
        __builtin_prefetch(&*(first + ((half) >> 1)), 0, 1);

        if (comp(*middle, value)) {
            first = ++middle;
            len -= half + 1;
        } else if (comp(value, *middle)) {
            last = middle;
            len = half;
        } else {
            auto mp1 = middle;
            return {
                std::lower_bound(first, middle, value, comp),
                std::upper_bound(++mp1, last, value, comp)
            };
        }
    }

    return {first, first};
}
#endif


gx::CIndex::nodes_view_t gx::CIndex::at(minword38_t minword) const
{
    ASSERT(m_finalized);

    const key38_t key{ minword.w };

    VERIFY((m_subcounts32 == nullptr) == (m_nodes_ptr == nullptr));
    VERIFY(!m_buckets.empty() ^ (m_subcounts32 != nullptr));

    const auto nodes_view = // depending on whether dealing with owning storage or memory-mapped (m_nodes_ptr)
        m_buckets.empty() ? nodes_view_t{ m_nodes_ptr + (key.key30 == 0 ? 0 : x_subcount(key.key30 - 1)),
                                          m_nodes_ptr + x_subcount(key.key30) }

                          : nodes_view_t{ m_buckets[key.key30].data(),
                                          m_buckets[key.key30].data() + m_buckets[key.key30].size() };

    auto begin_end =
            std::equal_range( // my_equal_range or std::equal_range
                 nodes_view.begin(), nodes_view.end(),
                 node_t{ seq_oid_t{}, pos1_t{}, key.key8 },
                 BY(_.subkey8));

    // If no hits, return close neighbors (ignoring subkey8 entirely, or lower bits thereof)

    if (begin_end.first != begin_end.second) {
        ; // have some hits
    } else if (nodes_view.end() - nodes_view.begin() <= 8) {

        // have few hits
        begin_end.first = nodes_view.begin();
        begin_end.second = nodes_view.end();

    } else {
        // within the subbucket, collect the hits where we have partial match on upper bits of subkey8
        // by extending in both directions.
        while (begin_end.second < nodes_view.end() && (begin_end.second->subkey8 >> 4) == (key.key8 >> 4)) {
            ++begin_end.second;
        }

        --begin_end.first; // convert start to rend
        while (begin_end.first >= nodes_view.begin() && (begin_end.first->subkey8 >> 4) == (key.key8 >> 4)) {
            --begin_end.first;
        }
        ++begin_end.first; // convert back to start
    }

    return { begin_end.first, begin_end.second };
}

#if 0
// Experimental implementation for: GP-33268

gx::CIndex::nodes_view_t gx::CIndex::at(minword38_t minword) const
{
    ASSERT(m_finalized);

    const auto key = key38_t{ minword.w };

    VERIFY((m_subcounts32 == nullptr) == (m_nodes_ptr == nullptr));
    VERIFY(!m_buckets.empty() ^ (m_subcounts32 != nullptr));

    auto hits = // depending on whether dealing with owning storage or memory-mapped (m_nodes_ptr)
        m_buckets.empty() ?
            nodes_view_t{ (key.key30 == 0 ? 0 : x_subcount(key.key30 - 1)) + m_nodes_ptr,
                                                x_subcount(key.key30)      + m_nodes_ptr }
          : nodes_view_t{ m_buckets[key.key30].data(),
                          m_buckets[key.key30].data() + m_buckets[key.key30].size() };

    const auto eqr = my::equal_range(
        hits.begin(), hits.end(),
        node_t{ seq_oid_t{}, pos1_t{}, key.key8 },
        [](const node_t& a, const node_t& b)
        {
            return a.subkey8 < b.subkey8;
        });
#if 1
    // If did not find any hits for the query, return the whole bucket
    // of approximately matching seeds (just on the 30-bit main key).
    return eqr.first != eqr.second ? nodes_view_t{ eqr.first, eqr.second } : hits;
#else
    // TODO: investigate whether it is preferable to instead further
    // narrow-down hits by progressively matching on the bits of subkey8
    // until have few hits.
    //
    // NB: test changes with "weird beasties" assemblies that don't have good coverage in gxdb.
    // e.g.
    // GCA_000762945.2__6973__Blattella_germanica__German_cockroach
    // GCA_014607495.2__252295__Carposina_sasakii__moths
    // GCF_016086655.2__72036__Lepeophtheirus_salmonis__salmon_louse
    // GCA_001574975.1__1344416__Gonapodya_prolifera_JEL478__monoblepharidomycetes

    int target_bit_bitmask = 0b10000000;
    int upper_bits_bitmask = 0b00000000;

    static const auto seeding_sensitivity_level = get_env("GX_SEEDING_SENSITIVITY_LEVEL", 16ul);
    while(hits.size() > seeding_sensitivity_level) {

        VERIFY((  hits.begin()->subkey8 & upper_bits_bitmask) == (key.key8 & upper_bits_bitmask));
        VERIFY(((hits.end()-1)->subkey8 & upper_bits_bitmask) == (key.key8 & upper_bits_bitmask));

        search_node.subkey8 = uint8_t((key.key8 & upper_bits_bitmask) | target_bit_bitmask);
        const auto it_mid = std::lower_bound(hits.begin(), hits.end(), search_node, by_subkey);

        const auto b = uint8_t(key.key8 & target_bit_bitmask);
        const auto subhits = b ? nodes_view_t{ it_mid, hits.end() }
                               : nodes_view_t{ hits.begin(), it_mid };
        if(subhits.size() == 0) {
            break;
        }

        upper_bits_bitmask |= target_bit_bitmask; // 0b1000000, 0b1100000, ...
        target_bit_bitmask >>= 1;
        hits = subhits;
    }
    return hits;
#endif
}
#endif

/////////////////////////////////////////////////////////////////////////

using db_version_t = uint32_t;
static const db_version_t k_exec_db_version = 100;
// lower two digits is the minor-version.

static const uint64_t k_file_magic_constant = 0xdad9083c0340076c; // random

void CIndex::to_stream(std::ostream& ostr) const
{
    VERIFY(m_finalized);

    ser::to_stream(ostr, &k_file_magic_constant);
    ser::to_stream(ostr, &k_exec_db_version);
    ser::write_padding(ostr);

    {
        std::vector<uint64_t> counts64(1 + m_buckets.size() / k_subcounts64_stride);
        std::vector<uint32_t> counts32(m_buckets.size());

        uint64_t total = 0;
        for (const auto i : irange{ m_buckets.size() }) {
            total += m_buckets[i].size();

            if (i % k_subcounts64_stride == 0) {
                counts64[ i/k_subcounts64_stride ] = total;
            }
            counts32[i] = uint32_t(total - counts64[ i/k_subcounts64_stride ]);
        }

        ser::to_stream(ostr, counts64);
        ser::to_stream(ostr, counts32);
    }

    for (const auto& nodes : m_buckets) {
        // Not as vectors (which would serialize the size at front of each),
        // but as simple concatenatino of nodes.
        ser::to_stream(ostr, nodes.data(), nodes.size());
    }
}

////////////////////////////////////////////////////////////////////////////

void CIndex::from_stream(std::istream& istr)
{
    const uint64_t magic_file_constant = ser::from_stream(istr);
    if (magic_file_constant != k_file_magic_constant) {
        GX_THROW("Unrecognized file content.");
    }

    const db_version_t db_version = ser::from_stream(istr);
    if (db_version/100 != k_exec_db_version/100) {
        GX_THROW("The version of the index is incompatible with the executable.");
    }

    ser::skip_padding(istr);

    auto memistr_p = dynamic_cast<ser::memistream*>(&istr);

    if (!memistr_p) { // deserialize in owning-mode

        const std::vector<uint64_t> counts64 = ser::from_stream(istr);
        VERIFY(counts64.size() == 1 + k_num_buckets / k_subcounts64_stride);

        const std::vector<uint64_t> counts32 = ser::from_stream(istr);
        VERIFY(counts32.size() == k_num_buckets);

        const auto get_subcount = [&](size_t i)
        {
            return counts64[i / k_subcounts64_stride] + counts32[i];
        };

        // will be using our own storage for nodes
        m_buckets.resize(k_num_buckets);
        for (const auto i : irange{ m_buckets.size() }) {
            auto& nodes = m_buckets[i];
            nodes.resize(get_subcount(i) - (i == 0 ? 0 : get_subcount(i - 1)));
            ser::from_stream(istr).to(nodes.data(), nodes.size());
        }

    } else { // deserialize in memory-mapped mode

        // subcounts64
        {
            uint64_t size = ser::from_stream(*memistr_p);
            VERIFY(size == 1 + k_num_buckets/k_subcounts64_stride);
            m_subcounts64 = ser::from_memistream<uint64_t>(*memistr_p, size);
        }

        // subcounts32
        {
            uint64_t size = ser::from_stream(*memistr_p);
            VERIFY(size == k_num_buckets);
            m_subcounts32 = ser::from_memistream<uint32_t>(*memistr_p, size);
        }

        m_nodes_ptr = reinterpret_cast<const node_t*>(memistr_p->begin() + memistr_p->tellg());

        const auto num_nodes = x_subcount(k_num_buckets - 1);

        memistr_p->seekg((const char*)(m_nodes_ptr + num_nodes) - memistr_p->begin());
    }

    m_finalized = true;
}

/////////////////////////////////////////////////////////////////////////
