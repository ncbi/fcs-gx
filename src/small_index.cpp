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
#include "small_index.hpp"
#include "types.hpp"
#include "seq_info.hpp"
#include <mutex>
#include <math.h>
#include <algorithm>

using namespace gx;

using rangeless::fn::operators::operator%=;
using rangeless::fn::operators::operator%;


// 20-bit key
struct key20_t
{
    // subkeys:
    uint16_t key12; // index into m_buckets
    uint8_t  key8;  // attached to the value

    key20_t(uint32_t k = 0)
      : key12(k & Ob1x(12))
      , key8((k >> 12) & Ob1x(8))
    {
        VERIFY(k >> 20 == 0);
    }
};


void gx::CSmallIndex::insert(const minword20_t minword, pos1_t pos)
{
    ASSERT(!m_finalized);

    VERIFY(minword.is_flipped == (pos < 0));

    const key20_t key{ minword.w };

    auto& nodes = m_buckets.at(key.key12).nodes;

    if (nodes.empty()) {
        nodes.reserve(32);
    }

    nodes.push_back(node_t{ pos, key.key8 });
    m_finalized = false;
}

/////////////////////////////////////////////////////////////////////////////


void gx::CSmallIndex::finalize()
{
    ASSERT(!m_finalized);

    for (auto& bucket : m_buckets) {
        offsets_t& offsets = bucket.offsets;
        nodes_t&   nodes   = bucket.nodes;

        nodes %= fn::unstable_sort();
        nodes %= fn::unique_adjacent();

        // fill offsets
        VERIFY(offsets.back() == 0); // expect cleared.

        // collect counts per sub-bucket
        for (const auto& node : nodes)  {
            ++offsets.at(node.subkey8);
        }

        // make cumulative
        for (const auto i : irange{ 1, offsets.size() }) {
            offsets[i] += offsets[i - 1];
        }

        VERIFY(offsets.back() == nodes.size());
    };

    m_finalized = true;
}

/////////////////////////////////////////////////////////////////////////


gx::CSmallIndex::nodes_view_t gx::CSmallIndex::at(minword20_t minword) const
{
    ASSERT(m_finalized);

    const key20_t key{ minword.w };

    const auto& bucket = m_buckets.at(key.key12);

    const auto& offsets = bucket.offsets;
    const auto& nodes = bucket.nodes;

    auto sb_beg = nodes.begin() + (key.key8 == 0 ? 0 : offsets.at(uint8_t(key.key8 - 1)));
    auto sb_end = nodes.begin() +                      offsets.at(key.key8);
    return { sb_beg, sb_end };
}

/////////////////////////////////////////////////////////////////////////
