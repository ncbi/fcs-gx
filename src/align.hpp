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
#include "segment.hpp"
#include "seq_info.hpp"
#include "small_index.hpp"

namespace gx
{


/////////////////////////////////////////////////////////////////////////////
// defined in refine_segs.cpp
void UngappedExtendSegsInPlace(const fasta_seq_t& qry_seq, const sbj_seq_t& sbj_seq, segs_view_t segs);

segment_t FindDownstreamContinuation(
                   const fasta_seq_t& qry_seq,
                     const sbj_seq_t& sbj_seq,
                            segment_t inp_seg,
                              uint8_t kmer_len = 24);

void GappedExtend(
                           const fasta_seq_t& qry_seq,
    const std::function<sbj_seq_t(seq_oid_t)> get_sbj_seq,
                                  segments_t& segs,
                                  const len_t min_len = 10);

/////////////////////////////////////////////////////////////////////////////
// defined in round2.cpp
CSmallIndex MakeQueryIndex(const fasta_seq_t& qry, CSmallIndex index = {});


std::vector<tax_id_t> SelectTaxaForRound2(
                         const segments_t& segs,
                        const seq_infos_t& sbj_infos,
                          const tax_map_t& tax_map);

segments_t SeedRound2( const CSmallIndex& index,
                       const fasta_seq_t& qry,
                         const sbj_seq_t& sbj_seq,
                          const seq_oid_t sbj_oid,
                       const segs_view_t& segs1,
                               segments_t ret = {});
}
