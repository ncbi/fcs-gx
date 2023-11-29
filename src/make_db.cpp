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
#define RANGELESS_ENABLE_TSV 1

#include "types.hpp"
#include "serial_util.hpp"
#include "index.hpp"
#include "seq_info.hpp"
#include "subbyte_array.hpp"
#include "segment.hpp"

#include <numeric>
#include <unordered_set>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;

/////////////////////////////////////////////////////////////////////////////

namespace gx
{
    ivls_t get_low_complexity_mask_narrow(const iupacna_seq_t&); // declared in align.cpp
}


/////////////////////////////////////////////////////////////////////////////
static ivls_t get_N_runs(const fasta_seq_t& inp_chunk, const len_t min_len = 10)
{
    ivls_t ivls{};
    auto n_run = ivl_t{};

    for (const auto i : irange{ inp_chunk.seq.size() }) {
        const auto c = inp_chunk.seq[i];
        auto pos = pos1_t(inp_chunk.offset + 1 + i);

        if (c != 'N' && c != 'n') {
            continue;
        } else if (n_run.endpos() == pos) {
            n_run.len++; // extend n-run
            continue;
        }

        if (n_run.len >= min_len) {
            ivls.push_back(n_run);
        }

        n_run.pos = pos;
        n_run.len = 1;
    }

    if (n_run.len >= min_len) {
        ivls.push_back(n_run);
    }

    return ivls;
}


/////////////////////////////////////////////////////////////////////////
// Add N-runs, low-complexity, hardmasks, softmasks; subtract those from inp-chunk's range
static ivls_t get_indexable_intervals(
                    const fasta_seq_t& inp_chunk,
                     const locs_map_t& softmask_map,
                     const locs_map_t& hardmask_map)
{
    const auto whole_chunk_ivl = ivl_t{ (int32_t)inp_chunk.offset + 1, (int32_t)inp_chunk.seq.size() };

    ivls_t ivls = get_low_complexity_mask_narrow(inp_chunk.seq)
         % fn::transform L((_.pos += (pos1_t)inp_chunk.offset, std::move(_))) // from chunk's to query coords
         % fn::append(get_N_runs(inp_chunk))
         % fn::append(
                 ivl_t::get_overlapping_v(
                     at_or_default(softmask_map, inp_chunk.seq_id),
                     whole_chunk_ivl))
         % fn::append(
                 ivl_t::get_overlapping_v(
                     at_or_default(hardmask_map, inp_chunk.seq_id),
                     whole_chunk_ivl));
#if 0
    // ivls.clear(); // index the whole chunk
    // NB: the resulting diffs with bacteroides-test looks peculiar:
    // start-coordinats match but the seg length sometimes don't,
    // indicating possible bug related to extensions.
#endif

    if (ivls.empty()) {
        ivls.push_back(whole_chunk_ivl);
        return ivls;
    }

    // convert masked-intervals into indexable-intervals
    ivl_t::sort_and_merge(ivls);
    ivls = ivl_t::invert(std::move(ivls));

    // front and back are 0-length boundaries; extend them all the way.
    {
        VERIFY(ivls.front().len == 0 && ivls.back().len == 0);
        ivls.front().len = ivls.front().pos - 1;
        ivls.front().pos = 1;
        const auto chunk_endpos = pos1_t(inp_chunk.offset + 1 + inp_chunk.seq.size());
        ivls.back().len = std::max(chunk_endpos, ivls.back().endpos()) - ivls.back().pos;
    }

    // constrain to inp-chunk boundaries
    for (auto& ivl : ivls) {
        ivl = ivl_t::intersect(ivl, whole_chunk_ivl);
    }

    // Filter to those that are at least k_word_tlen
    ivls %= fn::where L(_.len >= CIndex::k_word_tlen);

    return ivls;
}

// Extract gene symbol [gene=...] from fasta deflines as formatted in from *cds_from_genomic.fna.gz files.
// Return empty string if can't parse the gene symbol.
static std::string get_gene_symbol_from_defline(const std::string& defline)
{
    auto b = defline.find(" [gene=");
    b = b == std::string::npos ? defline.size() : b + 7;
    const auto e = defline.find(']', b);
    auto gene_symbol = e == std::string::npos ? "" : defline.substr(b, e - b);

#if 0
    // e.g. cds_XP_001338671.2_15793 [gene=si: ch211-223l2.4] [db_xref=GeneID:100002400]
    if (gene_symbol.find(' ') != std::string::npos) {
        std::cerr << "Warning: defline contains a gene symbol with a space:" << defline << "\n\n";
    }
#endif

    return gene_symbol;
}
static const bool test_get_gene_symbol = []
{
    VERIFY("LOC100385748" == get_gene_symbol_from_defline(
        ">lcl|NC_048383.1_cds_XP_035147155.1_1 [gene=LOC100385748] [db_xref=GeneID:100385748] etc."));
    return true;
}();


/////////////////////////////////////////////////////////////////////////////
void gx::MakeDb(     std::istream& fasta_istr,
                     std::istream& seq_id2tax_id_istr,
                     std::istream* taxa_istr_ptr,
                     std::istream* hardmask_istr_ptr,
                     std::istream* softmask_istr_ptr,
                const std::string& out_path) // /path/to/out_db.gxi
{
    auto get_num_locs = [](const locs_map_t& loc_map)
    {
        return fn::cfrom(loc_map)
             % fn::transform L(_.second.size()) 
             % fn::foldl_d([](size_t ret, size_t n) { return ret + n; });
    };

    const auto hardmask_map = hardmask_istr_ptr ? LoadLocsMap(*hardmask_istr_ptr) : locs_map_t{};
    if (hardmask_istr_ptr) {
        std::cerr << "Loaded hardmask map: " << get_num_locs(hardmask_map) << " locs.\n";
    }

    const auto softmask_map = softmask_istr_ptr ? LoadLocsMap(*softmask_istr_ptr) : locs_map_t{};
    if (softmask_istr_ptr) {
        std::cerr << "Loaded softmask map: " << get_num_locs(softmask_map) << " locs.\n";
    }

    /////////////////////////////////////////////////////////////////////////

    using id_tax_vec_t = std::vector<std::pair<seq_id_str_t, tax_id_t>>;
    const auto id_tax_vec = id_tax_vec_t{ [&]
    {
        auto ret = id_tax_vec_t{};
        namespace tsv = rangeless::tsv;
        ConsumeMetalineHeader(seq_id2tax_id_istr, GX_TSV_HEADER__SEQ_ID_MAPPING);
        for (const tsv::row_t& row : tsv::from(seq_id2tax_id_istr)) {
            VERIFY(row.size() == 2);
            ret.emplace_back( seq_id_str_t{             row[0] },
                                  tax_id_t{ tsv::to_num(row[1])});
        }
        std::cerr << "Loaded seq-ids: " << ret.size() << "\n";

        ret %= fn::unique_all(); // also sorts
        return ret;
    }()};

    const tax_map_t tax_map = LoadTaxa(taxa_istr_ptr);
    /////////////////////////////////////////////////////////////////////////

    seq_infos_t seq_infos{};
    {
        // The first element with oid=0 is an empty sentinel
        // (used in conjunction with repeat-markers, so we
        // don't confuse it with a real subject)
        seq_infos.reserve(1 + id_tax_vec.size());
        seq_infos.push_back(seq_info_t{});

        for (const auto& id_and_tax : id_tax_vec) {
            seq_infos.emplace_back(
                    seq_oid_t{},        // oid, will assign below
                    id_and_tax.second,  // tax-id
                    0UL,                // len
                    id_and_tax.first,   // seq-id
                    0UL);               // offset-in-db
        }

        // sort by tax-id, such that all seq-oids that we're about to assign
        // are collated by taxa, as we assume later; then by length(descending).
        seq_infos %= fn::unstable_sort_by L(std::make_pair(_.tax_id, int64_t(_.length)*-1));

        VERIFY(seq_infos.front().tax_id == tax_id_t{});
        for (const auto i : irange{ seq_infos.size() }) {
            seq_infos[i].seq_oid = seq_oid_t{ (uint32_t)i };
        }
    }

    /////////////////////////////////////////////////////////////////////////

    const auto seq_id2oid = [&]
    {
        auto ret = std::map<seq_id_str_t, seq_oid_t>{};
        for (const auto& x : seq_infos) {
            if (!ret.emplace(seq_id_str_t(x.get_seq_id()), x.seq_oid).second) {
                GX_THROW(std::string{} + "Duplicate seq_id unexpected: " + x.get_seq_id());
            }
        }
        return ret;
    }();

    /////////////////////////////////////////////////////////////////////////

    auto in_scope = [&](const fasta_seq_t& inp_chunk) -> bool
    {
        VERIFY(!inp_chunk.seq_id.empty());

        if (!seq_id2oid.count(inp_chunk.seq_id)) {

            std::cerr << "Skipping not-in-scope: " << inp_chunk.seq_id << "\n";
            return false;

        } else if (str::contains(inp_chunk.defline, "alternate locus")) {

            std::cerr << "Skipping alt-locus: " << inp_chunk.defline << "\n";
            return false;

        } else if (str::contains(inp_chunk.defline, "UNVERIFIED:")) {

            std::cerr << "Skipping UNVERIFIED: " << inp_chunk.defline << "\n";
            return false;

        } else if (    hardmask_map.count(inp_chunk.seq_id)
                   && !hardmask_map.at(inp_chunk.seq_id).empty()
                   &&  hardmask_map.at(inp_chunk.seq_id).front() == ivl_t::s_whole_ivl)
        {
            std::cerr << "Skipping whole-hardmasked: " << inp_chunk.seq_id << "\n";
            return false;
        }

        return true;
    };

    /////////////////////////////////////////////////////////////////////////

    auto update_seq_infos = [ &, last_seq_id = seq_id_str_t() ]
                            (fasta_seq_t inp_chunk) mutable -> fasta_seq_t
    {
        if (inp_chunk.seq_id != last_seq_id) { // reached next query

            const auto oid = seq_id2oid.at(inp_chunk.seq_id);
            auto& si = seq_infos[+oid];

            VERIFY(si.get_seq_id() == inp_chunk.seq_id);

            //std::cerr << "\rIndexing tax_id: " << +si.tax_id << "; seq-id: " << inp_chunk.seq_id << "; oid:" << +oid << "...";

            if (!last_seq_id.empty()) {
                const auto& prev_si = seq_infos.at(seq_id2oid.at(last_seq_id));
                si.offset_in_seq_db = prev_si.offset_in_seq_db
                                    + sbj_seq_t::subbyte_array_t::s_size_in_bytes(prev_si.length);
            }
            last_seq_id = inp_chunk.seq_id;
        }

        inp_chunk.seq_oid = seq_id2oid.at(inp_chunk.seq_id);
        seq_infos[+inp_chunk.seq_oid].length = uint32_t(inp_chunk.offset + inp_chunk.seq.size());

        return inp_chunk;
    };


    /////////////////////////////////////////////////////////////////////////

    CIndex index{};

    std::atomic_size_t num_bases_in_scope = {};

    // will be doing this part in parallel
    auto update_index_and_translate_to2bit =
        [&](fasta_seq_t inp_chunk) -> std::pair<fasta_seq_t, sbj_seq_t>
    {
        const bool is_in_frame_CDS =  // GP-34257
            str::contains(inp_chunk.seq_id, "cds_")
         && (   str::startswith(inp_chunk.seq, "ATG")
             || (   inp_chunk.seq.size() % 3 == 0
                 && (   str::endswith(inp_chunk.seq, "TAA")
                     || str::endswith(inp_chunk.seq, "TAG")
                     || str::endswith(inp_chunk.seq, "TGA"))));

        // Will use small k_stride for proks, and 2*k_stride for euks
        // to give more sensitivity to proks.
        const auto tax_id = seq_infos.at(inp_chunk.seq_oid).tax_id;
        const auto stride =
            is_in_frame_CDS ? 18  // multiple of 3 - frame-preserving; otherwise frame-rotating.
          : CIndex::k_stride * (tax_map.count(tax_id) && tax_map.at(tax_id).is_prok_or_virus() ? 1 : 2);


        // Apply hardmasking.
        const auto whole_chunk_ivl = ivl_t{ (int32_t)inp_chunk.offset + 1, (int32_t)inp_chunk.seq.size() };
        for (const auto& ivl : at_or_default(hardmask_map, inp_chunk.seq_id))
            if (const auto sub_ivl = ivl_t::intersect(ivl, whole_chunk_ivl); sub_ivl.len > 0)
                for (const auto i : irange{ sub_ivl.pos, sub_ivl.endpos() })
        {
            inp_chunk.seq[i - 1 - inp_chunk.offset] = 'N';
        }

        const auto ivls = get_indexable_intervals(inp_chunk, softmask_map, hardmask_map);
        num_bases_in_scope += sum_by(ivls, L(_.len));

        for (const auto& ivl : ivls) {
            const size_t start = size_t(ivl.pos      - 1) - inp_chunk.offset; // 0-based chunk-local coordinate
            const size_t end   = size_t(ivl.endpos() - 1) - inp_chunk.offset; // 0-based chunk-local coordinate
            VERIFY(start < inp_chunk.seq.size());

            for(auto kmer = kmer_ci_t(inp_chunk.seq, CIndex::k_word_tlen, start); kmer && kmer.i < end; ++kmer) {
                size_t i_pos = inp_chunk.offset + kmer.i;

                if ((i_pos + 1 - CIndex::k_word_tlen) % stride != 0) {
                    // Choosing positions where the word-start is in proper coding phase,
                    // in case we are indexing a CDS with a frame-preserving %3==0 stride.
                    continue;
                }

                const auto hmer = CIndex::hmer38_t{ kmer.buf };
                const auto pos1 = as_ivl(i_pos, CIndex::k_word_tlen, hmer.is_flipped).pos;
                index.insert(hmer, inp_chunk.seq_oid, pos1);
            }
        }

        auto nuc2_seq = sbj_seq_t{ inp_chunk.seq };
        return std::make_pair(std::move(inp_chunk), std::move(nuc2_seq));
    };

    /////////////////////////////////////////////////////////////////////////

    static const size_t k_chunk_stride = 100000UL;
    static const size_t k_chunk_overlap = 100;

    /////////////////////////////////////////////////////////////////////////
    // We will be rolling our own super-simple sequence local-storage for
    // aligner-stages that require access to sequence. I tried using blastdb
    // and object-manager for that, but turns out it doesn't scale (access
    // is extremely slow for TB-sized databases. Related SB-3063.
    //
    // Instead, we will encode each query as 2-bit-coding seq, and dump them
    // all to one file, recording the offset for each seq-id.
    // When using the index, we will memory-map the file and then can obtain
    // the view into the desired seq-interval with zero-overhead.

    VERIFY(str::endswith(out_path, ".gxi"));
    std::ofstream o_seq{ str::replace_suffix(out_path, ".gxi", ".gxs") };
    VERIFY(o_seq);

    auto dump_seq = [&, prev_oid = seq_oid_t()](std::pair<fasta_seq_t, sbj_seq_t> inp) mutable
    {
        // Dump the 2-bit-encoded sequence chunk to o_seq,
        // adjusting it to remove the chunk-overlaps, such that
        // the output file contains concatenation of all input seqs.

        const fasta_seq_t& chunk = inp.first;
        sbj_seq_t&     nuc2_seq = inp.second;

        if (chunk.seq.empty()) {
            return;
        }

        VERIFY(chunk.seq.size() == nuc2_seq.size());

        // To ensure that truncated-to-stride nuc2_seq internal storage
        // is aligned to byte boundary (one byte fits 4 bps)
        static_assert(k_chunk_stride % 4 == 0, "");

        const auto si = seq_infos.at(chunk.seq_oid);

        const size_t nonoverlapping_size = // excluding overlap-tail
              chunk.seq.size() == k_chunk_stride + k_chunk_overlap
           && chunk.offset + chunk.seq.size() <= si.length ? k_chunk_stride
                                                           : nuc2_seq.size();
        nuc2_seq.resize(nonoverlapping_size);

        // std::cerr << chunk.seq_id << "; oid=" << chunk.seq_oid
        //           << "@" << chunk.offset << "; chunk-size:" << chunk.seq.size() << "\n";

        // Upon reaching the next sequence, verify that tellp
        // is consistent with the computed si.offset_in_seq_db
        // (i.e. concatenation of truncated chunks is consistent
        // with sum of query lengths).
        VERIFY(chunk.seq_oid == prev_oid || si.offset_in_seq_db == size_t(o_seq.tellp()));
        prev_oid = chunk.seq_oid;

        VERIFY(o_seq.write(reinterpret_cast<const char*>(nuc2_seq.data().data()), (std::streamsize)nuc2_seq.data().size_bytes()));
    };

    /////////////////////////////////////////////////////////////////////////

    auto t = timer{};

    static const auto num_cores = get_env("GX_NUM_CORES", std::min(32U, std::thread::hardware_concurrency() - 1));
    auto fasta_reader = MakeFastaReader(fasta_istr, k_chunk_stride, k_chunk_overlap);
    if (num_cores == 0) {
        // In num_cores==0 will do the entire processing using the simple loop, in-this-thread,
        // without setting up the pipeline to get simpl(er) stack traces when debugging.
        for (fasta_seq_t qry_chunk : fasta_reader) {
            if (in_scope(qry_chunk)) {
                dump_seq(
                    update_index_and_translate_to2bit(
                        update_seq_infos(
                            std::move(qry_chunk)))); // TODO: missing CDS-group-and-filter logic here (see below)
            }
        }
    } else {
        std::move(fasta_reader)
# if 1
        // GP-34257. For CDS, if can parse gene from defline, group by gene
        // and pick one longest per gene (expecting them to be adjacent in the input stream)
        // NB: may be split into multiple chunk if large, e.g. TTN
      % fn::group_adjacent_if([](const fasta_seq_t& a, const fasta_seq_t& b)
        {
            const auto sym_a = get_gene_symbol_from_defline(a.defline);
            const auto sym_b = get_gene_symbol_from_defline(b.defline);
            return !sym_a.empty() && !sym_b.empty() && sym_a == sym_b;
        })
      % fn::transform([](std::vector<fasta_seq_t> group)
        {
            VERIFY(!group.empty());
            const auto id_of_longest = std::max_element(group.begin(), group.end(), BY(_.seq.size()))->seq_id;
            return std::move(group) % fn::where L(_.seq_id == id_of_longest);
        })
      % fn::concat()
#endif
      % fn::where(in_scope)
      % fn::transform(update_seq_infos)
      % fn::transform_in_parallel(update_index_and_translate_to2bit).queue_capacity(num_cores)
      % fn::for_each(dump_seq);
    }

    std::ofstream o_seq_info{ str::replace_suffix(out_path, ".gxi", ".seq_info.tsv") };
    VERIFY(o_seq_info);
    size_t num_seqs = 0;
    size_t num_bases = 0;
    o_seq_info << "#seq-id\tlength\ttax-id\tspecies\tcommon-name\tblast-div\tgx-div" << std::endl;
    for (const auto& si : seq_infos)
        if (+si.tax_id && si.length > 0) // skip seqs that were out of scope (did not compute length for them)
    {
        const auto ti = tax_map.at(si.tax_id);
        o_seq_info  << si.get_seq_id()
            << "\t" << si.length
            << "\t" << +si.tax_id
            << "\t" << ti.species
            << "\t" << ti.common_name
            << "\t" << ti.blast_taxdiv
            << "\t" << ti.gx_taxdiv
            << "\n";

        num_seqs++;
        num_bases += si.length;
    };
    o_seq_info.close();

    /////////////////////////////////////////////////////////////////////////

    std::cerr << "Indexed "
              << float(num_bases)/1e9f          << " Gbp; "
              << float(num_bases_in_scope)/1e9f << " Gbp in scope; "
              << num_seqs                       << " seqs in "
              << float(t)/60                    << " minutes; "
              << float(num_bases)/1e6/float(t)  << " Mbp/s.\n\n";

    {
        std::time_t time_now = std::time(0);
        char buf[128];
        strftime(buf, sizeof(buf), "%Y-%m-%d", std::localtime(&time_now));

        VERIFY(str::endswith(out_path, ".gxi"));
        std::ofstream o_meta{ str::replace_suffix(out_path, ".gxi", ".meta.jsonl") };
        o_meta << "{\"build-date\":\""  << std::string{ buf }
               << "\", \"seqs\":"       << num_seqs
               << ", \"Gbp\":"          << float(num_bases)/1e9f
               << "}\n";
    }

    t = timer{};
    std::cerr << "Finalizing..." << std::endl;
    index.finalize(seq_infos);
    std::cerr << "Finalized index in " << float(t)/60 << " minutes.\n\n";

    t = timer{};
    std::cerr << "Serializing..." << std::endl;
    std::ofstream ostr(out_path);
    VERIFY(ostr);
    ser::to_stream(ostr, seq_infos);
    index.to_stream(ostr);
    std::cerr << "Serialized index in " << float(t)/60 << " minutes.\n\n";

    t = timer{};
    std::cerr << "Destroying..." << std::endl;
    index = CIndex{};
    std::cerr << "Destroyed index in " << float(t)/60 << " minutes.\n\n";
}


#if 0
The orf-based sampling for indexing business does not seem to improve sensitivity

//////////////////////////////////////////////////////////////////////////
//
namespace orfs
{
    struct ivl_t
    {
        pos1_t start;
        pos1_t len;
        pos1_t endpos() const
        {
            return start + len;
        }

        std::string to_string() const
        {
            return "[" + std::to_string(start) + ".." + std::to_string(start + len - 1) + "]";
        }
    };
    using ivls_t = std::vector<ivl_t>;

    static inline
    void GetOrfsInFrame(const fasta_seq_t& seq,
                                    pos1_t start_pos,
                                    pos1_t end_pos,
                                       int min_orf_len,
                                   ivls_t& dest)
    {
        std::string codon{"NNN"};

        pos1_t last_stop_codon_pos = 1;
        pos1_t start_codon_pos = 1; // when equal to last_stop_codon_pos means not yet encountered
        for (pos1_t pos = start_pos; pos + 2 < end_pos; pos += 3) {
            codon[0] = seq.at1_local(pos);
            codon[1] = seq.at1_local(pos + 1);
            codon[2] = seq.at1_local(pos + 2);

            const bool is_stop_codon = codon == "TAG"
                                    || codon == "TAA"
                                    || codon == "TGA";

            const bool is_start_codon = codon == "ATG";

            if (is_start_codon && start_codon_pos == last_stop_codon_pos) {
                start_codon_pos = pos;
            }

            if (!is_stop_codon) {
                continue;
            }

            auto start_pos = start_codon_pos;
            while ((pos - start_pos) % 3 != 0) { // adjust to be in-frame
                ++start_pos;
            }
            auto ivl = ivl_t{ start_pos, pos + 3 - start_pos }; // +3 to convert to endpos of this stop-codon
            if (ivl.len >= min_orf_len) {
                dest.push_back(ivl);
            }
            last_stop_codon_pos = pos;
            start_codon_pos = pos;
        }
    }

    static inline
    ivls_t GetOrfsInSixFrames(const fasta_seq_t& seq, int min_orf_len)
    {
        ivls_t ret;
        auto stop_pos = pos1_t{ int32_t(seq.seq.size()) };
        GetOrfsInFrame(seq, pos1_t{1}, stop_pos+1, min_orf_len, ret);
        GetOrfsInFrame(seq, pos1_t{2}, stop_pos+1, min_orf_len, ret);
        GetOrfsInFrame(seq, pos1_t{3}, stop_pos+1, min_orf_len, ret);

        GetOrfsInFrame(seq, stop_pos * -1 + 0, pos1_t{0}, min_orf_len, ret);
        GetOrfsInFrame(seq, stop_pos * -1 + 1, pos1_t{0}, min_orf_len, ret);
        GetOrfsInFrame(seq, stop_pos * -1 + 2, pos1_t{0}, min_orf_len, ret);
        return ret;
    }

    // vector will have values set to 1 at absolute 0-based positions of word stop-positions.
    //
    // In every ORF                ...TAG
    // we will take the word ...xx-xx-xx (last dinucleotide is the first two bases of the stop codon);
    //
    // Then we will sample the words going upstream at stride (in-frame).
    static inline
    std::vector<bool> SamplePointsForIndexing(const fasta_seq_t& seq, int coding_stride = 18, int noncoding_stride=19)
    {
        VERIFY(coding_stride % 3 == 0);
        VERIFY(noncoding_stride % 3 != 0);
        VERIFY(noncoding_stride > coding_stride);

        std::vector<bool> ret(seq.seq.size(), false);
        auto ivls = GetOrfsInSixFrames(seq, 300);

        ivls %= fn::sort_by L(_.start);

        for (ivl_t& ivl : ivls) {
            VERIFY(ivl.len % 3 == 0);
            VERIFY(ivl.start > 0 || (ivl.start < 0 && ivl.start + ivl.len <= 0));

            // Drop the last nucleotide of the stop-codon (last nucleootide of a codon is not part of the word template).
            ivl.len--;

            // Truncate the start such that when we sample words of k_word_tlen at coding_stride,
            // the last word's end corresponds to interval's end, i.e. tiled sampled words in this
            // interval will fit exactly in the truncated interval, regardless of orientation.
            {
                const auto prefix_overhang = (ivl.len - CIndex::k_word_tlen) % coding_stride;
                ivl.start += prefix_overhang;
                ivl.len -= prefix_overhang;
                VERIFY((ivl.len - CIndex::k_word_tlen) % coding_stride == 0);
            }

            // make strand-agnostic
            if (ivl.start < 0) {
                ivl.start = flip_pos1(ivl.start, ivl.len);
            }
        }

        // drop shadowed
        {
            using fn::operators::operator%=;
            ivls %= fn::sort_by L(_.start);
            ivl_t* prev = nullptr;
            for (ivl_t& ivl : ivls) {
                if (!prev) {
                    prev = &ivl;
                } else if (prev->endpos() >= ivl.endpos()) {
                    ivl.len = 0;
                } else {
                    prev = &ivl;
                }
            }
            ivls %= fn::where L(_.len > 0);
        }

        for (const ivl_t& ivl : ivls)
            for (auto pos = ivl.start + ivl.len; pos >= ivl.start; pos -= coding_stride)
        {
            size_t i = pos > 0 ? pos - 2 // -1 to convert endpos to stop-pos, and -1 to make 0-based
                               : ((pos - CIndex::k_word_tlen) * -1) - 1; // convert to start, flip, make 0-based
            ret.at(i) = 1;
        }

#if 1
        size_t last_set_i = -1;
        for (const auto i : irange{ ret.size() })
            if (ret[i] || i + i == ret.size())
        {
            if (last_set_i + noncoding_stride < i) {
                auto offset = ((i - last_set_i) % noncoding_stride) / 2;
                for (size_t j = last_set_i + offset; j < i; j += noncoding_stride) {
                    ret[j] = 1;
                }
            }

            last_set_i = i;
        }
#endif

        return ret;
    }
    // Sample points from every orf, starting with stop, going upstream with stride = 15 (or 12 or 18?)
    //
    //
    //

}
//////////////////////////////////////////////////////////////////////////
#endif
