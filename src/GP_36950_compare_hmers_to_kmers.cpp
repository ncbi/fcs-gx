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

#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// GP-36950 - compare sensitivity of traditional k-mers vs. GX h-mers

#include "types.hpp"
#include "segment.hpp"
#include "serial_util.hpp"
#include <unordered_set>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;
namespace tsv = rangeless::tsv;

namespace GP_36950
{

// return all seqs in fasta, concatenated, delimited by "N"
static std::string read_fasta(const std::string& filename)
{
    auto fasta_istr = ser::open_istream(filename);
    auto ret = std::string{};
    for (const auto& seq : MakeFastaReader(*fasta_istr)) {
        if (   str::contains(seq.defline, "plastid")
            || str::contains(seq.defline, "plasmid")
            || str::contains(seq.defline, "mito")
            || str::contains(seq.defline, "organelle"))
        {
            continue;
        }

        ret += seq.seq;
        ret += "N";
    }
    return ret;
}

// return 38-bit orientation-invariant "traditonal" 19-kmer hash
static uint64_t get_regular_kmer(kmer_bufs_t bufs)
{
    auto w = bufs.twobit & Ob1x(38); // take the 19-mer in lower 38 bits of 2-bit-coding buf.
    return std::min(w, revcomp_twobits(w, 19)); // choose consistent orientation.
}

// return 38-bit orientation-invariant 38-hmer hash, using 1-bit coding alphabet
static uint64_t get_ungapped_hmer(kmer_bufs_t bufs)
{
    auto w = bufs.onebit & Ob1x(38);
    return std::min(w, revcomp_bits(w, 38));
}

// return 38-bit orientation-invariant 56-hmer hash, using 1-bit coding alphabet, and gapped template.
static uint64_t get_gapped_hmer(kmer_bufs_t bufs)
{
    auto w = drop_every_3rd_bit(bufs.onebit) & Ob1x(38);
    return std::min(w, revcomp_bits(w, 38));
}

template<typename GetMerHash>
static void run_one(
    const char*        label, 
    size_t             mer_len,
    GetMerHash         get_mer_hash, // one of the three functions defined above
    const std::string& seq1,
    const std::string& fasta_path1,
    const std::string& seq2,
    const std::string& fasta_path2)
{
    const auto extract_mers = [&](const auto& seq)
    {
        auto ret = std::unordered_set<uint64_t>{};
        size_t n = 0;
        process_kmers(seq, mer_len, [&](size_t pos, kmer_bufs_t bufs)
        {
            ret.insert(get_mer_hash(bufs));
            if (n++ % 1000000 == 0) {
                std::cerr << "."; // progress-bar
            }
        });
        std::cerr << "|";
        return ret;
    };

    const auto& mers1 = extract_mers(seq1);
    const auto& mers2 = extract_mers(seq2);

    auto common_mers = std::unordered_set<uint64_t>{};

    auto get_coverage_and_count = [&](const auto& seq, const auto& other_mers) 
    -> std::pair<size_t, int64_t>
    {
        // pos1_t maxpos goes up to 2G, whereas we are dealing with whole-genome seq
        // that can be > 3G for euks, so we'll split the intervals into 8 buckets,
        // the first one for positions 0-1G, second one for 1G-2G, etc. See k & j below.
        auto ivlss = std::vector<ivls_t>(16);
        size_t n = 0;
        process_kmers(seq, mer_len, [&](size_t pos, kmer_bufs_t bufs)
        {
            if (n++ % 1000000 == 0) {
                std::cerr << "."; // progress-bar
            }

            const auto h = get_mer_hash(bufs);
            if (other_mers.count(h)) {
                const auto k = pos / 1000000000ul;
                const auto j = pos % 1000000000ul;

                const auto ivl = ivl_t{ pos1_t(j + 1), len_t(mer_len) };
                try {
                    ivl_t::push_or_merge(ivlss[k], ivl);
                } catch (...) {
                    std::cerr << pos << " " << k << " " << j << ivlss[k].back().to_string() << " " << ivl.to_string() << "\n";               
                    throw;
                }
                common_mers.insert(h);
            }
        });
        std::cerr << "|";

        return { 
            gx::sum_by(ivlss, L(ivl_t::sum_lens(_))), 
            gx::sum_by(ivlss, L(_.size()))
        };
    };
    const auto cvg1 = get_coverage_and_count(seq1, mers2);
    const auto cvg2 = get_coverage_and_count(seq2, mers1);
    std::cerr << std::endl;

    // round to two decimal places.
    auto round = [](double x) { return double(uint64_t(x * 100 + 0.5)) / 100.0; };

    const double jaccard_similarity_pct = 
        round(    double(common_mers.size()) * 100.0
                / double(mers1.size() + mers2.size() - common_mers.size()));

    std::cout << "\n" << label
              << "\t" << mer_len 

              << "\t" << fasta_path1
              << "\t" << seq1.size()
              << "\t" << mers1.size()
              << "\t" << round(double(cvg1.first) * 100 / double(seq1.size()))
              << "\t" << cvg1.second

              << "\t" << fasta_path2
              << "\t" << seq2.size()
              << "\t" << mers2.size()
              << "\t" << round(double(cvg2.first) * 100 / double(seq2.size()))
              << "\t" << cvg2.second

              << "\t" << common_mers.size()
              << "\t" << jaccard_similarity_pct
              << "\t" << round(double(cvg1.first + cvg2.first) * 100 / double(seq1.size() + seq2.size()))
              << "\n";
}

void run();
void run()
{
    const auto fasta_path1 = get_env("GP_36950_GENOME1", std::string{});
    const auto fasta_path2 = get_env("GP_36950_GENOME2", std::string{});
    const auto seq1 = read_fasta(fasta_path1);
    const auto seq2 = read_fasta(fasta_path2);

    std::cout << "#label\tmer_len\tseq1\tlen1\tmers_set1\tpct_cvg1\tivl_count1\tseq2\tlen2\tmers_set2\tpct_cvg2\tivl_count2\tmers_isectn\tmers_pct_Jaccard\tpct_cvg\n";

    run_one("Traditional k-mers",                 19, get_regular_kmer,  seq1, fasta_path1, seq2, fasta_path2);
    run_one("Reduced-alphabet h-mers (ungapped)", 38, get_ungapped_hmer, seq1, fasta_path1, seq2, fasta_path2);
    run_one("Reduced-alphabet h-mers (gapped)",   56, get_gapped_hmer,   seq1, fasta_path1, seq2, fasta_path2);
}

}
