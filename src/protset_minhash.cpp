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
#define RANGELESS_ENABLE_TSV 1
#include "util.hpp"
#include "segment.hpp"
#include <set>
#include <queue>
#include <smmintrin.h>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;

namespace tsv = rangeless::tsv;

////////////////////////////////////////////////////////////////////////////////////////
// IUPAC-AA -> [0..4) or [0..16)
static std::string make_reduced_alphabet(const std::string& alphabet_grouping, uint8_t dflt)
{
    auto table = std::string(256, dflt);
    int8_t group_ord = 0;
    for (const auto& c : alphabet_grouping) {
        if (c == ' ' ) {
            group_ord++;
        } else {
            table[c] = group_ord;
            table[std::toupper(c)] = group_ord;
        }
    }
    return table;
}

////////////////////////////////////////////////////////////////////////////////////////
// Input is protetin fasta for a taxon.
// Output is 1kb base-64 min-hash signature.
std::string gx::MakeProtsetMinhash(std::istream& fasta_istr)
{
    constexpr auto mer_bitwidth = 30;
    constexpr auto minhash_size = 1024ul;

    static const auto important_aas = []
    {
        auto tbl = std::vector<bool>(256, 0);
        // https://pubmed.ncbi.nlm.nih.gov/8601843 PDEKGC
        // Alternatively BLOSUM high-weight AAs on diagonal: CHPW,
        // but seems to produce lower Jaccard concordances.
        static const std::string s = "PDEK";

        //std::cerr << "Anchor AAs: " << s << "\n";
        for (const auto aa : s) {
            tbl[aa] = true;
        }
        return tbl;
    }();

#if 1
    constexpr auto alphabet_bitwidth = 2;
    static const auto reduced_alphabet = make_reduced_alphabet("CFYW* MLIV GPATSN EDHQRK", 3);
#else
    constexpr auto alphabet_bitwidth = 4;
    static const auto reduced_alphabet = make_reduced_alphabet("C F Y W* M L IV G P ATSN E D H Q R K", 9);
#endif
    constexpr auto min_i = (mer_bitwidth + alphabet_bitwidth - 1) / alphabet_bitwidth - 1;

    // https://www.researchgate.net/publication/221596285_Bottom-k_sketches_Better_and_more_efficient_estimation_of_aggregates
    auto unique_mers_count = 0ul;
    auto bottom_sketch = std::priority_queue<uint64_t>{}; // capped at 100k
    {
        auto seen_seqs = std::set<uint64_t>();
        auto seen_mers = std::vector<bool>(1 << mer_bitwidth);

        for (fasta_seq_t seq : MakeFastaReader(fasta_istr, -1ul, 0ul, validate_iupacna_t::no))
            if (seen_seqs.insert(uint64_hash(seq.seq.data(), seq.seq.size())).second)
        {
            seq.seq.push_back('*');
            auto buf = uint64_t{ 0 };

            for (const auto i : irange{ seq.seq.size() }) {
                const auto aa = seq.seq[i];

                buf = buf << alphabet_bitwidth | reduced_alphabet[aa];
                const auto mer = buf & Ob1x(mer_bitwidth);

                if (i < min_i || !important_aas[seq.seq[i - (min_i / 2)]] || seen_mers[mer]) {
                    continue;
                }

                ++unique_mers_count;
                seen_mers[mer] = true;
                bottom_sketch.push(uint64_hash(mer));

                if (bottom_sketch.size() > minhash_size*100) {
                    bottom_sketch.pop();
                }
            }
        }
    }

    // Note: we are computing the minhashes over the bottom sketch as optimization.
    //
    // NB: bottom-sketch-size is 100x of the minhash to minimize collisions
    // (same bottom-sketch element ending up minimal in different hash-families).
    //
    // https://en.wikipedia.org/wiki/MinHash

    auto minhash_sig = std::vector<uint64_t>(minhash_size, uint64_t(-1));
    while (!bottom_sketch.empty()) {
        const auto h = bottom_sketch.top();
        bottom_sketch.pop();
        for (const auto i : irange{ minhash_sig.size() }) {
            minhash_sig[i] = std::min(minhash_sig[i], uint64_hash(h ^ (i << 32)));
        }
    }

    // Truncating each hash to 6 bits, so we can nicely encode it in base64.
    // Probability of a collision is 1/64, so this will boost the Jaccard index estimate by (1-J)/64.
    static const char* base64_tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                    "abcdefghijklmnopqrstuvwxyz"
                                    "0123456789+/";
    (void)unique_mers_count;
    //std::cerr << "Unique mers count:" << unique_mers_count << "\n";

    return unique_mers_count == 0 ? "NULL"
         : fn::cfrom(minhash_sig)
         % fn::transform L(base64_tbl[_ & 63])
         % fn::to(std::string{});
}

// https://stackoverflow.com/questions/15313646/fast-counting-the-number-of-equal-bytes-between-two-arrays
static size_t count_equal_bytes(const std::string& s, const std::string& t)
{
    VERIFY(s.size() == t.size());
    VERIFY(s.size() % 16 == 0);

    __m128i vsum = _mm_set1_epi32(0);
    for (size_t i = 0; i < s.size(); i += 16) {
        __m128i vs, vt, v, vh, vl, vtemp;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
        // NB: _mm_loadu_si128 does not require alignment, so cast-align warning is irrelevant.
        vs    = _mm_loadu_si128((const __m128i *)&s[i]); // load 16 chars from input
        vt    = _mm_loadu_si128((const __m128i *)&t[i]);
#pragma clang diagnostic pop
        v     = _mm_cmpeq_epi8(vs, vt);         // compare
        vh    = _mm_unpackhi_epi8(v, v);        // unpack compare result into 2 x 8 x 16 bit vectors
        vl    = _mm_unpacklo_epi8(v, v);
        vtemp = _mm_madd_epi16(vh, vh);         // accumulate 16 bit vectors into 4 x 32 bit partial sums
        vsum  = _mm_add_epi32(vsum, vtemp);
        vtemp = _mm_madd_epi16(vl, vl);
        vsum  = _mm_add_epi32(vsum, vtemp);
    }

    // get sum of 4 x 32 bit partial sums
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));
    return _mm_cvtsi128_si32(vsum);
}


void gx::PairwiseCompareMinHashes(std::istream& istr, std::ostream& ostr)
{
    static const auto min_reportable_threshold = get_env("GX_PROT_MINHASH_MIN_THRESHOLD", 0.3f);

    auto inps = std::vector<std::pair<std::string, std::string>>{};

    for (const tsv::row_t& row : tsv::from(istr)) {
        VERIFY(row.size() == 2);
        if (row[1] != "NULL") {
            VERIFY(row[1].size() == 1024);
            inps.emplace_back(std::make_pair(row[0], row[1]));
        }
    }

    auto t = timer{};

    // TODO: can parallelize, but for now fast-enough.
    for (const auto i : irange{ inps.size() })
        for (const auto j : irange{ i + 1, inps.size() })
    {
        auto f = (float)count_equal_bytes(inps[i].second, inps[j].second) / 1024.0f;
        f = std::max(f - (1 - f) / 64, 0.0f); // adjust for expected number of collisions
        if (f >= min_reportable_threshold) {
            ostr << inps[i].first << "\t" << inps[j].first << "\t" << f << "\n";
        }
    }

    std::cerr << "Processed " << inps.size()*(inps.size()-1)/2
              << " pairwise comparisons of " << inps.size()
              << " minhashes in " << t << "s.\n";
}
