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

#include "util.hpp"
#include "types.hpp"
#include "ext/unordered_dense.h"
//#include "ext/robin_hood.h"
#include "serial_util.hpp"

#include <iomanip>
#include <queue>
#include <bitset>

using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;
using fn::operators::operator<<=;

/////////////////////////////////////////////////////////////////////////////

static const size_t k_fasta_chunk_stride = 1000000UL;
static const size_t k_fasta_chunk_overlap = 50;

using word_t = uint64_t;
using words_t = std::deque<word_t>;

#if 0
template<typename Xs, typename Comp>
void par_sort(Xs& xs, const Comp& comp, size_t num_threads, size_t min_size_for_async)
{
	const auto n = xs.end() - xs.begin();

	if (num_threads <= 1 || n < min_size_for_async) {
        std::sort(xs.begin(), xs.end(), std::ref(comp));
	} else {
		const auto half_num_threads = num_threads / 2;
		const auto mid = xs.begin() + n / 2;
		auto fut1 = std::async(std::launch::async, [&]
		{
			par_sort(xs.begin(), mid, std::ref(comp), half_num_threads);
		});
		par_sort(mid, xs.end(), std::ref(comp), num_threads - half_num_threads);
		fut1.get();
		std::inplace_merge(xs.begin(), mid, xs.end());
	}
}
// TODO: implement misra-gries summary as vector-of-items and hash-map of counts.
// If an incoming item is in hash-map, increment the count, otherwise push into vec.
//
// If vec grew beyond capacity, sort it; put non-unique items into hash-map; drop singletons
#endif

void gx::ExtractConsensusRepeats(const std::string& genome_fasta_path, std::ostream& ostr)
{
    static const size_t kmer_len = 64; // one-bit-coding
    static const size_t bitwidth = kmer_len;
    static const uint64_t kmer_mask = Ob1x(bitwidth);

    const auto revcomp = L(revcomp_bits(_, kmer_len));
    const auto get_minword = L(std::min(_, revcomp(_)));

    using counts_map_t = ankerl::unordered_dense::map<word_t, uint32_t>;
    //using counts_map_t = robin_hood::unordered_map<word_t, uint32_t>;


    const auto misra_gries_prune_step = [&](counts_map_t& m, size_t capacity) -> int
    {
        const auto orig_size = m.size();
        if (m.size() < capacity) {
            return 0;
        }

        auto m2 = counts_map_t{};
        m2.reserve(capacity / 20);
        for (auto& kv : m) {
            if (kv.second > 1) {
                m2[kv.first] = kv.second - 1;
            }
        }
        m = std::move(m2);
        std::cerr << "Misra-Gries prune-step: " << orig_size << " -> " << m.size() << "      \n";
        return 1;
    };

    const counts_map_t counts = [&]
    {
        auto m = counts_map_t{};
        uint32_t num_misra_gries_prunings = 0;

        static const auto cached_kmers_path = get_env("GX_CACHED_KMERS_FILE", std::string{});
        if (!cached_kmers_path.empty()) {
            std::ifstream ifstr{ cached_kmers_path };
            errno = 0;
            if (ifstr) {
                namespace tsv = rangeless::tsv;
                for (const tsv::row_t& row : tsv::from(ifstr)) {
                    m[tsv::to_num(row[0])] = tsv::to_num(row[1]);
                }
                return m;
            }
        }

        static const auto capacity = get_env("GX_KMERS_MISRA_GRIES_SUMMARY_CAPACITY", 200000000ul);
        m.reserve(capacity);

        /* TODO: transform in parallel: seq -> counts_map_t
         * Fold: Merge each counts_map_t into master-map; while size > capacity, apply prune-step.
         */

        size_t num_low_complexity_words = 0;
        auto fasta_istr_ptr = ser::open_istream(genome_fasta_path);
        for (auto seq : MakeFastaReader(*fasta_istr_ptr, k_fasta_chunk_stride, k_fasta_chunk_overlap)) {
            if (str::contains(seq.defline, "alternate locus")) {
                // for GRCh there are many alt-loci, which makes the counts over-represented.
                continue;
            }

            process_kmers(seq.seq, kmer_len, [&](size_t, kmer_bufs_t bufs)
            {
                num_misra_gries_prunings += misra_gries_prune_step(m, capacity);
                const auto mw = get_minword(bufs.onebit & kmer_mask);


                // Will skip words that are rich in 00* 11* 10*, i.e. have high bitcount skew
                const auto k = std::max(
                        __builtin_popcountll(mw),
                        __builtin_popcountll(mw ^ 0xAAAAAAAAAAAAAAAAul)); // flip every other bit

                if (k >= 5 && k <= 59) { // ignoring words with high bit-skew
                    ++m[mw];
                } else {
                    ++num_low_complexity_words;
                }
            });
        }

        std::cerr << "Final prune step...\n";
        num_misra_gries_prunings += misra_gries_prune_step(m, m.size());
        const auto max_count = m.empty() ? 0 
                             : std::max_element(
                                     m.begin(), m.end(), 
                                     fn::by::make_comp L(_.second)
                               )->second;
        std::cerr << "Loaded; map-size:" << m.size() 
                  << "; prunings: "      << num_misra_gries_prunings
                  << "; max-count: "     << max_count
                  << "; lc-words: "      << num_low_complexity_words
                  << "\n";

        for (auto& kv : m) {
            kv.second += num_misra_gries_prunings;
        }

        if (!cached_kmers_path.empty()) {
            std::ofstream ofstr{ cached_kmers_path };
            for (const auto& kv : m) {
                ofstr << kv.first << "\t" << kv.second << "\n";
            }
        }

        return m;
    }();

    const auto get_count = [&](word_t w) -> counts_map_t::mapped_type
    {
        return at_or_default(counts, get_minword(w));
    };

    // Get most frequent word that continues w (overhangs by 1bp).
    const auto get_next = [&](word_t w) -> word_t
    {
#if 0
        const word_t w0 = (w << 1) & kmer_mask; // LSBs frame contains 0, i.e. A
        const word_t w1 = w0 | 1u;
        return std::max(
            std::make_pair(get_count(w0), w0),
            std::make_pair(get_count(w1), w1)
        ).second; 
#else
        // Instead of just looking at the next word position ahead,
        // (two candidates with 1 and 0 in the LSB),
        // we get better extensions if we look several bits ahead
        // (e.g. 4), and select the best continuation among the 2^4 candidates.

        const word_t w0 = (w << 4) & kmer_mask;
        auto best = std::make_pair(uint32_t(0), word_t{});
        for (size_t i = 0; i < 16; i++) {
            const auto w_i = w0 | i;
            best = std::max(best, std::make_pair(get_count(w_i), w_i));
        }

        return ((w << 1) & kmer_mask) | ((best.second >> 3) & 1);
#endif
    };

    auto seen_set = ankerl::unordered_dense::set<word_t>{};
    const auto is_seen = [&](word_t w) -> bool
    {
        const auto mw = get_minword(w);
        if (seen_set.count(mw)) {
            return true;
        }

        // also check every neighbor within 1 mismatch.
        for (size_t i = 0; i < kmer_len; i++) {
            if (seen_set.count(get_minword(mw ^ (1ul << i)))) {
                return true;
            }
        }
        return false;
    };

    const auto mark_seen = [&](word_t w) -> void
    {
        const auto mw = get_minword(w);
        seen_set.insert(mw);

#if 0
        // also mark every neighbor within 1 mismach
        for (size_t i = 0; i < kmer_len; i++) {
            seen_set.insert(get_minword(mw ^ (1ul << i)));
        }
#endif
    };

    static const auto min_seed_count = get_env("GX_CONSENSUS_REPEAT_MIN_SEED_SUP", 25ul);
    static const auto min_ext_count  = get_env("GX_CONSENSUS_REPEAT_MIN_EXT_SUP", 10ul);
    static const auto min_len        = get_env("GX_CONSENSUS_REPEAT_MIN_LEN", kmer_len * 2);

    struct repeat_t
    {
        word_t seed;
        words_t kmers;
        size_t weight;

        using residue_counts_t = std::array<uint32_t, 4>;

        // size == this->len(); counts of A/C/G/T at every pos.
        std::vector<residue_counts_t> residue_counts;

        size_t len() const
        {
            return kmers.size() + kmer_len - 1;
        }

        void validate() const
        {
            // Verify that all kmers overlap by len(kmer)-1 bp.
            for (size_t i = 1; i < kmers.size(); i++) {
                VERIFY(( kmers[i] == (((kmers[i - 1] << 1) | (kmers[i] & 1)) & kmer_mask) ));
            }
        }
    };
    using repeats_t = std::vector<repeat_t>;
    auto repeats = repeats_t{};

    auto seeds = // non-const because will pop-back from this.
        fn::cfrom(counts)
      % fn::transform L(std::make_pair(_.first, get_count(_.first))) // key-and-count
      % fn::where L(_.second >= min_seed_count)
      % fn::sort_by L(_.second); // sort by count, most-frequent last.

    std::cerr << "Seeds: " << seeds.size() << "\n";
    const auto t = timer{};

    // Assemble repeats: extend from every seed (starting from most-frequent) in both directions.
    for(; !seeds.empty(); seeds.pop_back()) {

        const auto seed = seeds.back().first;
        VERIFY(get_minword(seed) == seed);

        if (is_seen(seed)) {
            continue;
        } else {
            mark_seen(seed);
        }

        auto repeat = repeat_t{ seed, words_t(1, seed) , 0ul, {} };

        // Greedily extend repeat.kmers chain both ways.

        // To properly extend terminals in palindromic cases, we can't simply
        // rely on orientation-invariant seen_set, since a repeat can contain
        // a k-mer in both orientations. Instead, will keep repeat-local
        // seen-set to keep track whether we've seen a k-mer in specific orientation
        // in this putative repeat.
        for (auto seen_lcl = ankerl::unordered_dense::set<word_t>{}; true; ) {

            const word_t next_w = get_next(repeat.kmers.back());
            const word_t prev_w = revcomp(get_next(revcomp(repeat.kmers.front())));

            const auto next_count = is_seen(next_w) || seen_lcl.count(next_w) ? 0ul : get_count(next_w);
            const auto prev_count = is_seen(prev_w) || seen_lcl.count(prev_w) ? 0ul : get_count(prev_w);

            //const auto next_count_alt = next_count == 0 ? 0 : get_count(next ^ 1ul);
            //const auto prev_count_alt = prev_count == 0 ? 0 : get_count(prev ^ (1ul << 63));

            if (next_count >= prev_count && next_count >= min_ext_count) {
                repeat.kmers.push_back(next_w);
                seen_lcl.insert(next_w); // NB: not applying minword()

            } else if (prev_count >= next_count && prev_count >= min_ext_count) {
                repeat.kmers.push_front(prev_w);
                seen_lcl.insert(prev_w); // NB: not applying minword()

            } else {
                break;
            }
        }

        // TODO: insert all neighbors within 1bp
        for (const auto& w : repeat.kmers) {
            mark_seen(w);
        }

        if (repeat.len() >= min_len) {
            //std::cerr << repeat.seed << "\t" << repeat.len() << "\n";
            repeat.validate();
            repeat.weight = sum_by(repeat.kmers, L(get_count(_)));
            repeat.residue_counts.resize(repeat.len());
            repeats.push_back(std::move(repeat));
        }
    }

    // filter weak predictions - take top 90% by weight
    {
        const auto total_weight = sum_by(repeats, L(_.weight));
        const auto target_weight = size_t(double(total_weight) * 0.9);
        auto cumulative_weight = 0ul;

        std::cerr << "count before filtering: " << repeats.size() 
                  << "; total len:" << sum_by(repeats, L(_.len())) << "\n";

        repeats = std::move(repeats)
                % fn::sort_by L(_.weight)
                % fn::reverse() // best-first by weight
                % fn::take_while L((cumulative_weight += _.weight) < target_weight)
                % fn::to_vector();
    }

    std::cerr << "Extended seeds in " << float(t) << "s.\n";
    std::cerr << "count after filtering: " << repeats.size() 
              << "; total len:" << sum_by(repeats, L(_.len())) << "\n";

    // The putative repeats are in { AG, CT } 1-bit alphabet.
    // 
    // We need to convert it to { A, C, G T }.
    // 
    // For that, we'll index k-mers in the putative repeats,
    // and will do another pass over the input fasta.
    //
    // If we find an indexed repeat-specific 1-bit-coding k-mer, we increment
    // the residue counts at the corresponding positions in the corresponding
    // putative repeat-sequence.
    //
    // Then, when outputting the repeat-sequences, we will take
    // the most frequent residue at every position.

    struct ord_and_pos_t
    {
        size_t ord = 0ul;
        size_t pos = 0ul;
    };

    size_t unexpected_conflicts_count = 0;
    const auto index = [&] // kmer -> (repeat_ord, pos)
    {
        auto ret = ankerl::unordered_dense::map<word_t, ord_and_pos_t>{}; 
        for (size_t i = 0; i < repeats.size(); i++) {
            for (size_t pos = 0; pos < repeats[i].kmers.size(); pos++) {
                const auto kmer = repeats[i].kmers[pos];

                if (ret.count(kmer) == 0) {
                    ret[kmer] = ord_and_pos_t{ i, pos };
                } else {
                    std::cerr
                        << "Key conflict:" << std::bitset<64>(kmer) << ":" 
                        << ret.at(kmer).ord << "@" << ret.at(kmer).pos << " vs. " 
                        << i                << "@" << pos << "\n";
                    ++unexpected_conflicts_count;
                }
            };
        }
        return ret;
    }();

    VERIFY(unexpected_conflicts_count == 0);

    /////////////////////////////////////////////////////////////////////////
    // Second pass - collect the position-specific residue counts.
    
    //auto istr = ser::open_istream("chr23.fa.gz");
    auto istr = ser::open_istream(genome_fasta_path);
    for (auto seq : MakeFastaReader(*istr, k_fasta_chunk_stride, k_fasta_chunk_overlap)) {
        process_kmers(seq.seq, kmer_len, [&](size_t pos, kmer_bufs_t bufs) {

            const auto key    = bufs.onebit & kmer_mask;
            const auto key_rc = revcomp(key);

            const int8_t ori  = index.count(key)    ? 1 
                              : index.count(key_rc) ? -1 
                              :                       0;

            // TODO: process separately in one orientation only, and then compare the results.
            if (!ori) {
                return;
            }

            const auto& v = ori == 1 ? index.at(key) : index.at(key_rc);

            for (size_t i = 0; i < kmer_len; i++) {
                char na2 = seq.seq[pos + i];
                na2 = (char)std::toupper((unsigned char)na2);
                na2 = na2 == 'A' ? 0 
                    : na2 == 'C' ? 1
                    : na2 == 'G' ? 2
                    : na2 == 'T' ? 3
                    :              0;

                na2 = (ori == 1) ? na2 : 3 - na2; // complement;

                auto tgt_pos = (ori == 1) ? v.pos + i
                                          : v.pos + kmer_len - 1 - i;

                ++repeats[v.ord].residue_counts.at(tgt_pos).at(na2);
            }
        });
    }

    const char* alphabet = "ACGT";
    for (size_t i = 0; i < repeats.size(); i++) {
        const auto& r = repeats[i];

        auto seq = std::string{};
        for (const auto& cts : r.residue_counts) {
            const auto it = std::max_element(cts.begin(), cts.end());
            seq.push_back(*it == 0 ? 'N' : alphabet[it - cts.begin()]);
        }

        // Drop low-complexity suffixes and prefixes.
        for (bool truncated = true; truncated; ) {
            truncated = false;
            for (const auto s : {
                "AAAAAA", "TTTTTT", "CCCCCC", "GGGGGG",
                "AGAGAG", "CTCTCT", "GAGAGA", "TCTCTC",
                "GTGTGT", "ACACAC", "TGTGTG", "CACACA",
                "ATATAT", "TATATA" })
            {
                if (str::endswith(seq, s)) {
                    truncated = true;
                    seq.resize(seq.size() - 2);
                }

                if (str::startswith(seq, s)) {
                    truncated = true;
                    seq = seq.substr(2);
                }
            }
        }

        if (seq.size() >= min_len) {
            ostr << ">lcl|repeat." << std::setfill('0') << std::setw(4) << i + 1
                 << " weight=" << r.weight 
                 << "\n" << seq << "\n";
        }
    }
}
