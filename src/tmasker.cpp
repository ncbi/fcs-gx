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

#include "util.hpp"
#include "tmasker.hpp"
#include "segment.hpp"
using namespace gx;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;

/////////////////////////////////////////////////////////////////////////////

// NB: as function so we don't get noise on stderr from get_env if it isn't called.
static size_t s_get_num_cores()
{
    static const auto s_num_cores = get_env("GX_NUM_CORES", std::min(32U, std::thread::hardware_concurrency() - 1));
    return s_num_cores;
}

static const size_t k_fasta_chunk_stride = 1000000UL;
static const size_t k_fasta_chunk_overlap = 50;

gx::CTmasker::CTmasker(std::istream& fasta_istr1)
{
    m_counts.resize(1<<30);
    const auto t = timer{};

    size_t total_len = 0;
    size_t short_seqs_len = 0;

    const auto fill_index_from_chunk = [&](fasta_seq_t inp_chunk)
    {
        const auto bits = bool_na_view_t{ inp_chunk.seq };
        auto buf = skip_mod3_buf_t::init_from(bits, CTmasker::k_word_tlen, 0);
        size_t last_N_pos = -1;

        for (const auto i : irange{ CTmasker::k_word_tlen, bits.size()}) {
            buf <<= bits[i];

            if (inp_chunk.seq[i] == 'N' || inp_chunk.seq[i] == 'n') {
                last_N_pos = i;
            } else if (i > last_N_pos + CTmasker::k_word_tlen) {
                this->insert(CTmasker::minword30_t{ buf });
            }
        }
        return inp_chunk;
    };

    std::cerr << "Collecting masking statistics..." << std::endl;

    MakeFastaReader(fasta_istr1, k_fasta_chunk_stride, k_fasta_chunk_overlap)

    // Ignoring alt-loci - if there are many of them for a particular locus,
    // it may make the counts for it look repeat-like.
    //
    // Ignore very short seqs - they are more likely to be garbage, and we
    // don't want to use them for the purpose of stats gathering.
    // There could also be a very large number of them (millions, e.g. GCA_900067735.1),
    // so we need some kind of minimum length threshold.
    //
    // Accumulate short_seqs_len and total_len while we're at it,
    // which we'll use downstream for determining N80 threshold.
  % fn::where([&](const fasta_seq_t& seq)
    {
        static const auto short_seq_len_thr = get_env("GX_TMASKER_SHORT_SEQ_LEN", 1000ul);
        const bool is_short = seq.seq.size() < short_seq_len_thr;
        short_seqs_len += seq.seq.size() * is_short;
        total_len      += seq.seq.size();

        return !is_short && !str::contains(seq.defline, "alternate locus");
    })

    // Fill in the index from the longer queries that are unlikely contaminant-repeats.
    // Pass the shorter ones downstream.
  % fn::transform_in_parallel([&](fasta_seq_t inp)
    {
        static const auto long_seq_len_thr = get_env("GX_TMASKER_LONG_SEQ_LEN", 100000ul);
        if (inp.seq.size() >= long_seq_len_thr) {
            fill_index_from_chunk(std::move(inp));
            inp.seq.clear(); // mark as used.
        }
        return inp;
    }).queue_capacity(s_get_num_cores())

  % fn::where L(!_.seq.empty()) // drop the ones that were used.

    // GP-31698, GP-32789.
    // From the remaining shorter seqs, take those least N80 in length,
    // because the abundance of short seqs could be a case of egregious
    // contamination that looks like transposons - we don't want to mask it out.
    //
    // E.g. phiX in GCA_900166895.1, or iridovirus in cricket (GP-32789)
  % fn::sort_by L(_.seq.size())

    // Drop shortest seqs while partial sum is below 20% of total_len
    //
    // NB: at this point total_len is done, because sort_by
    // forces eager-evaluation, so all stages above are done.
  % fn::drop_while L((short_seqs_len += _.seq.size()) * 5 < total_len)

    // Fill-in the index from the remainder of short queries longer than N80
  % fn::transform_in_parallel(fill_index_from_chunk).queue_capacity(s_get_num_cores())

  % fn::for_each L(void(_)); // iterate over lazy-seq to actually do the work.

    /////////////////////////////////////////////////////////////////////////
    // compute and set m_baseline
    {
        auto histogram = std::array<size_t, 256>{};
        for (const auto& x : m_counts) {
            ++histogram[x];
        }

        size_t sum1 = 0, sum2 = 0;
        for (const size_t i : irange{ 1, 11 }) {
            sum1 += histogram[i] * i;
            sum2 += histogram[i];
        }

        m_baseline = float(sum1)/float(sum2+1);
    }

    std::cerr << "Collected masking stats:  "
              << float(total_len)/1e9f         << " Gbp; "
              << float(t)                      << "s; "
              << float(total_len)/1e6/float(t) << " Mbp/s. "
              << "Baseline: " << m_baseline << "\n\n";
}

gx::CTmasker::CTmasker(const std::string& fasta_path)
{
    if (fasta_path.empty()) {
        return; // construct empty tmasker if path is empty
    }
    auto istream = open_ifstream(fasta_path);
    *this = CTmasker{ istream };
}

gx::ivls_t gx::CTmasker::find_repeats(const iupacna_seq_t& seq) const
{
    if (!this->collected_stats()) {
        return ivls_t{};
    }

    /*
        When making changes to the parameters or algorithm, evaluate results by
        running megablast, e.g.

        echo -e "NC_000023\nNC_000024" | getfasta | zstd -c > chr23_24.fa.zst
        zstdcat chr23_24.fa.zst | gx find-repeats --repeats-basis-fa=<(zcat data/GCF_000001405.39_GRCh38.p13_genomic.fna.gz| pv -brat) \
                                                  --outfmt=lcase | zstd -c > chr23_24.mfa.zst
        zstdcat chr23_24.mfa.zst | head -n2 > chr23.mfa
        zstdcat chr23_24.mfa.zst | tail -n2 > chr24.mfa
        time blastn -task megablast -outfmt 7 -query chr23.mfa -subject chr24.mfa -lcase_masking > chr23_chr24.tm.mbhits
    */

    // minimum support - minimum count of occurences of a repeat-specific LSH-word
    static const float  min_repeat_sup_factor = get_env("GX_TMASKER_MIN_REPEAT_SUP_FACTOR", 30.0f);
    static const size_t min_repeat_sup_cutoff = get_env("GX_TMASKER_MIN_REPEAT_SUP_CUTOFF", 50ul);
    const size_t min_repeat_sup = std::min(min_repeat_sup_cutoff, size_t(m_baseline * min_repeat_sup_factor));
    VERIFY(min_repeat_sup <= 250); // because m_counts are capped

    static const auto min_repeat_len = get_env("GX_TMASKER_MIN_REPEAT_LEN", 100);
    VERIFY(min_repeat_len >= CTmasker::k_word_tlen);

    const auto bits = bool_na_view_t{ seq };
    auto buf = skip_mod3_buf_t::init_from(bits, CTmasker::k_word_tlen, 0);
    auto ivls = ivls_t{};

    for (const auto i : irange{ CTmasker::k_word_tlen, bits.size() }) {
        buf <<= bits[i];

        if (this->at(buf) < min_repeat_sup) {
            continue;
        }

        const auto ivl = ivl_t{ pos1_t(i + 2 - CTmasker::k_word_tlen), CTmasker::k_word_tlen };
        // +1 to convert stop-pos to end-pos;
        // +1 to make 1-based;
        // -k_word_tlen to convert to start-pos

        // merge if high-overlap (i.e. avoid merging suprious at this stage)
        if (!ivls.empty() && ivls.back().endpos() + 3 >= ivl.endpos()) {
            ivls.back().len = ivl.endpos() - ivls.back().pos;
        } else {
            ivls.push_back(ivl);
        }
    }

    // Drop suprious.
    ivls %= fn::where L(_.len > CTmasker::k_word_tlen + 10);

    // Dilate by few bp on both ends.
    static const int k = 5;
    for (auto& ivl : ivls) {
        if (ivl.pos > k) {
            ivl.pos -= k;
            ivl.len += k;
        }

        if (ivl.endpos()-1 + k <= (int64_t)seq.size()) {
            ivl.len += k;
        }
    }

    // Merge neighbors within min_repeat_len distance.
    for (const auto i : irange{ 1, ivls.size()})
        if (ivls[i-1].endpos() + min_repeat_len > ivls[i].pos)
    {
        ivls[i] = ivl_t::collapse(ivls[i-1], ivls[i]);
        ivls[i-1] = ivl_t{};
    }

    ivls %= fn::where L(_.len > min_repeat_len);

    return ivls;
}


void gx::CTmasker::process_fasta(std::istream& fasta_istr, std::ostream& ostr) const
{
    enum class outfmt_t { unset = 0, locs, lcase_fa, hardmask_fa};
    static const auto outfmt = (outfmt_t)get_env("GX_REPEATS_OUTFMT", 1);
    VERIFY(1 <= (int)outfmt && (int)outfmt <= 3);

    if (outfmt == outfmt_t::locs) {
         std::string header = GX_TSV_HEADER__LOCS;
         VERIFY(header.back() == ']');
         header.pop_back();
         header += ", \"repeat-regions\"]";
         ostr << header << "\n#seq-id\\tfrom1\\tto1\n";
    }

    auto t = timer{};

    MakeFastaReader(fasta_istr,
                    outfmt == outfmt_t::locs ? k_fasta_chunk_stride : 10000000000,
                    k_fasta_chunk_overlap)

  % fn::transform_in_parallel([&](fasta_seq_t inp_chunk)
    {
        auto ivls = this->find_repeats(inp_chunk.seq);

        // not expecting chunking in lcase or hardmasking modes
        VERIFY(outfmt == outfmt_t::locs || inp_chunk.offset == 0);

        if (outfmt == outfmt_t::hardmask_fa) {
            for (const auto& ivl : ivls)
                for (const auto pos : irange{ ivl.pos, ivl.endpos() })
            {
                inp_chunk.seq[pos - 1] = 'N';
            }

        } else if (outfmt == outfmt_t::lcase_fa) {
            // uppercase fasta
            for (auto& c : inp_chunk.seq) {
                c = (c >= 'a' && c <= 'z') ? char(c + ('A' - 'a')) : c;
            }

            // lowercase within intervals
            for (const auto& ivl : ivls)
                for (const auto pos : irange{ ivl.pos, ivl.endpos() })
            {
                char& c = inp_chunk.seq[pos - 1];
                c = char(c + ('a' - 'A'));
            }
        }

        return std::make_pair(std::move(inp_chunk), std::move(ivls));
    }).queue_capacity(s_get_num_cores())

  % fn::for_each([&](std::pair<fasta_seq_t, ivls_t> result)
    {
        auto& inp_chunk = result.first;
        const auto& ivls = result.second;

        if (outfmt != outfmt_t::locs) { // output is lcase or hardmasked fasta
            VERIFY(!inp_chunk.defline.empty() && inp_chunk.defline.front() == '>');
            ostr << inp_chunk.defline << "\n" << inp_chunk.seq << std::endl;
            return;
        }

        // TODO: merge intervals in input-chunk-overlaps

        for (const auto& ivl : ivls) {
            ostr <<         inp_chunk.seq_id
                 << "\t" << inp_chunk.offset + ivl.pos
                 << "\t" << inp_chunk.offset + ivl.endpos() - 1
                 << "\n";
        }

        // If empty, output a zero-length sentinel interval just to report the seq-id.
        if (ivls.empty()) {
            ostr <<         inp_chunk.seq_id
                 << "\t" << 0
                 << "\t" << -1
                 << "\n";
        }
    });

    std::cerr << "Processed repeats in " << float(t) << "s.\n";
}
