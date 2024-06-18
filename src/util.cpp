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
#include <map>
#include <cctype>

#include "types.hpp"
#include "segment.hpp"
#include "serial_util.hpp"
#include "ext/json5.hpp"

// for rusage
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

// for mmap
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

using namespace gx;

namespace fn = rangeless::fn;

using fn::operators::operator%; // see fn.hpp
using fn::operators::operator%=;

/////////////////////////////////////////////////////////////////////////////

uint64_t gx::revcomp_bits(uint64_t w, uint8_t num_bits)
{
    static const uint64_t k1  = 0x5555555555555555UL;
    static const uint64_t k2  = 0x3333333333333333UL;
    static const uint64_t k4  = 0x0F0F0F0F0F0F0F0FUL;

    w = ((w >>  1) & k1)  | ((w & k1)  <<  1); // swap adjacent bits
    w = ((w >>  2) & k2)  | ((w & k2)  <<  2); // swap bit-pairs
    w = ((w >>  4) & k4)  | ((w & k4)  <<  4); // swap nibbles

#if defined(__GNUC__)
    w = __builtin_bswap64(w);
#elif defined(_MSC_VER)
    w = _byteswap_uint64(w);
#else
    static const uint64_t k8  = 0x00FF00FF00FF00FFUL;
    static const uint64_t k16 = 0x0000FFFF0000FFFFUL;
    w = (( w >>  8 ) & k8  ) | (( w & k8  ) <<  8 ); // swap single bytes
    w = (( w >> 16 ) & k16 ) | (( w & k16 ) << 16 ); // swap 2-byte frames
    w = (( w >> 32 )       ) | (( w       ) << 32 ); // swap 4-byte frames
#endif

    return ~w >> (64 - num_bits); // flip the bits and shift into lsbs.
}


// Reverse-complement 2-bit-encoded word ending in LSB
uint64_t gx::revcomp_twobits(uint64_t w, uint8_t word_len)
{
    ASSERT(word_len <= 32);

    // First, reverse the order of letters.
    static const uint64_t k2  = 0x3333333333333333UL;
    static const uint64_t k4  = 0x0F0F0F0F0F0F0F0FUL;
    w = (( w >>  2 ) & k2  ) | (( w & k2  ) <<  2 ); // swap bit-pairs
    w = (( w >>  4 ) & k4  ) | (( w & k4  ) <<  4 ); // swap nibbles

#if defined(__GNUC__)
    w = __builtin_bswap64(w);
#elif defined(_MSC_VER)
    w = _byteswap_uint64(w);
#else
    static const uint64_t k8  = 0x00FF00FF00FF00FFUL;
    static const uint64_t k16 = 0x0000FFFF0000FFFFUL;
    w = (( w >>  8 ) & k8  ) | (( w & k8  ) <<  8 ); // swap single bytes
    w = (( w >> 16 ) & k16 ) | (( w & k16 ) << 16 ); // swap 2-byte frames
    w = (( w >> 32 )       ) | (( w       ) << 32 ); // swap 4-byte frames
#endif
    // The word is in MSBs after reversal, so shift it back into LSBs,
    // and flip the bits to complement the 2-bit bases.
    return ~w >> (2*(32 - word_len));
}

static const bool test_revcomp_twobits = []
{
                 //  T G G C C C A A A A
    uint64_t x1  = 0b11101001010100000000;

                 //  T T T T G G G C C A
    uint64_t y1e = 0b11111111101010010100;
    uint64_t y1a = gx::revcomp_twobits(x1, 10);
    VERIFY(y1a == y1e);
    return true;
}();


uint64_t gx::drop_every_3rd_bit(uint64_t w)
{
    static const uint64_t k1 = 0b0011000011000011000011000011000011000011000011000011000011000011;
    static const uint64_t k2 = 0b1111000000001111000000001111000000001111000000001111000000001111;
    static const uint64_t k3 = 0b0000000011111111000000000000000011111111000000000000000011111111;
    static const uint64_t k4 = 0b1111111111111111000000000000000000000000000000001111111111111111;
    static const uint64_t k5 = 0b0000000000000000000000000000000011111111111111111111111111111111;

    w = (w & k1) | ((w & (k1 <<  3)) >>  1);
    w = (w & k2) | ((w & (k2 <<  6)) >>  2);
    w = (w & k3) | ((w & (k3 << 12)) >>  4);
    w = (w & k4) | ((w & (k4 << 24)) >>  8);
    w = (w & k5) | ((w & (k5 << 48)) >> 16);
    return w;
}

static const bool test_drop_every_3rd_bit = []
{
    uint64_t x1  = 0b1011010110100010111000110101010111101001001111000110100001010010;
    uint64_t y1e = 0b1'11'10'10'00'10'11'00'10'01'10'11'01'01'01'11'00'10'00'01'10'10;
    uint64_t y1a = drop_every_3rd_bit(x1);
    VERIFY(y1a == y1e);
    return true;
}();

/////////////////////////////////////////////////////////////////////////////

using page_residency_t = std::vector<unsigned char>;

static int get_pct_pages_in_core(const page_residency_t& is_resident)
{
    const auto num_in_core = std::count_if(is_resident.begin(), is_resident.end(), L(_ != 0));
    return int(num_in_core * 100 / is_resident.size());
};

// Wrapper around mincore(2) https://man7.org/linux/man-pages/man2/mincore.2.html
// that avoids "problematic" invocations that otherwise cause performance issues.
static page_residency_t get_page_residency(const std::string_view sv, page_residency_t ret = {})
{
    static const auto page_size = (size_t)sysconf(_SC_PAGESIZE);
    const auto num_pages = (sv.size() + page_size - 1) / page_size;

    ret.clear();
    ret.resize(num_pages);

    // mincore call can be slow (up to 15 seconds even if all pages are in-core),
    // so first do a quick estimate accessing a random sample of pages and 
    // checking if the process had any page-faults.        
    const auto num_pagefaults_orig = ser::get_pagefault_count();
    static const auto num_iterations = 100ul;
    
    {
        const auto rand_seed = static_cast<uint64_t>(std::time(nullptr));
        volatile size_t pos = 0; // prevent from optimizing-out
        for (size_t i = 0; i < num_iterations; i++) {
            pos = uint64_hash(rand_seed ^ i ^ (uint64_t(sv[pos]) << 32)) % sv.size();
        }
    }

    const auto num_pagefaults = ser::get_pagefault_count() - num_pagefaults_orig;

    //std::cerr << "Num-pagefaults: " << num_pagefaults << std::endl;

    if (num_pagefaults <= 1) {
        // Almost no page-faults - report as-if all pages are in-core.
        std::fill(ret.begin(), ret.end(), 1);

    } else if (float(num_pagefaults) / float(num_iterations) > 0.8f) {
        // More than 80% page-faults - report as-if 0% in-core.
        //
        // NB: there's a strange performance regression that emerged:
        // When prefetching a file from disk that is 0% in-core,
        // (e.g. as-if after vmtouch -e, or from cold storage)
        // the prefetching loop in prefetch_mmapped_pages is extremely slow
        // (the program appears stuck and requires CTRL-C).
        //
        // Re-running a second time, when the file is even 0.1% in core,
        // the rest of the file is prefetched at expected rate (~1GiB/s from VAST).
        //
        // Avoiding the mincore call below when the file is (nearly) 0% in-core
        // seems to fix the issue (this else-clause), although I don't know how -
        // it may be related to the following:
        //
        // https://lwn.net/Articles/778437/
        // https://lwn.net/Articles/776801/
        //
        // I think might be a manifestation of a low-level OS issue,
        // rather than a bug in this code.
        std::fill(ret.begin(), ret.end(), 0);

    } else if (0 != mincore(const_cast<char*>(sv.data()), sv.size(), ret.data())) {
        
        [[ maybe_unused ]] static const bool printed_once = []
        {
            std::cerr << "Note: mincore() failed.\n";
            return true;
        }();

        std::fill(ret.begin(), ret.end(), 0);

    } else if (get_pct_pages_in_core(ret) == 100) {
        // Mincore "succeeded", but "lied" for security reasons (see 778437 above).
        //
        // "Interestingly, in the cases where mincore() 
        // does not return actual page-cache residency information,
        // it reports all pages as being present."
        //
        // In this case report all pages as not-present in order to
        // do full-prefetech, as not to erroneously skip prefetching everything.
        std::fill(ret.begin(), ret.end(), 0);
    }

    return ret;
}

size_t ser::get_pagefault_count()
{
    rusage r{};
    getrusage(RUSAGE_SELF, &r);
    return r.ru_majflt;
}


/////////////////////////////////////////////////////////////////////////////
// if `force`, then touch pages even though they are resident, to also prevent minor-page-faults
void ser::prefetch_mmapped_pages(const std::string& filename, std::string_view sv, bool force)
{
    static const auto num_phys_pages = (size_t)sysconf(_SC_PHYS_PAGES);
    static const auto page_size      = (size_t)sysconf(_SC_PAGESIZE);
    VERIFY(page_size >= 512);

    // mincore requires the starting address to be a multiple of page-size
    while (sv.size() > 0 && (size_t)sv.data() % page_size != 0) {
        sv = sv.substr(1);
    }

    if (sv.empty()) {
        return;
    }

    auto is_resident = get_page_residency(sv);
    const auto num_pages = is_resident.size();
    const auto pct_pages_in_core = get_pct_pages_in_core(is_resident);

    for (static bool printed_once = false;
         !printed_once && num_pages > num_phys_pages;
         printed_once = true)
    {
        std::cerr << "\033[91m" // red
                  << R"(
    Warning: The host does not have enough physical memory for the gx-database.
    The execution is likely to be extremely slow due to disk thrashing, and suitable only for tiny genomes (e.g. bacteria).
    Virtual memory paging does not provide adequate performance.
    See https://github.com/ncbi/fcs/wiki/FCS-GX for details.

        )" << "\033[0m";
    }

    if (pct_pages_in_core == 100 && !force) {
        return;
    }

    for (static bool printed_once = false;
         !printed_once && pct_pages_in_core < 100;
         printed_once = true)
    {
        std::cerr << R"(
    GX requires the database to be entirely in RAM to avoid thrashing.
    Consider placing the database files in a non-swappable tmpfs or ramfs.
    See https://github.com/ncbi/fcs/wiki/FCS-GX for details.
    Will prefetch (vmtouch) the database pages to have the OS cache them in main memory.

        )";
    }


    std::cerr << "\n\n" << filename << " is " << pct_pages_in_core << "% in RAM.\n";

    // Prefetching-loop:
    const auto elapsed = timer{};
    auto last_pct_processed = 0UL;
    for (const auto i : irange{ num_pages }) {
        if (force || !is_resident.at(i)) {
            const volatile char c = sv.at(i * page_size); // touching non-resident pages
            (void)c;
        }

        // The rate of this is about the same as vmtouch -t, which does essentially the same thing.
        // However, cat file_on_disk > /dev/null, which also warms the cache, is 50% faster. HOW??
        // NB: __builtin_prefetch, which prefetches from RAM into CPU-cache, is of no use here.

        // Update progress-message.
        const auto pct_processed = i * 100 / num_pages;
        if (pct_processed != last_pct_processed) {
            std::cerr << "Prefetching " << filename << " " 
                      << pct_processed << "%...                         \r";
            last_pct_processed = pct_processed;
        }
    }

    // update residency info for reporting.
    is_resident = get_page_residency(sv, std::move(is_resident));

    std::cerr << "\nPrefetched " << filename << " in " << float(elapsed) << "s; "
              << float(sv.size())/1e9f/float(elapsed) << " GB/s. "
              << "The file is " << get_pct_pages_in_core(is_resident) << "% in RAM.\n";
}



tax_map_t gx::LoadTaxa(std::istream* istr)
{
    auto ret = tax_map_t{};
    namespace tsv = rangeless::tsv;

    ret[tax_id_t(0)] = taxon_t{ "NULL", "NULL", "NULL", "NULL", taxdiv_oid_t{} };

    // Columns:
    // #tax_id   species-name   common-name         BLAST-name  GX-div-name
    // ----------------------------------------------------------------------
    // 102107    Prunus mume    Japanese apricot    eudicots    plants
    if (!istr) {
        return ret;
    }

    ConsumeMetalineHeader(*istr, GX_TSV_HEADER__TAXA);

    auto gx_taxdivs = std::vector<std::string>{};
    for (const tsv::row_t& row : tsv::from(*istr)) {
        const auto tax_id = (tax_id_t)tsv::to_num(row[0]);
        VERIFY(row.size() == 5);
        auto taxon = taxon_t{ row[1],     // species
                              row[2],     // common
                              row[3],     // blastgroup-name
                              row[4],     // div
                              taxdiv_oid_t{} };

        if (taxon.gx_taxdiv == "synthetic") {
            taxon.gx_taxdiv = "synt:synthetic"; // temporary work-around to support older gxdbs. GP-34646
        }

        // GP-31064
        // Strip "our" suffixes that are sequence-specific, not taxon-specific.
        static const auto suffixes = std::vector<std::string>{{", mitochondrion", ", plastid", ", plasmid"}};
        for (const auto& suffix : suffixes) {
            taxon.species = str::replace_suffix(std::move(taxon.species), suffix, "");
        }

        ret[tax_id] = taxon;
        gx_taxdivs.push_back(taxon.gx_taxdiv);
        VERIFY(gx_taxdivs.back() != "NULL");
    }

    // Assign taxdiv-ids based on taxdiv-name
    {
        gx_taxdivs %= fn::unique_all();
        VERIFY(std::is_sorted(gx_taxdivs.begin(), gx_taxdivs.end()));

        static_assert(sizeof(taxdiv_oid_t) == 1);
        VERIFY(gx_taxdivs.size() <= 250);

        for (auto& kv : ret)
            if (+kv.first) // oid == 0 for tax_id == 0
        {
            const auto it = std::lower_bound(gx_taxdivs.begin(), gx_taxdivs.end(), kv.second.gx_taxdiv);
            VERIFY(it < gx_taxdivs.end());
            kv.second.taxdiv_oid = taxdiv_oid_t(it - gx_taxdivs.begin() + 1); // 0 is reserved for tax_id=0
        }
    }

    return ret;
}


std::string gx::ConsumeMetalineHeader(std::istream& istr, std::string header)
{
    if (istr.peek() != '#') {
        // For now don't make it a hard error to allow working easier with grepped (filtered) inputs.
        std::cerr << "Warning: missing header '" << header << "'\n";
        return "";
    }

    VERIFY(istr);

    VERIFY(str::endswith(header, "]]")); // e.g. "[[\"GX taxonomic divisions\",1,1]]"
    VERIFY(std::count(header.begin(), header.end(), ',') == 2ul);

    // truncate the header to contain only the format name and the major-version,
    // e.g. "[[\"GX taxonomic divisions\",1,
    {
         const auto pos = header.find_last_of(',');
         VERIFY(pos + 3 < header.size());
         header.resize(pos + 1);
    }

    std::string line;
    std::getline(istr, line);
    header = str::replace(header, ", ", ",");
    line   = str::replace(line, ", ", ",");

    if (!str::startswith(line, header)) {
        GX_THROW("Expected the first line of the input file to begin with header: \n" + header + "\nfound: \n" + line + "\n");
    }
    VERIFY(istr);
    return line;
}

extern const char* const g_git_revision; // defined in main.cpp

std::string gx::MakeMetaLine(std::string header)
{
    VERIFY(str::startswith(header, "##"));

    const std::string time_now_str = [&]
    {
       auto const now = std::chrono::system_clock::to_time_t(
               std::chrono::system_clock::now());
       std::string s{ std::ctime(&now) };
       while (!s.empty() && s.back() == '\n') {
           s.pop_back();
       }
       return s;
    }();

    static const std::string s_extra = []
    {
        auto p = std::getenv("GX_METALINE_JSON_EXTRA");
        auto ret = !p ? "" : std::string(p);
        return ret;
    }();

    auto info = json5::value_t{};
    info["git-rev"] = g_git_revision;
    info["run-date"] = time_now_str;
    if (!s_extra.empty()) {
        info["extra"] = json5::parse(s_extra);
    }

    auto j = json5::parse(header.substr(2));
    j.get<json5::array_t>().push_back(info);
    header = "##" + json5::to_string(j);
    return header;
}


void ser::validate_eof(const std::string& filename, uint64_t magic_constant_eof)
{
    const size_t size = std::ifstream(filename, std::ifstream::ate | std::ifstream::binary).tellg();
    VERIFY(size > 8);

    std::ifstream ifstr{ filename };
    VERIFY(ifstr);
    ifstr.seekg(size - sizeof(uint64_t));

    const uint64_t actual_magic_constant = ser::from_stream(ifstr);
    if (actual_magic_constant != magic_constant_eof) {
        GX_THROW("File " + filename + " corrputed - missing expected file footer.");
    }
}


std::string_view ser::mmap(const std::string& path)
{
    const size_t size        = std::ifstream(path, std::ifstream::ate | std::ifstream::binary).tellg();
    const int    fd          = open(path.c_str(), O_RDONLY);
    const auto   close_guard = fn::make_scope_guard([&]{ close(fd); });

    if (fd <= 0) {
        GX_THROW("Could not open " + path);
    }

    // https://man7.org/linux/man-pages/man2/munmap.2.html
    // Additional flags of interest:
    // MAP_POPULATE // prevents TLB misses, but takes 40 seconds for RAM-cached db
    // MAP_LOCKED   // takes ~100 seconds for RAM-cached db
    //
    // MAP_HUGETLB  // only works in conjunction with MAP_ANONYMOUS;
    //              // must be backed by hugetlbfs (?)
    //              // requires privileged access (?)
    //              // https://www.kernel.org/doc/Documentation/vm/hugetlbpage.txt
    // MAP_HUGE_1GB // might use in conjuction with MAP_HUGETLB
    errno = 0;
    void* ptr = ::mmap(nullptr/*starting-address*/, size, PROT_READ, MAP_SHARED, fd, 0/*offset*/);

    if (ptr == MAP_FAILED) {
        GX_THROW("mmap failed; errno=" + std::to_string(errno));
    }

    return std::string_view{ (const char*)ptr, size };
}


std::unique_ptr<std::istream> ser::open_istream(std::string path)
{
    VERIFY(path != "-");
    VERIFY(path != "");

    const bool is_manifest_file = gx::str::endswith(path, ".mft");
    const auto file_ext = // excluding .mft
        get_file_extension(path.substr(0, path.size() - (is_manifest_file ? 4 : 0)));

    const auto orig_errno = errno;
    errno = 0;

    auto ifstr = std::make_unique<std::ifstream>(path, std::ifstream::in);

    if (errno || !*ifstr) {
        GX_THROW("Failed to open file: " + path + " - " + strerror(errno));
    }

    errno = orig_errno;

    const std::string decompressor_cmd = 
          file_ext == ".gz" && system("command -v minigzip >/dev/null") == 0
                             ? "minigzip -c -d"
        : file_ext == ".gz"  ?     "gzip -c -d"
        : file_ext == ".zstd"?     "zstd -c -d"
        : file_ext == ".lz4" ?      "lz4 -c -d"
        : file_ext == ".bz2" ?    "bzip2 -c -d"
        : file_ext == ".xz"  ?       "xz -c -d"
        : file_ext == ".lzma"?       "xz -c -d"
        :                                    "";

    static const bool enable_pv = 
           get_env("GX_ENABLE_PV", false) && (system("command -v pv >/dev/null") == 0);
    static const auto pv_cmd = std::string{ enable_pv ? " | pv -Wbrat " : "" };

    // NB: do not use zcat, bzcat, xzcat, zstdcat, lz4cat, because
    // in the context of xargs when processing a manifest,
    // when one fails, xargs continues to plow through the rest
    // before finally erroring-out the entire pipeline, so instead
    // we `xargs cat` the files, and pipe the output into the decompressor.

    // single-enquote path for passing to shell.
    path = "'" + str::replace(std::move(path), "'", "'\\''") + "'";

    if (is_manifest_file) {
        // Using awk to filter out empty and #-lines from manifest
        // because grep returns non-zero retcode if it doesn't find matching lines.
        return std::make_unique<pipe_istream>(
            "set -eo pipefail; cat " + path 
          + " | awk '!/^($|#)/' | xargs -n1 cat "
          + (decompressor_cmd.empty() ? "" : " | " + decompressor_cmd + pv_cmd)
        );
    } else if (decompressor_cmd != "") { // single compressed file
        return std::make_unique<pipe_istream>(decompressor_cmd + " < " + path + pv_cmd);
    } else {
        return ifstr;
    }
}
