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

    const auto num_pages = (sv.size() + page_size - 1) / page_size;

    for (static bool printed_once = false; !printed_once && num_pages > num_phys_pages; printed_once = true) {
        std::cerr << "\033[91m" // red
                  << R"(
    Warning: The host does not have enough physical memory for the gx-database.
    The execution is likely to be extremely slow due to disk thrashing, and suitable only for tiny genomes (e.g. bacteria).
    Virtual memory paging does not provide adequate performance.
    See https://github.com/ncbi/fcs/wiki/FCS-GX for details.

        )" << "\033[0m";
    }

    // mincore will set values to non-0 for corresponding resident pages.
    auto is_resident = std::vector<unsigned char>(num_pages, 0);

    auto get_pct_pages_in_core = [&]
    {
        const auto get_rusage = []
        {
            rusage r{};
            getrusage(RUSAGE_SELF, &r);
            return r;
        };

        // mincore call is slow (up to 15 seconds even if all pages are in-core),
        // so first do a quick estimate accessing a random sample of pages and 
        // checking if the process had any page-faults.        
        const auto ru_before = get_rusage();
        const auto rand_seed = static_cast<uint64_t>( time(NULL) );
        volatile size_t pos = 0; // prevent from optimizing-out
        for(size_t i = 0; i < 100; i++) {
            pos = uint64_hash(rand_seed ^ i ^ (uint64_t(sv[pos]) << 32) ) % sv.size();
        }

        if (get_rusage().ru_majflt == ru_before.ru_majflt) {
            return 100ul; // assume 100% in-core
        } else {
            VERIFY(0 == mincore(const_cast<char*>(sv.data()), sv.size(), is_resident.data()));
            return std::count_if(is_resident.begin(), is_resident.end(), L(_ != 0)) * 100 / num_pages;
        }
    };

    auto pct_pages_in_core = get_pct_pages_in_core();

    if (pct_pages_in_core == 100 && !force) {
        return;
    }

    for (static bool printed_once = false; !printed_once && pct_pages_in_core < 100; printed_once = true) {
        std::cerr << R"(
    GX requires the database to be entirely in RAM to avoid thrashing.
    Consider placing the database files in a non-swappable tmpfs or ramfs.
    See https://github.com/ncbi/fcs/wiki/FCS-GX for details.
    Will prefetch (vmtouch) the database pages to have the OS cache them in main memory.

        )";
    }

    const auto elapsed = timer{};
    auto last_pct_processed = 0UL;
    volatile char c = 0;  // volatile to prevent optimizing-out the sv access.

    for (const auto i : irange{ num_pages }) {
        c = (force || !is_resident.at(i)) && sv.at(i * page_size); // touching non-resident pages

        // The rate of this is about the same as vmtouch -t, which does essentially the same thing.
        // However, cat file_on_disk > /dev/null, which also warms the cache, is 50% faster. HOW??
        // NB: __builtin_prefetch, which prefetches from RAM into CPU-cache, is of no use here.

        // update progress message
        if (const auto pct_processed = i * 100 / num_pages; pct_processed != last_pct_processed) {
            std::cerr << "Prefetching " << filename << " " << pct_processed << "%...                         \r";
            last_pct_processed = pct_processed;
        }
    }
    (void)c;

    pct_pages_in_core = get_pct_pages_in_core();
    std::cerr << "\nPrefetched " << filename << " in " << float(elapsed) << "s; "
              << float(sv.size())/1e9f/float(elapsed) << " GB/s. "
              << "The file is " << pct_pages_in_core << "% in RAM.\n";
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
       if (!s.empty() && s.back() == '\n') {
           s.pop_back();
       }
       return s;
    }();

    VERIFY(header.back() == ']');
    header.back() = ',';

    header += " {\"git-rev\":\"";
    header += g_git_revision;
    header += "\", \"run-date\":\"";
    header += time_now_str;
    header += "\"";

    static const std::string s_extra = []
    {
        auto p = std::getenv("GX_METALINE_JSON_EXTRA");
        auto ret = !p ? "" : std::string(p);

        // If does not look like json-array, object, or string, wrap as string.
        if (!ret.empty() && ret.find_first_of("[{\"") != 0) {
            VERIFY(ret.find('\"') == std::string::npos);
            VERIFY(ret.front() != ' ' && ret.back() != ' ');
            ret = "\"" + ret + "\"";
        }
        return ret;
    }();

    if (!s_extra.empty()) {
        header += ", \"extra\":";
        header += s_extra;
    }

    header += "}]";
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
