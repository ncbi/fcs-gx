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


// Capitalization conventions.
// Defaults: std/boost.
// Macros: ALL_CAPS.
// Classes that "do things": PascalCase, prefixed with C, e.g CMyClass.
// Fields in above classes prefixed with "m_"; static fields prefixed with "s_".
// Structs that represent product-types (tuple-like): snake_case, suffixed with "_t", e.g. my_type_t.
// typedefs: same as above: e.g. using things_t = std::vector<thing_t>;
// Free-functions in gx:: namespace having external linkage: PascalCase.
// internal-linkage free-functions: snake_case.

#include <iostream>
#include <fstream>
#include <chrono>
#include <utility>

#define RANGELESS_ENABLE_TSV 1
#include <ext/fn.hpp>

// In order to facilitate versioning and the exchange of metadata via TSV files
// between a producer and a consumer, the first line in a TSV-file shall begin
// with a single line of metadata, prefixed with '##'.
//
// The second line shall contain column-names, prefixed with '#'.
//
// In general, the content of meta-line is unspecified, e.g. could be just an
// unstructured comment, a magic-value, or could contain structured information
// in a secret encoding between the producer and the consumer.
//
// For GX tabular inputs and outputs, The '##'-line shall contain a JSON-array,
// with the fist element specifying the format name and version as:
// ["$FORMAT_NAME",$FORMAT_MAJOR_VER,$FORMAT_MINOR_VER]. The array may contain
// additional 0 or more elements as pertains to the particular format.
// E.g.
//
// ##[["GX taxonomic divisions",1,1], { "comment" : "some additional metadata" }]
// #tax-id  species       common-name     BLAST-div       div
// 9606     Homo sapiens  human           primates        primates
// ...
//
// The consumer may, but does not have to verify the column names,
// as the format specifier shall be bound to a specific number of columns, their types, and semantic meaning.
//
// E.g. if the producer code is changing the meaning of column-2, chromosome-position,
// to mean 1-based instead of 0-based, they will correspondingly update the major-version,
// indicating not-backward-compatible change, whereas simply changing column name
// from "chr" to "Chromosome" or adding an extra-column without changing the order or
// meaning of other columns will only require a minor-version update.
//
#define GX_TSV_HEADER__TAXA             "##[[\"GX taxonomic divisions\",1,1]]"           // input for gx-make-db and gx-taxify
#define GX_TSV_HEADER__SEQ_ID_MAPPING   "##[[\"GX seq-id to tax-id mapping\",1,1]]"      // input for gx-make-db
#define GX_TSV_HEADER__HITS             "##[[\"GX hits\",2,1]]"                          // output of gx-align align, input of gx-taxify
// ver-2: _transposon and _low_complexity reporting

#define GX_TSV_HEADER__PRE_TAXONOMY_RPT "##[[\"GX taxonomy pre-analysis report\",3,1]]"  // output of gx-taxify, input to classify_taxonomy.py
#define GX_TSV_HEADER__FCS_GENOME_RPT   "##[[\"FCS genome report\",2,1]]"
// ver-2: divs with tax-kingdom prefixes
// ver-3: col-3 is comma-delimited, containig low-complexity and transposon lengths

#define GX_TSV_HEADER__LOCS             "##[[\"GX locs\",1,1]]"                          // 3-column TSV of seq-id, start1, stop1(inclusive)
                                                                                         // start1==stop1==0 represents whole-loc.

static inline std::string make_source_location_str(const char* file, const char* function, size_t line)
{
    auto ret = std::string{ file };
    if (const auto pos = ret.rfind('/'); pos != std::string::npos) {
        ret = ret.substr(pos+1); // keep just the filename from the full path.
    }
    return std::move(ret) + ":" + std::to_string(line) + " in " + function + "(...)";
}
#define GX_SOURCE_LOCATION_STR make_source_location_str(__FILE__, __func__, __LINE__)

#define GX_THROW(msg) throw std::runtime_error(GX_SOURCE_LOCATION_STR + ": " + (msg))

#define VERIFY(x) do { if (__builtin_expect(!(x), 0)) GX_THROW("Assertion failed: "#x); } while(0);
#define ASSERT(x) assert((x))

namespace fn = rangeless::fn;

// Abbreviated unary lambda, e.g. std::find_if(cats.begin(), cats.end(), L(_.name=="Spot"));
#define L(expr) ([&](auto&& _){ return expr; })

// e.g. std::max_element(ivls.begin(), ivls.end(), BY(_.len));
#define BY fn::by::make_comp L


// returns std::string_view into a byte-range of an embedded blob NAME linked with
// -Wl,--format=binary -Wl,NAME -Wl,--format=default
#define GET_EMBEDDED_BLOB(NAME)                       \
[]{                                                   \
    extern const char _binary_ ## NAME ## _start[];   \
    extern const char _binary_ ## NAME ## _end[];     \
    const char* b =   _binary_ ## NAME ## _start;     \
    const char* e =   _binary_ ## NAME ## _end;       \
    return std::string_view(b, size_t(e-b));          \
}()

namespace gx
{

struct str
{
    static inline bool startswith(const std::string& s, const std::string& prefix)
    {
        return s.length() >= prefix.length()
            && 0 == s.compare(0, prefix.length(), prefix);
    }

    static inline bool endswith(const std::string& s, const std::string& suffix)
    {
        return s.length() >= suffix.length()
            && 0 == s.compare(s.length() - suffix.length(), suffix.length(), suffix);
    }

    static inline bool contains(const std::string& s, const std::string& substr)
    {
        return s.find(substr) != std::string::npos;
    }

    [[ nodiscard ]]
    static inline std::string replace_suffix(
        std::string s,
        const std::string& old_sfx,
        const std::string& new_sfx)
    {
        if (str::endswith(s, old_sfx)) {
            s.resize(s.size() - old_sfx.size());
            s += new_sfx;
        }
        return s;
    }

    [[ nodiscard ]]
    static inline std::string replace(
        std::string s, // taking by value, returning via NRVO
        const std::string& substr,
        const std::string& new_substr,
        size_t count = size_t(-1))
    {
        VERIFY(!substr.empty());

        size_t n = 0;

        for (auto pos = s.find(substr);
             pos != std::string::npos && n < count;
             pos = s.find(substr, pos + new_substr.size()))
        {
            s.replace(pos, substr.length(), new_substr);
            ++n;
        }
        return s;
    }
};

/////////////////////////////////////////////////////////////////////////////
// e.g. static const auto my_var = get_env("MY_VAR", 42);
template<typename T> // T is either numeric or a std::string
T get_env(const char* name, T default_value)
{
    const char* p = std::getenv(name);
    if (!p) {
        return default_value;
    }
    const bool is_empty = *p == 0;

    auto ret = T{};
    if constexpr (std::is_arithmetic<T>::value) {
        VERIFY(!errno);
        ret = is_empty ? default_value : (T)rangeless::tsv::to_num(p);
    } else {
        static_assert(std::is_same<T, std::string>::value);
        ret = is_empty ? std::move(default_value) : p;
    }

    std::cerr << "Using " << name << "=" << ret << (is_empty ? " (default)" : "") << std::endl;
    return ret;
}

/////////////////////////////////////////////////////////////////////////////
/// \brief A simple timer.
struct timer
{
    /// returns the time elapsed since timer's instantiation, in seconds.
    /// To reset: `my_timer = timer_t{}`.
    operator float() const
    {
        return (float)std::chrono::duration_cast<std::chrono::nanoseconds>(
            clock_t::now() - start_timepoint).count() * 1e-9f;
    }
private:
    using clock_t = std::chrono::steady_clock;
    clock_t::time_point start_timepoint = clock_t::now();
};


/////////////////////////////////////////////////////////////////////////////

constexpr uint64_t Ob1x(uint64_t n)
{
    return ~(~uint64_t(0) << n);
}
static_assert(Ob1x(5) == 0b11111, "");

/////////////////////////////////////////////////////////////////////////////

static inline uint64_t uint64_hash(uint64_t x)
{
    // Warning: replacing the implementation of this function with
    // something else will silently break compatibility with existing
    // serialized indexes containing bloom-filters created with
    // this function, so don't.
    x += 0x9E3779B97F4A7C15ULL; // 2^64 / phi
    x ^= x >> 30; // note: this op is bijective, even with shift amt < 32
    x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27;
    x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;

    return x;
}

// fold the hash function above over some POD array
template<typename T>
static uint64_t uint64_hash(const T* arr, size_t n)
{
    uint64_t h = 0;
    const char* beg = reinterpret_cast<const char*>(arr);
    const char* end = reinterpret_cast<const char*>(arr + n);

    while (beg + 8 < end) {
        uint64_t x = *reinterpret_cast<const uint64_t*>(beg);
        h = uint64_hash(h ^ x);
        beg += 8;
    }
    uint64_t x = 0;
    memcpy(&x, beg, size_t(end - beg));
    return uint64_hash(h ^ x);
}


/////////////////////////////////////////////////////////////////////////////
// for(const auto i : irange{0, 42ul}) { // or irange{ 42ul }
//       std::cerr << i << "\n";
// }
// NB: type of i is that of end-value (i.e. unsigned long in this case).
template<typename T>
struct irange
{
    struct it_t
    {
        T i;
        it_t& operator++() { ++i; return *this; }
        bool operator!=(const it_t& other) const { return i < other.i; }
        const T& operator*() const { return i; }
    };
    it_t m_it;
    it_t m_end;

    explicit irange(T e) : m_it{0}, m_end{e} {}

#if 0
    // NB: type deduction will fail when b and e are of different types,
    // e.g. irange{ 1, vec.size() }
    explicit irange(T b, T e) : m_it{ b }, m_end{ e }
# else
    template<typename Tb>
    explicit irange(Tb b, T e) : m_it{ (T)b }, m_end{ e }
    {
        static_assert(std::is_same<Tb, T>::value || sizeof(Tb) < sizeof(T), "");
        if constexpr (std::is_signed_v<Tb> && std::is_unsigned_v<T>) {
            VERIFY(b >= 0);
        }
    }
#endif

    it_t begin() const { return m_it; }
    it_t end() const { return m_end; }
};


template<typename Container>
auto make_vec_of_iterators(Container&& cont)
{
    auto ret = std::vector<decltype(cont.begin())>{};
    ret.reserve(cont.size());
    for (auto it = cont.begin(); it != cont.end(); ++it) {
        ret.push_back(it);
    }
    return ret;
}

/////////////////////////////////////////////////////////////////////////////
// reverse-complement lower num_bits.
uint64_t revcomp_bits(uint64_t buf, uint8_t num_bits);

// drop every third bit starting with 3'rd lsb
uint64_t drop_every_3rd_bit(uint64_t buf);

// for each transitively adjacent-equal-by-get-key group of elements, invoke fn on the group.
template<typename Iterable, typename KeyFn, typename DoFn>
void for_each_group_by(
              Iterable& iterable,
                  KeyFn get_key,
                   DoFn fn)
{
    auto get_range = [&](auto beg)
    {
        if (beg == iterable.end()) {
            return fn::from(beg, iterable.end());
        }

        const auto key = get_key(*beg); // By const auto& here instead?

        auto it_end = std::find_if(std::next(beg), iterable.end(), [&](const auto& x)
        {
            return !(key == get_key(x));
        });

        return fn::from(beg, it_end);
    };

    for (auto r = get_range(iterable.begin());
              r.begin() != r.end();
              r = get_range(r.end()))
    {
        fn(r);
    }
}

// similar to above, except using binary CompareEqual evaluated on pairs of adjacent elements
template<typename Iterable, typename CompareEqual, typename DoFn>
void for_each_group_if(
              Iterable& iterable,
           CompareEqual eq,
                   DoFn fn)
{
    auto get_range = [&](auto beg)
    {
        if (beg == iterable.end()) {
            return fn::from(beg, iterable.end());
        }

        // adjacent_find returns iterator ito the first element of the adjacent pair that
        // compares equal, while we are looking the end of the run of equal elements, i.e.
        // the first pair of elements that does not compare equal.
        auto it_end = std::adjacent_find(beg, iterable.end(), [&](const auto& a, const auto& b)
        {
            return !eq(a, b);
        });

        if (it_end != iterable.end()) {
            ++it_end;
        }

        return fn::from(beg, it_end);
    };

    for (auto r = get_range(iterable.begin());
              r.begin() != r.end();
              r = get_range(r.end()))
    {
        fn(r);
    }
}

/////////////////////////////////////////////////////////////////////////////
template<typename Map>
auto at_or_default(const Map& m, const typename Map::key_type& key) -> const typename Map::mapped_type&
{
    static const typename Map::mapped_type default_value{};
    const auto it = m.find(key);
    return it == m.end() ? default_value : it->second;
}


/////////////////////////////////////////////////////////////////////////////
// e.g. sum_by(segs, L(_.len))
template<typename Xs, typename F>
auto sum_by(const Xs& xs, const F& fn) -> decltype(fn(*xs.begin()))
{
    decltype(fn(*xs.begin())) ret = 0;
    for (const auto& x : xs) {
        ret += fn(x);
    }
    return ret;
}


/////////////////////////////////////////////////////////////////////////////
/// For executing code on scope-exit (via exception or normally).
/// See https://en.cppreference.com/w/cpp/error/uncaught_exception
template<typename UnaryInvokable>
auto make_scope_guard(UnaryInvokable fn/*fn(bool is_exception){ ... }*/)
{
    struct guard_t
    {
        UnaryInvokable fn;
        const int orig_uncaught_ex_count;
        // If the user code needs to be able to dismiss the guard,
        // they can conditionally do it in fn's body.

        ~guard_t() noexcept(false) // allow fn to throw without aborting.
        {
            const bool is_exception = orig_uncaught_ex_count < std::uncaught_exceptions();
            fn(is_exception);
        }
    };
    return guard_t{ std::move(fn), std::uncaught_exceptions() };
}

/// For executing code on scope-exit during stack-unwinding.
/// See https://en.cppreference.com/w/cpp/error/uncaught_exception
template<typename NullaryInvokable>
auto make_exception_scope_guard(NullaryInvokable fn)
{
    return make_scope_guard([fn = std::move(fn)](bool is_exception)
    {
        if (is_exception) {
            fn();
        }
    });
}

/////////////////////////////////////////////////////////////////////////////

static inline std::ifstream open_ifstream(const std::string& path)
{
    VERIFY(!errno);
    auto ret = std::ifstream{ path };
    if (!ret) {
        errno = 0;
        GX_THROW("Could not open file: " + path);
    }
    return ret;
}

static inline std::unique_ptr<std::ifstream> open_ifstream_opt(const std::string& path)
{
    VERIFY(path != "-");
    return path.empty() ? std::unique_ptr<std::ifstream>{}
                        : std::make_unique<std::ifstream>(open_ifstream(path));
}

/////////////////////////////////////////////////////////////////////////////


#if RANGELESS_FN_ENABLE_PARALLEL
template<typename Container, typename F>
auto for_each_in_parallel( Container& xs,
                             unsigned num_threads,
                                    F transform_in_place)
  -> decltype(transform_in_place(*xs.begin()), void())
{
    if (num_threads <= 1) {
        for (auto& x : xs) {
            transform_in_place(x);
        }
        return;
    }

    namespace fn = rangeless::fn;
    using fn::operators::operator%; // see fn.hpp

    fn::seq([&, i = 0UL]() mutable
    {
        return i < num_threads ? i++ : fn::end_seq();
    })

  % fn::transform_in_parallel([&](size_t k)
    {
        // each thread will process elements where i % num_threads == k, where 0<k<num_threads
        size_t i = 0;

        for (auto& x : xs) {
            if (i % num_threads == k) {
                transform_in_place(x);
            }
            ++i;
        }
        return 0;
    }).queue_capacity(num_threads)

  % fn::for_each([](int){});
}

#endif

// Process GX-hits and create .taxonomy.rpt output.
void Taxify( std::istream& istr,
             std::istream& taxa_istr,
        const std::string& hardmask_locs_path,
        const std::string& asserted_div,
        const std::string& db_path,
        const std::string& out_path);

// Create gxdb.
void MakeDb( std::istream& fasta_istr,
             std::istream& seq_id2tax_id_istr,
             std::istream* taxa_istr,     // can be nullptr
             std::istream* hardmask_istr, // can be nullptr
             std::istream* softmask_istr, // can be nullptr
        const std::string& out_path);     // /path/to/out_db.gxi

// Align queries against gxdb
void ProcessQueries(const std::string& db_path,
                    const std::string& taxa_path,               // can be empty
                    const std::string& fasta_for_repeats_path,  // can be empty
                         std::istream& fasta_istr,
                         std::ostream& ostr);

// Retrieve sequences from gxdb as fasta.
void GetFasta(const std::string& db_path, std::istream& istr, std::ostream& ostr);

// Apply edits from action-report to the input fasta
void ApplyActionReport( std::istream& fasta_istr, 
                        std::istream& action_report_istr,
                        std::ostream& ostr, 
                        std::ostream* contam_fasta_out_ofstr_ptr,
                        size_t min_seq_len,
                        bool silent = false);

// Split fasta on N-runs
void SplitFasta(std::istream& fasta_istr, std::ostream& ostr);

// Will verify that the stream starts with expected header (e.g. GX_TSV_HEADER__TAXA)
// (ignoring minor version) and consume it from the stream (throw otherwise).
std::string ConsumeMetalineHeader(std::istream& istr, std::string header);

// Construct a string for a specified header, e.g.
// ##[["GX hits",1,1], {"git-rev":"5d84581", "run-date":"Mon Sep 20 10:13:01 2021", "extra":"contents of $GX_METALINE_JSON_EXTRA"}]
std::string MakeMetaLine(std::string header);

// Returns queries chunked at every multiple of chunk_stride of length chunk_stride+overlap_size
struct fasta_seq_t;
enum class validate_iupacna_t : bool { no, yes };
fn::any_seq_t<fasta_seq_t> MakeFastaReader(
        std::istream& istr,
        size_t chunk_stride = 4000000000UL,
        size_t overlap_size = 100UL,
        validate_iupacna_t  = validate_iupacna_t::yes);

std::string MakeProtsetMinhash(std::istream& prot_fasta);

void PairwiseCompareMinHashes(std::istream& istr, std::ostream& ostr);

} // namespace gx
