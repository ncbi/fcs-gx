// See PUBLIC DOMAIN NOTICE at the bottom.
#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <variant>
#include <sstream>
#include <charconv>
#include <cmath>

namespace gx
{

/// JSON5-parser / JSON-writer https://json5.org/
class json5
{
public:
    struct value_t;

    // Object as vector-of-(key, val) instead of map:key->val to control and preserve the order of keys.
    using object_t  = std::vector<std::pair<std::string, value_t>>; // {"foo": 42, "bar": [1, 2, 3]}
    using array_t   = std::vector<value_t>;                         // [null, false, 0, 0.0, "", [], {}]
    using int_t     = int64_t;                                      // /^[+-]?(0x)?\d+$/
    using real_t    = long double;  // NB: write with ".0" suffix as necessary to preserve floatness.
    using variant_t = std::variant<std::nullptr_t, bool, int_t, real_t, std::string, array_t, object_t>;

    struct value_t : public variant_t
    {
        using variant_t::variant;
        using variant_t::operator=;

        /// e.g `if (j.is<json5::real_t>()) ...`
        template<typename T> auto is() const -> bool { return std::holds_alternative<T>(*this); }

        /// e.g. `j.get<json5::real_t>() > 42.0`
        template<typename T> auto get() const -> const T& { return std::get<T>(*this); }
        template<typename T> auto get()       ->       T& { return std::get<T>(*this); }

        /// access element of a json-array; throw if not an array or out of bounds
        auto at(size_t i) const -> const value_t&
        {
            return std::get<array_t>(*this).at(i);
        }

        /// access element of a json-object; throw if not an object, or key is missing
        auto at(std::string_view key) const -> const value_t&
        {
            for (auto& kv : std::get<object_t>(*this))
                if (kv.first == key)
            {
                return kv.second;
            }
            throw std::runtime_error("No such key in the object: " + std::string(key));
        }

        /// access element of a json-array (implicitly convert to array; append if i == size())
        auto operator[](size_t i) -> value_t&
        {
            if (!std::holds_alternative<array_t>(*this)) {
                *this = array_t{};
            }
            array_t& arr = std::get<array_t>(*this);
            
            if (i == arr.size()) { // implicit append.
                arr.resize(i + 1);
            }
            return arr.at(i);
        }

        /// access first element of a json-object matching key (implicitly convert to object; insert if missing)
        auto operator[](std::string_view key) -> value_t&
        {
            if (!std::holds_alternative<object_t>(*this)) {
                *this = object_t{};
            }
            object_t& obj = std::get<object_t>(*this);

            for (auto& kv : obj)
                if (kv.first == key)
            {
                return kv.second;
            }

            obj.emplace_back(std::string(key), value_t{});
            return obj.back().second;
        }

        /// return count of elements matching key (normally 0 or 1)
        auto count(std::string_view key) const -> size_t
        {
            size_t n = 0;
            for (auto& kv : std::get<object_t>(*this)) {
                n += kv.first == key;
            }
            return n;
        }
    };

    static_assert(sizeof(variant_t) == sizeof(value_t), "Can't have additional fields");

    /////////////////////////////////////////////////////////////////////////

    static auto parse(std::string_view json_str, size_t* last_pos_ptr = nullptr) -> value_t
    {
        auto parser = json5(json_str);
        auto ret = parser.parse_value();

        if (last_pos_ptr) {
            *last_pos_ptr = parser.i_;
        } else {
            // expecting single JSON record in json_str
            parser.skipws();
            parser.verify(parser.i_ == json_str.size(), "unexpected trailing non-whitespace characters");
        }

        return ret;
    }

    /// read whole istr into a string and parse, expecting single record.
    static auto parse(std::istream& istr)
    {
        return parse(std::string{
            std::istreambuf_iterator<char>(istr),
            std::istreambuf_iterator<char>()
        });
    }

    /// Serialize as JSON (without JSON5 extensions; may emit NaN, Infinity, -Infinity for corresponding values of real_t)
    inline friend std::ostream& operator<<(std::ostream& os, const value_t& v)
    {
        s_to_stream(os, v);
        return os;
    }

    static auto to_string(const value_t& value) -> std::string
    {
        static thread_local std::ostringstream oss{};
        oss.clear();
        oss.str("");

        s_to_stream(oss, value);
        return oss.str();
    }

#if GX_JSON5_ENABLE_RUN_TESTS
    static void s_run_tests();
#endif

private:
    /////////////////////////////////////////////////////////////////////////
    std::string_view data_ = {};
    size_t           i_    = 0;  // pos into data_; 0 <= i_ <= data_.size()

    json5(std::string_view json_str)
        : data_{ json_str }, i_{ 0 }
    {}

    char peek() const
    {
        verify(i_ < data_.size(), "unexpected end of JSON");
        return data_[i_];
    }

    bool skip(char c)
    {
        skipws();
        return i_ < data_.size() && peek() == c && (i_ += 1);
    }

    bool skip(std::string_view sv, bool do_skipws = true)
    {
        if (do_skipws) {
            skipws();
        }
        return data_.substr(i_, sv.size()) == sv && (i_ += sv.size());
    }

    // skip whitespace and comments.
    void skipws()
    {
        // NB: skip(..., false) in this function to avoid mutual infinite recursion
        while (i_ < data_.size()) {
            while (i_ < data_.size() && std::isspace(peek())) {
                ++i_;
            }
            if (skip("//", false) || skip("#", false)) { // or # -comment
                while (i_ < data_.size() && peek() != '\n') {
                    ++i_;
                }
            } else if (skip("/*", false)) { /* comment */
                while (i_ < data_.size() && !skip("*/", false)) {
                    ++i_;
                }
            } else {
                break;
            }
        }
    }

    void verify(bool cond, const char* msg) const
    {
        if (!cond) {
            throw std::runtime_error(
                "Invalid JSON at pos " + std::to_string(i_)
              + " - " + msg
              + ": " + std::string(data_.substr(i_, 20))
              + "..."
            );
        }
    }

    value_t parse_value()
    {
             if (skip("null"))    return nullptr;
        else if (skip("true"))    return true;
        else if (skip("false"))   return false;
        else if (peek() == '\'')  return parse_string(); // ' ... '
        else if (peek() == '"')   return parse_string(); // " ... "
        else if (peek() == '[')   return parse_array();  // [ ... ]
        else if (peek() == '{')   return parse_object(); // { ... }
        else                      return parse_number();
    }

    static bool s_is_allowed_in_names(char c)
    {
        // NB: '$' allowed because in JSON5 object keys may be an ECMAScript 5.1 IdentifierName.
        return std::isalnum((unsigned char)c) || c == '$' || c == '_' || c == '-';
    };

    std::string parse_string()
    {
        auto str = std::string{};
        const char quote = skip('"') ? '"' : skip('\'') ? '\'' : 0;

        while (i_ < data_.size() && (quote ? peek() != quote : s_is_allowed_in_names(peek()))) {
            if (!quote || peek() != '\\') {
                str += data_.at(i_++);
            } else if ( skip("\\\r\n"   , false)
                     || skip("\\\n"     , false)
                     || skip("\\\r"     , false)
                     || skip("\\\\u2028", false)
                     || skip("\\\\u2029", false))
            {
                continue; // skipping newline in json5 multiline string
            } else if (!skip("\\u")) { 
                str += parse_escaped_char();
            } else if (const auto cp1 = (unsigned)parse_hex(4); 0xD800 <= cp1 && cp1 <= 0xDBFF) {
                // cp1 is high-surrogate in \uXXXX\uXXXX pair
                verify(skip("\\u"), "expected UT8-8 low-surrogate");
                const auto cp2 = (unsigned)parse_hex(4);
                verify(0xDC00 <= cp2 && cp2 <= 0xDFFF, "invalid UTF-8 low-surrogate value");
                s_codepoint_to_utf8(0x10000 + (((cp1 - 0xD800) << 10) | (cp2 - 0xDC00)), str);
            } else {
                s_codepoint_to_utf8(cp1, str);
            }
        }

        verify(!quote || skip(quote), "unterminated string");
        return str;
    }

    char parse_escaped_char()
    {
        verify(skip('\\'), "expected \\");
        verify(peek() != '0' || i_ == data_.size() || !std::isdigit(data_.at(i_ + 1)), "\\0 cannot be followed by a digit");

        switch (peek()) { // https://en.wikipedia.org/wiki/Escape_character
            case '\'': ++i_; return '\'';
            case '"' : ++i_; return '"' ;
            case '\\': ++i_; return '\\';
            case 'b' : ++i_; return '\b';
            case 'f' : ++i_; return '\f';
            case 'n' : ++i_; return '\n';
            case 'r' : ++i_; return '\r';
            case 't' : ++i_; return '\t';
            case '0' : ++i_; return '\0';               // JSON5-specific
            case 'v' : ++i_; return '\v';               // JSON5-specific
            case 'x' : ++i_; return (char)parse_hex(2); // JSON5-specific \x00..\xFF
            default  :       return data_.at(i_++);     // JSON5-specific
        }
    }

    unsigned long parse_hex(size_t len = std::string::npos) // len is 2 (e.g \xFF) or 4 (e.g \uFFFF)
    {
        unsigned short value = 0;
        const auto s = data_.substr(0 + i_, len);
        const auto result = std::from_chars(s.data(), s.data() + s.size(), value, 16);
        verify(!int(result.ec) && (len == std::string::npos || result.ptr == s.data() + len), "invalid number");
        i_ += result.ptr - s.data();
        return value;
    }

    value_t parse_number() // -> int_t | real_t; capture 123 as int_t; 123.0 or 1.23e+2 as real_t
    {
        const int sign = skip('-') ? -1 : (skip('+'), 1);
        if (skip("Infinity") || skip("Inf") || skip("inf") || skip("INF")) {
            // * JSON5 specifies "Infinity" and "NaN", but std::to_chars emits "inf" and "nan", so we support both.
            // * std::numeric_limits<long double>::infinity() == HUGE_VALL, (not INFINITY)
            return HUGE_VALL * sign;
        } else if (skip("NaN") || skip("nan") || skip("NAN")) {
            return std::nanl("");
        } else if (skip("0x") || skip("0X")) {
            const auto x = int_t(parse_hex());
            verify(x >= 0, "can't represent value as int_t");
            return x * sign;
        }

        auto real_value   = real_t(0);
        const auto s      = data_.substr(i_);
        const auto result = std::from_chars(s.data(), s.data() + s.size(), real_value);
        verify(!int(result.ec), "invalid number");
        real_value *= sign;

        // Check the prefix that was parsed as number for floating-point-related chars.
        const bool is_int = s.substr(0, result.ptr - s.data()).find_first_of(".eE") == std::string::npos;
        const size_t parsed_len = result.ptr - s.data();

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
        verify(!is_int || s.at(0) != '0' || parsed_len == 1, "integers can't start with '0'");
        verify(!is_int || real_t(int_t(real_value)) == real_value, "can't represent value as int_t");
#pragma GCC diagnostic pop

        i_ += parsed_len;
        return is_int ? value_t(int_t(real_value)) : value_t(real_value);
    }

    array_t parse_array()
    {
        auto array = array_t{};
        verify(skip('['), "expected '['");

        while (!skip(']')) {
            array.push_back(parse_value());
            verify(skip(',') || peek() == ']', "expected ',' or ']'"); // JSON5: allowing trailing comma before ']'
        }
        return array;
    }

    object_t parse_object()
    {
        auto object = object_t{};
        verify(skip('{'), "expected '{'");

        while (!skip('}')) {
            std::string key = parse_string();
            verify(skip(':'), "expected ':' after key");
            object.emplace_back(std::move(key), parse_value());
            verify(skip(',') || peek() == '}', "expected ',' or '}'"); // JSON5: allowing trailing comma before '}'
        }
        return object;
    }

    static void s_write_escaped(std::ostream& oss, const std::string& arg)
    {
        oss.put('"');

        for (size_t i = 0; i < arg.size(); i++) {
            const char c = arg[i];

            if (' ' <= c && c <= '~' && c != '\\' && c != '"') { // printable
                oss.put(c);
                continue;
            }

            oss.put('\\');

            // Handle UTF8 newlines: \u2028 and \u2029 https://www.compart.com/en/unicode/U+2028
            if (   (unsigned char)c == 0xE2
                && i + 2 < arg.size()
                && (   (unsigned char)arg[i + 1] == 0x80)
                && (   (unsigned char)arg[i + 2] == 0xA8 
                    || (unsigned char)arg[i + 2] == 0xA9))
            {
                oss << ((unsigned char)arg[i + 2] == 0xA8 ? "u2028" : "u2029");
                i += 2;
                continue;
            }

            // \0 and \v not allowed in json - must encode as \uXXXX.
            // Do not need to escape single-quote.
            // Same encoding as by python's print(json.dumps("".join([chr(i) for i in range(0,256)])))
            switch(c) {
            case '"':  oss.put('"');  break;
            case '\\': oss.put('\\'); break;
            case '\b': oss.put('b');  break;
            case '\f': oss.put('f');  break;
            case '\n': oss.put('n');  break;
            case '\r': oss.put('r');  break;
            case '\t': oss.put('t');  break;
            default: // as UTF-8 escape-sequence
                const char* const digits = "0123456789ABCDEF";
                oss << "u00";
                oss.put(digits[(unsigned char)(c) >> 4]);
                oss.put(digits[(unsigned char)(c) & 0xF]);
            }
        }

        oss.put('"');
    }

    static void s_write_floating(std::ostream& oss, const real_t& arg)
    {
        static thread_local auto buf = std::string{};
        buf.clear();
        buf.resize(256);

        auto res = std::to_chars(
            buf.data(),
            buf.data() + buf.size(),
            arg,
            std::chars_format::general,
            (int)oss.precision()
        );
        if (int(res.ec)) {
            throw std::runtime_error("can't print " + std::to_string(arg));
        }
        buf.resize(res.ptr - buf.data());

        if (buf.find_first_not_of("-0123456789") == std::string::npos) {
            buf += ".0"; // i.e as "42.0" rather than as "42" to preserve type-info
        } else if (std::isnan(arg)) { // json5 requires NaN and Infinity instead of nan and inf
            buf = "NaN";
        } else if (std::isinf(arg)) {
            buf = arg < 0 ? "-Infinity" : "Infinity";
        }
        oss << buf;
    }

    static void s_to_stream(std::ostream& oss, const value_t& value)
    {
        std::visit([&](auto&& arg)
        {
            using T = std::decay_t<decltype(arg)>;

            constexpr bool is_object     = std::is_same_v<T, object_t>;
            constexpr bool is_collection = std::is_same_v<T, array_t> || is_object;

                 if constexpr (std::is_same_v<T, std::nullptr_t>) oss << "null";
            else if constexpr (std::is_same_v<T, bool>)           oss << (arg ? "true" : "false");
            else if constexpr (std::is_same_v<T, std::string>)    s_write_escaped(oss, arg);
            else if constexpr (std::is_floating_point_v<T>)       s_write_floating(oss, arg);
            else if constexpr (!is_collection)                    oss << arg;
            else if constexpr (is_collection)
            {
                oss << (is_object ? "{" : "[");
                bool is_first = true;
                for (const auto& item : arg) {
                    oss << (is_first ? "" : ", ");
                    is_first = false;
                    if constexpr (is_object) {
                        s_write_escaped(oss, item.first);
                        oss << ": ";
                        s_to_stream(oss, item.second);
                    } else {
                        s_to_stream(oss, item);
                    }
                }
                oss << (is_object ? '}' : ']');
            }
        }, value);
    }

    // related: https://github.com/dropbox/json11/blob/master/json11.cpp
    //          https://github.com/nlohmann/json/blob/develop/single_include/nlohmann/json.hpp
    static void s_codepoint_to_utf8(unsigned codepoint, std::string& utf8)
    {
        if (codepoint <= 0xFF) { // NB: not 0x7F - encoding extended ascii 0x80..0xFF as single byte like python does
            utf8 += char(codepoint);
        } else if (codepoint <= 0x7FF) {
            utf8 += char(0xC0 | (0x1F & (codepoint >> 6 )));
            utf8 += char(0x80 | (0x3F & (codepoint      )));
        } else if (codepoint <= 0xFFFF) {
            utf8 += char(0xE0 | (0x0F & (codepoint >> 12)));
            utf8 += char(0x80 | (0x3F & (codepoint >> 6 )));
            utf8 += char(0x80 | (0x3F & (codepoint      )));
        } else if (codepoint <= 0x10FFFF) {
            utf8 += char(0xF0 | (0x07 & (codepoint >> 18)));
            utf8 += char(0x80 | (0x3F & (codepoint >> 12)));
            utf8 += char(0x80 | (0x3F & (codepoint >> 6 )));
            utf8 += char(0x80 | (0x3F & (codepoint      )));
        } else {
            throw std::invalid_argument("Invalid codepoint: " + std::to_string(codepoint));
        }
    }
};
} // namespace gx


#if GX_JSON5_ENABLE_RUN_TESTS
#include <iostream>

#define JSON5_CHECK(x) do { if (__builtin_expect(!(x), 0)) throw std::runtime_error("Assertion failed on line " + std::to_string(__LINE__) + ": "#x); } while(0);
void gx::json5::s_run_tests()
{
    {
        const std::string json_str = R"(
        {
            // comment
            # yaml-style comment
            "foo": 42.0, // trailing-comment
            /*
             * multiline comment
             */
            bar: [true, false, null, {}, [[]]],   // unquoted-key
            "baz": [0, -0.0, 0.0, +0.0, -1, 1, +1, -2.0, .2e+1, +2.0, .3, 3., -.3, +.3, 0xFF, inf, Inf, INF, -Infinity, nan, NaN],

            "fred": "a multiline \
string /* with sneaky */ // commentses",

            // trailing comma in array and parent object below
            "corge": ["unescaped 'single-quotes' here", 'unescaped "double-quotes" here',],
        }   // trailing comment)";

        auto j = json5::parse(json_str);

        JSON5_CHECK(j.count("baz") == 1);
        JSON5_CHECK(j.at("baz").is<json5::array_t>());
        JSON5_CHECK(j.at("baz").at(6).get<json5::int_t>() == 1);
        JSON5_CHECK(j.at("baz").at(7).get<json5::real_t>() == -2.0);

        j["qux"][0] = 456.0; // convert to array; set first element to 456.0 (implicit append)

        const std::string expected = R"({"foo": 42.0, "bar": [true, false, null, {}, [[]]], "baz": [0, -0.0, 0.0, 0.0, -1, 1, 1, -2.0, 2.0, 2.0, 0.3, 3.0, -0.3, 0.3, 255, Infinity, Infinity, Infinity, -Infinity, NaN, NaN], "fred": "a multiline string /* with sneaky */ // commentses", "corge": ["unescaped 'single-quotes' here", "unescaped \"double-quotes\" here"], "qux": [456.0]})";
        const std::string out = json5::to_string(j);
        //std::cerr << expected << "\n" << out << "\n";

        JSON5_CHECK(out == expected);
    }

    {
        // from https://github.com/miloyip/nativejson-benchmark/tree/master/data/roundtrip
        const std::string json_str = R"([null, true, false, 0, "foo", {}, [0, 1], {"foo": "bar"}, {"a": null, "foo": "bar"}, -1, -2147483648, -1234567890123456789, -9223372036854775808, 1, 2147483647, 4294967295, 1234567890123456789, 9223372036854775807, 0.0, -0.0, 1.2345, -1.2345, 5e-324, 2.225073858507201e-308, 2.2250738585072014e-308, 1.7976931348623157e+308])";
        
        std::stringstream oss{};
        oss.precision(17); // round-trip high-precision double values in json_str
        oss << json5::parse(json_str);

        //std::cerr << "\n\n" << json_str << "\n" << oss.str() << "\n";
        JSON5_CHECK(oss.str() == json_str);
    }

    // exercise encoding and decoding non-printable characters.
    {
        std::string s; // will round-trip
        for (size_t i = 0; i < 256; i++) {
            s.push_back(i);
        }
        const auto json_str = json5::to_string(s);
        JSON5_CHECK(json5::parse(json_str).get<std::string>() == s);
        JSON5_CHECK(json5::to_string(json5::parse(R"("\x00\x01\x02\x03")")) == R"("\u0000\u0001\u0002\u0003")");

        // check surrogate-pair decoding, e.g. https://www.compart.com/en/unicode/U+1F3BC
        JSON5_CHECK(json5::parse(R"("\uD83C\uDFBC")") == json5::parse(R"("\xf0\x9f\x8e\xbc")"));

        s = R"("\u2028\b\n\r\f\t\u2029")"; // check round-trip of escapes an utf-8 newlines
        JSON5_CHECK(json5::to_string(json5::parse(s)) == s);
    }
}
#undef JSON5_CHECK
#endif


/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*/
