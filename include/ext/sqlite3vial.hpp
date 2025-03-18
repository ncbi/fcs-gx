// sqlite3-vial - single-header c++17 convenience wrapper for sqlite3.h C API

#pragma once
#ifndef SQLITE3VIAL_HPP
#define SQLITE3VIAL_HPP

#include <sqlite3.h>
#include <string>
#include <string_view>
#include <memory>
#include <utility>
#include <tuple>
#include <optional>
#include <type_traits>
#include <stdexcept>
#include <cassert>
#include <cstdio>

namespace sqlite3vial
{

/// -------------------------------------------------------------------------
/// Type for binding and retrieving blobs.
struct blob_t
{
    const void* ptr;
    std::size_t size;
};

/// -------------------------------------------------------------------------
/// Will be thrown if the underlying C api communicates an error via errcode.
struct exception_t : public std::exception 
{
    explicit exception_t(std::string message, int ext_erc = 0) noexcept
        : message_{ std::move(message) }
        , extended_errcode_{ ext_erc }
    {}
    const char* what() const noexcept override { return message_.c_str(); }

    std::string message_ = {};
    int extended_errcode_ = 0;  // output of sqlite3_extended_errcode
};

/// -------------------------------------------------------------------------
using destructor_t = void(*)(void*); // https://sqlite.org/c3ref/bind_blob.html

/// -------------------------------------------------------------------------
namespace detail
{
    template<typename T> struct is_optional                      : std::false_type{};
    template<typename T> struct is_optional<std::optional<T>>    : std::true_type{};
    // NB: if need more of these, implement via is_specialization_of
    
    template<typename T> inline constexpr bool is_optional_v     = is_optional<T>::value;
    template<typename T> inline constexpr bool is_always_false_v = false;

    /// ---------------------------------------------------------------------
    /// Function Traits
    struct function_traits_impl
    {
        // metafunction to get I'th Arg-type from Args (can also do it with std::tuple, but that's extra compile-overhead)
        template <std::size_t I,             typename... Args> struct get_type;
        template <std::size_t I, typename T, typename... Args> struct get_type<I, T, Args...> : get_type<I - 1, Args...> {};
        template <               typename T, typename... Args> struct get_type<0, T, Args...> { using type = T; };

        template<typename Ret, typename... Args>
        struct base
        {
            static constexpr std::size_t arity     = sizeof...(Args);
            using return_type                      = Ret;
            template<std::size_t N> using arg_type = typename get_type<N, Args...>::type;
        };

        template<typename T>
        struct traits : public traits<decltype(&T::operator())> {};

        template<typename Class, typename Ret, typename... Args> // for non-mutable lambdas
        struct traits<Ret(Class::*)(Args...) const> : public base<Ret, Args...> {};

        template<typename Class, typename Ret, typename... Args> // for mutable lambdas
        struct traits<Ret(Class::*)(Args...)>       : public base<Ret, Args...> {};

        template<                typename Ret, typename... Args> // for free functions
        struct traits<Ret(*)(Args...)>              : public base<Ret, Args...> {};
    };

    template<typename F> using function_traits = function_traits_impl::traits<F>;

    /// ----------------------------------------------------------------------
    /// replace nullptr str with ???" for safe appending
    static inline auto make_safe(const char* str) -> const char*
    {
        return str ? str : "???";
    }


    /// ----------------------------------------------------------------------
    /// We need this just for numercis, so no need for full-fledged abi::__cxa_demangle
    template<typename T>
    auto get_typename() -> const char*
    {
        static const char* typenames[] = { "int8_t", "int16_t", "int32_t", "int64_t", "uint8_t", "uint16_t", "uint32_t", "uint64_t" };
        constexpr auto sz = sizeof(T);
        constexpr auto k = sz == 1 ? 0 
                         : sz == 2 ? 1 
                         : sz == 4 ? 2 
                         :           3;

             if constexpr (std::is_same_v<T, bool>)         return "bool";
        else if constexpr (std::is_same_v<T, float>)        return "float";
        else if constexpr (std::is_same_v<T, double>)       return "double";
        else if constexpr (std::is_same_v<T, long double>)  return "long double";
        else if (std::is_integral_v<T> && sz <= 8)          return typenames[std::is_unsigned_v<T> * 4 + k];
        else                                                return typeid(T).name();
    }
    
    /// ----------------------------------------------------------------------
    /// For safe integral conversions.
    template<typename TOut, typename TIn> 
    auto value_preserving_cast(TIn x) -> TOut
    {
        if constexpr (std::is_same_v<TOut, TIn>) {
            return x;
        }

        const auto out = static_cast<TOut>(x);
        if (x != static_cast<TIn>(out) || (x < 0 && std::is_unsigned_v<TOut>)) {
            throw exception_t{
                "Cast of value " + std::to_string(x)    + " from type '" + get_typename<TIn>()
                + "' to type '"  + get_typename<TOut>() + "' is not value-preserving."
            };
        }
        return out;
    }

    /// ----------------------------------------------------------------------
    /// convert the value in val_ptr to T, which may be integral, double, std::string_view, std::string, const char*, blob_t, or std::optional thereof
    template<typename T>
    auto as(::sqlite3_value* val_ptr) -> T
    {
        assert(val_ptr);

        constexpr bool is_constructible_from_c_str =
               std::is_same_v<T, std::string_view> 
            || std::is_same_v<T, std::string>
            || std::is_same_v<T, const char*>;

        const auto db_type = sqlite3_value_type(val_ptr);

        const bool is_type_compatible =
               std::is_same_v<T, ::sqlite3_value*> 
            || detail::is_optional_v<T>   // will recurse on T::value_type
            || (db_type == SQLITE_INTEGER && std::is_integral_v<T>)
            || (db_type == SQLITE_FLOAT   && std::is_same_v<T, double>)
            || (db_type == SQLITE_TEXT    && is_constructible_from_c_str)
            || (db_type == SQLITE_BLOB    && std::is_same_v<T, blob_t>);
 
        if (!is_type_compatible) {
            // https://sqlite.org/c3ref/c_blob.html
            static const char* sql_typenames[] = { "???", "INTEGER", "FLOAT", "TEXT", "BLOB", "NULL" };

            throw exception_t{
                std::string{ "Can't convert db-value is of db-type " } + sql_typenames[db_type % (int)sizeof(sql_typenames)]
                + " to requested incompatbile callback's arg-type '" + get_typename<T>() + "'"
            };
        }
             if constexpr (std::is_same_v<T, sqlite3_value*>)   return val_ptr;
        else if constexpr (std::is_integral_v<T>)               return value_preserving_cast<T>(sqlite3_value_int64(val_ptr));
        else if constexpr (std::is_same_v<T, double>)           return sqlite3_value_double(val_ptr);
        else if constexpr (std::is_same_v<T, blob_t>)           return blob_t{ sqlite3_value_blob(val_ptr), (std::size_t)sqlite3_value_bytes(val_ptr) };
        else if constexpr (is_constructible_from_c_str)         return T{ (const char*)sqlite3_value_text(val_ptr) };
        else if constexpr (detail::is_optional_v<T>)            return db_type == SQLITE_NULL ? T{} : as<typename T::value_type>(val_ptr);
        else static_assert(detail::is_always_false_v<T>, "Unsupported arg type - expecting integral, double, std::string_view, std::string, blob_t, or std::optional thereof.");
    }

    /// ----------------------------------------------------------------------
    /// Invoke appropriate version of sqlite3_result() depending on the type of T
    template<typename T>
    auto set_sqlite3_result(::sqlite3_context* ctx, const T& value) -> void
    {
        assert(ctx);

        // NB: similar dispatch as in bind_arg(...) below
             if constexpr (std::is_same_v<T, std::nullptr_t>)   sqlite3_result_null(        ctx);
        else if constexpr (std::is_integral_v<T>)               sqlite3_result_int64(       ctx, value_preserving_cast<int64_t>(value));
        else if constexpr (std::is_floating_point_v<T>)         sqlite3_result_double(      ctx, value);
        else if constexpr (std::is_same_v<T, blob_t>)           sqlite3_result_blob64(      ctx, value.ptr, value.size, SQLITE_TRANSIENT);
        else if constexpr (detail::is_optional_v<T>)            value ? set_sqlite3_result( ctx, *value) : sqlite3_result_null(ctx);
        else if constexpr (std::is_convertible_v<T, std::string_view>) {
            const auto sv = std::string_view{ value };
            sqlite3_result_text64(ctx, sv.data(), sv.size(), SQLITE_TRANSIENT, SQLITE_UTF8);
        } else {
            static_assert(detail::is_always_false_v<T>, "Unsupported type - expecting nullptr_t, numeric, convertible-to-std::string_view, blob_t, or std::optional thereof.");
        }
    }

    /// ----------------------------------------------------------------------
    /// Invoke appropriate version of sqlite3_bind() depending on the type of T
    template<typename T>
    auto bind_arg(::sqlite3_stmt* stmt, int index, const T& value, destructor_t destructor = SQLITE_TRANSIENT) -> int
    {
        assert(stmt);
        auto rc = SQLITE_OK;

             if constexpr (std::is_same_v<T, std::nullptr_t>)   rc = sqlite3_bind_null(   stmt, index);
        else if constexpr (std::is_integral_v<T>)               rc = sqlite3_bind_int64(  stmt, index, value_preserving_cast<int64_t>(value));
        else if constexpr (std::is_floating_point_v<T>)         rc = sqlite3_bind_double( stmt, index, value);
        else if constexpr (std::is_same_v<T, blob_t>)           rc = sqlite3_bind_blob64( stmt, index, value.ptr, value.size, destructor);
        else if constexpr (detail::is_optional_v<T>)            return value ? bind_arg(  stmt, index, *value, destructor) 
                                                                             : bind_arg(  stmt, index, nullptr);
        else if constexpr (std::is_convertible_v<T, std::string_view>) {
            const auto sv = std::string_view{ value };
            rc = sqlite3_bind_text64(stmt, index, sv.data(), sv.size(), destructor, SQLITE_UTF8);
        } else {
            static_assert(detail::is_always_false_v<T>, "Unsupported type - expecting nullptr_t, numeric, convertible-to-std::string_view, blob_t, or std::optional thereof.");
        }

        return rc;
    }

    /// ----------------------------------------------------------------------
    /// return fn(as<Arg0>(argv[0]), as<Arg1>(argv[1]), ...) with fn's Arg-types
    template<typename F, std::size_t... I>
    auto invoke(F& fn, sqlite3_value** argv, std::index_sequence<I...>)
    {
        assert(argv);
        return fn(as<typename function_traits<F>::template arg_type<I>>(argv[I])...);
    }

    /// ----------------------------------------------------------------------
    [[noreturn]]
    inline auto throw_ex(::sqlite3& db, std::string(msg)) -> void
    {
        const int extended_errcode = sqlite3_extended_errcode(&db);
        throw exception_t{
            std::move(msg) 
                + "\n(sqlite3_extended_errcode: " + std::to_string(extended_errcode)
                + ", sqlite3_errmsg: " + make_safe(sqlite3_errmsg(&db)) + ")"
          , extended_errcode
        };
    }
}

/// -------------------------------------------------------------------------
/// RAII-wrapper for ::sqlite3_stmt pointer that calls sqlite3_finalize on destruction.
class stmt_t 
{
public:

    /// ---------------------------------------------------------------------
    stmt_t(::sqlite3& db, std::string_view sql, unsigned int prep_flags = 0)
    {
        auto stmt_ptr = (::sqlite3_stmt*)nullptr;
        
        const int rc = 
            sqlite3_prepare_v3(&db, sql.data(), (int)sql.size(), prep_flags, &stmt_ptr, nullptr);

        stmt_ = stmt_ptr_t{ stmt_ptr, &sqlite3_finalize };

        if (rc != SQLITE_OK) {
            // NB: not via x_throw, because it can't get the sql unless prepare succeeded
            detail::throw_ex(db, std::string{ "Failed to prepare statement: "} + sql.data());
        }

        bound_param_flags_.resize((std::size_t)sqlite3_bind_parameter_count(stmt_.get()));
    }

    /// ---------------------------------------------------------------------
    auto get_ptr() -> ::sqlite3_stmt*
    {
        return stmt_.get();
    }

    /// ---------------------------------------------------------------------    
    template<typename T>
    auto bind_arg(const int index, const T& value, destructor_t destructor = SQLITE_TRANSIENT) -> stmt_t&
    {
        assert(index > 0 && index <= (int)bound_param_flags_.size());

        // must reset before rebinding, in case fetch() was invoked. (NB: this does not reset the parameters)
        reset(false);

        if (const auto rc = detail::bind_arg(stmt_.get(), index, value, destructor); rc != SQLITE_OK) {
            x_throw("Failed to bind parameter " + std::to_string(index));
        }

        bound_param_flags_.at((std::size_t)index - 1) = 1; // mark as bound
        return *this;
    }

    /// ---------------------------------------------------------------------
    template<typename T>
    auto bind_arg(const char* col_name, const T& value, destructor_t destructor = SQLITE_TRANSIENT) -> stmt_t&
    {
        assert(col_name);
        const int col_idx = sqlite3_bind_parameter_index(stmt_.get(), col_name);
        if (col_idx == 0) {
            x_throw(std::string{"No such column-name: "} + col_name);
        }
        return bind_arg(col_idx, value, destructor);
    }

    /// ---------------------------------------------------------------------
    template<typename... Args>
    auto bind_args(Args&&... args) -> stmt_t&
    {
        reset(true);   // Reset and clear the current bindings.
        int index = 1; // NB: in C-api the column indexes are 1-based when binding, and 0-based when retrieving.
        (bind_arg(index++, std::forward<Args>(args)), ...);
        return *this;
    }

    /// ---------------------------------------------------------------------
    /// Invoke `fn` on rows in resultset; break-out and return as soon as it returns something truthy.
    /// Return default-constructed Ret{} if fn never returns something truthy.
    /// On subsequent invocation(s), continue where left off (without resetting the statement).
    ///
    /// F's args can be int, int64_t, double, string_view, blob_t, or std::optional thereof, or sqlite3_value*
    template<typename F, typename Ret = typename detail::function_traits<F>::return_type>
    auto fetch(F fn) -> Ret
    {
        const auto orig_rc = sqlite3_errcode(sqlite3_db_handle(stmt_.get()));

        // if this is not a subsequent invocation, check for unbound args.
        if (orig_rc != SQLITE_ROW && orig_rc != SQLITE_DONE)
            if (const auto pos = bound_param_flags_.find((char)0, 0ull); pos != std::string::npos)
        {
            x_throw(
                "Statement has unbound param #" + std::to_string(pos + 1) 
                + ", name:" + detail::make_safe(sqlite3_bind_parameter_name(stmt_.get(), (int)pos + 1))
            );
        }

        constexpr auto arity = detail::function_traits<F>::arity;
        constexpr auto ix_seq = std::make_index_sequence<arity>{};

        auto row_count = 0ul;
        if (orig_rc != SQLITE_DONE) // in case fetchone() is called more than once after reaching SQLITE_DONE
            for (int rc = sqlite3_step(stmt_.get()); rc != SQLITE_DONE; rc = sqlite3_step(stmt_.get()))
        {
            ++row_count;
            if (rc != SQLITE_ROW) {
                x_throw("Failed to execute statement; error on row " + std::to_string(row_count));
            } else if (const auto num_cols = (std::size_t)sqlite3_column_count(stmt_.get()); num_cols < arity) {
                x_throw("Resultset contains "                    + std::to_string(num_cols)
                        + " columns; callback expects at least " + std::to_string(arity));
            }

            try {
                if constexpr (std::is_void_v<Ret>) {
                    x_invoke(fn, ix_seq);
                } else if (auto x = x_invoke(fn, ix_seq); x) {
                    return x; // callback returned something truthy - pass it out.
                }
            } catch (...) {
                std::fprintf(stderr, "Uncaught exception in callback while processing row #%ld: ", row_count);
                x_print_current_row(stderr);
                throw;
            }
        }

        if constexpr (!std::is_void_v<Ret>) {
            return Ret{};
        }   
    }

    /// ---------------------------------------------------------------------
    auto exec() -> void
    {
        fetch([](){});
    }

    /// ---------------------------------------------------------------------
    /// Return next row from the query as std::optional<Arg> or std::optional<std::tuple<Args...>>
    /// (empty-optional if no more rows)
    /// NB: naming similar to python.sqlite3
    template<typename... Args>
    auto fetchone() // -> std::optional<Arg> if single arg, else std::optional<std::tuple<Args...>>
    {
        static_assert(!(... || (std::is_same_v<Args, blob_t> || std::is_same_v<Args, std::string_view> || std::is_same_v<Args, const char*>)), "Can't have blob_t, string_view, or const char* in output type(s) as they will become dangling. Get as std::string or use stmt_t::fetch(fn) instead.");

        if constexpr ((sizeof...(Args)) == 1) {
            return fetch([](Args... args) { return std::make_optional(std::move(args)...); });
        } else {
            return fetch([](Args... args) { return std::make_optional(std::make_tuple(std::move(args)...)); });
        }
    }

    auto reset(bool clear_args = false) -> void
    {
        if (const auto rc = sqlite3_reset(stmt_.get()); rc != SQLITE_OK) {
            x_throw("Failed to reset statement.");
        }

        if (const auto rc = !clear_args? SQLITE_OK : sqlite3_clear_bindings(stmt_.get()); rc != SQLITE_OK) {
            x_throw("Failed to clear bindings.");
        }
    }

private:
    using stmt_ptr_t = std::unique_ptr<::sqlite3_stmt, decltype(&sqlite3_finalize)>;

    stmt_ptr_t  stmt_              = stmt_ptr_t{ nullptr, &sqlite3_finalize };
    std::string bound_param_flags_ = {};    // contains 0 at positions corresponding to unbound params

    /// ---------------------------------------------------------------------
    [[noreturn]]
    auto x_throw(std::string msg) -> void
    {
        detail::throw_ex(
            *sqlite3_db_handle(stmt_.get()),
            std::move(msg) + "\nSQL:\n----\n" + detail::make_safe(sqlite3_sql(stmt_.get())) + "\n----"
        );
    }


    /// ---------------------------------------------------------------------
    void x_print_current_row(FILE* fp)
    {
        const auto stmt = stmt_.get();

        for (int i = 0, col_count = sqlite3_column_count(stmt); i < col_count; ++i) {
            std::fprintf(fp, "%s %s:", i == 0 ? "{" : ",", detail::make_safe(sqlite3_column_name(stmt, i)));

            switch (const int col_type = sqlite3_column_type(stmt, i); col_type) {
                case SQLITE_BLOB:
                    std::fprintf(fp, "BLOB[%d]", sqlite3_column_bytes(stmt, i));
                    break;
                case SQLITE_NULL:
                    std::fprintf(fp, "(null)");
                    break;
                default:
                    const char* qw = col_type == SQLITE_TEXT ? "\"" : "";
                    std::fprintf(fp, "%s%s%s", qw, (const char*)sqlite3_column_text(stmt, i), qw);
                    break;
            }
        }

        std::fprintf(fp, " }\n");
    }

    /// ---------------------------------------------------------------------
    /// Used from x_invoke() to extract a value of a column as type expected by
    /// the corresponding (index'th) arg of the user-provided callback function.
    template<typename T>
    auto x_get(const int index) -> T
    {
        try {
            return detail::as<T>(sqlite3_column_value(stmt_.get(), index));
        } catch (std::exception& e) {
            x_throw("Can't convert column #" + std::to_string(index + 1)
                    + ", a.k.a '"            + detail::make_safe(sqlite3_column_name(stmt_.get(), index))
                    + "': "                  + e.what());
        }
    }

    /// ---------------------------------------------------------------------
    /// Invoke fn(x_get<Arg0>(0), x_get<Arg1>(1), ...) with fn's arg-types.
    template<typename F, std::size_t... I>
    auto x_invoke(F& fn, std::index_sequence<I...>)
    {
        return fn(x_get<typename detail::function_traits<F>::template arg_type<I>>(I)...);
    }
};

///--------------------------------------------------------------------------
class db_t
{
public:
    db_t(
        const std::string_view path,
        const int              flags  = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
        const char*            zVfs   = nullptr
    )
    {
        auto ptr      = (::sqlite3*)nullptr;
        const auto rc = sqlite3_open_v2(path.data(), &ptr, flags, zVfs);
        db_ptr_       = db_ptr_t(ptr, &sqlite3_close);
        
        if (rc != SQLITE_OK) {
            detail::throw_ex(*ptr, std::string{ "Failed to open db " } + path.data());
        }
    }

    /// Create db_t from an existing sqlite3 instance (non-owning; will not call sqlite3_close on destruction)
    db_t(::sqlite3& db)
        : db_ptr_{ db_ptr_t{ &db, [](::sqlite3*){ return 0; }}} // with no-op deleter
    {}

    /// ---------------------------------------------------------------------
    auto get_ptr() -> ::sqlite3*
    {
        return db_ptr_.get();
    }

    /// ---------------------------------------------------------------------
    auto make_stmt(std::string_view sql) -> stmt_t
    {
        return stmt_t{ *db_ptr_.get(), sql };
    }

    /// ---------------------------------------------------------------------
    template<typename F>
    auto transaction(F fn, const std::string_view begin_sql = "BEGIN TRANSACTION") -> void
    {
        make_stmt(begin_sql).exec(); // BEGIN [EXCLUSIVE|IMMEDIATE|DEFERRED] TRANSACTION
        try {
            fn();
        } catch (...) {
            make_stmt("ROLLBACK TRANSACTION").exec();
            throw;
        }
        make_stmt("COMMIT TRANSACTION").exec();
    }

    /// ---------------------------------------------------------------------
    // For aggregate funcions, fn1 is `xStep` and fn2 is `xFinal` https://www.sqlite.org/c3ref/create_function.html
    template<typename F1, typename F2 = std::nullptr_t>
    auto register_function(std::string name, F1 fn1, F2 fn2 = nullptr) -> void
    {
        assert(!name.empty());

        // We'll need to call fn1 and fn2 from the lambda-wrappers below,
        // which must be decayable to function-pointer, so it must
        // be non-capturing, so fn1 and fn2 must be made static to be visible. 
        //
        // We'll access them via s_fn1_ptr,and s_fn2_ptr, and defer setting
        // these until the very end, after sqlite3_create_function succeeds below.
        //
        // NB: obviously, the calling code must ensure that any 
        // reference captures in `fn`s must outlive their invocations.
        static auto s_fn1_ptr = (decltype(&fn1))nullptr;
        static auto s_fn2_ptr = (decltype(&fn2))nullptr;
        (void)s_fn2_ptr; // silence 'unused' compiler warning
        
        if (s_fn1_ptr) {
            throw exception_t{
                std::string{ "Failed to register function " } + name
                + " - already registered function of type " + typeid(F1).name()
                + " - only can do it once per F's type."
            };
        }

        // NB: these must be static constexpr in order to be visible from non-capturing lambda in clang
        static constexpr auto arity  = detail::function_traits<F1>::arity;
        static constexpr auto ix_seq = std::make_index_sequence<arity>{};
        using ret_t = typename detail::function_traits<F1>::return_type;

        static_assert(std::is_null_pointer_v<F2> != std::is_void_v<ret_t>, "F1 must return a value iff F2 is not specified");

        // unary+ decays the non-capturing lambda into registereable function-pointer.
        static const auto fn1_wr = +[](sqlite3_context* context, int argc, sqlite3_value** argv)
        {
            (void)argc;
            assert(argc == arity); // this is also ensured internaly in sqlite3.

            try {
                // detail::invoke converts db-values from argv to the arg-types expected by s_fn.
                if constexpr (std::is_void_v<ret_t>) { // lambda is a wrapper for xStep in aggregate
                    detail::invoke(*s_fn1_ptr, argv, ix_seq);
                } else {                               // lambda is a wrapper for xFunc
                    auto result = detail::invoke(*s_fn1_ptr, argv, ix_seq);
                    detail::set_sqlite3_result(context, result);
                }
            } catch (std::exception& e) {
                sqlite3_result_error(context, e.what(), -1);
            }
        };

        auto rc = SQLITE_OK;

        if constexpr (std::is_same_v<F2, std::nullptr_t>) { // registering a regular xFunc function
            rc = sqlite3_create_function(db_ptr_.get(), name.c_str(), arity, SQLITE_UTF8, nullptr, fn1_wr, nullptr, nullptr);
        } else { // registering an aggregate function - xStep & xFinal pair; fn2_wr is sqlite3-xFinal wrapper
            static const auto fn2_wr = +[](sqlite3_context* context)
            {
                try {
                    auto result = (*s_fn2_ptr)();
                    detail::set_sqlite3_result(context, result);
                } catch (std::exception& e) {
                    sqlite3_result_error(context, e.what(), -1);
                }
            };
            rc = sqlite3_create_function(db_ptr_.get(), name.c_str(), arity, SQLITE_UTF8, nullptr, nullptr, fn1_wr, fn2_wr);
        }

        if (rc != SQLITE_OK) {
            detail::throw_ex(*db_ptr_, "Failed to register funciton " + name);
        }

        static auto s_fn1 = std::move(fn1);
        static auto s_fn2 = std::move(fn2);
        s_fn1_ptr = &s_fn1;
        s_fn2_ptr = &s_fn2;
    }

private:
    using db_ptr_t   = std::unique_ptr<::sqlite3, decltype(&sqlite3_close)>;
    db_ptr_t db_ptr_ = db_ptr_t{ nullptr, &sqlite3_close };
};



#if SQLITE3VIAL_ENABLE_TESTS
/// -------------------------------------------------------------------------
static inline void test() 
{
    /*  
        // Basic usage:
        auto db = sqlite3vial::db_t(db_path);

        db.make_stmt(sql)
            .bind_args(arg1, arg2, ...)
            .fetch([&](
                int64_t               col1, 
                std::string_view      col2,
                std::optional<double> nullable_col3,
                sqlite3vial::blob_t   col4,
                ...
            )
        {
            // Consume row
        });
    */

    auto db = sqlite3vial::db_t(":memory:"); // closes on destruction

    // NB: STRICT enforces the strict-typing during insertions. https://sqlite.org/stricttables.html
    // Disabled, requires v3.37.0
    db.make_stmt(R"(
        CREATE TABLE Foos (
            id           INTEGER PRIMARY KEY,
            name         TEXT NOT NULL,
            real         REAL NOT NULL,
            data         BLOB,
            nullable_int INTEGER
        ) -- STRICT
    )").exec();

    using opt_int64_t = std::optional<int64_t>;

    // Example of inserting some rows 
    std::fprintf(stderr, "\nUsing sqlite3 v%s.\nInserting three rows...\n", sqlite3_libversion());
    {
        const uint16_t blob_data[] = {1, 2, 3, 4, 5};
        const auto blob = sqlite3vial::blob_t{ blob_data, sizeof(blob_data) };

        auto stmt = db.make_stmt("INSERT INTO Foos (name, real, data, nullable_int) VALUES (?, ?, ?, ?)");
        stmt.bind_args("First row",  1.1111, blob,    1111          ).exec();
        stmt.bind_args("",           0.0f,   nullptr, opt_int64_t{} ).exec();
    }
    db.make_stmt("INSERT INTO Foos VALUES (3, 'Third row', 3.333, NULL, 3333)").exec();
    std::fprintf(stderr, "%d rows affected.\n", sqlite3_total_changes(db.get_ptr()));

    // Examples of returning scalar values from query as std::optional<T> or std::optional<std::tuple<Ts...>>
    {
        assert(db.make_stmt("SELECT COUNT(*) FROM Foos").fetchone<int64_t>() == std::make_optional(3));

        // If the db-type is nullable: stmt.fetchone<std::optional<T>>(),
        // and the return type is std::optional<std::optional<T>>:
        //      the outer optional is not empty iff we have a row, 
        //      the inner optional is not empty iff the db-value is not null.
        assert(db.make_stmt("SELECT nullable_int FROM Foos WHERE id =  1").fetchone<opt_int64_t>() == std::make_optional(opt_int64_t{1111}));
        assert(db.make_stmt("SELECT nullable_int FROM Foos WHERE id =  2").fetchone<opt_int64_t>() == std::make_optional(opt_int64_t{}));
        assert(db.make_stmt("SELECT nullable_int FROM Foos WHERE id = -1").fetchone<opt_int64_t>() == std::optional<opt_int64_t>{});

        // Example of binding named args.
        assert(db.make_stmt("SELECT COUNT(*) FROM Foos WHERE id >= @id AND real > @real")
                 .bind_args(3, 0.0)                         // Bind all args by position.
                 .bind_arg(    1, 3).bind_arg(      2, 0.0) // Equivalent to the above: bind args by 1-based column-index.
                 .bind_arg("@id", 3).bind_arg("@real", 0.0) // Equivalent to the above: bind args by name.
                 .fetchone<int64_t>() == std::make_optional(1));

        // Get multiple columns as std::optional<std::tuple<...>>,
        // consuming all rows until encountering empty-tuple:
        auto stmt = db.make_stmt("SELECT id, name FROM Foos");
        auto concat_result = std::string{};
        for (auto opt_tpl = stmt.fetchone<uint64_t, std::string>();
                  opt_tpl;
                  opt_tpl = stmt.fetchone<uint64_t, std::string>())
        {
            concat_result += std::to_string(std::get<0>(*opt_tpl)) + ":'" + std::get<1>(*opt_tpl) + "';  ";
        }
        assert(concat_result == "1:'First row';  2:'';  3:'Third row';  ");
        assert(!stmt.fetchone()); // expecting empty-optional after processed all rows
        assert(!stmt.fetchone()); // this subsequent call must be tested too
    }

    // Example of registering a user-function; will call it from the SELECT below.
    // Args and return-type must be type-compatible with db-types;
    // std::optional for nullable values and sqlite3_value* for heterogeneous types.
    auto exec_count = 0ul;
    db.register_function("my_length", [&](opt_int64_t arg)
    {
        ++exec_count;
        return !arg ? opt_int64_t{} : (int64_t)std::to_string(*arg).size();
    });

    // Example of consuming the results of a query with a lambda passed to fetch()
    {
        auto rowcount = 0ull;
        const auto max_id = 3;
        db.make_stmt("SELECT *, my_length(nullable_int) FROM Foos where id <= @max_id")
            .bind_args(max_id)
            .fetch([&](
                int64_t                id,
                std::string_view       name,           // Can also be std::string or const char*
            //  int64_t                real,           // Will throw at runtime because db-type is SQLITE_FLOAT.
                double                 real,
                std::optional<blob_t>  data,           // using std::optional<> for nullable columns.
            //  int64_t                nullable_int,   // Will throw at runtime because db-type is NULL and type is not optional.
            //  std::optional<uint8_t> nullable_int,   // will throw at runtime because db-value does not fit in uint8_t.
                opt_int64_t            nullable_int,
                sqlite3_value*         value           // As sqlite3_value* if column contains heterogeneous types.
            // ,int64_t                extra_col       // Will throw at runtime - too many columns.
            // ,float                  flt_col         // Will not compile - not a supported type.
            // ,auto                   generic_col     // Will not compile - types must be explicit.
            ) mutable                                  // Lambda can be mutable.
        {
            ++rowcount;
            std::fprintf(stderr,
                "Row %lld | id: %ld, name:'%s', real:%f, blob-size:%s, nullable_int:%s, value:%s\n",
                rowcount, id, name.data(), real,
                (data         ? std::to_string(data->size).c_str()    : "(null)"),
                (nullable_int ? std::to_string(*nullable_int).c_str() : "(null)"),
                sqlite3_value_text(value)
            );

            // throw std::runtime_error("NB: If the callback throws, current row will be printed to stderr during the stack-unwinding");

            // return 42; // returning any truthy value will break-out and return this value from fetch()
            // subsequent invocation of exec will continue on the fetchone row (unless reset() is called)
        });
        assert(rowcount == 3);
        assert(exec_count == 3);
    }

    // Example of using a user-defined aggregate function that computes average over non-null values.
    {
        auto sum = 0.0;
        auto n = 0ul;
        db.register_function(
            "my_avg",
            // sqlite3 xStep function that updates the state for the aggregate.
            [&](std::optional<double> x) { n += x ? 1u : 0u; sum += x ? *x : 0; },
            // sqlite3 xFinal nullary function that computes the aggregate from state.
            [&]{ return n == 0 ? 0.0 : sum / (double)n; }
        );

        const auto opt_diff = db.make_stmt("SELECT ABS(AVG(real) - my_avg(real)) FROM Foos").fetchone<double>();
        (void)opt_diff; // avoiding false "unsued" warning
        assert(opt_diff && 0 <= *opt_diff && *opt_diff <= 0); // checking equality with <=, avoiding float-equal warning
    }

    // Transaction example.
    try {
        db.transaction([&]
        {
            // implicit BEGIN TRANSACTION
            db.make_stmt("UPDATE Foos SET name = 'Oops' WHERE id = 1" ).exec();
            throw std::runtime_error("Exiting by exception will roll back the transaction.");
            // implicit COMMIT TRANSACTION
        });
    } catch (const std::runtime_error&) {
        // Demonstrate that the Foos.name did not change after the rollback.
        const auto opt_name = db.make_stmt("SELECT name From Foos WHERE id = 1").fetchone<std::string>();
        assert(opt_name.value() == "First row");
    }

    std::fprintf(stderr, "\nAll tests passed successfully.\n");
}
#endif // GX_SQLITE3_ENABLE_TESTS

} // namespace sqlite3

#endif // GX_SQLITE3_HPP
