set(CMAKE_CXX_STANDARD 17)

add_compile_options(
    -pthread 
    -g                              # Always include debug symbols enable stack-traces with gdb or libbackward.
    -msse4.2                        # grep for <smmintrin.h> for the code using SSE intrinsics.
    -fno-strict-aliasing            # Allow UB-free reinterpret_cast-ing of PODs from memory-mapped memory.
    -fwrapv                         # Assume 2s-complement arithemtic for signed-int overflow. (consider -ftrapv in debug mode).
    -ftemplate-backtrace-limit=1    # Make template-related error-messages less verbose.
    -Wall -Wextra                   # Everything below is warnings-related.
    -Wconversion -Wdeprecated-declarations -Wstrict-aliasing
    -Wlogical-op -Wmissing-declarations -Wpedantic -Wshadow
    -Wsuggest-attribute=format -Wswitch -Wpointer-arith
    -Wmissing-include-dirs -Winvalid-pch -Wmissing-format-attribute
    -Wfloat-equal -Winit-self
    # -Woverloaded-virtual -Wsign-promo can't conditionally disable these when building sqlite3
    -Wvla -Wredundant-decls
    -Wcast-align -Wcast-qual
    -Wno-misleading-indentation
    -Wno-sign-conversion            # Temporary
    -Wno-unknown-pragmas            # Ignore clang-specific pragmas in GCC
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(
        -Wno-unknown-warning-option     # To silence clang for gcc-specific flags
        -Wno-unused-local-typedef       # Noise from clang
        -Wno-unused-function            # Noise from clang
    )
endif()

include_directories(${PROJECT_SOURCE_DIR}/include/)

set_source_files_properties(main.cpp PROPERTIES COMPILE_FLAGS "-DGIT_REVISION=`../../build-support/git_describe.sh`")

file(GLOB_RECURSE SOURCES *.cpp)
add_executable(gx ${SOURCES})

#############################################################################
# SQLITE3

# for dynamically linking with system sqlite3
# find_package(SQLite3 REQUIRED)
# target_link_libraries(gx SQLite::SQLite3)

# For building and statically linking sqlite3
include(FetchContent)
FetchContent_Declare(
  sqlite3
  URL      https://www.sqlite.org/2024/sqlite-amalgamation-3460100.zip
  URL_HASH SHA3_256=af6aae8d3eccc608857c63cf56efbadc70da48b5c719446b353ed88dded1e288
)
FetchContent_MakeAvailable(sqlite3)
add_library(sqlite3 STATIC ${sqlite3_SOURCE_DIR}/sqlite3.c)
target_compile_options(sqlite3 PRIVATE -w)  # disable warnings

target_compile_definitions(sqlite3 PRIVATE
    SQLITE_ENABLE_API_ARMOR=1
    SQLITE_THREADSAFE=2
    SQLITE_ENABLE_FTS5=1
    SQLITE_ENABLE_FTS3_PARENTHESIS=1
    SQLITE_ENABLE_RTREE=1   # may need for loc-indexing
    SQLITE_ENABLE_COLUMN_METADATA=1
    SQLITE_USE_URI=1
    SQLITE_TEMP_STORE=1 # in memory
    SQLITE_OMIT_LOAD_EXTENSION=1
)

target_include_directories(sqlite3 PUBLIC ${sqlite3_SOURCE_DIR})
#############################################################################

# link the LICENSE blob
set(LINK_LICENSE "-Wl,--format=binary -Wl,../../LICENSE -Wl,--format=default")

# NB: for statically linking in Galaxy/bioconda, may need to link rt and pthread as follows:
#-Wl,--whole-archive -lrt -lpthread -Wl,--no-whole-archive
target_link_libraries(gx PRIVATE sqlite3 -flto -lrt -lpthread ${LINK_LICENSE})

