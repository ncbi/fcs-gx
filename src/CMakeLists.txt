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
    -Woverloaded-virtual -Wsign-promo -Wvla -Wredundant-decls
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

link_libraries("-flto -pthread -Wl,--format=binary -Wl,../../LICENSE -Wl,--format=default") # see show-license command in main.cpp

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)

file(GLOB_RECURSE SOURCES *.cpp)

add_executable(gx ${SOURCES})
