cmake_minimum_required(VERSION 3.1)
project(gx)

set(gx_VERSION_MAJOR 1)
set(gx_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -O0 -ggdb3 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -static -O2 -DNDEBUG")  # TODO: investigate whether can/should enable -march=x86-64-v4

add_subdirectory(src)
