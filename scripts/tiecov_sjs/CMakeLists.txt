cmake_minimum_required(VERSION 3.10)
project(tiebrush)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

set(CMAKE_CXX_STANDARD 11)

include("cmake/zlib.cmake")
#include("cmake/libdeflate.cmake")
include("cmake/htslib.cmake")

include(FetchContent)

FetchContent_Declare(
        gclib
        GIT_REPOSITORY https://github.com/gpertea/gclib.git
)

FetchContent_GetProperties(gclib)
if(NOT gclib_POPULATED)
    FetchContent_Populate(gclib)
endif()

message(STATUS "gclib source dir: ${gclib_SOURCE_DIR}")
message(STATUS "gclib binary dir: ${gclib_BINARY_DIR}")

include_directories(${gclib_SOURCE_DIR})

add_executable(tiebrush src/tiebrush.cpp
                        src/GSam.cpp
                        src/GSam.h
                        src/tmerge.cpp
                        src/tmerge.h
                        ${gclib_SOURCE_DIR}/GStr.h
                        ${gclib_SOURCE_DIR}/GStr.cpp
                        ${gclib_SOURCE_DIR}/GArgs.h
                        ${gclib_SOURCE_DIR}/GArgs.cpp
                        ${gclib_SOURCE_DIR}/GBase.h
                        ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(tiebrush htslib zlib)
target_link_libraries(tiebrush ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread)

add_executable(tiecov src/tiecov.cpp
                    src/GSam.cpp
                    src/GSam.h
                    src/tmerge.cpp
                    src/tmerge.h
                    ${gclib_SOURCE_DIR}/GStr.h
                    ${gclib_SOURCE_DIR}/GStr.cpp
                    ${gclib_SOURCE_DIR}/GArgs.h
                    ${gclib_SOURCE_DIR}/GArgs.cpp
                    ${gclib_SOURCE_DIR}/GBase.h
                    ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(tiecov htslib zlib)
target_link_libraries(tiecov ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread)

add_executable(tiecov_sjs src/tiecov_sjs.cpp
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${gclib_SOURCE_DIR}/GStr.h
        ${gclib_SOURCE_DIR}/GStr.cpp
        ${gclib_SOURCE_DIR}/GArgs.h
        ${gclib_SOURCE_DIR}/GArgs.cpp
        ${gclib_SOURCE_DIR}/GBase.h
        ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(tiecov_sjs htslib zlib)
target_link_libraries(tiecov_sjs ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread)

add_executable(untie src/untie.cpp
                src/GSam.cpp
                src/GSam.h
                src/tmerge.cpp
                src/tmerge.h
                ${gclib_SOURCE_DIR}/GStr.h
                ${gclib_SOURCE_DIR}/GStr.cpp
                ${gclib_SOURCE_DIR}/GArgs.h
                ${gclib_SOURCE_DIR}/GArgs.cpp
                ${gclib_SOURCE_DIR}/GBase.h
                ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(untie htslib zlib)
target_link_libraries(untie ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread)

add_executable(extract_pairs src/extract_pairs.cpp
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${gclib_SOURCE_DIR}/GStr.h
        ${gclib_SOURCE_DIR}/GStr.cpp
        ${gclib_SOURCE_DIR}/GArgs.h
        ${gclib_SOURCE_DIR}/GArgs.cpp
        ${gclib_SOURCE_DIR}/GBase.h
        ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(extract_pairs htslib zlib)
target_link_libraries(extract_pairs ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread)

add_executable(tieclean src/tieclean.cpp
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${gclib_SOURCE_DIR}/GStr.h
        ${gclib_SOURCE_DIR}/GStr.cpp
        ${gclib_SOURCE_DIR}/GArgs.h
        ${gclib_SOURCE_DIR}/GArgs.cpp
        ${gclib_SOURCE_DIR}/GBase.h
        ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(tieclean htslib zlib)
target_link_libraries(tieclean ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread)
