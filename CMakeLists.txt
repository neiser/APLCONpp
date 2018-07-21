cmake_minimum_required (VERSION 3.0)
project(APLCONpp CXX C)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

# set also default CMAKE_BUILD_TYPE to Release
include(ConfigSafeGuards)

if((CMAKE_CXX_COMPILER_ID MATCHES GNU) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(DEBUG_FLAGS "-O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FLAGS}")
else()
  message(FATAL_ERROR "Only GNU/Clang CXX and GNU Fortran compiler supported")
endif()

# figure out compile flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(DEFAULT_CXX_COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${BUILD_TYPE}})
set(DEFAULT_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}})

message(STATUS "*** Build Type:                " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags (CXX):      " ${DEFAULT_CXX_COMPILE_FLAGS})
message(STATUS "*** Compiler Flags (Fortran):  " ${DEFAULT_Fortran_COMPILE_FLAGS})

add_subdirectory(APLCON)

add_library(aplcon++ SHARED
    src/APLCON.hpp src/detail/APLCON_hpp.hpp
    src/wrapper/APLCON.cpp src/wrapper/APLCON.h
    $<TARGET_OBJECTS:aplcon>)

# Benchmarking business
option(EnableBenchmark "Build benchmark binaries" OFF)
if(EnableBenchmark)
  set(BENCHMARK_ENABLE_TESTING OFF)
  include(ConfigGBench)
  add_subdirectory(bench)
endif()

if(COVERAGE)
  message(STATUS "Coverage build, enforce Debug build type")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
endif()


# enable testing at top-level to make it available to CTest
enable_testing()

# use some concurrency for tests
if(NOT CTEST_PARALLEL_JOBS)
  set(CTEST_PARALLEL_JOBS 2)
endif()

add_subdirectory(test)
