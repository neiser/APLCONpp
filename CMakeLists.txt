cmake_minimum_required (VERSION 3.0)
project(APLCONpp CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

include(ConfigSafeGuards)

if((CMAKE_CXX_COMPILER_ID MATCHES GNU) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    AND (CMAKE_Fortran_COMPILER_ID MATCHES GNU))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(DEBUG_FLAGS "-O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  # suppress warnings in release mode (just not to confuse people)
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w")
else()
  message(FATAL_ERROR "Only GNU/Clang CXX and GNU Fortran compiler supported")
endif()

add_subdirectory(APLCON)
add_subdirectory(src)
add_subdirectory(src/old)

option(EnableBenchmark "Build benchmark binaries" OFF)

if(EnableBenchmark)
  set(BENCHMARK_ENABLE_TESTING OFF)
  include(ConfigGBench)
  add_subdirectory(bench)
endif()


# enable testing at top-level to make it available to CTest
enable_testing()

# use some concurrency for tests
if(NOT CTEST_PARALLEL_JOBS)
  set(CTEST_PARALLEL_JOBS 2)
endif()

add_subdirectory(test)
