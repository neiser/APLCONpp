cmake_minimum_required (VERSION 3.0)
project(APLCONpp CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

include(ConfigSafeGuards)

if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

set(BENCHMARK_ENABLE_TESTING OFF)

if((CMAKE_CXX_COMPILER_ID MATCHES GNU) OR (CMAKE_CXX_COMPILER_ID MATCHES Clang))
   set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11")
   set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
   set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

include(ConfigGBench)

add_subdirectory(APLCON)
add_subdirectory(src)
add_subdirectory(src/old)
add_subdirectory(bench)

# enable testing here to make it available to CTest
enable_testing()
# use some concurrency for tests
if(NOT CTEST_PARALLEL_JOBS)
  set(CTEST_PARALLEL_JOBS 2)
endif()

add_subdirectory(test)
