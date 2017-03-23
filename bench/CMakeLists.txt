find_package(Threads REQUIRED)
include_directories(../src)

add_executable(bench_APLCON bench_APLCON.cc)

target_link_libraries(bench_APLCON aplcon++
    benchmark ${CMAKE_THREAD_LIBS_INIT})
