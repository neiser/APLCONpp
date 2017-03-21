find_package(Threads REQUIRED)

add_executable(bench_APLCON bench_APLCON.cc)

target_link_libraries(bench_APLCON aplcon++
    benchmark ${CMAKE_THREAD_LIBS_INIT})