find_package(Threads REQUIRED)

add_executable(bench_APLCON bench_APLCON.cc)
add_executable(bench_APLCON_old bench_APLCON_old.cc)

target_link_libraries(bench_APLCON aplcon++
    benchmark ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bench_APLCON_old old_aplcon++
    benchmark ${CMAKE_THREAD_LIBS_INIT})
