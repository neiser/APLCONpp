find_package(Threads REQUIRED)

add_executable(main main.cc wrapper/APLCON.f90 wrapper/APLCON.h)
target_link_libraries(main
  aplcon
  benchmark ${CMAKE_THREAD_LIBS_INIT}
)
