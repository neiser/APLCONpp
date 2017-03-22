add_library(aplcon++
    APLCON.hpp detail/APLCON_hpp.hpp
    wrapper/APLCON.f90 wrapper/APLCON.h)

add_library(aplcon++_interface INTERFACE)
target_include_directories(aplcon++_interface INTERFACE .)

target_link_libraries(aplcon++
    aplcon aplcon++_interface
)

add_executable(main main.cc)
target_link_libraries(main aplcon++)