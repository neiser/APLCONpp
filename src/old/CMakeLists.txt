INCLUDE_DIRECTORIES(. ..)
add_library(old_aplcon++
  APLCON.cc
  APLCON.hpp
  detail/APLCON_hpp.hpp
  detail/APLCON_cc.hpp
  detail/APLCON_ostream.hpp
  ../wrapper/APLCON.f90
  ../wrapper/APLCON.h
  )
target_link_libraries(old_aplcon++ aplcon)
