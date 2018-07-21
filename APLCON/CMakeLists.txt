# those flags (from original Makefile)
# are crucial to make the program run properly
# otherwise you end up in an endless loop...
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic -fno-backslash")

# build the library
set(SRCS_LIB
  aploop.f90.cpp
  chprob.f90.h
  condutil.f90.h
  helper.h
  )

add_library(aplcon OBJECT ${SRCS_LIB})
set_target_properties(aplcon PROPERTIES COMPILE_FLAGS "-fPIC")
