if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  message(FATAL_ERROR "Only GNU Fortran compiler supported")
endif()

# those flags (from original Makefile)
# are crucial to make the program run properly
# otherwise you end up in an endless loop...
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic -fno-backslash")

# build the library
set(SRCS_LIB
  a12prof.F
  aplist.F
  aploop.F
  aplprint.F
  aplusopt.F
  chprob.F
  condutil.F
  dcsplines.F
  fpl.F
  )

add_library(aplcon OBJECT ${SRCS_LIB})
