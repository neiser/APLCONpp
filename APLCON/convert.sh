#!/bin/bash -e
gfortran -c condutil.F -save-temps
gfortran -c aploop.F -save-temps
gfortran -c chprob.F -save-temps
gfortran -c aplusopt.F -save-temps
gfortran -c aplist.F -save-temps
rm *.s *.o
mv condutil.f90 condutil.f90.F
mv aploop.f90 aploop.f90.F
mv chprob.f90 chprob.f90.F
mv aplusopt.f90 aplusopt.f90.F
mv aplist.f90 aplist.f90.F
f2c *.f90.F
rm *.f90.F
clang-format -i *.f90.c

