#!/bin/bash -e
gfortran -c condutil.F -save-temps
gfortran -c aploop.F -save-temps
gfortran -c chprob.F -save-temps
rm *.s *.o
mv condutil.f90 condutil.f90.F
mv aploop.f90 aploop.f90.F
mv chprob.f90 chprob.f90.F
f2c *.f90.F
rm *.f90.F
