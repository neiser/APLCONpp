#!/bin/bash -e
gfortran -c condutil.F -save-temps
gfortran -c aploop.F -save-temps
gfortran -c chprob.F -save-temps
rm *.s *.o
fable.cout aploop.f90 chprob.f90 condutil.f90 --namespace aplcon > aploop.cpp
rm *.f90
