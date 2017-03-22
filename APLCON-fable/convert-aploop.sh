#!/bin/sh
gfortran -w -c aploop.F -save-temps -o /dev/null
rm aploop.s
fable.cout aploop.f90 > aploop.cpp
rm aploop.f90
