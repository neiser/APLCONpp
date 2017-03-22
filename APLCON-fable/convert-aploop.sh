#!/bin/bash
function make_f90 {
    gfortran -w -c $1.F -save-temps -o /dev/null
    rm $1.s
}
LIST="condutil.f90 aploop.f90 chprob.f90"
for i in $LIST; do
    make_f90 ${i%.*}
done
fable.cout --namespace aplcon_fortran $LIST > APLCON-fable.cpp
rm $LIST
