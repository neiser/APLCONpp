#!/bin/bash
LIST="aploop"
for i in $LIST; do
    gcc -E -traditional -x c $i.F | f2c > $i.F.c
done
