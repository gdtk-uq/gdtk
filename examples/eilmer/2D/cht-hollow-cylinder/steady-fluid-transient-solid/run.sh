#!/bin/bash
# run.sh
var=`date +"%F"`

# run steady-state solver calculation
/usr/bin/time -f "RSS=%M elapsed=%E" mpirun -np 4 e4-cht-dist --job=cylinder --snapshot-start=0 --verbosity=1 | tee log_${var}.dat

#e4-cht-shared --job=cylinder --snapshot-start=0 --verbosity=1
