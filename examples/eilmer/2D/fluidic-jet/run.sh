#!/bin/bash
# run.sh
# e4shared --run --job=fj --verbosity=1 --max-cpus=4
mpirun -np 4 e4mpi --run --job=fjl --verbosity=1
e4shared --post --job=fjl --tindx-plot=all --vtk-xml \
         --add-vars="mach,pitot,total-p,total-h"

