#!/bin/bash
# run-with-staged-prep.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua

e4shared --prep-grids --job=ffs
# e4shared --prep-flow --job=ffs
e4-prep-parallel --job=ffs --n-workers=7 --nb-per-task=3

# e4shared --run --job=ffs --verbosity=1 --max-cpus=3
mpirun -np 3 e4mpi --run --job=ffs --verbosity=1
e4shared --post --job=ffs --tindx-plot=all --vtk-xml --add-vars="mach"
