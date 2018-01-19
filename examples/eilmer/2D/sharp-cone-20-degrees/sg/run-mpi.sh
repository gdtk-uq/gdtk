#!/bin/bash
# run-mpi.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cone20
mpirun -np 2 e4mpi --run --job=cone20 --verbosity=1 --max-cpus=1
# e4shared --post --job=cone20 --vtk-xml --add-vars="mach,pitot,total-p,total-h"
