#!/bin/bash
# run-mpi.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cone20
mpirun -np 8 --oversubscribe e4mpi --run --job=cone20 --verbosity=1
e4shared --post --job=cone20 --vtk-xml --add-vars="mach,pitot,total-p,total-h"
