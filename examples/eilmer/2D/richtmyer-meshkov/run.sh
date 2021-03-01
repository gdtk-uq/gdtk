#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=rmi
mpirun -np 6 e4mpi --run --job=rmi --verbosity=1
e4shared --post --job=rmi --tindx-plot=all --vtk-xml --add-vars="mach"
