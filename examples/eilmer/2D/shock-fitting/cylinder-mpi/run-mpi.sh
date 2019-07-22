#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cyl-sf-mpi
mpirun -np 8 e4mpi --run --job=cyl-sf-mpi --verbosity=1
e4shared --post --job=cyl-sf --tindx-plot=all --vtk-xml
