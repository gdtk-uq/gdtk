#!/bin/bash

prep-gas ideal-air.inp ideal-air.lua
e4shared --job=test-import --prep
mpirun -np 5 e4mpi --job=test-import --run
e4shared --job=test-import --post --vtk-xml --tindx-plot=all

echo "Grid may be viewed in paraview. Expected output is 5 blocks."
echo "> paraview plot/test-import.pvd"
