#!/bin/bash
# run.sh
ugrid_partition cone20.su2 mapped_cells 4 2 
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cone20
e4shared --post --job=cone20 --tindx-plot=last --vtk-xml
mpirun -np 4 e4mpi --run --job=cone20 --verbosity=1 --threads-per-mpi-task=1
e4shared --post --job=cone20 --add-vars="mach,pitot,total-p,total-h" --vtk-xml \
         --extract-streamline="0.1,0.11,0.0" --output-file="streamline.dat"
