#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --custom-script --script-file=makegrid.lua
ugrid_partition grid.su2 mapped_cells 4 2
mkdir -p su2grid
mv block_* su2grid
e4shared --prep --job=mabey
mpirun -np 4 e4-nk-dist --job=mabey --verbosity=1
#e4-nk-shared --job=mabey --verbosity=1 --max-cpus=1
e4shared --post --job=mabey --vtk-xml --add-vars="mach,pitot,total-p,total-h"
