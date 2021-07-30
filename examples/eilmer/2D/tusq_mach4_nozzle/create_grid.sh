#!/bin/bash
# create_grid.sh
python3 gmsh_Mach4_nozzle.py
ugrid_partition mach4nozzleStage1.su2 mapped_cells 100 2
rm -r GRID
mkdir GRID
mv block* GRID/.