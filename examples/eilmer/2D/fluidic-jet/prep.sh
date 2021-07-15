#!/bin/bash
# prep.sh
# ugrid_partition fjl-gmsh.su2 mapped_cells 4 2 
# prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=fjl
e4shared --post --job=fjl --tindx-plot=last --vtk-xml
