#!/bin/bash
# prep.sh
ugrid_partition icon-gmsh.su2 mapped_cells 4 2 
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=icon
e4shared --post --job=icon --tindx-plot=last --vtk-xml
