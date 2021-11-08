#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=sphere_cone
e4shared --post --job=sphere_cone --tindx-plot=all --vtk-xml
e4shared --run --job=sphere_cone --verbosity=1 --max-cpus=4
e4shared --post --job=sphere_cone --tindx-plot=all --vtk-xml
