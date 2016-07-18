#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cyl50
e4shared --run --job=cyl50 --verbosity=1 --max-cpus=4
e4shared --post --job=cyl50 --tindx-plot=all --vtk-xml --add-vars="mach,pitot"
