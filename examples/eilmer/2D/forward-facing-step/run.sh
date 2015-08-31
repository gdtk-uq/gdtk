#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=ffs
e4shared --run --job=ffs --verbosity=1 --max-cpus=3
e4shared --post --job=ffs --tindx-plot=all --vtk-xml --add-vars="mach"
