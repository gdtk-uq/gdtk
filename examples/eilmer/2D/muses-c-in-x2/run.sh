#!/bin/bash
# run.sh
prep-gas ideal-air-13.inp ideal-air-13-gas-model.lua
e4shared --prep --job=cap
time e4shared --run --job=cap --verbosity=1 --max-cpus=3
e4shared --post --job=cap --vtk-xml --add-vars="mach,pitot,total-p,total-h" \
         --tindx-plot=all
