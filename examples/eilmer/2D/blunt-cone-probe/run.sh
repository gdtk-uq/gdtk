#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=bcp
e4shared --run --job=bcp --verbosity=1 --max-cpus=6
e4shared --post --job=bcp --vtk-xml --tindx-plot=all --add-vars="mach,pitot,total-p,total-h"
