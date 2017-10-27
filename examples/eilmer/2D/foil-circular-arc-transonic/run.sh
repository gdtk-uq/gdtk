#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=arc
e4shared --run --job=arc --max-cpus=4 --verbosity=1
e4shared --post --job=arc --vtk-xml --add-vars="mach,total-p" --tindx-plot=all
