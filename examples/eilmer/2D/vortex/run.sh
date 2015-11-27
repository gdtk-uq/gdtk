#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=vtx
e4shared --run --job=vtx --verbosity=1 --max-cpus=1
e4shared --post --job=vtx --vtk-xml --tindx-plot=all --add-vars="mach,pitot,total-p,total-h"
