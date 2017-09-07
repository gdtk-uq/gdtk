#!/bin/sh
#  run.sh
#  Run the Navier-Stokes solver for the nozzle simulation. 
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=nozzle
e4shared --run --job=nozzle --verbosity=1
e4shared --post --job=nozzle --tindx-plot=all --vtk-xml --add-vars="mach,pitot,total-p,total-h"
