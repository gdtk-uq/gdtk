#!/bin/sh
#  back_run.sh
#  Exercise the Navier-Stokes solver for the conical nozzle 
#  as used by Back, Massier and Gier (1965) AIAA J. 3(9):1606-1614.
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=back
e4shared --run --job=back --verbosity=1 --max-cpus=2
e4shared --post --job=back --tindx-plot=all --vtk-xml --add-vars="mach,pitot,total-p,total-h"
