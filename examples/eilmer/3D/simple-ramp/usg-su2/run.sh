#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=ramp
e4shared --run --job=ramp --verbosity=1
e4shared --post --job=ramp --tindx-plot=all --vtk-xml --add-vars="mach,pitot"
