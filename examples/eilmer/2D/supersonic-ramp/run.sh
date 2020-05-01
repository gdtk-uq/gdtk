#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=ramp
e4-nk-shared --job=ramp --verbosity=1
e4shared --post --job=ramp --vtk-xml --add-vars="mach,pitot,total-p,total-h"


