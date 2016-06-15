#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=coles
e4shared --run --job=coles --verbosity=1 --max-cpus=1
e4shared --post --job=coles --vtk-xml --add-vars="mach,pitot,total-p,total-h"
