#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=os
e4shared --run --job=os --verbosity=1 --max-cpus=1
e4shared --post --job=os --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
