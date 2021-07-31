#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cone
e4shared --run --job=cone --verbosity=1 --max-cpus=2
e4shared --post --job=cone --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
