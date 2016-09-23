#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=conep
e4shared --run --job=conep --verbosity=1 --max-cpus=2
e4shared --post --job=conep --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
