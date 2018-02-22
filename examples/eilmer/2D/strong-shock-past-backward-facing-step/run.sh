#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=corner
e4shared --run --job=corner --verbosity=1 --max-cpus=2
e4shared --post --job=corner --vtk-xml --tindx-plot=all \
	 --add-vars="mach,pitot,total-p,total-h"
