#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=sodsf
e4shared --run --job=sodsf --verbosity=1 --max-cpus=2
e4shared --post --job=sodsf --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
