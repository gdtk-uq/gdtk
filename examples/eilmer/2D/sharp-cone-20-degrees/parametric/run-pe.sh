#!/bin/bash
# run-pe.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=conepe
e4shared --run --job=conepe --verbosity=1 --max-cpus=3
e4shared --post --job=conepe --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
