#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cyl
e4shared --run --job=cyl --verbosity=1 --max-cpus=6 | tee run.log
e4shared --post --job=cyl --vtk-xml --tindx-plot=all \
	 --add-vars="mach,pitot,total-p,total-h"
