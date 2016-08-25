#!/bin/bash
# run-time.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=lam_flat_plate_time
e4shared --run --job=lam_flat_plate_time --verbosity=1 --max-cpus=4
e4shared --post --job=lam_flat_plate_time --tindx-plot=all --vtk-xml --add-vars="mach,pitot"
