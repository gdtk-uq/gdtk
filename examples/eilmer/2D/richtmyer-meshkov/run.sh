#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=rmi
e4shared --run --job=rmi --verbosity=1 --max-cpus=4 > LOGFILE_RUN
e4shared --post --job=rmi --tindx-plot=all --vtk-xml --add-vars="mach"
