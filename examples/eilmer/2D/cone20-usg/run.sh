#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cone20
e4shared --run --job=cone20 --verbosity=1 --max-cpus=2
e4shared --post --job=cone20 --extract-streamline="0.1,0.11,0.0" --output-file="streamline.dat" --vtk-xml --add-vars="mach,pitot,total-p,total-h"
