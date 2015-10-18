#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=swbli
e4shared --run --job=swbli --verbosity=1 --max-cpus=4
e4shared --post --job=swbli --vtk-xml --add-vars="mach,pitot,total-p,total-h"
