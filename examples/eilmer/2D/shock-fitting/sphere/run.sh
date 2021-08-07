#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=sphere
e4shared --run --job=sphere --verbosity=1 --max-cpus=4
e4shared --post --job=sphere --tindx-plot=all --vtk-xml
