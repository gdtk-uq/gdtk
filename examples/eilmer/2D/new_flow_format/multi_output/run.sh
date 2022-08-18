#!/bin/bash
# run.sh

prep-gas ideal-air.inp ideal-air-gas-model.lua

e4shared --prep --job=ffs
e4shared --run --job=ffs --verbosity=1 --max-cpus=3
e4shared --post --job=ffs --tindx-plot=all --vtk-xml

e4shared --post --job=ffs --tindx-plot="1,2,3,4,5" --vtk-xml --plotTag="average"
e4shared --post --job=ffs --tindx-plot=last --vtk-xml --plotTag="DFT"
e4shared --post --job=ffs --tindx-plot="1,2,3,4,5" --vtk-xml --plotTag="viscous_gradient"

