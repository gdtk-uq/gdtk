#!/bin/bash
# run.sh
DGD_REPO=${DGD_REPO:=${HOME}/dgd}
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua cea-lut-air.lua
e4shared --prep --job=sphere
e4shared --run --job=sphere --verbosity=1 --max-cpus=4
e4shared --post --job=sphere --tindx-plot=all --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
