#!/bin/bash
# prep.sh

# Start by copying the gas model files from the source-code repository.
DGD_REPO=${DGD_REPO:=${HOME}/dgd}
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .

# Do the preparation of grid and initial flow state, only.
e4shared --prep --job=t4m4
e4shared --post --job=t4m4 --tindx-plot=0 --vtk-xml
