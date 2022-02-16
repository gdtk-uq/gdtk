#!/bin/bash
# prep.sh

# Start by copying the gas model files from the source-code repository.
DGD_REPO=${DGD_REPO:=${HOME}/gdtk}
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .

# Do the preparation of grid and initial flow state, only.
# Note that we use the debug flavour of the code for the CEA gas model
# that gets used in the calculation of the throat conditions.
e4shared --prep --job=t4m4b

# Postprocessing uses only the LUT gas model.
e4shared --post --job=t4m4b --tindx-plot=0 --vtk-xml
