#! /bin/bash
# Prepare gas files for the air driving equilibrium-air shock tube.
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua .
prep-gas ideal-air.inp ideal-air-gas-model.lua
