#! /bin/bash
# Prepare gas files for the helium driving equilibrium-air shock tube.
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua .
prep-gas ideal-helium.inp ideal-helium-gas-model.lua
