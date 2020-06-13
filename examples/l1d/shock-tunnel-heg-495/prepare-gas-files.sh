#! /bin/bash
# Prepare gas files for the HEG shock tube simulation.
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
prep-gas ideal-air.inp ideal-air-gas-model.lua
prep-gas thermally-perfect-He-Ar.inp he-ar-gas-model.lua
