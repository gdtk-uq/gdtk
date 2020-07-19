#! /bin/bash
# Prepare gas files for the T4 shock tube simulation, argon driver.
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
prep-gas ideal-air.inp ideal-air-gas-model.lua
prep-gas ideal-argon.inp ideal-argon-gas-model.lua
