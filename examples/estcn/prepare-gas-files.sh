#! /bin/bash
# Prepare gas files for the estcn demonstrations.
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
cp ${DGD_REPO}/src/gas/sample-data/cea-air5species-gas-model.lua .
cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
cp ${DGD_REPO}/src/gas/sample-data/ideal-air-gas-model.lua .
cp ${DGD_REPO}/src/gas/sample-data/air-5sp-eq.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .
prep-gas he-n2-o2-mix.inp he-n2-o2-mix-gas-model.lua
prep-gas air-5sp-1T.inp air-5sp-1T.lua
