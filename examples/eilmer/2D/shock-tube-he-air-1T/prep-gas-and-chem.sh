#!/bin/bash
# prep-gas-and-chem.sh
cp ${DGD_REPO}/src/gas/sample-data/cea-air5species-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
prep-gas ideal-helium.inp ideal-helium-gas-model.lua
prep-gas he-air-5sp-1T.inp he-air-5sp-1T-gas-model.lua
prep-chem he-air-5sp-1T-gas-model.lua GuptaEtAl-air-reactions.lua he-air-5sp-1T-chem.lua
