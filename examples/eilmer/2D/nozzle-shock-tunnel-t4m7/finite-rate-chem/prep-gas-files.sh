#! /bin/bash
# prep-gas-files.sh
cp ${DGD_REPO}/src/gas/sample-data/cea-air5species-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
prep-gas air-5sp-1T.inp 5sp-gas-model.lua
prep-chem 5sp-gas-model.lua GuptaEtAl-air-reactions.lua 5sp-chem-model.lua
