#! /bin/bash
# prepare-gas-files.sh
cp ${DGD_REPO}/src/gas/sample-data/cea-air5species-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
prep-gas air-5sp-1T.inp air-5sp-1T.lua
prep-chem air-5sp-1T.lua GuptaEtAl-air-reactions.lua air-5sp-1T-reactions.lua
