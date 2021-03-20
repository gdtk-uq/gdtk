#!/bin/bash
# prep-chem.sh
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-7sp-1T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
sed -i '/^SUBSET_SELECTION =/s/5/7/' GuptaEtAl-air-reactions.lua
prep-gas air-7sp-1T.inp air-7sp-gas-model.lua
prep-chem air-7sp-gas-model.lua GuptaEtAl-air-reactions.lua air-7sp-chemistry.lua
