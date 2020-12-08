#!/bin/bash
# prep-gas-and-chem.sh
cp ${DGD_REPO}/src/gas/sample-data/cea-air11species-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-11sp-2T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-2T.lua .
sed -i '/^NO_SPECIES =/s/5/11/' GuptaEtAl-air-2T.lua
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
prep-gas ideal-helium.inp ideal-helium-gas-model.lua
prep-gas air-11sp-2T.inp air-11sp-2T.lua
prep-chem air-11sp-2T.lua GuptaEtAl-air-2T.lua air-11sp-2T-chem.lua
