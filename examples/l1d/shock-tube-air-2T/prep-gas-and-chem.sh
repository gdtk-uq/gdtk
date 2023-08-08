#!/bin/bash
# prep-gas-and-chem.sh
cp ${DGD_REPO}/src/gas/sample-data/cea-air11species-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-11sp-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua .
sed -i '/^SUBSET_SELECTION =/s/5/11/' GuptaEtAl-air-reactions-2T.lua
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
sed -i '/^WITH_ELECTRONS = /s/false/true/' air-energy-exchange.lua
prep-gas ideal-helium.inp ideal-helium-gas-model.lua
prep-gas air-11sp-gas-model.lua air-11sp-2T.gas
prep-chem air-11sp-2T.gas GuptaEtAl-air-reactions-2T.lua air-11sp-2T.chem
prep-kinetics air-11sp-2T.gas air-11sp-2T.chem air-energy-exchange.lua air-11sp-2T.exch
