#!/bin/bash
# prep-gas-and-chem.sh
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-reactions-2T.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
sed -i 's/11-species-air/5-species-air/' GuptaEtAl-air-reactions-2T.lua
prep-gas air-5sp-gas-model.lua gm.lua
prep-chem gm.lua GuptaEtAl-air-reactions-2T.lua rr.lua
prep-kinetics gm.lua air-energy-exchange.lua ee.lua
