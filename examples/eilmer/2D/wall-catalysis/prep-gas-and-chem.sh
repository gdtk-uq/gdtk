#!/bin/bash
# prep-gas-and-chem.sh
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-2T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-2T.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
prep-gas air-5sp-2T.inp air-5sp-2T.lua
prep-chem air-5sp-2T.lua GuptaEtAl-air-2T.lua air-5sp-2T-chem.lua
