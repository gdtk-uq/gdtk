#! /bin/bash
# prepare-gas-files-air.sh
#
# For equilibrium thermo.
cp ${DGD_REPO}/src/gas/sample-data/cea-air5species-gas-model.lua .
cp ${DGD_REPO}/src/gas/sample-data/air-5sp-eq.lua .
#
# For 1-temperature thermo with finite-rate chemistry.
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
prep-gas air-5sp-1T.inp air-5sp-1T.lua
prep-chem air-5sp-1T.lua GuptaEtAl-air-reactions.lua air-5sp-1T-reactions.lua
#
# For 2-temperature thermo with finite-rate chemistry.
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-5sp-2T.inp .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/GuptaEtAl-air-2T.lua .
cp ${DGD_REPO}/examples/kinetics/air-chemistry-2T/air-energy-exchange.lua .
prep-gas air-5sp-2T.inp air-5sp-2T.lua
prep-chem air-5sp-2T.lua GuptaEtAl-air-2T.lua air-5sp-6r-2T.lua
