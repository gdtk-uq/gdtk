#!/bin/bash
# prep-gas-and-chem.sh
cp ${DGD_REPO}/src/gas/sample-data/cea-air11species-gas-model.lua .
prep-gas ideal-helium.inp ideal-helium-gas-model.lua
prep-gas air-11sp-2T.inp air-11sp-2T.lua
prep-chem air-11sp-2T.lua GuptaEtAl-air-2T.lua air-11sp-2T-chem.lua
