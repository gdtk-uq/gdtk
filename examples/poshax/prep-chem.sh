#!/bin/bash
# prep-chem.sh
prep-gas air-7sp.inp air-7sp-gas-model.lua
prep-chem air-7sp-gas-model.lua GuptaEtAl-air-reactions.lua air-7sp-chemistry.lua
