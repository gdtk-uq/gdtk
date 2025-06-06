#!/bin/bash

# Exercise the build-uniform-lut script using an 13 species equilibrium air model
#
# Notes:
# - This script requires some extra build steps:
#    $ cd ~/gdtk/src/gas
#    $ make libgas.so install
#
# - Then ensure that the following is added to .bashrc
#    export PYTHONPATH=${PYTHONPATH}:${DGD}/lib
#    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DGD}/lib
#
# @author: Nick Gibbons (n.gibbons@uq.edu.au)

cp $DGD_REPO/src/gas/sample-data/air-13sp-1T-input.lua ./
cp $DGD_REPO/src/gas/sample-data/air-13sp-eq-gas-model.lua ./

lmr prep-gas -i air-13sp-1T-input.lua -o air-13sp-gas-model.lua
build-uniform-lut --gas-model=air-13sp-eq-gas-model.lua --table-name=air13
