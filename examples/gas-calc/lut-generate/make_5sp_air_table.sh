#!/bin/bash

# Exercise the build-uniform-lut script using a 5 species equilibrium air model
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

cp $DGD_REPO/src/gas/sample-data/air-5sp-1T-input.lua ./
cp $DGD_REPO/src/gas/sample-data/air-5sp-eq.lua ./

lmr prep-gas -i air-5sp-1T-input.lua -o air-5sp-1T.lua
build-uniform-lut --gas-model=air-5sp-eq.lua --table-name=air5
