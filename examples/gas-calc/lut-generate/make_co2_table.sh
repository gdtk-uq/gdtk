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

cp $DGD_REPO/src/gas/sample-data/co2-5sp-1T-input.lua ./
cp $DGD_REPO/src/gas/sample-data/air-5sp-eq.lua ./co2-5sp-eq.lua

lmr prep-gas -i co2-5sp-1T-input.lua -o co2-5sp-1T.lua
sed -i 's/air/co2/g' co2-5sp-eq.lua
sed -i 's/N2=0.79/CO2=1.0/g' co2-5sp-eq.lua
sed -i 's/O2=0.21//g' co2-5sp-eq.lua
build-uniform-lut --gas-model=co2-5sp-eq.lua --table-name=co2
