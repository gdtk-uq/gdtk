#!/bin/bash
# run.sh
prep-gas kiock-gas.inp kiock-gas-model.lua
e4shared --prep --job=kiock
time e4shared --run --job=kiock --verbosity=0 --max-cpus=4
e4shared --post --job=kiock --vtk-xml --tindx-plot=all --add-vars="mach,pitot,total-p,total-h"
paraview plot/kiock.pvd
