#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=mabey
e4zsss --job=mabey --verbosity=1
e4shared --post --job=mabey --vtk-xml --add-vars="mach,pitot,total-p,total-h"
