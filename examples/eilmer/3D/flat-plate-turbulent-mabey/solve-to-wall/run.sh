#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=mabey_3D
e4shared --run --job=mabey_3D --verbosity=1 --max-cpus=8
