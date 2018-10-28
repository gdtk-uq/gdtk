#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=eklund
e4shared --run --job=eklund --verbosity=1 --max-cpus=8
