#! /bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
l1d4-prep --job=piston-with-valve
l1d4 --run-simulation --job=piston-with-valve
l1d4 --piston-history --job=piston-with-valve --pindx=0
