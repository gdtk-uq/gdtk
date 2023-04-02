#! /usr/bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep-grid --job=cone20
e4shared --prep-flow --job=cone20
e4shared --job=cone20 --run
