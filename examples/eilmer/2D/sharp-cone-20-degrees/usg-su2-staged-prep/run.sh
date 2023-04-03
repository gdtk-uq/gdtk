#! /usr/bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep-grid --job=cone20
e4shared --prep-flow --job=cone20
e4shared --run --job=cone20
e4shared --post --job=cone20 --vtk-xml
