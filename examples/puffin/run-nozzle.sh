#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=nozzle
puffin --job=nozzle
puffin-post --job=nozzle --output=vtk

