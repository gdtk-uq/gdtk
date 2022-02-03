#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=cone10
puffin --job=cone10
puffin-post --job=cone10 --output=vtk

