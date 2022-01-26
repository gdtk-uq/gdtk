#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=sample
puffin --job=sample
puffin-post --job=sample --output=vtk

