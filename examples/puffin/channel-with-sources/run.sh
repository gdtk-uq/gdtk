#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=source_test
puffin --job=source_test
puffin-post --job=source_test --output=vtk
#puffin-post --job=case1 --output=stream --cell-index=$ --stream-index=0
