#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=nozzle
puffin --job=nozzle
puffin-post --job=nozzle --output=vtk
puffin-post --job=nozzle --output=stream --cell-index=$ --stream-index=0
gnuplot plot-p.gnuplot
