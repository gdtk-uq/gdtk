#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=salas-nozzle
puffin --job=salas-nozzle
puffin-post --job=salas-nozzle --output=vtk
puffin-post --job=salas-nozzle --output=stream --cell-index=$ --stream-index=0
