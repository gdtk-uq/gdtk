#! /bin/bash
prep-gas ideal-air.inp ideal-air-gas-model.lua
puffin-prep --job=salas-scramjet
puffin --job=salas-scramjet
puffin-post --job=salas-scramjet --output=vtk

