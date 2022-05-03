#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=fp
e4-nk-shared --job=fp --verbosity=1
e4shared --post --job=fp --vtk-xml
python3 wht.py
