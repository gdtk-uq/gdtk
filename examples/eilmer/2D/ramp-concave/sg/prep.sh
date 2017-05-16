#!/bin/bash
# prep.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cubic-ramp
e4shared --post --job=cubic-ramp --vtk-xml --tindx-plot=0
