#!/bin/bash
# prep-ram2.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=ram2
e4shared --post --job=ram2 --tindx-plot=0 --vtk-xml

echo "At this point, we should have a starting grid"
echo "Run run-ram2.sh next"

