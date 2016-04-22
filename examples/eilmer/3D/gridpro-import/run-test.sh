#!/bin/bash

prep-gas ideal-air.inp ideal-air.lua
e4shared --job=test-import --prep
e4shared --job=test-import --post --vtk-xml

echo "Grid may be viewed in paraview. Expected output is 3 blocks."
echo "> paraview plot/test-import.pvd"

