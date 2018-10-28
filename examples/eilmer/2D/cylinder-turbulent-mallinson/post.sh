#!/bin/bash

# Extracts Paraview files for flowfield visualisation
e4shared --post --job=mallinson --vtk-xml --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mallinson --output-file="mallinson-y-wall.dat" \
  --slice-list=":,:,0,0" --add-vars="mach,pitot,total-p,total-h"
