#!/bin/bash
# prep.sh

# prep gas model file
prep-gas gm-ideal-air.inp ideal-air.gas

# prep flow solver files
e4shared --prep --job=cylinder

# generate a plot of the inflow conditions
e4shared --post --job=cylinder --vtk-xml --tindx-plot=0
