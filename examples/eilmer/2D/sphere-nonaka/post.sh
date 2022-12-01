#!/bin/bash

# 0. Make normal vtk files for visual inspection
e4shared --post --job=nonaka --vtk-xml --tindx-plot=last

# 1. residual extraction and plotting
gnuplot plot-residuals.gplot

# 2. shock location and comparison to experiment
e4shared --custom-script --script-file=shock-shape.lua
gnuplot plot-shock-shape.gplot
