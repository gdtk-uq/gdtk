#!/bin/bash

# 1. residual extraction and plotting
gnuplot plot-residuals.gplot

# 2. shock location and comparison to experiment
e4shared --custom-script --script-file=shock-shape.lua
gnuplot plot-shock-shape.gplot
