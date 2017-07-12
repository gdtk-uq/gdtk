#!/bin/bash
# make-temperature-plot.sh

# Extract the stagnation line data from the steady flow field.
e4shared --post --job=n90 --output-file=n90_100_iy1.data --tindx-plot=5 \
    --slice-list="0,:,1,0"
gnuplot plot_comparison.gnu


