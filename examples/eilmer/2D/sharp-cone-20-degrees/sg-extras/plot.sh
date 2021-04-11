#!/bin/bash
# plot.sh
# Compute coefficient of pressure for the history point
# and plot it against the previously computed high-res data.
#
# PJ, 2016-09-22
#
awk -f cp.awk hist/cone20-blk-1-cell-20.dat > cone20_cp.dat
gnuplot plot_cp.gnuplot

