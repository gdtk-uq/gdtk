#!/bin/bash
# plot.sh
# Compute coefficient of pressure for the history point
# and plot it against the previously computed high-res data.
#
# PJ, 2016-09-22, 2018-01-28 updated for 6-block MPI simulation
#
awk -f cp.awk hist/cone20-blk-4-cell-4.dat > cone20_cp.dat
gnuplot plot_cp.gnuplot

