#!/bin/bash
# plot.sh
# Compute coefficient of pressure for the history point
# and plot it against the previously computed high-res data.
#
# PJ, 2016-09-22, 2018-01-28 updated for 6-block MPI simulation
#     2021-10-19, updated for 8-block MPI simulation and new history names
#
awk -f cp.awk hist/cone20-blk-4-cell-9.dat.0 > cone20_cp.dat
gnuplot plot_cp.gnuplot

