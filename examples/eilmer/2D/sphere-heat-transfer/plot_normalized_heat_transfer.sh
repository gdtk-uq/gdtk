#!/bin/bash
# plot_normalized_heat_transfer.sh
# Plot the heat transfer that has been extracted from the loads files.
# PJ 2021-08-01

gnuplot <<EOF
set term postscript eps enhanced 20
set output "normalized_heat_transfer.eps"
set style line 1 linetype 1 linewidth 3.0
set xlabel "angle from stagnation point, degrees"
set ylabel "q/q_s"
set logscale y
# set yrange [0.1:2.0]
set yrange [0.1:1.2]
set title "Normalized heat transfer to R6.6mm sphere with Ms=8"
# set key top right
set key bottom left
set label "source-flow setting" at 60,1.02
set label "r=0.1 x0=-0.1066m" at 60,0.9
plot "normalized_heat_transfer.data" using 1:3 \
          title "Eilmer4 simulation" with lines ls 1, \
     "kemp_theory.dat" using 1:2 title "Kemp-Rose-Detra theory" \
          with linespoints, \
     "kemp_experiment.dat" using 1:2 title "Kemp-Rose-Detra experiment" \
          with points
EOF
