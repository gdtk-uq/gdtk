#!/bin/bash
# plot_heat_transfer.sh

# Get the flow data around the surface of the sphere.
# This runs up the east boundary face of all blocks.
e4shared --post --job=sphere --tindx-plot=last --slice-list=":,$,:,:" \
         --output-file=sphere_surface_conditions.dat

awk -f compute_heat_transfer.awk sphere_surface_conditions.dat > sphere_heat_transfer.dat

gnuplot <<EOF
set term postscript eps enhanced 20
set output "sphere_norm_heat_transfer.eps"
set style line 1 linetype 1 linewidth 3.0
set xlabel "angle from stagnation point, degrees"
set ylabel "q/q_s"
set logscale y
# set yrange [0.1:2.0]
set yrange [0.1:1.2]
set title "Normalised heat transfer to R6.6mm sphere with Ms=8"
# set key top right
set key bottom left
plot "sphere_heat_transfer.dat" using 1:3 \
          title "Eilmer4 simulation" with lines ls 1, \
     "kemp_theory.dat" using 1:2 title "Kemp-Rose-Detra theory" \
          with linespoints, \
     "kemp_experiment.dat" using 1:2 title "Kemp-Rose-Detra experiment" \
          with points
EOF
