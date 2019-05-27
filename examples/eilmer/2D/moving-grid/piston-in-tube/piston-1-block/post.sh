#!/bin/bash

# We'll ask for all of the flow field snapshots so
# that we can make an animation in Paraview.
e4shared --post --job=pit1 --vtk-xml --tindx-plot=all

# Run the custom post script to check our energy and mass balances
e4shared --custom-script --script-file="balanceCheck.lua"

# Plot the trajectory
gnuplot <<EOF
set term postscript eps enhanced 20
set output "velocity-time.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "Piston trajectory, 1 FluidBlock"
set xlabel "time, ms"
set ylabel "velocity, m/s"
set xtic 5.0
set ytic 50.0
set yrange [0:300]
set key bottom right
plot "piston.data" using (\$1*1000):3 title "Eilmer", \
     "../analytic/piston.data" using (\$1*1000):3 title "analytic" with lines
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "velocity-distance.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "Piston trajectory, 1 FluidBlock"
set xlabel "position, m"
set ylabel "velocity, m/s"
set xtic 1.0
set ytic 50.0
set xrange [0:8]
set yrange [0:300]
set key bottom right
plot "piston.data" using (\$2):3 title "Eilmer", \
     "../analytic/piston.data" using (\$2):3 title "analytic" with lines
EOF
