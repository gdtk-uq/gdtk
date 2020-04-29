#!/bin/bash
# plot-piston-data.sh

# Generate reference data
pushd ./analytic; ruby piston.rb > piston.data
popd

# Plot the trajectory
gnuplot <<EOF
set term postscript eps enhanced 20
set output "velocity-time.eps"
set style line 1 linetype 1 linewidth 3.0
set title "Piston trajectory, L1d4 simulation"
set xlabel "time, ms"
set ylabel "velocity, m/s"
set xtic 5.0
set ytic 50.0
set yrange [0:300]
set key bottom right
plot "piston-0000-history.data" using (\$2*1000):4 title "L1d4", \
     "./analytic/piston.data" using (\$1*1000):3 title "analytic" with lines
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "velocity-distance.eps"
set style line 1 linetype 1 linewidth 3.0
set title "Piston trajectory, L1d4 simulation"
set xlabel "position, m"
set ylabel "velocity, m/s"
set xtic 1.0
set ytic 50.0
set xrange [0:8]
set yrange [0:300]
set key bottom right
plot "piston-0000-history.data" using (\$3):4 title "L1d4", \
     "./analytic/piston.data" using (\$2):3 title "analytic" with lines
EOF
