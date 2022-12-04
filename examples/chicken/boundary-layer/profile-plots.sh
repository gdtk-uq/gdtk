#!/bin/bash
# profile-plots.sh
# PJ 2022-11-24, 2022-12-04 Short plate
#
awk -f integral-thicknesses.awk slice-t0005.data > thicknesses.txt

gnuplot<<EOF
set term postscript eps enhanced 20
set output "velx-velocity-comparison.eps"
set title "Short plate velocity near wall at x = 1.0 m"
set xlabel "vx/vx_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./slice-t0005.data" using ((\$11)/1390.0):((\$2)*1000) \
     title "Chicken" with points pt 1, \
     "boundary-layer-profile-clbl.data" using (\$3):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF


gnuplot<<EOF
set term postscript eps enhanced 20
set output "density-comparison.eps"
set title "Short plate density near wall at x = 1.0 m"
set xlabel "rho/rho_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./slice-t0005.data" using ((\$7)/0.01176):((\$2)*1000) \
     title "Chicken" with points pt 1, \
     "boundary-layer-profile-clbl.data" using (\$5):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF

