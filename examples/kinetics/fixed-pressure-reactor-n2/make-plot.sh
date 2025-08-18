#!/bin/bash
# make-plot.sh
# Display history for the fixed-pressure-reactor simulation.
#
# PJ & RG, 2018-04-21, 2019-11-26, 2025-08-16

# $ python3 fpreactor.py

gnuplot <<EOF
set terminal postscript eps 20
set output "reactor-4000K-1bar-T.eps"
#
set title "Nitrogen in a fixed-pressure reactor."
set mxtics 5
set xtics 0,100,600
set xlabel 't, microseconds'
#
set mytics 5
set ytics 4000,500,7000
set yrange[4000:7000]
set ylabel 'T, K'
#
set nokey
# set key bottom right
#
plot 'fpreactor.data' using (\$1*1.0e6):2 with lines linetype 1
EOF


gnuplot <<EOF
set terminal postscript eps 20
set output "reactor-4000K-1bar-rho.eps"
#
set title "Nitrogen in a fixed-pressure reactor."
set mxtics 5
set xtics 0,100,600
set xlabel 't, microseconds'
#
set mytics 5
set ytics 20,20,80
set yrange[20:80]
set ylabel 'density, g/m**3'
#
set nokey
# set key bottom right
#
plot 'fpreactor.data' using (\$1*1.0e6):(\$8*1.0e3) with lines linetype 1
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "reactor-4000K-1bar-massf-N.eps"
#
set title "Nitrogen in a fixed-pressure reactor."
set mxtics 5
set xtics 0,100,600
set xlabel 't, microseconds'
#
set mytics 5
set yrange[0:0.2]
set ytics 0,0.05,0.2
set ylabel 'mass fraction N'
#
set nokey
# set key bottom right
#
plot 'fpreactor.data' using (\$1*1.0e6):5 with lines linetype 1
EOF
