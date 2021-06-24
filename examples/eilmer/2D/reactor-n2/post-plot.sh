#!/bin/bash
# post-plot.sh
# Display history for the fixed-volume-reactor simulation. 
#
# PJ & RG, 2018-04-21

gas-calc fvreactor.lua

gnuplot <<EOF
set terminal postscript eps 20
set output "reactor-4000K-1bar-T.eps"
#
set title "Nitrogen in a fixed-volume reactor."
set mxtics 5
set xtics 0,100,400
set xlabel 't, microseconds'
#
set mytics 5
set ytics 4000,500,7000
set yrange[4000:7000]
set ylabel 'T, K'
#
set key bottom right
#
plot 'fvreactor.data' using (\$1*1.0e6):2 with lines linetype 1 title 'fvreactor script', \
     'hist/reactor-blk-0-cell-0.dat.0' using (\$1*1.0e6):21 with points pointtype 4 title 'Eilmer'
EOF


gnuplot <<EOF
set terminal postscript eps 20
set output "reactor-4000K-1bar-p.eps"
#
set title "Nitrogen in a fixed-volume reactor."
set mxtics 5
set xtics 0,100,400
set xlabel 't, microseconds'
#
set mytics 5
set ytics 1,0.2,2
set yrange[1:2]
set ylabel 'p, bar'
#
set key bottom right
#
plot 'fvreactor.data' using (\$1*1.0e6):(\$3/1.0e5) with lines linetype 1 title 'fvreactor script', \
     'hist/reactor-blk-0-cell-0.dat.0' using (\$1*1.0e6):(\$10/1.0e5) with points pointtype 4 title 'Eilmer'
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "reactor-4000K-1bar-massf-N.eps"
#
set title "Nitrogen in a fixed-volume reactor."
set mxtics 5
set xtics 0,100,400
set xlabel 't, microseconds'
#
set mytics 5
set yrange[0:0.2]
set ytics 0,0.05,0.2
set ylabel 'mass fraction N'
#
set key bottom right
#
plot 'fvreactor.data' using (\$1*1.0e6):5 with lines linetype 1 title 'fvreactor script', \
     'hist/reactor-blk-0-cell-0.dat.0' using (\$1*1.0e6):18 with points pointtype 4 title 'Eilmer'
EOF
