#!/bin/bash
# run.sh
chkn-prep --job=sod --binary
chkn-run --job=sod --binary
chkn-post --job=sod --binary --vtk-xml --tindx=all
chkn-post --job=sod --binary --slice="0,0,0,:,0,0;1,0,0,:,0,0" --tindx=$ --output-file=tube-data

gnuplot<<EOF
set output "tube-density.eps"
set term postscript eps enhanced 20
set xlabel "x, m"
set ylabel "density, kg/m^3"
set yrange [0.0:1.2]
plot "tube-data-t0001.data" using 1:7 title "chicken run"
EOF
