#! /bin/bash
# plot-surface-p.sh
# PJ, 2016-05-26
e4shared --post --job=bcp --tindx-plot=last --add-vars="mach,pitot,total-p,total-h" \
         --slice-list="2:11,:,0,0" --output-file="surface.data"

gnuplot <<EOF
set term postscript eps 20
set output 'surface_p.eps'

set title 'Static pressure along cone surface'
set key bottom right
# set xrange [0.0:20.0]
set yrange [0.0:18.0]
set xlabel 'x, mm'
set ylabel 'p, kPa'

set style line 1 linetype 1 linewidth 3.0 
set arrow from 15.0,11.1 to 17.0,11.1 nohead linestyle 1
set label "Value from\nNACA 1135\nChart 6" at 14.0,16.0 right
set arrow from 14.0,16.0 to 16.0,11.1 head

plot 'surface.data' using (\$1*1000):(\$9/1000) title 'inviscid, t=100us' with lines 
EOF
