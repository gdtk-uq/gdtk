#!/bin/bash
# plot-surface-mach.sh
# Plot the Mach number around the blade surface.

e4shared --post --job=kiock --tindx-plot=last --output-file=surface-flow.data \
         --add-vars="mach,pitot,total-p,total-h" \
         --slice-list="0,:,0,0;1,:,0,0;2,:,0,0;3,:,0,0;4,$,:,0;5,:,$,0;6,:,$,0"
#                      0south  1south  2south  3south  4east   5north  6north

awk -f compute-surface-mach.awk surface-flow.data > surface-mach.data

gnuplot<<EOF
set term postscript eps 20
set output "surface-mach.eps"
set xlabel "chord location"
set ylabel "surface Mach number"
set xrange [-0.2:1.2]
set yrange [-0.2:1.2]
set key bottom right

plot "kiock_experimental_RG.dat" using 1:2 title "Rhode-St.-Genese" with points pt 1, \
     "kiock_experimental_GO.dat" using 1:2 title "Goettingen" with points pt 2, \
     "kiock_experimental_BS.dat" using 1:2 title "Braunschweig" with points pt 3, \
     "kiock_experimental_OX.dat" using 1:2 title "Oxford" with points pt 4, \
     "surface-mach.data" using 1:2 title "Eilmer4" with points pt 5

EOF
