# plot-t-along-axis.gnuplot
# Plot the transrotational and vibrational temperatures
set term postscript eps 20
set output 'ttv-along-axis.eps'
set title 'Normalized temperatures along nozzle axis'
set xlabel 'x, cm'
set ylabel 'T/T0'
set xrange [-2.0:6.0]
set yrange [0.3:1.1]
set key top right
plot 'normalized-tvib-experiment.data' using 1:2:3 title 'Tvib/T0 experiment' with errorbars pt 7 ps 1 lw 2, \
     'axis-Blackman.data' using ($1*100.0):($19/5727.0) title 'Ttrans/T0' with lines ls 2 lw 2, \
     'axis-Blackman.data' using ($1*100.0):($21/5727.0) title 'Tvib/T0 (Blackman)' with lines ls 1 lw 2

