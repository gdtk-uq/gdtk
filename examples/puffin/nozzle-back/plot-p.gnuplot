# plot-p.gnuplot
set term postscript eps 20
set output 'normalized-pressure.eps'

set title 'Pressure along the nozzle wall'
set xlabel 'distance from throat (inches)'
set ylabel 'p/pt'

set xrange [0.0:3.0]
set yrange [0:0.6]
plot 'nozzle/nozzle-streamline-0-cell-59.data' using ($1/0.0254):($7/500.0e3) \
        title "Puffin calculation" with lines lw 2, \
     'back-experiment.data' using 1:2 title "experiment" with points pt 4 ps 1.5
