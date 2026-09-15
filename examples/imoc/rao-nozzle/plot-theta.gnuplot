# plot-theta.gnuplot
set term postscript eps 20
set output 'rao-nozzle-theta.eps'

set title 'Rao thrust-optimised nozzle A, wall angle'
set xlabel 'x/r_A'
set ylabel 'theta, degrees'
set key right top

set xrange [0:9]
set yrange [0:40]
plot 'nozzle-wall.data' using 1:4 title "gdtk.imoc" with lines lw 2, \
     'kunze-wall-theta.data' using 1:2 title "Kunze (2020) IMOC, digitised" \
        with lines dt 2 lw 1, \
     'rao-wall-theta.data' using 1:2 title "Rao (1958), digitised" \
        with points pt 4 ps 1.5
