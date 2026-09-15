# plot-mach.gnuplot
set term postscript eps 20
set output 'rao-nozzle-mach.eps'

set title 'Rao thrust-optimised nozzle A, wall Mach number'
set xlabel 'x/r_A'
set ylabel 'M'
set key right bottom

set xrange [0:9]
set yrange [1:4]
plot 'nozzle-wall.data' using 1:3 title "gdtk.imoc" with lines lw 2, \
     'kunze-wall-mach.data' using 1:2 title "Kunze (2020) IMOC, digitised" \
        with lines dt 2 lw 1, \
     'rao-wall-mach.data' using 1:2 title "Rao (1958), digitised" \
        with points pt 4 ps 1.5
