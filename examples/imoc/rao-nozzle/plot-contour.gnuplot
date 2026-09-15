# plot-contour.gnuplot
set term postscript eps 20
set output 'rao-nozzle-contour.eps'

set title 'Rao thrust-optimised nozzle A, contour'
set xlabel 'x/r_A'
set ylabel 'r/r_A'
set key left top

set xrange [0:9]
set yrange [0:5]
plot 'nozzle-wall.data' using 1:2 title "gdtk.imoc" with lines lw 2, \
     'kunze-contour.data' using 1:2 title "Kunze (2020) IMOC, digitised" \
        with lines dt 2 lw 1, \
     'rao-contour.data' using 1:2 title "Rao (1958), digitised" \
        with points pt 4 ps 1.5
