# plot-trajectory.gnuplot
#
set term pdfcairo font ",16"
set output 'piston-position.pdf'
set title 'Piston position.'
set ylabel 'x, m'
set xlabel 't, ms'
set key top left
plot 'piston.data' using ($1*1000):($2) title 'simulation' with points pointtype 4 pointsize 1, \
     'analytic-trajectory.data' using ($1*1000):($2) title 'analytic' with lines linewidth 2

set output 'piston-velocity.pdf'
set title 'Piston velocity.'
set ylabel 'v, m/s'
set xlabel 't, ms'
set key top left
plot 'piston.data' using ($1*1000):($3) title 'simulation' with points pointtype 4 pointsize 1, \
     'analytic-trajectory.data' using ($1*1000):($3) title 'analytic' with lines linewidth 2
