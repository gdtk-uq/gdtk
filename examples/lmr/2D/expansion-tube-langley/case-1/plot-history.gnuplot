# plot-history.gnuplot
set term pdfcairo font ",12"
set output 'pressure-history-at-exit-centreline-1r.pdf'
set title 'Pressure history at exit centerline case-1r'
set xrange [4.0:5.4]
set yrange [0:3]
set xlabel 't, ms'
set ylabel 'p, kPa'
set key off
plot 'lmrsim/hist/blk-0020-cell-122.dat.1' using ($1*1000):($8*0.001) with lines lw 2

