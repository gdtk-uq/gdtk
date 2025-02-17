# plot-p-history.gnuplot
set term pdfcairo font ",12"
set output 'pressure-history-at-exit-centreline.pdf'
set title 'Pressure history at exit centerline for condition 3'
set xrange [3.0:5.0]
set yrange [0:30]
set xlabel 't, ms'
set ylabel 'p, kPa'
set key off
plot 'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):(column('p')*0.001) \
    with lines lw 2

