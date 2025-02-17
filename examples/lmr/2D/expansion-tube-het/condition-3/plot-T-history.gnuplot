# plot-T-history.gnuplot
set term pdfcairo font ",12"
set output 'temperature-history-at-exit-centreline.pdf'
set title 'Temperature history at exit centerline for condition 3'
set xrange [3.0:5.0]
set yrange [0:4000]
set xlabel 't, ms'
set ylabel 'T, K'
set key off
plot 'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):'T' \
    with lines lw 2

