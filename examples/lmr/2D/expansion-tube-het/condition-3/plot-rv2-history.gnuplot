# plot-rv2-history.gnuplot
set term pdfcairo font ",12"
set output 'pitot-history-at-exit-centreline.pdf'
set title 'Approx Pitot pressure history at exit centerline for condition 3'
set xrange [3.0:5.0]
set yrange [0:1200]
set xlabel 't, ms'
set ylabel 'rho*v^2, kPa'
set key off
plot 'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):(0.001*column('rho')*column('vel.x')**2) \
    with lines lw 2

