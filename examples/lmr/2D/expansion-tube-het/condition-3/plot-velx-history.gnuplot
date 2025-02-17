# plot-velx-history.gnuplot
set term pdfcairo font ",12"
set output 'velocity-history-at-exit-centreline.pdf'
set title 'Velocity history at exit centerline for condition 3'
set xrange [3.0:5.0]
set yrange [0:4000]
set xlabel 't, ms'
set ylabel 'vel.x, m/s'
set key off
plot 'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):'vel.x' \
    with lines lw 2

