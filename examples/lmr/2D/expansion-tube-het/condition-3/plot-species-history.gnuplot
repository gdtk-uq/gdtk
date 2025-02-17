# plot-species-history.gnuplot
set term pdfcairo font ",12"
set output 'species-history-at-exit-centreline.pdf'
set title 'Species mass-fractions at exit centerline for condition 3'
set xrange [3.0:5.0]
set yrange [-0.10:1.4]
set xlabel 't, ms'
set ylabel 'mass fraction'
set key top right
plot 'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):'massf-CO2' title 'CO2' \
    with lines lw 2, \
    'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):'massf-CO' title 'CO' \
    with lines lw 2, \
    'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):'massf-O' title 'O' \
    with lines lw 2, \
    'lmrsim/hist/hc-01-blk-0013-cell-347.dat.0' \
    using (column('time')*1000):'massf-He' title 'He' \
    with lines lw 2
