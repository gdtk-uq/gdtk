# contour-p-zoom-into-shock.gnuplot
# Make a contour map of the xt-datasets from all 4 gas slugs,
# but zoom into the shock-processing part of the simulation.
# Usage:
# $ gnuplot < contour-p-zoom-into-shock.gnuplot
#
# PJ 2020-05-23, 2020-07-16
#
set terminal pdfcairo
set output "contour-p-zoom-into-shock.pdf"
set title "x,t-diagram of log10(p) in T4 Shot 11316"
set xlabel "x, m"
set ylabel "t, ms"
set view map
set contour base
set surface
set pm3d map
unset key
unset clabel
set cntrparam linear
set cntrparam levels incremental 5.5,0.02,7.30
set yrange [310.0:340.0]
set xrange [24:36]
splot "slug-0000-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0000-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0003-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0003-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black'
