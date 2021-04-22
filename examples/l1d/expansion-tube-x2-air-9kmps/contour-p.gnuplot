# contour-p.gnuplot
# Make a contour map of the xt-datasets from all 3 gas slugs.
# Usage:
# $ gnuplot < contour-p.gnuplot
#
# PJ 2020-05-23, 2020-07-22
#
set terminal pdfcairo
set output "contour-p.pdf"
set title "x,t-diagram of log10(p) in Drummond Tunnel"
set xlabel "x, m"
set ylabel "t, ms"
set view map
set contour base
set surface
set pm3d map
unset key
unset clabel
set cntrparam linear
set cntrparam levels incremental 2.05,0.05,6.50
set yrange [23.0:26.0]
set xrange [4.0:15.0]
splot "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0003-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0003-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black'
