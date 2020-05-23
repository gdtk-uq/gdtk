# contour-p.gnuplot
# Make a contour map of the xt-datasets from all 4 gas slugs.
# Usage:
# $ gnuplot < contour-p.gnuplot
#
# PJ 2020-05-23
#
set terminal pdfcairo
set output "contour-p.pdf"
set title "x,t-diagram of log10(p) in T4 Shot 9945"
set xlabel "x, m"
set ylabel "t, ms"
set view map
set contour base
set surface
set pm3d map
unset key
unset clabel
set cntrparam linear
set cntrparam levels incremental 2.05,0.05,7.50
# set yrange [0:1.0]
# set xrange [0:3.0]
splot "slug-0000-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0000-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0003-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0003-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black'
