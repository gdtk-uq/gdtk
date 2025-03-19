# contour.gnuplot
# Make a contour map of the xt-datasets from all 3 gas slugs.
# Usage:
# $ gnuplot < contour.gnuplot
#
# PJ 2020-05-20
#
set terminal pdfcairo
set output "contour.pdf"
set title "x,t-diagram of log10(p) in expansion tube with piston diaphragm"
set xlabel "x, m"
set ylabel "t, ms"
set view map
set contour base
set surface
set pm3d map
unset key
unset clabel
set cntrparam linear
set cntrparam levels incremental 3.0,0.10,5.9
set yrange [0:0.9]
set xrange [0:2.5]
splot "slug-0000-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0000-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0001-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0002-xtdata-p.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black'
