# contour-T-zoom-into-shock.gnuplot
# Make a contour map of the xt-datasets from all 4 gas slugs,
# but zoom into the shock-processing part of the simulation.
# Usage:
# $ gnuplot < contour-T-zoom-into-shock.gnuplot
#
# PJ 2020-05-23
#
set terminal pdfcairo
set output "contour-T-zoom-into-shock.pdf"
set title "x,t-diagram of T in HEG Shot 495"
set xlabel "x, m"
set ylabel "t, ms"
set view map
set contour base
set surface
set pm3d map
unset key
unset clabel
set cntrparam linear
set cntrparam levels incremental 301,300,10000
set yrange [154.0:170.0]
set xrange [32.0:51.5]
splot "slug-0000-xtdata-T.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0000-xtdata-T.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0001-xtdata-T.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0001-xtdata-T.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0002-xtdata-T.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0002-xtdata-T.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black', \
      "slug-0003-xtdata-T.data" using ($1):($2*1000):($3) with pm3d at b, \
      "slug-0003-xtdata-T.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black'
