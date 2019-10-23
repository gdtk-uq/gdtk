# Make a contour map of the xt-data.
# Usage:
# $ gnuplot < contour.gnuplot
#
# PJ 2019-10-23
#
set terminal pdfcairo
set output "contour.pdf"
set title "x,t-diagram of log10(p) in expansion tube demo"
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
set yrange [0:1.0]
set xrange [0:2.5]
splot "xtlogp.data" using ($1):($2*1000):($3) with pm3d at b, \
      "xtlogp.data" using ($1):($2*1000):($3) nosurface with lines lc rgb 'black'
