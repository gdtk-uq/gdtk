# Make a contour map of the xt-data.
# Usage:
# $ gnuplot < contour.gnuplot
#
# PJ 2020-02-13
#
set terminal pdfcairo
set output "contour.pdf"
set title "x,t-diagram of temperature"
set xlabel "x, m"
set ylabel "t, ms"
set view map
set contour base
set surface
set pm3d map
unset key
unset clabel
set cntrparam linear
set cntrparam levels incremental 250,20,390
set yrange [0:0.6]
set xrange [0:1.0]
splot "xt.data" using ($2):($3*1000):($5) with pm3d at b, \
      "xt.data" using ($2):($3*1000):($5) nosurface with lines lc rgb 'black'
