# pressure.gnuplot
set term postscript eps 20
set output "pressure.eps"
set title "Static pressure along the plate pf/p0=1.4"
set ylabel "p, kPa"
set yrange [0:10]
set xlabel "x, mm"
set key left top
plot "./bl.data" using ($1*1000.0):($9/1000.0) title "Eilmer" with lines, \
     "../notes/fig6b-pressure.data" using ($1*25.4-1.0):($2/0.125*6.3) \
     title "Hakkinen Fig.6b" with points pt 4
