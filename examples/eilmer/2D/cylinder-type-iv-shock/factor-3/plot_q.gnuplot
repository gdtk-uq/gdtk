set term postscript eps enhanced 20
set output "surface_heat_transfer.eps"
set style line 1 linetype 1 linewidth 3.0
set title "Surface heat transfer"
set xlabel "angle, degrees"
set xrange [-90:90]
set ylabel "q, W/cm**2"
set yrange [0:60]
set key top right
plot "surface.data" using 1:3 title "180x360" with lines
