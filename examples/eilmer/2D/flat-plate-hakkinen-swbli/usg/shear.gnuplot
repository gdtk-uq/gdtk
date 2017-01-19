# shear.gnuplot
set term postscript eps 20
set output "shear.eps"
set title "Shear-stress coefficient along the plate pf/p0=1.4"
set ylabel "Cf"
set yrange [-0.001:0.004]
set xlabel "x, mm"
set key right top
plot "./shear.data" using ($1*1000.0):($3) title "Eilmer" with lines, \
     "./shear.data" using ($1*1000.0):($4) title "Blasius" with lines, \
     "../notes/fig6b-shear.data" using ($1*25.4-1.0):($2) \
     title "Hakkinen Fig.6b" with points pt 4
