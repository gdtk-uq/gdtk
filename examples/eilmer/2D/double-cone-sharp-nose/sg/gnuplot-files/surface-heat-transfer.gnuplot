# surface-heat-transfer.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer.eps'
set title 'Double-cone sharp-nose, heat-flux along surface'
set ylabel 'q, kW/m**2'
set xlabel 'x, mm'
set yrange [0:1500]
set key top right
plot '../eilmer4-surface.data' using ($1*1000):($5/1000) with lines \
     lw 3.0 title 'Eilmer4 k*dT/dy', \
     './indented-cone-heat-transfer.data' \
     using ($2*74.69+34.76):($9*11.377) \
     title 'CUBRC Run 35 x-affine' with points pt 4
