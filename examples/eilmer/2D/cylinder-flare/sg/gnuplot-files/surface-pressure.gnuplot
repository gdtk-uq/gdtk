# surface-pressure.gnuplot
set term postscript eps 20
set output 'surface-pressure.eps'
set title 'Cylinder with extended flare, pressure along surface'
set ylabel 'p, Pa'
set xlabel 'x, mm'
set key below
plot '../eilmer4-surface.data' using ($1*1000):($3) with lines \
     lw 3.0 title 'Eilmer4', \
     './cylinder-extended-flare-pressure.data' \
     using ($2*101.7):($10*6894.8) \
     title 'CUBRC Run 14' with points pt 4
