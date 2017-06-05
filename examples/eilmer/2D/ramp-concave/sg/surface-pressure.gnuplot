# surface-pressure.gnuplot
set term postscript eps 20
set output 'surface-pressure.eps'
set title 'Cubic ramp, pressure along surface'
set ylabel 'p, Pa'
set xlabel 'x, mm'
set key top left
plot './loads/t10-loads.dat' using ($1*1000):($10) with lines \
     lw 3.0 title 'Eilmer4', \
     '../notes/mohammadian-figure-9-p_p_inf.data' \
     using ($1*25.4):($2*66.43) \
     title 'Mohammadian (1972)' with points pt 4
