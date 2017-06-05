# surface-heat-transfer.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer.eps'
set title 'Cubic ramp, heat-flux along surface'
set ylabel 'St.Re^(1/2)'
set xlabel 'x, mm'
set yrange [0:15]
set key top left
plot './stanton.data' using ($1*1000):($2) with lines \
     lw 3.0 title 'Eilmer4', \
     '../notes/mohammadian-figure-10-stanton.data' \
     using ($1*25.4):($2) \
     title 'Mohammadian (1972)' with points pt 4
