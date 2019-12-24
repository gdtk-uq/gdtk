set term postscript eps 20
set output "massf-H2.eps"
set title "Supersonic diffuser with hydrogen combustion"
set ylabel "massf-H2"
set xlabel "x, m
set xrange [0:2.1]
set yrange [0:0.01]
set key top right
plot "drummond-flow.data" using 1:8 title "steady flow analysis" with lines ls 1 lw 2, \
     "notes/drummond-spectral-17-H2.tsv" using "x":"massf-H2" title "Drummond Spectral" with points pt 6 ps 2
