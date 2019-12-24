set term postscript eps 20
set output "massf-H2O.eps"
set title "Supersonic diffuser with hydrogen combustion"
set ylabel "massf-H2O"
set xlabel "x, m
set xrange [0:2.1]
set yrange [0:0.07]
set key bottom right
plot "drummond-flow.data" using 1:11 title "steady flow analysis" with lines ls 1 lw 2, \
     "notes/drummond-spectral-17-H2O.tsv" using "x":"massf-H2O" title "Drummond Spectral" with points pt 6 ps 2
