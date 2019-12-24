set term postscript eps 20
set output "massf-OH.eps"
set title "Supersonic diffuser with hydrogen combustion"
set ylabel "massf-OH"
set xlabel "x, m
set xrange [0:2.1]
plot "drummond-flow.data" using 1:10 title "steady flow analysis" with lines ls 1 lw 2, \
     "notes/drummond-spectral-17-OH.tsv" using "x":"massf-OH" title "Drummond Spectral" with points pt 6 ps 2
