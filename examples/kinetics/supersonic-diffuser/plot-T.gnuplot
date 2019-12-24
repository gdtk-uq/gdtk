set term postscript eps 20
set output "temperature.eps"
set title "Supersonic diffuser with hydrogen combustion"
set ylabel "T, K"
set xlabel "x, m
set xrange [0:2.1]
set yrange [1200:2400]
plot "drummond-flow.data" using 1:5 title "steady flow analysis" with lines ls 1 lw 2, \
     "notes/drummond-spectral-17-T.tsv" using "x":"T" title "Drummond Spectral" with points pt 6 ps 2
