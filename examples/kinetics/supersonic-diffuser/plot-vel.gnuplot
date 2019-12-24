set term postscript eps 20
set output "velocity.eps"
set title "Supersonic diffuser with hydrogen combustion"
set ylabel "v, km/s"
set xlabel "x, m
set xrange [0:2.1]
set yrange [0.6:2.4]
set key bottom right
plot "drummond-flow.data" using 1:($7/1000) title "steady flow analysis" with lines ls 1 lw 2, \
     "notes/drummond-spectral-17-vel.tsv" using "x":($2/1000) title "Drummond Spectral" with points pt 6 ps 2
