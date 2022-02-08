# plot-parameter-evolution.gnuplot
set term png
set output 'parameter-evolution.png'
set ylabel 'Parameter value'
set xlabel 'Objective evaluations'
set xrange [0:80]
set yrange [-20:40]
set key bottom right
plot 'progress.txt' using ($1):($3) title 'initial, degrees' with points pointtype 1 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($4) title 'cone, degrees' with points pointtype 2 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($5*1000) title 'deltas[0], mm' with points pointtype 3 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($6*1000) title 'deltas[1], mm' with points pointtype 4 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($7*1000) title 'deltas[2], mm' with points pointtype 5 pointsize 1 linewidth 2
