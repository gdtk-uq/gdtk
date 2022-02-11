# plot-parameter-evolution.gnuplot
set term png
set output 'parameter-evolution.png'
set ylabel 'Angle, degree'
set xlabel 'Objective evaluations'
set xrange [0:80]
set yrange [-10:45]
set key top right
plot 'progress.txt' using ($1):($3) title 'initial' with points pointtype 1 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($4) title 'alpha' with points pointtype 2 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($5) title 'beta' with points pointtype 3 pointsize 1 linewidth 2, \
     'progress.txt' using ($1):($6) title 'cone' with points pointtype 4 pointsize 1 linewidth 2
