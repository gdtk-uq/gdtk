# plot-thrust-evolution.gnuplot
set term png
set output 'thrust-evolution.png'
set ylabel 'Thrust N'
set xlabel 'Objective evaluations'
set xrange [0:80]
set yrange [1800:2300]
set nokey
plot 'progress.txt' using ($1):($7*-1.0) with points pointtype 6 pointsize 2 linewidth 2
