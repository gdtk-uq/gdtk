# plot-poshax.sh

gnuplot <<EOF
set term postscript eps 20
set output "stagnation-line-T.eps"
set xlabel "x, m"
set ylabel "T/T_{/Symbol \245}"
set yrange [0:12]
set xrange [-0.1:0.5]

plot "stagnation-line-poshax.data" using 1:(\$2/300) title "T" with lines lt 1, \
     "stagnation-line-poshax.data" using 1:(\$9/300) title "Tvib" with lines lt 2

EOF
