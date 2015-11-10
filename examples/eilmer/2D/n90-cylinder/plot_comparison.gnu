set term postscript eps 20
set output "n90_compare_T_stag_line.eps"
set title "Stagnation line to 90mm Cylinder in M10 nitrogen flow"
set xlabel "x, m"
set ylabel "temperature, K"
set xrange [-0.020:0.0]
set xtics 0.005
set yrange [0:20000]
set ytics 5000.0
set key left top
plot "stag_sebo.dat" using 1:7 title "CEVCATS-N Reference" with lines, \
     "n90_100_iy1.data" using 1:22 title "Nonequilibrium"

