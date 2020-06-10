# plot-history.sh

echo "Compression tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "x3s561-pcomp-two-bump.eps"
set title "X3 Shot 561: Compression Tube Pressure"
set xlabel "t, s"
set ylabel "p, MPa"
set xrange [-0.20:0.50]
set yrange [0:20.0]
set key top right
plot "x3r-561/history-loc-0001.data" using (\$1-0.1705):(\$5/1.0e6) title "L1d4 cells 600:300:50" with lines linestyle 1, \
     "exp-data/x3s561/x3s561.lvm" using (\$1-1.0522):(\$2/1000.0) title "x3s561 Experiment" with lines linestyle 2
EOF

gnuplot <<EOF
set term postscript eps 20
set output "x3s561-pcomp-one-bump.eps"
set title "X3 Shot 561: Compression Tube Pressure"
set xlabel "t, s"
set ylabel "p, MPa"
set xrange [-0.10:0.10]
set yrange [0:20.0]
set key top right
plot "x3r-561/history-loc-0001.data" using (\$1-0.1705):(\$5/1.0e6) title "L1d4 cells 600:300:50" with lines linestyle 1, \
     "exp-data/x3s561/x3s561.lvm" using (\$1-1.0522):(\$2/1000.0) title "x3s561 Experiment" with lines linestyle 2
EOF
