# plot-history-simple.sh

echo "Compression tube pressure"

gnuplot <<EOF
set term postscript eps 20
set output "x3s561-pcomp-simple-loc-0000.eps"
# set title "X3 Shot 561: Compression Tube Pressure"
set xlabel "t, s"
set ylabel "p, MPa"
set xrange [140.0:185.0]
set yrange [0:18.0]
set key top right
plot "x3r-561/history-loc-0000.data" using (\$1*1000):(\$5/1.0e6) title "L1d4 loc 0 cells 600:300:50" with lines linestyle 1
EOF

gnuplot <<EOF
set term postscript eps 20
set output "x3s561-pcomp-simple-loc-0001.eps"
# set title "X3 Shot 561: Compression Tube Pressure"
set xlabel "t, s"
set ylabel "p, MPa"
set xrange [140.0:185.0]
set yrange [0:18.0]
set key top right
plot "x3r-561/history-loc-0001.data" using (\$1*1000):(\$5/1.0e6) title "L1d4 loc 1 cells 600:300:50" with lines linestyle 1
EOF
