# plot.sh
# Sod's 1-D shock tube exercise as a 3D simulation overkill.
# PJ 2024-03-03
#
gnuplot<<EOF
set term postscript eps
set output "sod_p.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "sod_new.dat" using "pos.x":"p" title "lmr" with points ps 1 pt 1, \
     "../analytic/analytic.data" using 1:3 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_rho.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:1.2]
plot "sod_new.dat" using "pos.x":"rho" title "lmr" with points ps 1 pt 1, \
     "../analytic/analytic.data" using 1:2 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_u.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "sod_new.dat" using "pos.x":"vel.x" title "lmr" with points ps 1 pt 1, \
     "../analytic/analytic.data" using 1:5 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_T.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "sod_new.dat" using "pos.x":"T" title "lmr" with points ps 1 pt 1, \
     "../analytic/analytic.data" using 1:4 title "analytic" with lines lt 1
EOF
