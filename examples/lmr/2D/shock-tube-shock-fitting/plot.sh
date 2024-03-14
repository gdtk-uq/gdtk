# plot.sh
# Sod's 1-D shock tube with fitted-shock boundary.
# Assume that the slice has been extracted from the flow solution
# and that the analytic calculation has already been run.

gnuplot<<EOF
set term postscript eps
set output "sodsf_p.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "eilmer.data" using 1:4 with points ps 1 pt 1, \
     "analytic.data" using 1:3 with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sodsf_rho.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:1.2]
plot "eilmer.data" using 1:3 with points ps 1 pt 1, \
     "analytic.data" using 1:2 with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sodsf_vel.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "eilmer.data" using 1:(-1.0*\$6) with points ps 1 pt 1, \
     "analytic.data" using 1:5 with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sodsf_T.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "eilmer.data" using 1:5 with points ps 1 pt 1, \
     "analytic.data" using 1:4 with lines lt 1
EOF
