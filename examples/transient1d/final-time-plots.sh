# final-time-plots.sh
#
# Plot the gas properties along the duct at the final time
# of the transient-1D flow simulation.
# Compare with the analytic data from the classic-shock-tube
# state-to-state analysis script.
#

gnuplot<<EOF
set term postscript eps 20
set output "final_p.eps"
set title "One-D Shock Tube (ideal air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "final.data" using 2:4 title "t1d" with points ps 1 pt 1, \
     "analytic.data" using 1:3 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "final_rho.eps"
set title "One-D Shock Tube (ideal air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:1.2]
plot "final.data" using 2:6 title "t1d" with points ps 1 pt 1, \
     "analytic.data" using 1:2 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "final_vx.eps"
set title "One-D Shock Tube (ideal air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "final.data" using 2:7 title "t1d" with points ps 1 pt 1, \
     "analytic.data" using 1:5 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "final_T.eps"
set title "One-D Shock Tube (ideal air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set key bottom right
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "final.data" using 2:5 title "t1d" with points ps 1 pt 1, \
     "analytic.data" using 1:4 title "analytic" with lines lt 1
EOF
