# sod_run_and_plot.sh
# Sod's 1-D shock tube exercise with a LUT+ideal gas mix.
#
e4shared --prep --job=sod
e4shared --run --job=sod --verbosity=1
e4shared --post --job=sod --tindx-plot=last --vtk-xml 
e4shared --post --job=sod --tindx-plot=last \
         --output-file=sod_new.dat --slice-list="0,:,0,0"

gnuplot<<EOF
set term postscript eps
set output "sod_p.eps"
set title "One-D Shock Tube (ideal driving eq-air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "sod_new.dat" using 1:9 title "air driving eq-air" with points ps 1 pt 1, \
     "../sod_old.dat" using 1:7 title "air driving air" with points ps 1 pt 2, \
     "../analytic/analytic.data" using 1:3 with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_rho.eps"
set title "One-D Shock Tube (ideal air driving eq-air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:1.2]
plot "sod_new.dat" using 1:5 title "air driving eq-air" with points ps 1 pt 1, \
     "../sod_old.dat" using 1:3 title "air driving air" with points ps 1 pt 2, \
     "../analytic/analytic.data" using 1:2 with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_vx.eps"
set title "One-D Shock Tube (ideal air driving eq-air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "sod_new.dat" using 1:6 title "air driving eq-air" with points ps 1 pt 1, \
     "../sod_old.dat" using 1:4 title "air driving air" with points ps 1 pt 2, \
     "../analytic/analytic.data" using 1:5 with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_T.eps"
set title "One-D Shock Tube (ideal air driving eq-air) at t = 0.6ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set key bottom right
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "sod_new.dat" using 1:22 title "air driving eq-air" with points ps 1 pt 1, \
     "../sod_old.dat" using 1:10 title "air driving air" with points ps 1 pt 2, \
     "../analytic/analytic.data" using 1:4 with lines lt 1
EOF
