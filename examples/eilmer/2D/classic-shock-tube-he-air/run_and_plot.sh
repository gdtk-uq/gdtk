# run_and_plot.sh
# Classic 1-D shock tube exercise with a LUT+ideal gas mix.
#
e4shared --prep --job=cst
mpirun -np 2 e4mpi --run --job=cst --verbosity=1
e4shared --post --job=cst --tindx-plot=last --vtk-xml
e4shared --post --job=cst --tindx-plot=last \
         --output-file=profile.data --slice-list="0,:,0,0;1,:,0,0"

gnuplot<<EOF
set term postscript eps
set output "cst-p.eps"
set title "One-D Shock Tube (ideal-He driving eq-air) at t = 100us"
set xlabel "x, m"
set ylabel "Pressure, MPa"
set xrange [0.0:1.0]
set yrange [0.0:32.0]
set key left bottom
plot "profile.data" using 1:(\$9/1.0e6) title "Eilmer4" with points ps 1 pt 6, \
     "analytic.data" using 1:(\$3/1.0e6) title "analytic" with lines lt 1 lw 3
EOF

gnuplot<<EOF
set term postscript eps
set output "cst-rho.eps"
set title "One-D Shock Tube (ideal-He driving eq-air) at t = 100us"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:5.0]
set key left bottom
plot "profile.data" using 1:5 title "Eilmer4" with points ps 1 pt 6, \
     "analytic.data" using 1:2 title "analytic" with lines lt 1 lw 3
EOF

gnuplot<<EOF
set term postscript eps
set output "cst-vx.eps"
set title "One-D Shock Tube (ideal-He driving eq-air) at t = 100us"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0.0:1.0]
set yrange [0.0:3500.0]
set key left top
plot "profile.data" using 1:6 title "Eilmer4" with points ps 1 pt 6, \
     "analytic.data" using 1:5 title "analytic" with lines lt 1 lw 3
EOF

gnuplot<<EOF
set term postscript eps
set output "cst-T.eps"
set title "One-D Shock Tube (ideal-He driving eq-air) at t = 100us"
set xlabel "x, m"
set ylabel "Temperature, K"
set key bottom right
set xrange [0.0:1.0]
set yrange [0.0:5000.0]
set key left bottom
plot "profile.data" using 1:20 title "Eilmer4" with points ps 1 pt 6, \
     "analytic.data" using 1:4 title "analytic" with lines lt 1 lw 3
EOF
