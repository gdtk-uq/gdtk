# run_and_plot.sh
# Sod's 1-D shock tube exercise, helium driving equilibrium air.
#
l1d4-prep --job=sod-he-air
l1d4 --run-simulation --job=sod-he-air
l1d4 --time-slice --job=sod-he-air --tindx=40

gnuplot<<EOF
set term postscript eps 20
set output "sod_p.eps"
set title "One-D Shock Tube, He driving air, at t = 0.4ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "slug-0000-tindx-0040-cells.data" using 1:6 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:6 title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic.data" using 1:3 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "sod_rho.eps"
set title "One-D Shock Tube, He driving air, at t = 0.4ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:0.5]
plot "slug-0000-tindx-0040-cells.data" using 1:5 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:5 title "L1d4 slug 1" with points ps 1 pt 1, \
     "analytic.data" using 1:2 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "sod_vel.eps"
set title "One-D Shock Tube, He driving air, at t = 0.4ms"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set key top left
set xrange [0.0:1.0]
set yrange [0.0:550.0]
plot "slug-0000-tindx-0040-cells.data" using 1:3 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:3 title "L1d4 slug 1" with points ps 1 pt 1, \
     "analytic.data" using 1:5 title "analytic" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "sod_T.eps"
set title "One-D Shock Tube, He driving air, at t = 0.4ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set key bottom left
set xrange [0.0:1.0]
set yrange [0.0:550.0]
plot "slug-0000-tindx-0040-cells.data" using 1:7 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:7 title "L1d4 slug 1" with points ps 1 pt 1, \
     "analytic.data" using 1:4 title "analytic" with lines lt 1
EOF
