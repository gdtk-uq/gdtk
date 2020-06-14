# run_and_plot.sh
# Sod's 1-D shock tube exercise, helium driving reacting nitrogen.
#
l1d4-prep --job=shock_he_n2
l1d4 --run-simulation --job=shock_he_n2
l1d4 --time-slice --job=shock_he_n2 --tindx=40

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_p.eps"
set title "One-D Shock Tube, He driving N2, at t = 80us"
set xlabel "x, m"
set ylabel "Pressure, kPa"
set xrange [0.0:1.0]
set yrange [0.0:1500.0]
plot "slug-0000-tindx-0040-cells.data" using 1:(\$6/1000) title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:(\$6/1000) title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic.data" using 1:(\$3/1000) title "analytic eq-N2" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_rho.eps"
set title "One-D Shock Tube, He driving N2, at t = 80us"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
#set yrange [0.0:0.5]
plot "slug-0000-tindx-0040-cells.data" using 1:5 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:5 title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic.data" using 1:2 title "analytic eq-N2" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_vel.eps"
set title "One-D Shock Tube, He driving N2, at t = 80us"
set xlabel "x, m"
set ylabel "Velocity, km/s"
set key top left
set xrange [0.0:1.0]
set yrange [0.0:6.0]
plot "slug-0000-tindx-0040-cells.data" using 1:(\$3/1000) title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:(\$3/1000) title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic.data" using 1:(\$5/1000) title "analytic eq-N2" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_T.eps"
set title "One-D Shock Tube, He driving N2, at t = 80us"
set xlabel "x, m"
set ylabel "Temperature, K"
set key top left
set xrange [0.0:1.0]
set yrange [0.0:16000.0]
plot "slug-0000-tindx-0040-cells.data" using 1:7 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0040-cells.data" using 1:7 title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic.data" using 1:4 title "analytic eq-N2" with lines lt 1
EOF
