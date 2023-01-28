# plot.sh
# Tamara's 1-D shock tube exercise, helium driving 1T air.
#
l1d4 --time-slice --job=he-air-1T --tindx=140
python3 analytic_he_air_eq.py
python3 analytic_he_air_frz.py

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_p.eps"
set title "One-D Shock Tube, He driving 1T air, at t = 2.8ms"
set xlabel "x, m"
set ylabel "Pressure, kPa"
set xrange [-2.0:9.0]
set yrange [0.0:200.0]
plot "slug-0000-tindx-0140-cells.data" using 1:(\$6/1000) title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0140-cells.data" using 1:(\$6/1000) title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic_eq.data" using 1:(\$3/1000) title "analytic eq-air" with lines lt 1, \
     "analytic_frz.data" using 1:(\$3/1000) title "analytic frz-air" with lines lt 2
EOF

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_rho.eps"
set title "One-D Shock Tube, He driving 1T air, at t = 2.8ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [-2.0:9.0]
set yrange [0.0:0.2]
plot "slug-0000-tindx-0140-cells.data" using 1:5 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0140-cells.data" using 1:5 title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic_eq.data" using 1:2 title "analytic eq-air" with lines lt 1, \
     "analytic_frz.data" using 1:2 title "analytic frz-air" with lines lt 2
EOF

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_vel.eps"
set title "One-D Shock Tube, He driving 1T air, at t = 2.8us"
set xlabel "x, m"
set ylabel "Velocity, km/s"
set key top left
set xrange [-2.0:9.0]
set yrange [0.0:4.0]
plot "slug-0000-tindx-0140-cells.data" using 1:(\$3/1000) title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0140-cells.data" using 1:(\$3/1000) title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic_eq.data" using 1:(\$5/1000) title "analytic eq-air" with lines lt 1, \
     "analytic_frz.data" using 1:(\$5/1000) title "analytic frz-air" with lines lt 2
EOF

gnuplot<<EOF
set term postscript eps 20
set output "shock_tube_T.eps"
set title "One-D Shock Tube, He driving 1T air, at t = 2.8ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set key top left
set xrange [-2.0:9.0]
set yrange [0.0:4000.0]
plot "slug-0000-tindx-0140-cells.data" using 1:7 title "L1d4 slug 0" with points ps 1 pt 1, \
     "slug-0001-tindx-0140-cells.data" using 1:7 title "L1d4 slug 1" with points ps 1 pt 2, \
     "analytic_eq.data" using 1:4 title "analytic eq-air" with lines lt 1, \
     "analytic_frz.data" using 1:4 title "analytic frz-air" with lines lt 2
EOF
