# plot.sh
# Sod's 1-D shock tube with fitted-shock boundary.
#
python3 sod-driving-left.py > sod-driving-left.transcript
e4shared --post --job=sodsf --tindx-plot=last \
         --output-file=eilmer.data --slice-list=":,:,0,0"

gnuplot<<EOF
set term postscript eps
set output "sodsf_p.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "eilmer.data" using 1:9 with points ps 1 pt 1, \
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
plot "eilmer.data" using 1:5 with points ps 1 pt 1, \
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
plot "eilmer.data" using 1:18 with points ps 1 pt 1, \
     "analytic.data" using 1:4 with lines lt 1
EOF
