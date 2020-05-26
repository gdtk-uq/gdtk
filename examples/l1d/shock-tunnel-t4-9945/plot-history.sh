# plot-history.sh

echo "Rescale the pressures and times to get convenient units and alignment."
awk '$1 != "#" {print $1*1.0e3 -259.58, $5/1.0e6}' t4-9945/history-loc-0004.data > pstag_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $6}' t4-9945/history-loc-0004.data > Tstag_sim.degreeK
awk -f add-time-column.awk 9945A.200 > T4-9945-spa.us-kPa
awk '$1 != "#" {print $1/1.0e3 -2.700, $2/1.0e3}' T4-9945-spa.us-kPa > pstag_exp.MPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $5/1.0e6}' t4-9945/history-loc-0001.data > pshock1_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $5/1.0e6}' t4-9945/history-loc-0002.data > pshock2_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $5/1.0e6}' t4-9945/history-loc-0003.data > pshock3_sim.MPa
awk -f add-time-column.awk 9945A.220 > T4-9945-ss.us-kPa
awk '$1 != "#" {print $1/1.0e3 -2.700, -0.5 + 10*$2}' T4-9945-ss.us-kPa > shock-sensor_exp.MPa
awk '$1 != "#" {print $1*1.0e3 -256.0, $5/1.0e6}' t4-9945/history-loc-0000.data > pcomp_sim.MPa

echo "Nozzle supply pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4-9945-pstag.eps"
set title "T4 Shot 9945: Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-1.0:6.0]
set yrange [0:40.0]
plot "pstag_sim.MPa" using 1:2 title "L1d4 Simulation" with lines linestyle 1, \
     "pstag_exp.MPa" using 1:2 title "Measured" with lines linestyle 2
EOF

echo "Nozzle supply temperature"
gnuplot <<EOF
set term postscript eps 20
set output "t4-9945-Tstag.eps"
set title "T4 Shot 9945: Nozzle Supply Temperature"
set xlabel "t, ms"
set ylabel "T, K"
set xrange [-1.0:6.0]
set yrange [0:10000.0]
plot "Tstag_sim.degreeK" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF

echo "Compression tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4-9945-pcomp.eps"
set title "T4 Shot 9945: Compression Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-15.0:10.0]
set yrange [0:40.0]
set key top left
plot "pcomp_sim.MPa" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF

echo "Shock tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4-9945-pshock.eps"
set title "T4 Shot 9945: Upstream Shock Tube Pressures"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-2.5:1.0]
set yrange [-4:10.0]
set key top right
plot "pshock1_sim.MPa" using 1:2 title "L1d4 x=30m" with lines linestyle 1, \
     "pshock2_sim.MPa" using 1:2 title "L1d4 x=32m" with lines linestyle 2, \
     "pshock3_sim.MPa" using 1:2 title "L1d4 x=34m" with lines linestyle 3, \
     "shock-sensor_exp.MPa" using 1:2 title "T4 shock sensors" with lines linestyle 4
EOF
