# plot-history.sh

echo "Rescale the pressures and times to get convenient units and alignment."
awk '$1 != "#" {print $1*1.0e3 -327.4, $5/1.0e6}' t4-11316/history-loc-0004.data > pstag_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -327.4, $6}' t4-11316/history-loc-0004.data > Tstag_sim.degreeK
# awk -f add-time-column.awk 11316A.200 > T4-11316-spa.us-kPa
awk -f add-time-column.awk 11316A.210 > T4-11316-spb.us-kPa
awk '$1 != "#" {print $1/1.0e3 -4.800, $2/1.0e3}' T4-11316-spb.us-kPa > pstag_exp.MPa
awk '$1 != "#" {print $1*1.0e3 -327.4, $5/1.0e6}' t4-11316/history-loc-0001.data > pshock1_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -327.4, $5/1.0e6}' t4-11316/history-loc-0002.data > pshock2_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -327.4, $5/1.0e6}' t4-11316/history-loc-0003.data > pshock3_sim.MPa
awk -f add-time-column.awk 11316A.220 > T4-11316-ss.us-kPa
awk '$1 != "#" {print $1/1.0e3 -5.200, -0.5 + 10*$2}' T4-11316-ss.us-kPa > shock-sensor_exp.MPa
awk '$1 != "#" {print $1*1.0e3 -320.0, $5/1.0e6}' t4-11316/history-loc-0000.data > pcomp_sim.MPa

echo "Nozzle supply pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4-11316-pstag.eps"
set title "T4 Shot 11316: Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-1.0:15.0]
set yrange [-5:45.0]
plot "pstag_sim.MPa" using 1:2 title "L1d4 Simulation" with lines linestyle 1, \
     "pstag_exp.MPa" using 1:2 title "Measured" with lines linestyle 2
EOF

echo "Nozzle supply temperature"
gnuplot <<EOF
set term postscript eps 20
set output "t4-11316-Tstag.eps"
set title "T4 Shot 11316: Nozzle Supply Temperature"
set xlabel "t, ms"
set ylabel "T, K"
set xrange [-1.0:15.0]
set yrange [0:2000.0]
plot "Tstag_sim.degreeK" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF

echo "Compression tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4-11316-pcomp.eps"
set title "T4 Shot 11316: Compression Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-15.0:10.0]
set yrange [0:30.0]
set key top left
plot "pcomp_sim.MPa" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF

echo "Shock tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4-11316-pshock.eps"
set title "T4 Shot 11316: Upstream Shock Tube Pressures"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-10:10]
set yrange [-20:20.0]
set key top left
plot "pshock1_sim.MPa" using 1:2 title "L1d4 x=30m" with lines linestyle 1, \
     "pshock2_sim.MPa" using 1:2 title "L1d4 x=32m" with lines linestyle 2, \
     "pshock3_sim.MPa" using 1:2 title "L1d4 x=34m" with lines linestyle 3, \
     "shock-sensor_exp.MPa" using 1:2 title "T4 shock sensors" with lines linestyle 4
EOF
