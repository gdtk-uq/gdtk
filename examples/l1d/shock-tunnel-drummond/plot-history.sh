#! /bin/sh
# plot-history.sh
#

# Recale the pressures and heat-transfer to get more convenient units.
awk '$1 != "#" {print $1*1000.0, $5/1.0e6}' dn2/history-loc-0001.data > pstag.Mpa
awk '$1 != "#" {print $1*1000.0, $6}' dn2/history-loc-0001.data > Tstag.Mpa
awk '$1 != "#" {print $1*1000.0, $10/-1.0e6}' dn2/history-loc-0002.data > heatflux.Mwatt
awk '$1 != "#" {print $1*1000.0, 0.92*$4*$2*$2/1.0e3}' dn2/history-loc-0004.data > pitot.kpa

gnuplot <<EOF
set term postscript eps 20
set output "dn2_ppitot.eps"
set title "Drummond Tunnel with N2 Driver, Pitot Pressure"
set xlabel "t, ms"
set ylabel "0.92*rho*vel^2, kPa"
set xrange [3.0:8.0]
set xtic 1.0
set key top left
plot "pitot.kpa" using 1:2 title "Simulation" with lines linestyle 1
EOF

gnuplot <<EOF
# was file heat_trans.gnu
set term postscript eps 20
set output "dn2_heatflux.eps"
set title "Drummond Tunnel with N2 Driver, Heat Transfer"
set xlabel "t, ms"
set ylabel "Q, MW/m**2"
set xrange [3.0:8.0]
set xtic 1.0
set key top left
plot "heatflux.Mwatt" using 1:2 title "Simulation" with lines linestyle 1
EOF

gnuplot <<EOF
# was file pstag.gnu
set term postscript eps 20
set output "dn2_pstag.eps"
set title "Drummond Tunnel with N2 Driver, Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [3.0:8.0]
set xtic 1.0
set yrange [0:2.0]
set key top left
plot "pstag.Mpa" using 1:2 title "Simulation" with lines linestyle 1
EOF

gnuplot <<EOF
# was file pstag.gnu
set term postscript eps 20
set output "dn2_Tstag.eps"
set title "Drummond Tunnel with N2 Driver, Nozzle Supply Temperature"
set xlabel "t, ms"
set ylabel "T, K"
set xrange [3.0:8.0]
set xtic 1.0
set yrange [0:1200.0]
set key top left
plot "Tstag.Mpa" using 1:2 title "Simulation" with lines linestyle 1
EOF

# Now set up a couple of files that describe the tube cross-section.
awk '$1 != "#" {print $1, $2/2.0}' dn2/tube.data > upper.profile
awk '$1 != "#" {print $1, -$2/2.0}' dn2/tube.data > lower.profile

gnuplot <<EOF
# was file d_tube.gnu
set term postscript eps 20
set output "dn2_tube.eps"
set title "Drummond Tunnel Facility Profile"
set xlabel "x, m"
set ylabel "y, m"
set xrange [-4:0.5]
set yrange [-0.1:0.1]
plot "upper.profile" using 1:2 title "" with lines linestyle 1, \
     "lower.profile" using 1:2 title "" with lines linestyle 1
EOF

gnuplot <<EOF
# was file d_noz.gnu
set term postscript eps 20
set output "dn2_noz.eps"
set title "Drummond Tunnel Nozzle Profile"
set xlabel "x, m"
set ylabel "y, m"
set xrange [0:0.3]
set yrange [-0.1:0.1]
plot "upper.profile" using 1:2 title "" with lines linestyle 1, \
     "lower.profile" using 1:2 title "" with lines linestyle 1
EOF
