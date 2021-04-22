#! /bin/sh
# plot-data.sh

# Data from the history files are in SI units

# Plotting the pressures along the X2 facility
# --------------------------------------------

# Rescale the pressure data to get more convenient units
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0002.data > psd1.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0003.data > psd2.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0004.data > psd3.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0005.data > pst1.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0006.data > pst2.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0007.data > pst3.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0008.data > pat1.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0009.data > pat2.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0010.data > pat3.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0011.data > pat4.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0012.data > pat5.kpa
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0013.data > pat6.kpa

# Plot shock and acceleration tube pressure data
gnuplot << EOF
set term png
set output "x2_condition1_shocktube.png"
set title "PCB pressure data along tunnel - X2 Nonequilibrium DP Condition 1 (air)"
set xlabel "t, ms"
set ylabel "p, kPa"
set xtic 100.0e-3
set xrange [24.3:24.7]
set yrange [0:100]
set key top left
plot "pst1.kpa" using 1:2 title "st1, grid2, inv" with lines linestyle 1, \
     "pst2.kpa" using 1:2 title "st2, grid2, inv" with lines linestyle 2, \
     "pst3.kpa" using 1:2 title "st3, grid2, inv" with lines linestyle 3
EOF

gnuplot << EOF
set term png
set output "x2_condition1_acceltube.png"
set title "PCB pressure data along tunnel - X2 Nonequilibrium DP Condition 1 (air)"
set xlabel "t, ms"
set ylabel "p, kPa"
set xtic 100.0e-3
set xrange [24.6:25.0]
set yrange [0:25]
set key top left
plot "pat1.kpa" using 1:2 title "at1, grid2, inv" with lines linestyle 1, \
     "pat2.kpa" using 1:2 title "at2, grid2, inv" with lines linestyle 2, \
     "pat3.kpa" using 1:2 title "at3, grid2, inv" with lines linestyle 3, \
     "pat4.kpa" using 1:2 title "at4, grid2, inv" with lines linestyle 4, \
     "pat5.kpa" using 1:2 title "at5, grid2, inv" with lines linestyle 5, \
     "pat6.kpa" using 1:2 title "at6, grid2, inv" with lines linestyle 6
EOF


# Plotting data for the nozzle exit
# ---------------------------------

# Rescale the pressure data to get more convenient units
awk '$1 != "#" {print $1*1000.0, $5/1.0e3}' x2_condition1/history-loc-0015.data > pressure.kpa
awk '$1 != "#" {print $1*1000.0, $6}' x2_condition1/history-loc-0015.data > temperature.K
awk '$1 != "#" {print $1*1000.0, 0.92*$4*$2*$2/1.0e3}' x2_condition1/history-loc-0015.data > pitot.kpa

# Plot nozzle supply pressure data
gnuplot << EOF
set term png
set output "x2_condition1_nozzleP.png"
set title "Nozzle supply pressure - X2 Nonequilibrium DP Condition 1 (air)"
set xlabel "t, ms"
set ylabel "p, kPa"
set xtic 100.0e-3
set xrange [25.0:26.0]
set key top left
plot "pressure.kpa" using 1:2 title "grid2, inv" with lines linestyle 1
EOF

# Plot nozzle supply temperature data
gnuplot << EOF
set term png
set output "x2_condition1_nozzleT.png"
set title "Nozzle supply temperature - X2 Nonequilibrium DP Condition 1 (air)"
set xlabel "t, ms"
set ylabel "T, K"
set xtic 100.0e-3
set xrange [25.0:26.0]
set key top left
plot "temperature.K" using 1:2 title "grid2, inv" with lines linestyle 1
EOF

# Plot nozzle pitot pressure data
gnuplot << EOF
set term png 
set output "x2_condition1_pitotP.png"
set title "Nozzle pitot pressure - X2 Nonequilibrium DP Condition 1 (air)"
set xlabel "t, ms"
set ylabel "0.92*rho*vel^2, kPa"
set xtic 100.0e-3
set xrange [25.0:26.0]
set yrange [0:500]
set key top left
plot "pitot.kpa" using 1:2 title "grid2, inv" with lines linestyle 1
EOF

