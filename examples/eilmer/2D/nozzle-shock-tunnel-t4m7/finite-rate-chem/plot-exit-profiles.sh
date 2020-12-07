#!/bin/sh
# plot-exit-profiles.sh
#

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Mach.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Mach number"
set xrange [0:0.14]
set yrange [0:8]
set parametric
const=7.0
set trange [0:0.2] #BL_thickness=6.42mm
set arrow from 0.1,0.0 to 0.1,8 nohead dt 3 lc 8 lw 1
set key at graph 0.6, graph 0.3
plot "nozzle-exit-chem.data" u (\$2):(\$26) \
      every::1 w lines lt 2 lc 8 lw 2 notitle,\
      t, const lc 8 dt 3 t "Target Mach number"   

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Pitot.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot pressure, kPa"
set xrange [0:0.14]
set yrange [0:270]
set key at graph 0.6, graph 0.3
set arrow from 0.1,0.0 to 0.1,270 nohead dt 3 lc 8 lw 1
plot "nozzle-exit-chem.data" u (\$2):((\$27)/1000) \
      every::1 w lines lt 2 lc 8 lw 2 notitle,\

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-staticP.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Static pressure, kPa"
set xrange [0:0.14]
set yrange [0:5.0]
set key at graph 0.6, graph 0.3
set arrow from 0.1,0.0 to 0.1,5.0 nohead dt 3 lc 8 lw 1
plot "nozzle-exit-chem.data" u (\$2):((\$9)/1000) \
      every::1 w lines lt 2 lc 8 lw 2 notitle,\

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-flow-divergence.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Flow divergence, degrees"
set xrange [0:0.14]
set yrange [-2:2]
set key at graph 0.65, graph 0.3
set parametric
const=0.0
set trange [0:0.2] 
set arrow from 0.1,-2 to 0.1,2 nohead dt 3 lc 8 lw 1
plot "nozzle-exit-chem.data" u (\$2):(atan((\$7)/(\$6)*180/pi)) \
      every::1 w lines lt 2 lc 8 lw 2 notitle,\
      t, const lc 8 dt 3 t "Target flow angularity"

EOF

ps2pdf -dEPSCrop profile-Mach.eps
ps2pdf -dEPSCrop profile-Pitot.eps
ps2pdf -dEPSCrop profile-staticP.eps
ps2pdf -dEPSCrop profile-flow-divergence.eps
rm *eps
