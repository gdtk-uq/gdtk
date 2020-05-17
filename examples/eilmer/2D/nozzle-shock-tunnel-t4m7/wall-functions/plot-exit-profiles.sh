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
set key at graph 0.6, graph 0.3
plot "nozzle-exit.data" u (\$2):(\$21) \
      every::1 w lines lt 2 lc 1 lw 2 t "wf, x_{tr}=100mm" 
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Pitot.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot pressure, kPa"
set xrange [0:0.14]
set yrange [0:110]
set key at graph 0.6, graph 0.3
plot "nozzle-exit.data" u (\$2):((\$22)/1000) \
      every::1 w lines lt 2 lc 1 lw 2 t "wf, x_{tr}=100mm"
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-staticP.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Static pressure, kPa"
set xrange [0:0.14]
set yrange [0:2.0]
set key at graph 0.6, graph 0.3
plot "nozzle-exit.data" u (\$2):((\$9)/1000) \
      every::1 w lines lt 2 lc 1 lw 2 t "wf, x_{tr}=100mm"
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-flow-divergence.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Flow divergence, degrees"
set xrange [0:0.14]
set yrange [-2:2]
set key at graph 0.6, graph 0.3
plot "nozzle-exit.data" u (\$2):(atan((\$7)/(\$6)*180/pi)) \
      every::1 w lines lt 2 lc 1 lw 2 t "wf, x_{tr}=100mm"
EOF

ps2pdf -dEPSCrop profile-Mach.eps
ps2pdf -dEPSCrop profile-Pitot.eps
ps2pdf -dEPSCrop profile-staticP.eps
ps2pdf -dEPSCrop profile-flow-divergence.eps
rm *eps
