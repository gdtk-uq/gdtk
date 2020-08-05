#!/bin/sh
# plot-exit-profiles.sh
#

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Mach.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Mach number"
set xrange [0:0.05]
#set yrange [0:5]
set key at graph 0.6, graph 0.3
plot "../contour-22-fine-grid/nozzle-exit.data" u (\$2):(\$21) \
       t "Eilmer3" w lines lt 1 lc 1 lw 2, \
     "nozzle-exit-eilmer4.data" u (\$2):(\$19) \
       t "Eilmer4" w points pt 2 lc 3 lw 2
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Pitot.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot pressure, kPa"
set xrange [0:0.05]
set yrange [0:1000]
set key at graph 0.6, graph 0.3
plot "../contour-22-fine-grid/nozzle-exit.data" u (\$2):((\$22)/1000) \
       t "Eilmer3" w lines lt 1 lc 1 lw 2, \
     "nozzle-exit-eilmer4.data" u (\$2):((\$20)/1000) \
       t "Eilmer4" w points pt 2 lc 3 lw 2
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-staticP.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Static pressure, kPa"
set xrange [0:0.05]
set yrange [0:50]
set key at graph 0.6, graph 0.3
plot "../contour-22-fine-grid/nozzle-exit.data" u (\$2):((\$9)/1000) \
       t "Eilmer3" w lines lt 1 lc 1 lw 2, \
     "nozzle-exit-eilmer4.data" u (\$2):((\$9)/1000) \
       t "Eilmer4" w points pt 2 lc 3 lw 2
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-flow-divergence.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Flow divergence, degrees"
set xrange [0:0.05]
set yrange [-2:2]
set key at graph 0.6, graph 0.3
plot "../contour-22-fine-grid/nozzle-exit.data" u (\$2):(atan((\$7)/(\$6)*180/pi)) \
       t "Eilmer3" w lines lt 1 lc 1 lw 2, \
     "nozzle-exit-eilmer4.data" u (\$2):(atan((\$7)/(\$6)*180/pi)) \
       t "Eilmer4" w points pt 2 lc 3 lw 2
EOF

ps2pdf -dEPSCrop profile-Mach.eps
ps2pdf -dEPSCrop profile-Pitot.eps
ps2pdf -dEPSCrop profile-staticP.eps
ps2pdf -dEPSCrop profile-flow-divergence.eps
rm *eps
