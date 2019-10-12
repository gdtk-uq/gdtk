#!/bin/bash
# post-profile.sh
e4shared --post --job=vortex --tindx-plot=0 \
         --add-vars="mach,total-h,entropy" \
         --extract-line="0.0,5.0001,0.0,10.0,5.0001,0.0,300" \
         --output-file=profile-init.data
e4shared --post --job=vortex --tindx-plot=last \
         --add-vars="mach,total-h,entropy" \
         --extract-line="0.0,5.0001,0.0,10.0,5.0001,0.0,300" \
         --output-file=profile-final.data

gnuplot <<EOF
set term postscript eps enhanced 20
set output "vortex-T-profile.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "Profile through vortex"
set xlabel "x, m"
set ylabel "static temperature, K"
set yrange [100:350]
set key bottom right
plot "profile-init.data" using 1:20 title "initial" with lines, \
     "profile-final.data" using 1:20 title "final" with points pt 4
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "vortex-p-profile.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "Profile through vortex"
set xlabel "x, m"
set ylabel "pressure, kPa"
set yrange [0:120]
set key bottom right
plot "profile-init.data" using (\$1):(\$9/1000) title "initial" with lines, \
     "profile-final.data" using (\$1):(\$9/1000) title "final" with points pt 4
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "vortex-rho-profile.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "Profile through vortex"
set xlabel "x, m"
set ylabel "density, kg/m^3"
set yrange [0:1.5]
set key bottom right
plot "profile-init.data" using 1:5 title "initial" with lines, \
     "profile-final.data" using 1:5 title "final" with points pt 4
EOF
