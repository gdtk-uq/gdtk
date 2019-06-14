#!/bin/sh
# plot.sh
#
# Plot the profile data at approx. x=368mm so that we can compare
# the simulation results with data from Mabey's 1953 experiment.
# The dimensional data from the experiment is reconstructed from 
# the normalised data recorded in the AGARD report 223. 
#
# Wilson Chan, 6 Sep 2010
#

gnuplot<<EOF
set term postscript eps monochrome enhanced dashlength 2.0 linewidth 2.0 size 9cm,6cm
set output "mabey-cf.eps"
set title "c_f profile along wall"
set xlabel "x, m"
set ylabel "c_f"
set yrange [0:0.002]
set key left top
plot "./viscous_data_along_wall.dat" using (\$1):(\$6) title "Eilmer3" with lines lt 1 lw 2, \
     "./experimental_data_cf.txt" using (\$1):(\$2):(0.10*0.001) \
     title "Mabey" with yerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps monochrome enhanced dashlength 2.0 linewidth 2.0 size 9cm,6cm
set output "mabey-u.eps"
set title "Velocity profile near wall at x = 0.368 m"
set xlabel "u-velocity, m/s"
set ylabel "y, mm"
set key left top
set xrange [0:800]
set yrange [0:25]
plot "./mabey-x-368mm.dat" using (\$6):((0.16-(\$2))*1000) title "Eilmer3" with lines lt 1 lw 2, \
     "./experimental_data.txt" using (\$7*712.89):((\$2)*1000):(0.05*712.89) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps monochrome enhanced dashlength 2.0 linewidth 2.0 size 9cm,6cm
set output "mabey-temperature.eps"
set title "Temperature profile near wall at x = 0.368 m"
set xlabel "Temperature, K"
set ylabel "y, mm"
set key right top
set xrange [0:350]
set yrange [0:25]
plot "./mabey-x-368mm.dat" using (\$20):((0.16-(\$2))*1000) title "Eilmer3" with lines lt 1 lw 2, \
     "./experimental_data.txt" using (\$8*62.157):((\$2)*1000):(0.05*62.157) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps monochrome enhanced dashlength 2.0 linewidth 2.0 size 9cm,6cm
set output "mabey-pitot.eps"
set title "Pitot-pressure profile near wall at x = 0.368 m"
set xlabel "P_{pitot}, kPa"
set ylabel "y, mm"
set key left top
set xrange [0:100]
set yrange [0:25]
plot "./mabey-x-368mm.dat" using ((\$22)/1000):((0.16-(\$2))*1000) title "Eilmer3" with lines lt 1 lw 2, \
     "./experimental_data.txt" using ((\$3)*3.1634):(\$2*1000):(0.05*84.57) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps monochrome enhanced dashlength 2.0 linewidth 2.0 size 9cm,6cm
set output "mabey-mach.eps"
set title "Mach number profile near wall at x = 0.368 m"
set xlabel "Mach number"
set ylabel "y, mm"
set key left top
set xrange [0:5]
set yrange [0:25]
plot "./mabey-x-368mm.dat" using (\$21):((0.16-(\$2))*1000) title "Eilmer3" with lines lt 1 lw 2, \
     "./experimental_data.txt" using ((\$6)*4.5099):((\$2)*1000):(0.05*4.5099) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps monochrome enhanced dashlength 2.0 linewidth 2.0 size 9cm,6cm
set output "mabey-dimensionless.eps"
set title "Velocity profile near wall at x=0.368m"
set xlabel "y^+"
set ylabel "u^+"
set xrange [1:2000]
set yrange [0:35]
set logscale x
set key left top
plot "./dimensionless_data_CFD.dat" using (\$3):(\$5) title "Eilmer3"\
       with linespoint lt 1 lw 2 pt 1, \
     "./dimensionless_data_exp.dat" using (\$3):(\$5) title "Mabey"\
       with points pt 6
EOF

epstopdf mabey-cf.eps
epstopdf mabey-u.eps
epstopdf mabey-temperature.eps
epstopdf mabey-pitot.eps
epstopdf mabey-mach.eps
epstopdf mabey-dimensionless.eps

rm *.eps
