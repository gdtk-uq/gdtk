#!/bin/sh
# plot-exit-profiles.sh
#

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Mach.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Mach number"
#set xrange [0:0.14]
set yrange [0:8]
set key at graph 0.8, graph 0.3
#unset key
plot "nozzle-exit-eq.data" \
      every::1 u (\$2):(\$21) \
      w lines dt 2 lc 7 lw 2 t "eq.",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):(\$26) \
      w lines lt 2 lc 2 lw 2 t "chem. noneq.",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):(\$30) \
      w lines lt 2 lc 3 lw 2 t "thermochem. noneq.",\

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Pitot.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot pressure, kPa"
set xrange [0:0.14]
set yrange [0:400]
set key at graph 0.6, graph 0.3
#unset key
plot "nozzle-exit-eq.data" \
      every::1 u (\$2):((\$22)/1000) \
      w lines dt 2 lc 7 lw 2 t "eq.",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$27)/1000) \
      w lines lt 2 lc 2 lw 2 t "chem. noneq.",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$31)/1000) \
      w lines lt 2 lc 3 lw 2 t "thermochem. noneq."

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Pitot-on-supply.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot-on-supply pressure"
set xrange [0:0.14]
set yrange [0:0.016]
set key at graph 0.65, graph 0.4
plot "nozzle-exit-eq.data" \
      every::1 u (\$2):((\$22)/19.3253e6) \
      w lines dt 2 lc 7 lw 2 t "eq.",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$27)/19.3253e6) \
      w lines lt 2 lc 2 lw 2 t "chem. noneq.",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$31)/19.3253e6) \
      w lines lt 2 lc 3 lw 2 t "thermochem. noneq.",\
     
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-staticP.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Static pressure, kPa"
set xrange [0:0.14]
set yrange [0:7.5]
set key at graph 0.6, graph 0.3
#unset key
plot "nozzle-exit-eq.data" \
      every::1 u (\$2):((\$9)/1000) \
      w lines dt 3 lc 1 lw 2 t "eq.",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$9)/1000) \
      w lines lt 2 lc 2 lw 2 t "chem. noneq.",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$9)/1000) \
      w lines lt 2 lc 3 lw 2 t "thermochem. noneq.",\

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-Ttr-Tve.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Temperature, K"
set xrange [0:0.14]
set yrange [0:1400]
set key at graph 0.7, graph 0.5
#unset key
plot "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$26)) \
      w lines lt 2 lc 1 lw 2 t "themochem. noneq., T_{tr}",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$28)) \
      w lines lt 2 lc 2 lw 2 t "themochem. noneq., T_{ve}",\

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-T.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Temperature, K"
set xrange [0:0.14]
set yrange [0:800]
set key at graph 0.6, graph 0.7
#unset key
plot "nozzle-exit-eq.data" \
      every::1 u (\$2):((\$20)) \
      w lines dt 2 lc 7 lw 2 t "eq.",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$25)) \
      w lines lt 2 lc 2 lw 2 t "chem. noneq.",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$26)) \
      w lines lt 2 lc 3 lw 2 t "themochem. noneq.",\

EOF

gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-massf.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Mass fractions"
set xrange [0:0.14]
set yrange [-0.6:0.8]
set key at graph 0.98, graph 0.35
set key samplen -6 maxrows 5 
set key font ",14"
#unset key
plot "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$18)) \
      w points pt 65 lc 1 ps 0.8 t "chem. noneq., N_2",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$19)) \
      w points pt 65 lc 2 ps 0.8 t "chem. noneq., O_2",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$20)) \
      w points pt 65 lc 3 ps 0.8 t "chem. noneq., N",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$21)) \
      w points pt 65 lc 4 ps 0.8 t "chem. noneq., O",\
     "nozzle-exit-chem.data" \
      every::1 u (\$2):((\$22)) \
      w points pt 65 lc 5 ps 0.8 t "chem. noneq., NO",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$19)) \
      w points pt 64 lc 7 ps 0.8 t "thermoch. noneq., N_2",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$20)) \
      w points pt 64 lc 4 ps 0.8 t "thermoch. noneq., O_2",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$21)) \
      w points pt 64 lc 2 ps 0.8 t "thermoch. noneq., N",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$22)) \
      w points pt 64 lc 3 ps 0.8 t "thermoch. noneq., O",\
     "nozzle-exit-noneq.data" \
      every::1 u (\$2):((\$23)) \
      w points pt 64 lc 1 ps 0.8 t "thermoch. noneq., NO",\

EOF


gnuplot<<EOF
set term postscript eps color enhanced size 7.5cm,5.625cm
set output "profile-flow-divergence.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Flow divergence, degrees"
set xrange [0:0.14]
set yrange [-2:2]
#set key at graph 0.8, graph 0.3
#unset key
plot  "nozzle-exit-eq.data" \
       every::1 u (\$2):(atan((\$7)/(\$6)*180/pi)) \
       w lines dt 2 lc 7 lw 2 t "eq.",\
      "nozzle-exit-chem.data" \
       every::1 u (\$2):(atan((\$7)/(\$6)*180/pi)) \
       w lines lt 2 lc 2 lw 2 t "chem. noneq.",\
      "nozzle-exit-noneq.data" \
       every::1 u (\$2):(atan((\$7)/(\$6)*180/pi)) \
       w lines lt 2 lc 3 lw 2 t "thermochem. noneq.",\

EOF

ps2pdf -dEPSCrop profile-Mach.eps
ps2pdf -dEPSCrop profile-Pitot.eps
ps2pdf -dEPSCrop profile-Pitot-on-supply.eps
ps2pdf -dEPSCrop profile-staticP.eps
ps2pdf -dEPSCrop profile-Ttr-Tve.eps
ps2pdf -dEPSCrop profile-T.eps
ps2pdf -dEPSCrop profile-massf.eps
ps2pdf -dEPSCrop profile-flow-divergence.eps
rm *eps
