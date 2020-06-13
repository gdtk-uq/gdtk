# plot-history.sh

echo "Rescale the pressures and times to get convenient units and alignment."
awk '$1 != "#" {print $1*1.0e3, $5/1.0e6}' heg-495/history-loc-0007.data > pstag_sim.MPa
awk '$1 != "#" {print $1*1.0e3, $6}' heg-495/history-loc-0007.data > Tstag_sim.degreeK
awk '$1 != "#" {print $1*1.0e3, $5/1.0e6}' heg-495/history-loc-0000.data > pcomp_sim.MPa

echo "Nozzle supply pressure"
gnuplot <<EOF
set term postscript eps 20
set output "HEG-495-pstag.eps"
set title "HEG Shot 495: Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [162.5:167.5]
set yrange [0:90.0]
plot "pstag_sim.MPa" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF

echo "Nozzle supply temperature"
gnuplot <<EOF
set term postscript eps 20
set output "HEG-495-Tstag.eps"
set title "HEG Shot 495: Nozzle Supply Temperature"
set xlabel "t, ms"
set ylabel "T, K"
set xrange [162.5:167.5]
set yrange [0:12000.0]
plot "Tstag_sim.degreeK" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF

echo "Compression tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "HEG-495-pcomp.eps"
set title "HEG Shot 495: Compression Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [152:167]
set yrange [0:80.0]
set key top left
plot "pcomp_sim.MPa" using 1:2 title "L1d4 Simulation" with lines linestyle 1
EOF
