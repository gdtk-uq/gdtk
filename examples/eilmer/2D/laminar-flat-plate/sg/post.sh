#! /bin/sh
# post.sh

# Extracts slices at x = 1.0 m. 
# This will allow us to examine the boundary layer profile.
e4shared --post --job=lam_flat_plate_march --tindx-plot=last --add-vars="mach,pitot" \
	 --output-file=profile-at-x1m-march.dat --extract-line="1.0,0.40,0.0,1.0,0.44,0.0,201"

awk -f ../integral-thicknesses.awk profile-at-x1m-march.dat > thicknesses-march.txt

gnuplot<<EOF
set term postscript eps enhanced 20
set output "u-velocity-comparison-march.eps"
set title "Velocity profile near wall at x = 1.0 m"
set xlabel "u/u_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./profile-at-x1m-march.dat" using ((\$6)/1390.0):((0.44-\$2)*1000) \
     title "e3march.py" with points pt 1, \
     "../boundary-layer-profile-clbl.data" using (\$3):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20
set output "density-comparison-march.eps"
set title "Density profile near wall at x = 1.0 m"
set xlabel "rho/rho_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./profile-at-x1m-march.dat" using ((\$5)/0.011835):((0.44-\$2)*1000) \
     title "e3march.py" with points pt 1, \
     "../boundary-layer-profile-clbl.data" using (\$5):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF
