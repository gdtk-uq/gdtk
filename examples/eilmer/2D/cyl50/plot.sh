# cyl50_plot.sh
# Plot the profiles of temperature and velocity toward the end of the cylinder.

# Extract the solution data and reformat.
e4shared --post --job=cyl50 --tindx-plot=2  --vtk-xml --add-vars="mach,pitot,total-p,total-h"

# Extract the profile near the downstream end of the cylinder.
e4shared --post --job=cyl50 --tindx-plot=2 --slice-list="1,23,:,0;2,23,:,0" --output-file=profile-i23.data

gnuplot <<EOF
set term postscript eps enhanced 20
set output "cyl50_profile_T.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "cyl50: Profile at x=0.917m"
set ylabel "y-R, mm"
set key top right
set xlabel "Temperature, K"
set yrange [0:25]
set xrange [200:270]
set style line 1 linetype 1 linewidth 4.0
plot "profile-i23.data" using (\$20):(\$2-0.005)*1000 title "50x50 grid" with points pt 4, \
     "cyl50_dimensional.dat" using (\$2):(\$1-0.005)*1000.0 title "spectral" with lines ls 1
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "cyl50_profile_ux.eps"
set style line 1 linetype 1 linewidth 1.0 
set title "cyl50: Profile at x=0.917m"
set ylabel "y-R, mm"
set xlabel "ux, m/s"
set key top left
set yrange [0:25]
set xrange [0:700]
set style line 1 linetype 1 linewidth 4.0
plot "profile-i23.data" using (\$6):(\$2-0.005)*1000 title "50x50 grid" with points pt 4, \
     "cyl50_dimensional.dat" using (\$3):(\$1-0.005)*1000.0 title "spectral" with lines ls 1
EOF
