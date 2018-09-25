# downstream_shock_tip.gnuplot
#
# Zachary J. Denman
# 2017-12-19
 
set term pdfcairo enhanced
set output "wall-temperature.pdf"

set title "Surface temperature" font "Helvetica,16"
set xlabel "Distance from stagnation point relative to nose radius" font "Helvetica,16" 
set ylabel "Temperature, K" font "Helvetica,16"
#set ytics nomirror

# set yrange [1200:3000]

set key font "Helvetica,16"

plot "surface-properties.dat" using 1:3 title "" with lines