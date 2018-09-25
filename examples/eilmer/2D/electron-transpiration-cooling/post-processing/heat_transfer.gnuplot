# downstream_shock_tip.gnuplot
#
# Zachary J. Denman
# 2017-12-19

set term pdfcairo enhanced
set output "heat-transfer.pdf"

set title "Wall heat transfer" font "Helvetica,16"
set xlabel "s/R" font "Helvetica,16"
set ylabel "qw [W/cm^2]" font "Helvetica,16"
set ytics nomirror
# set yrange [0:500]

set key font "Helvetica,16"

plot "surface-properties.dat" using 1:3 title "" with lines
