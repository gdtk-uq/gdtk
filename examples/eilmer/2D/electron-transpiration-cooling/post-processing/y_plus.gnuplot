# downstream_shock_tip.gnuplot
#
# Zachary J. Denman
# 2017-12-19

set term pdfcairo enhanced
set output "y-plus.pdf"

set title "y+ along cylinder-wedge"
set xlabel "s/R"
set ylabel "y+"
set ytics nomirror
plot "surface-properties.dat" using 1:4 title "" with lines
