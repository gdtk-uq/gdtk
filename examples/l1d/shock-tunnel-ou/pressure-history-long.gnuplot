set term postscript eps 20
set output "pressure-history-long.eps"
set title "Exp D shock tube pressure history"
set xlabel "t, ms"
set ylabel "p, kPa"
set xrange [0.0:25.0]
set yrange [0:1800.0]
plot "he-air/history-loc-0000.data" using ($1*1000):($5/1000) title "667" with lines linestyle 1, \
     "he-air/history-loc-0001.data" using ($1*1000):($5/1000) title "668" with lines linestyle 2
