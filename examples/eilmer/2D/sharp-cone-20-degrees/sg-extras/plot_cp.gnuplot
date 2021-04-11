set term postscript eps enhanced 20
set output "cone20_cp.eps"
set style line 1 linetype 1 linewidth 3.0 
set title "History of pressure coefficient at 2/3 position on cone surface"
set xlabel "time, ms"
set ylabel "C_p"
set xtic 1.0
set ytic 0.1
set yrange [0:0.5]
set key bottom right
set arrow from 5.2,0.387 to 5.8,0.387 nohead linestyle 1
set label "Value from\nNACA 1135\nChart 6" at 5.0,0.3 right
set arrow from 5.0,0.3 to 5.5,0.387 head
plot "cone20_cp.dat" using 1:2 title "10x40+30x40", \
     "cone20_cp_hi-res.dat" using 1:2 title "80x320+240x320" with lines

