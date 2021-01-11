set term pdfcairo size 15cm,10cm font "courier,15"
set output "n90_compare_T_stag_line.pdf"
set title "Stagnation line to 90mm cylinder in 5km/s nitrogen flow"
set xlabel "x, m"
set ylabel "temperature, K"
set xrange [-0.015:0.0]
# set xtics 0.005
set yrange [0:12000]
# set ytics 5000.0
set key right top font ",13"
plot "stag_sebo.dat" using 1:7 title "CEVCATS-N Reference" with lines, \
     "n90-stag-prof-Park.data" using 1:21 title "Eilmer: Park (1993) rates" w linesp pt 1, \
     "n90-stag-prof-KH.data" using 1:21 title "Eilmer: Kewley Hornung (1974) rates" w linesp pt 4 ps 0.7
     

