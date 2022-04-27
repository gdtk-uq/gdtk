# lines-of-code.gnuplot
set term postscript eps 20
set output "lines-of-code.eps"
set xdata time
set timefmt "%d/%m/%Y"
set format x "%d/%m/%Y"
set xrange ["01/01/2009":"01/01/2023"]
set xtics rotate
set yrange [0:300000]
set key off
set xlabel "Date"
set ylabel "Lines of code, doc"
set title "Source code and documentation development"
set key top left
plot "./lines-of-code-eilmer3.data" using 1:2 title "e3code" with points pointtype 4 pointsize 2.0, \
"./lines-of-code-eilmer4.data" using 1:2 title "e4code" with points pointtype 5 pointsize 2.0, \
"./lines-of-doc-eilmer3.data" using 1:2 title "e3doc" with points pointtype 6 pointsize 2.0, \
"./lines-of-doc-eilmer4.data" using 1:2 title "e4doc" with points pointtype 7 pointsize 2.0

