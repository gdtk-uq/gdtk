#Gnuplot script file

clear
set terminal pngcairo
set output "divB_max_plot.png"
set xlabel "t/t_final"
set ylabel "divB_max"
set logscale y
set title "Maximum Magnetic Divergence of the Kelvin-Helmholtz Instability"
plot filename u 1:2 t "clean" ps 0 with linespoints, filename u 1:6 t "unclean" ps 0 with linespoints
set output "divB_RMS_plot.png"
set xlabel "t/t_final"
set ylabel "divB_RMS"
set logscale y
set title "RMS Magnetic Divergence of the Kelvin-Helmholtz Instability"
plot filename u 1:2 t "clean" ps 0 with linespoints, filename u 1:7 t "unclean" ps 0 with linespoints
