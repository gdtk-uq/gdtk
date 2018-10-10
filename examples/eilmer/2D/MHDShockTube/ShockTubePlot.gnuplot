set terminal pngcairo
set output "ShockTubeCapture_with_waves.png"
set xlabel "x"
set ylabel "{/Symbol r}" enhanced
set xrange [0: 1]
set title "Sod's Shock Tube"
plot 'ShockTubeNumerical.dat' u 1:2 notitle ps 0 with lines, 'ShockTubeAnalytical.dat' u 1:4 notitle ps 0 with lines

