# plot_data.gnu

set arrow 1 from 25,15.1597 to 45,15.1597 nohead lw 4
set arrow 2 from 25,10.5723 to 45,10.5723 nohead lw 4
set arrow 3 from 33,13.5 to 33,15.1597
set arrow 4 from 33,12.5 to 33,10.5723
set label 3 "CEA calculation" at 3.5,13
set term postscript eps enhanced color 20
set output "profile_T_rho.eps"
set title "Temperature and density distribution behind a normal shock\nM_{/Symbol \245}=12.28, T_{/Symbol \245}=300.0 K, p_{/Symbol \245}=133.32 Pa"
set xlabel "distance behind shock (cm)"
set ylabel "temperature, density ratio"
set xrange [0.01:100]
set yrange [6:24]
set logscale x
set mytics 2
set grid xtics ytics mxtics mytics
set label 1 "T/T_{{/Symbol \245}}" at 0.2,21
set label 2 "{/Symbol r}/{/Symbol r}_{{/Symbol \245}}" at 0.2,8
plot "air-Mach_12.3.data" u ($1*100):($2/300.0) t 'poshax' w l lt 1 lw 2, \
     "ref-data/marrone_fig4_T_ratio.g3data" u 1:2 t 'Marrone 1963' w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($4/1.539e-03) notitle w l lt 1 lw 2, \
     "ref-data/marrone_fig4_rho_ratio.g3data" u 1:2 notitle w l lt 2 lw 2


unset label 1
unset label 2
unset label 3
unset arrow 1
unset arrow 2

set arrow 1 from 25,1.083 to 45,1.083 nohead lw 4
set label 1 "N2" at 50,1.083
set arrow 2 from 25,2.91166e-3 to 45,2.91166e-3 nohead lw 4
set label 2 "O2" at 50,2.91166e-3
set arrow 3 from 25,0.022515 to 45,0.022515 nohead lw 4
set label 3 "NO" at 50,0.017
set arrow 4 from 25,3.66628e-5 to 45,3.66628e-5 nohead lw 4
set label 4 "NO+,e-" at 30,2.0e-5
set arrow 5 from 25,0.027282 to 45,0.027282 nohead lw 4
set label 5 "N" at 50,0.04
set arrow 6 from 25,0.51851 to 45,0.51851 nohead lw 4
set label 6 "O" at 50,0.51851

set output "profile_moles.eps"
set title "Species concentrations behind a normal shock in air\nM_{/Symbol \245}=12.28, T_{/Symbol \245}=300.0 K, p_{/Symbol \245}=133.32 Pa"
set xlabel "distance behind shock (cm)"
set ylabel "moles per original mole"
set logscale y
set mytics 10
set yrange [1.0e-7:2]
set xrange [0.01:100]
set grid nomxtics nomytics
set format y "10^{%L}"
set key bottom right
plot "ref-data/marrone_fig3_N2.g3data" using 1:2 title 'Marrone 1963' w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($8/4.092e-01)  title 'poshax' w l lt 1 lw 2, \
     "ref-data/marrone_fig3_O2.g3data" using 1:2 notitle w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($10/4.092e-01) notitle w l lt 1 lw 2, \
     "ref-data/marrone_fig3_NO.g3data" using 1:2 notitle w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($16/4.092e-01) notitle w l lt 1 lw 2, \
     "ref-data/marrone_fig3_O.g3data" using 1:2 notitle w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($14/4.092e-01) notitle w l lt 1 lw 2, \
     "ref-data/marrone_fig3_N.g3data" using 1:2 notitle w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($12/4.092e-01) notitle w l lt 1 lw 2, \
     "ref-data/marrone_fig3_NO+_e.g3data" using 1:2 notitle w l lt 2 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($18/4.092e-01) notitle w l lt 1 lw 2, \
     "air-Mach_12.3.data" u ($1*100):($20/4.092e-01) notitle w l lt 1 lw 2
