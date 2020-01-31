# plot_data.gnu

set term pngcairo enhanced size 2048,1024 # color 20
set output "profile_T_rho.png"
set title "Temperature distribution behind a normal shock\nv_{/Symbol \245}=800 m/s, T_{/Symbol \245}=210.0 K"
set xlabel "distance behind shock (cm)"
set ylabel "temperature, K"
set xrange [0:20]
set yrange [0:5e4]
# set logscale x
set mytics 2
set grid xtics ytics mxtics mytics
plot "co2-v_8.00.data" u ($1*100):($2) notitle w l lt 1 lw 2, \
     # "co2-v_8.00.data" u ($1*100):($4/1.539e-03) notitle w l lt 1 lw 2

X1 = 1.1368e-6
X2 = 5.03846e-5
X3 = 0.00171427
X4 = 1.17687e-5
X5 = 0.0558665
X6 = 0.148932
X7 = 0.0458555
X8 = 0.00562125
X9 = 2.48167e-5
X10 = 4.10967e-5
X11 = 2.19272e-5
X12 = 7.72231e-9
X13 = 2.16055e-6
X14 = 4.87942e-8
X15 = 1.56597e-5
X16 = 0.000103141
X17 = 0.000550328
X18 = 0.000671338
X_sum = X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18

# set arrow 1 from 25,1.083 to 45,1.083 nohead lw 4
set label 1 "CO2" at 20.05, X1/X_sum
# set arrow 2 from 25,2.91166e-3 to 45,2.91166e-3 nohead lw 4
set label 2 "N2" at 20.05, X2/X_sum
# set arrow 3 from 25,0.022515 to 45,0.022515 nohead lw 4
set label 3 "Ar" at 20.05, X3/X_sum
# set arrow 4 from 25,3.66628e-5 to 45,3.66628e-5 nohead lw 4
set label 4 "O2" at 20.05, X4/X_sum
# set arrow 5 from 25,0.027282 to 45,0.027282 nohead lw 4
set label 5 "CO" at 20.05, X5/X_sum
# set arrow 6 from 25,0.51851 to 45,0.51851 nohead lw 4
set label 6 "O" at 20.05, X6/X_sum
# set arrow 7 from 25,1.083 to 45,1.083 nohead lw 4
set label 7 "C" at 20.05, X7/X_sum
# set arrow 8 from 25,2.91166e-3 to 45,2.91166e-3 nohead lw 4
set label 8 "N" at 20.05, X8/X_sum
# set arrow 9 from 25,0.022515 to 45,0.022515 nohead lw 4
set label 9 "C2" at 20.05, X9/X_sum
# set arrow 10 from 25,3.66628e-5 to 45,3.66628e-5 nohead lw 4
set label 10 "CN" at 20.05, X10/X_sum
# set arrow 11 from 25,0.027282 to 45,0.027282 nohead lw 4
set label 11 "NO" at 20.05, X11/X_sum
# set arrow 12 from 25,0.51851 to 45,0.51851 nohead lw 4
set label 12 "NCO" at 1, 0.3e-6 #X12/X_sum
# set arrow 13 from 25,1.083 to 45,1.083 nohead lw 4
set label 13 "NO+" at 20.05, X13/X_sum
# set arrow 14 from 25,2.91166e-3 to 45,2.91166e-3 nohead lw 4
set label 14 "O2+" at 20.05, X14/X_sum
# set arrow 15 from 25,0.022515 to 45,0.022515 nohead lw 4
set label 15 "CO+" at 20.05, X15/X_sum
# set arrow 16 from 25,3.66628e-5 to 45,3.66628e-5 nohead lw 4
set label 16 "O+" at 20.05, X16/X_sum
# set arrow 17 from 25,0.027282 to 45,0.027282 nohead lw 4
set label 17 "C+" at 20.05, X17/X_sum
# set arrow 18 from 25,0.51851 to 45,0.51851 nohead lw 4
set label 18 "e-" at 20.05, X18/X_sum

set output "profile_moles.png"
set rmargin at screen 0.95
set title "Species concentrations behind a normal shock in air\nv_{/Symbol \245}=8000, T_{/Symbol \245}=210.0 K, p_{/Symbol \245}=133.32 Pa"
set xlabel "distance behind shock (cm)"
set ylabel "moles per original mole"
set logscale y
set mytics 10
set yrange [1e-7:1]
set xrange [0.0:20]
set grid nomxtics nomytics
set format y "10^{%L}"
set key bottom right
plot "co2-v_8.00.data" u ($1*100):($8/X_sum)  notitle w l lt 1 lw 2, \
     "co2-v_8.00.data" u ($1*100):($10/X_sum) notitle w l lt 2 lw 2, \
     "co2-v_8.00.data" u ($1*100):($12/X_sum) notitle w l lt 3 lw 2, \
     "co2-v_8.00.data" u ($1*100):($14/X_sum) notitle w l lt 4 lw 2, \
     "co2-v_8.00.data" u ($1*100):($16/X_sum) notitle w l lt 5 lw 2, \
     "co2-v_8.00.data" u ($1*100):($18/X_sum) notitle w l lt 6 lw 2, \
     "co2-v_8.00.data" u ($1*100):($20/X_sum) notitle w l lt 7 lw 2, \
     "co2-v_8.00.data" u ($1*100):($22/X_sum) notitle w l lt 8 lw 2, \
     "co2-v_8.00.data" u ($1*100):($24/X_sum) notitle w l lt 9 lw 2, \
     "co2-v_8.00.data" u ($1*100):($26/X_sum) notitle w l lt 10 lw 2, \
     "co2-v_8.00.data" u ($1*100):($28/X_sum) notitle w l lt 11 lw 2, \
     "co2-v_8.00.data" u ($1*100):($30/X_sum) notitle w l lt 12 lw 2, \
     "co2-v_8.00.data" u ($1*100):($32/X_sum) notitle w l lt 13 lw 2, \
     "co2-v_8.00.data" u ($1*100):($34/X_sum) notitle w l lt 14 lw 2, \
     "co2-v_8.00.data" u ($1*100):($36/X_sum) notitle w l lt 15 lw 2, \
     "co2-v_8.00.data" u ($1*100):($38/X_sum) notitle w l lt 16 lw 2, \
     "co2-v_8.00.data" u ($1*100):($40/X_sum) notitle w l lt 17 lw 2, \
     "co2-v_8.00.data" u ($1*100):($42/X_sum) notitle w l lt 18 lw 2
