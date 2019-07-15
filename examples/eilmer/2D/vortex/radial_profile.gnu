# radial_profile.gnu

set term postscript eps enhanced 20
set output "radial_profile_p.eps"
set title "Inviscid Vortex"
set xlabel "r/r_i"
set ylabel "p/p_i"
# set yrange [1.0:4.5]
set key bottom right
plot "radial_profile_0.dat" using 1:2 title "exact" with lines, \
     "radial_profile_45.dat" using 1:2 title "45 degrees", \
     "radial_profile_90.dat" using 1:2 title "exit plane"

set term postscript eps enhanced 20
set output "radial_profile_vel.eps"
set title "Inviscid Vortex"
set xlabel "r/r_i"
set ylabel "vel/vel_i"
# set yrange [0.7:1.0]
set key top right
plot "radial_profile_0.dat" using 1:3 title "exact" with lines, \
     "radial_profile_45.dat" using 1:3 title "45 degrees", \
     "radial_profile_90.dat" using 1:3 title "90 degrees"

set term postscript eps enhanced 20
set output "radial_profile_T.eps"
set title "Inviscid Vortex"
set xlabel "r/r_i"
set ylabel "T/T_i"
# set yrange [1.0:1.7]
set key bottom right
plot "radial_profile_0.dat" using 1:5 title "exact" with lines, \
     "radial_profile_45.dat" using 1:5 title "45 degrees", \
     "radial_profile_90.dat" using 1:5 title "90 degrees"
