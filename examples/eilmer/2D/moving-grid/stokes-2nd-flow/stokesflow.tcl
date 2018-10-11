# stokesflow.tcl
# Functions to compute velocity profile for stokes 2nd flow
#

proc velocity { angle y } {
    # Constants for the particular problem
    set u_max 50.   
    set omega 300.  
    set mu 1.8e-5   
    set rho 1.225   

    set nu    [expr $mu / $rho ]
    set a     [expr -1.0 * $y * sqrt($omega / (2.0*$nu) )]
    set u     [expr $u_max * exp($a) * cos( $angle + $a )]
}; # end proc velocity

# Now, compute a number of positions and velocities
set output_file [open ideal.dat w]
puts $output_file "# ideal velocity profile"
puts $output_file "# y (m), u1 (m/s), u2 (m/s), u3 (m/s), u4 (m/s)"
set pi [expr atan(1) * 4 ]
for {set y 0.0} {$y <= 1e-3} {set y [expr $y + 1.0e-5]} {
   set u1 [velocity [expr 0*$pi/2] $y]
   set u2 [velocity [expr 1*$pi/2] $y]
   set u3 [velocity [expr 2*$pi/2] $y]
   set u4 [velocity [expr 3*$pi/2] $y]
   puts $output_file "$y $u1 $u2 $u3 $u4"
}; # end for t
close $output_file
