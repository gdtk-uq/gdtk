# piston.tcl
# Functions to compute ideal piston speed and position
#
proc position { t } {
    # Constants for the particular problem
    set a_R 347.17
    set A   7.854e-5
    set P_R 1.0e5
    set g   1.4
    set m_p 0.001

    set t_bar [expr $P_R * $A * $t / ($m_p * $a_R)]
    set t1    [expr 2.0 / ($g + 1.0)]
    set t2    [expr 1.0 + $t_bar * ($g + 1.0) / 2.0]
    set t3    [expr 1.0 + $t_bar - pow($t2, $t1)]
    set x_bar [expr 2.0 / ($g - 1.0) * $t3]

    # Convert back to dimensional form
    set x [expr $x_bar * $m_p * $a_R * $a_R / ($P_R * $A)]
}; # end proc position

proc velocity { t } {
    # Constants for the particular problem
    set a_R 347.17
    set A   7.854e-5
    set P_R 1.0e5
    set g   1.4
    set m_p 0.001

    set t_bar [expr $P_R * $A * $t / ($m_p * $a_R)]
    set t1    [expr (1.0 - $g) / ($g + 1.0)]
    set t2    [expr 1.0 + $t_bar * ($g + 1.0) / 2.0]
    set t3    [expr 1.0 - pow($t2, $t1)]
    set V_bar [expr 2.0 / ($g - 1.0) * $t3]

    # Convert back to dimensional form.
    set V [expr $V_bar * $a_R]
}; # end proc velocity

# Now, compute a number of positions and velocities
set output_file [open ideal.dat w]
puts $output_file "# ideal piston acceleration"
puts $output_file "# t (sec), x (m), V (m/s)"
for {set t 0.0} {$t <= 50.01e-3} {set t [expr $t + 1.0e-3]} {
   set x [position $t]
   set V [velocity $t]
   puts $output_file "$t $x $V"
}; # end for
close $output_file
