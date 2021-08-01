# compute_heat_transfer.awk
# Calculate heat transfer from cell-centre data
# that has been extracted from a flow soution.
#
BEGIN {
    R_sphere = 6.6e-3; # Radius of hemisphere
    T_wall = 296.0
    k_wall = 0.02624; # assumed ideal gas value
    q_stag = -1.0; # dummy value, to be replaced by stagnation-point value
    print "# theta_degrees  q_W/m^^2  q/q_stag"
}

NR > 1 {
    x = $1; y = $2; r = sqrt(x*x + y*y);
    dr = r - R_sphere
    theta = atan2(y, -x); # hemisphere centred on (0,0) nose in -x direction
    theta_degrees = theta *180.0/3.14159
    p = $9; T = $18; mu = $11; k = $12;
    q = k_wall * (T - T_wall)/dr
    if (q_stag < 0.0) { q_stag = q }
    print theta_degrees, q, q/q_stag
}
