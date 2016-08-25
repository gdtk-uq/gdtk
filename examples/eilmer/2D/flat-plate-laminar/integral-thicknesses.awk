# integral-thicknesses.awk
# Compute integral thicknesses for the boundary layer as per AGARD report.
# PJ, October 2007 Coles turbulent case
#     October 2010 laminar flat plate case

BEGIN {
    Rgas = 287.1;     # J/kg/K
    g = 1.4;          # ratio of specific heats
    p_inf = 1.013e3;  # Pa
    T_inf = 300.0;    # degrees K
    # Subscript D indicates quantity at the outer edge of the boundary layer.
    rho_D = p_inf / (Rgas * T_inf);  # density, kg/m**3 
    u_D = 1390.0;     # velocity, m/s
    # Replace rho_D and u_D with values observed in the profile file.
    rho_D = 0.01183;  # kg/m**3
    u_D = 1389.3;     # m/s
    delta = 0.015;    # boundary-layer edge, distance from wall in metres
    y_offset = 0.44;  # y-position of plate, m
    rhou_D = rho_D * u_D;
    d_1 = 0.0;        # displacement thickness
    d_2 = 0.0;        # momentum-deficit
    integrating = 0; 
    print "rho_D=", rho_D, " u_D=", u_D
}

$1 != "#" {
    # We are approaching the wall from outside the boundary layer.
    y = y_offset - $2;
    rho = $5;
    u = $6;
    if ( integrating == 0 && y < delta ) integrating = 1;
    if ( integrating ) {
	print "y=", y, "rho=", rho, "u=", u;
	rhou = 0.5 * (rho * u + rho_old * u_old);
	dy = (y_old - y);
	d_1 += (1.0 - rhou/rhou_D) * dy;
	d_2 += rhou/rhou_D * (1.0 - rhou/rhou_D) * dy;
    }
    y_old = y;
    rho_old = rho;
    u_old = u;
}

END {
    print "d_1=", d_1;
    print "d_2=", d_2;
}
