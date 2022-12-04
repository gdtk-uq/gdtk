# integral-thicknesses.awk
# Compute integral thicknesses for the boundary layer as per AGARD report.
# PJ, October 2007 Coles turbulent case
#     October 2010 laminar flat plate case
#     November 2022 plate is now at the y=0 boundary for the chicken example.
#     2022-12-04 Tune for short-plate example.

BEGIN {
    Rgas = 287.1;     # J/kg/K
    g = 1.4;          # ratio of specific heats
    p_inf = 1.013e3;  # Pa
    T_inf = 300.0;    # degrees K
    # Subscript D indicates quantity at the outer edge of the boundary layer.
    rho_D = p_inf / (Rgas * T_inf);  # density, kg/m**3
    u_D = 1390.0;     # velocity, m/s
    # Replace rho_D and u_D with values observed in the profile file.
    rho_D = 0.01176;  # kg/m**3
    u_D = 1390.0;     # m/s
    delta = 0.015;    # boundary-layer edge, distance from wall in metres
    rhou_D = rho_D * u_D;
    d_1 = 0.0;        # displacement thickness
    d_2 = 0.0;        # momentum-deficit
    print "rho_D=", rho_D, " u_D=", u_D
}

NR > 1 {
    # Skip over the starting comment line.
    # We are progressing from the wall through the boundary layer.
    y = $2;
    rho = $7;
    u = $11;
    print "NR=", NR, "y=", y, "rho=", rho, "u=", u;
    if ( NR > 2 && y < delta ) {
	rhou = 0.5 * (rho * u + rho_old * u_old);
	dy = y - y_old;
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
