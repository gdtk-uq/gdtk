# cp.awk
# Scan a history file, picking out pressure and scaling it
# to compute coefficient of pressure.
#
# PJ, 2016-09-22
#
BEGIN {
    Rgas = 287.1; # J/kg.K
    p_inf = 95.84e3; # Pa
    T_inf = 1103; # K
    rho_inf = p_inf / (Rgas * T_inf)
    V_inf = 1000.0; # m/s
    q_inf = 0.5 * rho_inf * V_inf * V_inf
    print "# rho_inf=", rho_inf, " q_inf=", q_inf
    print "# t,ms cp"
}

$1 != "#" {
    t = $1; p = $10
    print t*1000.0, (p - p_inf)/q_inf
}

END {}
