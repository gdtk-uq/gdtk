# scale-heat-flux.awk
# PJ 2017-06-06
BEGIN {
    rho_inf = 5.521e-3; # kg/m^^3
    p_inf = 66.43; # Pa
    u_inf = 1589.8; # m/s
    T_inf = 41.92; # K
    T_wall = 296.0; # K
    T_0 = 1300.0; # K
    specific_heat = 1004.5; # J/kg.K
    mu_ref = 1.716e-5; # Pa.s
    T_ref = 273.0; # K
    S = 111.0; # K
    mu_inf = mu_ref * (T_inf/T_ref)*sqrt(T_inf/T_ref) * (T_ref+S)/(T_inf+S);
    print("# x(m) St.Re^1/2");
}

NR > 1 {
    x = $1; q = $5;
    Rex = rho_inf*u_inf*x/mu_inf
    St = q/(rho_inf*u_inf*specific_heat*(T_0 - T_wall))
    print(x, St*sqrt(Rex))
}
