# compute-shear.awk
# Invoke with the command line:
# $ awk -f compute-shear.awk bl.data > shear.data
#
# PJ, 2016-11-01
#
BEGIN {
    rho_inf = 0.1315 # kg/m**3
    velx_inf =  514.0 # m/s
    T_inf = 164.4 # K
    # Sutherland expression for viscosity
    mu_ref = 1.716e-5; T_ref = 273.0; S_mu = 111.0
    mu_inf = (T_inf/T_ref)*sqrt(T_inf/T_ref)*(T_ref+S_mu)/(T_inf+S_mu)*mu_ref
    print("# x(m)  tau_w(Pa)  Cf   y_plus")
}

NR > 1 {
    x = $1; y = $2; rho = $5; velx = $6; mu = $11; k = $12
    dvelxdy = (velx - 0.0) / y # Assuming that the wall is straight down at y=0
    tau_w = mu * dvelxdy    # wall shear stress
    Cf = tau_w / (0.5*rho_inf*velx_inf*velx_inf)
    if (tau_w > 0.0) abs_tau_w = tau_w; else abs_tau_w = -tau_w;
    vel_tau = sqrt(abs_tau_w / rho) # friction velocity
    y_plus = vel_tau * y * rho / mu
    Rex = rho_inf * velx_inf * x / mu_inf
    Cf_blasius = 0.664 / sqrt(Rex) 
    print(x, tau_w, Cf, Cf_blasius, y_plus)
}
