# ideal_gas_flow.py
"""
One-dimensional steady flow of an ideal gas.

Author:
   PA Jacobs
   Centre for Hypersonics, School of Engineering
   The University of Queensland

Versions:
   30-Sep-1994: Xplore version
   16-May-2004: Python equivalent adapted from the Xplore version.
   27-Feb-2012: Use relative import in cfpylib
   29-Dec-2019: Port to the Eilmer4 collection, update to Python3.
   17-jul-2022: Bring in Maciej's additions (theta_cone_flowfield).

Contents:

One-dimensional flows:
   * Isentropic flow relations.
     State zero (0) refers to the stagnation condition.
     State star is the sonic (throat) condition.
   * 1D (Normal) Shock Relations
     State 1 is before the shock and state 2 after the shock.
     Velocities are in a shock-stationary frame.
   * 1-D flow with heat addition (Rayleigh-line)
     State star is the (hypothetical) sonic condition.

Two-dimensional flows:
   * Prandtl-Meyer functions
   * Oblique-shock relations
   * Taylor-Maccoll conical flow
"""

from math import *
import numpy
from gdtk.numeric.zero_solvers import secant as solve
from gdtk.numeric.zero_solvers import newton as solve_newton

# ---------------------------------------------------------------
# Isentropic flow

def A_Astar(M, g=1.4):
    """
    Area ratio A/Astar for an isentropic, quasi-one-dimensional flow.

    M: Mach number at area A
    g: ratio of specific heats
    Returns: A/Astar
    """
    t1 = (g + 1.0) / (g - 1.0)
    m2 = M**2
    t2 = 1.0 / m2 * (2.0 / (g + 1.0) * (1.0 + (g - 1.0) * 0.5 * m2))**t1
    t2 = sqrt(t2)
    return t2

def T0_T(M, g=1.4):
    """
    Total to static temperature ratio for an adiabatic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: T0/T
    """
    return 1.0 + (g - 1.0) * 0.5 * M**2

def p0_p(M, g=1.4):
    """
    Total to static pressure ratio for an isentropic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: p0/p
    """
    return (T0_T(M, g))**( g / (g - 1.0) )

def r0_r(M, g=1.4):
    """
    Stagnation to free-stream density ratio for an isentropic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: r0/r
    """
    return (T0_T(M, g))**(1.0 / (g - 1.0))

# -----------------------------------------------------------------
# 1-D normal shock relations.

def m2_shock(M1, g=1.4):
    """
    Mach number M2 after a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: M2
    """
    numer = 1.0 + (g - 1.0) * 0.5 * M1**2
    denom = g * M1**2 - (g - 1.0) * 0.5
    return sqrt(numer / denom)

def r2_r1(M1, g=1.4):
    """
    Density ratio r2/r1 across a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: r2/r1
    """
    numer = (g + 1.0) * M1**2
    denom = 2.0 + (g - 1.0) *M1**2
    return numer / denom

def v2_v1(M1, g=1.4):
    """
    Velocity ratio v2/v1 across a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: v2/v1
    """
    return 1 / r2_r1(M1, g)

def p2_p1(M1, g=1.4):
    """
    Static pressure ratio p2/p1 across a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: p2/p1
    """
    return 1.0 + 2.0 * g / (g + 1.0) * (M1**2 - 1.0)

def T2_T1(M1, g=1.4):
    """
    Static temperature ratio T2/T1 across a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: T2/T1
    """
    return  p2_p1(M1, g) / r2_r1(M1, g)

def p02_p01(M1, g=1.4):
    """
    Stagnation pressure ratio p02/p01 across a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: p02/p01
    """
    t1 = (g + 1.0) / (2.0 * g * M1**2 - (g - 1.0))
    t2 = (g + 1.0) * M1**2 / (2.0 + (g - 1.0) * M1**2)
    return t1**(1.0/(g-1.0)) * t2**(g/(g-1.0))

def ds_Cv(M1, g=1.4):
    """
    Nondimensional entropy change ds/Cv across a normal shock.

    M1: Mach number of incoming flow
    g: ratio of specific heats Cp/Cv
    Returns: ds/Cv
    """
    t1 = p2_p1(M1, g)
    t2 = r2_r1(M1, g)
    return log(t1 * t2**g)

def pitot_p(p1, M1, g=1.4):
    """
    Pitot pressure for a specified Mach number free-stream flow.

    Will shock the gas if required.

    M1: Mach number of incoming flow
    g: ratio of specific heats
    Returns: Pitot pressure (absolute)
    """
    if M1 > 1.0:
        p2 = p2_p1(M1,g)*p1
        M2 = m2_shock(M1, g)
        return p0_p(M2, g)*p2
    else:
        return p0_p(M1, g)*p1


# -----------------------------------------------------------------
# 1-D flow with heat addition (Rayleigh-line)

def T0_T0star(M, g=1.4):
    """
    Total temperature ratio for flow with heat addition.

    M: initial Mach number
    g: ratio of specific heats
    Returns: T0/T0star where T0 is the total temperature of the initial flow
        and T0star is the total temperature that would be achieved
        if enough heat is added to get to sonic conditions.
    """
    term1 = (g + 1.0) * M**2
    term2 = (1.0 + g * M**2)**2
    term3 = 2.0 + (g - 1.0) * M**2
    return term1 / term2 * term3

def M_Rayleigh(T0T0star, g=1.4):
    """
    Computes M from Total Temperature ratio for Rayleigh-line flow.

    T0T0star: total temperature ratio (star indicating sonic conditions)
    g: ratio of specific heats
    Returns: initial Mach number of flow

    Note that supersonic flow is assumed for the initial guess.
    """
    def f_to_solve(m): return T0_T0star(m, g) - T0T0star
    return solve(f_to_solve, 2.5, 2.4)

def T_Tstar(M, g=1.4):
    """
    Static temperature ratio T/Tstar for Rayleigh-line flow.

    M: initial Mach number
    g: ratio of specific heats
    Returns: T/Tstar where T is the static temperature of the initial flow
      and Tstar is the static temperature that would be achieved
      if enough heat is added to get to sonic conditions.
    """
    return M**2 * ( (1.0 + g) / (1.0 + g * M**2) )**2

def p_pstar(M, g=1.4):
    """
    Static pressure ratio p/pstar for Rayleigh-line flow.

    M: initial Mach number
    g: ratio of specific heats
    Returns: p/pstar where p is the static pressure of the initial flow
      and pstar is the static pressure that would be achieved
      if enough heat is added to get to sonic conditions.
    """
    return (1.0 + g) / (1.0 + g * M**2)

def r_rstar(M, g=1.4):
    """
    Density ratio r/rstar for Rayleigh-line flow.

    M: initial Mach number
    g: ratio of specific heats
    Returns: r/rstar where r is the density of the initial flow
      and rstar is the density that would be achieved
      if enough heat is added to get to sonic conditions.
    """
    return 1.0 / M**2 / (1.0 + g) * (1.0 + g * M**2)

def p0_p0star(M, g=1.4):
    """
    Stagnation pressure ratio p0/p0star for Rayleigh-line flow.

    M: initial Mach number
    g: ratio of specific heats
    Returns: p0/p0star where p0 is the total pressure of the initial flow
      and p0star is the total pressure that would be achieved
      if enough heat is added to get to sonic conditions.
    """
    term1 = (2.0 + (g - 1.0) * M**2) / (g + 1.0)
    term2 = g / (g - 1.0)
    return (1.0 + g) / (1.0 + g * M**2) * term1**term2

# -----------------------------------------------------------------
# Prandtl-Meyer functions

def PM1(M, g=1.4):
    """
    Prandtl-Meyer function.

    M: Mach number
    g: ratio of specific heats
    Returns: Prandtl-Meyer function value (in radians)
    """
    if M > 1.0:
        t1 = M**2 - 1.0
        t2 = sqrt((g - 1.0) / (g + 1.0) * t1)
        t3 = sqrt(t1)
        t4 = sqrt((g + 1.0) / (g - 1.0))
        nu = t4 * atan(t2) - atan(t3)
    else:
        nu = 0.0
    return nu

def PM2(nu, g=1.4):
    """
    Inverse Prandtl-Meyer function.

    nu: Prandtl-Meyer function value (in radians)
    g: ratio of specific heats
    Returns: Mach number

    Solves the equation PM1(m, g) - nu = 0, assuming supersonic flow.
    """
    def f_to_solve(m): return PM1(m, g) - nu
    return solve(f_to_solve, 2.0, 2.1)

# -----------------------------------------------------------------
# Oblique shock relations
# beta is shock angle wrt on-coming stream direction (in radians)
# theta is flow deflection wrt on-coming stream (in radians)

def beta_obl(M1, theta, g=1.4, tol=1.0e-6):
    """
    Oblique shock wave angle.

    M1: upstream Mach number
    theta: flow deflection angle (radians)
    Returns: shock angle with respect to initial flow direction (radians)
    """
    if M1 < 1.0: raise Exception("M1 is subsonic")
    sign_beta = -1 if theta < 0.0 else 1
    theta = abs(theta)
    b1 = asin(1.0/M1);
    if theta < tol:
        # Small deflection will produce a very weak shock.
        return sign_beta * b1
    b2 = b1 * 1.05
    def f_to_solve(beta): return theta_obl(M1, beta, g) - theta
    f1 = f_to_solve(b1)
    if abs(f1) < tol: return sign_beta * b1
    f2 = f_to_solve(b2)
    if abs(f2) < tol: return sign_beta * b2
    return sign_beta * solve(f_to_solve, b1, b2, tol=tol)

def beta_obl_newt(M1, theta, g=1.4, tol=1.0e-6):
    """
    Oblique shock wave angle.
    calculation using Newton's methods

    M1: upstream Mach number
    theta: flow deflection angle (radians)
    Returns: shock angle with respect to initial flow direction (radians)
    """
    if M1 < 1.0: raise Exception("M1 is subsonic")
    sign_beta = -1 if theta < 0.0 else 1
    theta = abs(theta)
    b0 = asin(1.0/M1)
    if theta < tol:
        # Small deflection will produce a very weak shock.
        return sign_beta * b0
    def fun(beta):
        m1sb = M1 * abs(sin(beta))
        m1cb = M1 * abs(cos(beta))
        if m1cb < 1.0: raise Exception("Subsonic normal Mach number: %g" % m1cb)
        t1 = 2.0 / tan(beta) * (m1sb**2 - 1.0)
        t2 = 1/(M1**2 * (g + cos(2.0 * beta)) + 2.0)
        return atan(t1*t2) - theta
    def fun_dash(beta):
        m1sb = M1 * abs(sin(beta))
        t1 = 2.0 / tan(beta) * (m1sb**2 - 1.0)
        t2 = 1/(M1**2 * (g + cos(2.0 * beta)) + 2.0)
        x = t1*t2
        df_dx = 1/(1+x**2)
        dt1_db = 2*M1**2 * cos(2*beta) + 2/(sin(beta)**2)
        dt2_db = (M1**2 * (g + cos(2.0 * beta)) + 2.0)**-2 * M1**2 * 2 *sin(2*beta)
        dx_db = dt1_db*t2 + t1*dt2_db
        return df_dx*dx_db
    return sign_beta * solve_newton(fun, fun_dash, 0.9*b0, tol=tol)

def beta_obl2(M1, p2_p1, g=1.4):
    """
    Oblique shock wave angle.

    M1: upstream Mach number
    p2_p1: static pressure ratio p2/p1 across the oblique shock
    Returns: shock angle with respect to initial flow direction (radians)
    """
    if M1 < 1.0: raise Exception("M1 is subsonic: %g" % M1)
    if p2_p1 < 1.0: raise Exception("Invalid p2_p1: %g" % p2_p1)
    dum1 = sqrt(((g+1.0)*p2_p1+g-1.0)/2.0/g)
    return asin(dum1/M1)

def theta_obl(M1, beta, g=1.4):
    """
    Compute the deflection angle given the shock wave angle.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: theta, flow deflection angle (radians)
    """
    m1sb = M1 * abs(sin(beta))
    m1cb = M1 * abs(cos(beta))
    if m1sb < 1.0: raise Exception("Subsonic normal Mach number: %g" % m1sb)
    t1 = 2.0 / tan(beta) * (m1sb**2 - 1.0)
    t2 = M1**2 * (g + cos(2.0 * beta)) + 2.0
    theta = atan(t1/t2)
    return theta

def M2_obl(M1, beta, theta, g=1.4):
    """
    Mach number after an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: M2, Mach number in flow after the shock
    """
    m1sb = M1 * abs(sin(beta))
    m1cb = M1 * abs(cos(beta))
    if m1sb < 1.0: raise Exception("Subsonic normal Mach number: %g" % m1sb)
    numer = 1.0 + (g - 1.0) * 0.5 * m1sb**2
    denom = g * m1sb**2 - (g - 1.0) * 0.5
    m2 = sqrt(numer / denom / (sin(beta - theta))**2 )
    return m2

def r2_r1_obl(M1, beta, g=1.4):
    """
    Density ratio r2/r1 across an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: r2/r1
    """
    m1sb = M1 * abs(sin(beta))
    m1cb = M1 * abs(cos(beta))
    if m1sb < 1.0: raise Exception("Subsonic normal Mach number: %g" % m1sb)
    numer = (g + 1.0) * m1sb**2
    denom = 2.0 + (g - 1.0) * m1sb**2
    return numer / denom

def vn2_vn1_obl(M1, beta, g=1.4):
    """
    Normal velocity ratio vn1/vn2 across an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: v2/v1
    """
    return 1.0/r2_r1_obl(M1, beta, g=g)

def v2_v1_obl(M1, beta, g=1.4):
    """
    Flow-speed ratio v2/v1 across an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: v2/v1
    """
    return sqrt((sin(beta) / r2_r1_obl(M1, beta, g))**2 + (cos(beta))**2)

def p2_p1_obl(M1, beta, g=1.4):
    """
    Static pressure ratio p2/p1 across an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: p2/p1
    """
    m1sb = M1 * abs(sin(beta))
    m1cb = M1 * abs(cos(beta))
    if m1sb < 1.0: raise Exception("Subsonic normal Mach number: %g" % m1sb)
    return 1.0 + 2.0 * g / (g + 1.0) * (m1sb**2 - 1.0)

def T2_T1_obl(M1, beta, g=1.4):
    """
    Static temperature ratio T2/T1 across an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: T2/T1
    """
    return p2_p1_obl(M1, beta, g) / r2_r1_obl(M1, beta, g)

def p02_p01_obl(M1, beta, g=1.4):
    """
    Ratio of stagnation pressures p02/p01 across an oblique shock.

    M1: upstream Mach number
    beta: shock angle with respect to initial flow direction (radians)
    Returns: p02/p01
    """
    m1sb = M1 * abs(sin(beta))
    m1cb = M1 * abs(cos(beta))
    if m1sb < 1.0: raise Exception("Subsonic normal Mach number: %g" % m1sb)
    t1 = (g + 1.0) / (2.0 * g * m1sb**2 - (g - 1.0))
    t2 = (g + 1.0) * m1sb**2 / (2.0 + (g - 1.0) * m1sb**2)
    return t1**(1.0/(g-1.0)) * t2**(g/(g-1.0))

#------------------------------------------------------------------------
# Taylor-Maccoll cone flow.

def taylor_maccoll_odes(z, theta, g=1.4):
    """
    The ODEs from the Taylor-Maccoll formulation.

    See PJ's workbook for Feb 2012 for details.
    We've packaged them formally so that we might one day use
    a more sophisticated ODE integrator requiring fewer steps.
    """
    rho, V_r, V_theta, h, p = z
    # Assemble linear system for determining the derivatives wrt theta.
    A = numpy.zeros((5,5), float)
    b = numpy.zeros((5,), float)
    A[0,0] = V_theta; A[0,2] = rho; b[0] = -2.0*rho*V_r - rho*V_theta/tan(theta)
    A[1,1] = 1.0; b[1] = V_theta
    A[2,1] = rho*V_r; A[2,2] = rho*V_theta; A[2,4] = 1.0
    A[3,1] = V_r; A[3,2] = V_theta; A[3,3] = 1.0
    A[4,0] = h*(g-1)/g; A[4,3] = rho*(g-1)/g; A[4,4] = -1.0
    dzdtheta = numpy.linalg.solve(A,b)
    return dzdtheta

def theta_cone(V1, p1, T1, beta, R=287.1, g=1.4, dtheta=-1.0e-5):
    """
    Compute the cone-surface angle and conditions given the shock wave angle.

    V1: speed of gas into shock
    p1: free-stream pressure
    T1: free-stream static temperature
    beta: shock wave angle wrt stream direction (in radians)
    R: gas constant
    g: ratio of specific heats
    dtheta: angular increment for integration to the cone surface (in radians)

    Returns: tuple of theta_c, V_c, p_c, T_c:
      theta_c is stream deflection angle in radians
      V_c is the cone-surface speed of gas in m/s
      p_c is the cone-surface pressure
      T_c is the cone-surface static temperature

    The computation starts with the oblique-shock jump and then integrates
    across theta until V_theta goes through zero.
    The cone surface corresponds to V_theta == 0.

    Versions:
    08-Mar-2012: This ideal-gas version adapted from the cea2_gas_flow version.
    24-Jun-2012: RJG added checks to catch the limiting case when beta < mu and
      a linear interpolation when beta is only slightly larger than mu (1% larger)
    June 2022: Pass in the angular increment.
    """
    # When beta is only this fraction larger than mu,
    # we'll apply a linear interpolation
    LINEAR_INTERP_SWITCH = 1.01
    # Free-stream properties and gas model.
    a1 = sqrt(g*R*T1)
    M1 = V1 / a1
    C_p = R * g / (g-1)
    h1 = C_p * T1
    rho1 = p1 / (R * T1)
    # Test beta in relation to the Mach angle, mu
    mu = asin(1.0/M1)
    beta2 = LINEAR_INTERP_SWITCH*mu
    #print "beta= ", beta, "mu= ", mu, " beta2= ", beta2
    if beta <= mu:
        # An infinitely weak shock angle
        return 0.0, V1, p1, T1
    if beta < beta2:
        # It is difficult to integrate between the shock and cone body
        # when the shock angle is only slightly larger than the Mach
        # angle. In this instance, find the value at LINEAR_INTER_SWITCH*mu
        # and linearly interpolate to find the value at beta
        (theta2, V2, p2, T2) = theta_cone(V1, p1, T1, beta2, R, g)
        frac = (beta - mu)/(beta2 - mu)
        theta_c = frac*theta2
        V = (1.0 - frac)*V1 + frac*V2
        p = (1.0 - frac)*p1 + frac*p2
        T = (1.0 - frac)*T1 + frac*T2
        return theta_c, V, p, T
    #
    # Start at the point just downstream the oblique shock.
    theta_s = theta_obl(M1, beta, g)
    M2 = M2_obl(M1, beta, theta_s, g)
    assert M2 > 1.0
    rho2 = rho1 * r2_r1_obl(M1, beta, g)
    V2 = V1 * v2_v1_obl(M1, beta, g)
    p2 = p1 * p2_p1_obl(M1, beta, g)
    T2 = T1 * T2_T1_obl(M1, beta, g)
    h2 = T2 * C_p
    theta = beta
    V_r = V2 * cos(beta - theta_s)
    V_theta = -V2 * sin(beta - theta_s)
    #
    # For integrating across the shock layer, the state vector is:
    z = numpy.array([rho2, V_r, V_theta, h2, p2])
    #
    while V_theta < 0.0:
        # Keep a copy for linear interpolation at the end.
        z_old = z.copy(); theta_old = theta
        # Do the update using a low-order method (Euler) for the moment.
        dzdtheta = taylor_maccoll_odes(z, theta, g)
        z += dtheta * dzdtheta; theta += dtheta
        rho, V_r, V_theta, h, p = z
        if False: print("DEBUG theta=", theta, "V_r=", V_r, "V_theta=", V_theta)
    #
    # At this point, V_theta should have crossed zero so
    # we can linearly-interpolate the cone-surface conditions.
    V_theta_old = z_old[2]
    frac = (0.0 - V_theta_old)/(V_theta - V_theta_old)
    z_c = z_old*(1.0-frac) + z*frac
    theta_c = theta_old*(1.0-frac) + theta*frac
    # At the cone surface...
    rho, V_r, V_theta, h, p = z_c
    T = h / C_p
    assert abs(V_theta) < 1.0e-6
    #
    return theta_c, V_r, p, T

def beta_cone(V1, p1, T1, theta, R=287.1, g=1.4, tol=1.0e-8, dtheta=-1.0e-5):
    """
    Compute the conical shock wave angle given the cone-surface deflection angle.

    V1: speed of gas into shock
    p1: free-stream pressure
    T1: free-stream static temperature
    theta: stream deflection angle (in radians)
    R: gas constant
    g: ratio of specific heats
    tol: tolerance on the computed angle of the cone surface (in radians)
    dtheta: angular increment for integration to the cone surface (in radians)

    Returns: shock wave angle wrt incoming stream direction (in radians)

    This ideal-gas version adapted from the cea2_gas_flow version, 08-Mar-2012.
    """
    # Free-stream properties and gas model.
    a1 = sqrt(g*R*T1)
    M1 = V1 / a1
    C_p = R * g / (g-1)
    h1 = C_p * T1
    rho1 = p1 / (R * T1)
    # Initial guess
    M1 = V1 / a1
    b1 = asin(1.0 / M1) * 1.01 # to be stronger than a Mach wave
    b2 = b1 * 1.05
    def error_in_theta(beta_guess):
        theta_guess, V_c, p_c, T_c = theta_cone(V1, p1, T1, beta_guess, R, g, dtheta)
        return theta_guess - theta
    return solve(error_in_theta, b1, b2, tol=tol, limits=[asin(1.0/M1), pi/2.0])

def beta_cone2(M1, theta, R=287.1, g=1.4, tol=1.0e-8, dtheta=-1e-5):
    """
    Compute the conical shock wave angle given the cone-surface deflection angle and
    free stream Mach number.

    M1: free stream Mach number
    theta: stream deflection angle (in radians)
    R: gas constant
    g: ratio of specific heats
    tol: tolerance on the computed angle of the cone surface (in radians)
    dtheta: angular increment for integration to the cone surface (in radians)

    Returns: shock wave angle wrt incoming stream direction (in radians)

    This version basically delegates work to beta_cone().
    """
    # Compute free stream velocity assuming unit value temperature
    T1 = 1.0
    a1 = sqrt(g*R*T1)
    V1 = M1*a1
    # Set free stream pressure to unit value
    p1 = 1.0
    # Now ready to call beta_cone()
    return beta_cone(V1, p1, T1, theta, R, g, tol, dtheta)

def theta_cone_flowfield(V1, p1, T1, beta, theta_cone, rays_num,
                         R=287.1, g=1.4, dtheta=-1.0e-5):
    """
    Returns the flowfield properties for a collection of rays
    through the conical shock layer.

    Maciej Grybko, University of Southern Queensland, 2022
    """
    # Free-stream properties and gas model.
    a1 = sqrt(g*R*T1)
    M1 = V1 / a1
    C_p = R * g / (g-1)
    h1 = C_p * T1
    rho1 = p1 / (R * T1)
    #
    # Start at the point just downstream the oblique shock.
    theta_s = theta_obl(M1, beta, g)
    M2 = M2_obl(M1, beta, theta_s, g)
    assert M2 > 1.0
    rho2 = rho1 * r2_r1_obl(M1, beta, g)
    V2 = V1 * v2_v1_obl(M1, beta, g)
    p2 = p1 * p2_p1_obl(M1, beta, g)
    T2 = T1 * T2_T1_obl(M1, beta, g)
    h2 = T2 * C_p
    theta = beta
    V_r = V2 * cos(beta - theta_s)
    V_theta = -V2 * sin(beta - theta_s)
    #
    # For integrating across the shock layer, the state vector is:
    z = numpy.array([rho2, V_r, V_theta, h2, p2])
    #
    # Save Mach number and flow direction for a number of rays
    M = [M2]
    flow_dir = [beta + atan(V_theta/V_r)] # flow direction
    theta_vec = [beta]                    # polar coordinate
    mu = [asin(1/M2)]                     # mach wave angle
    #
    S = beta - theta_cone                  # sum of all theta increments
    n = rays_num - 2
    q = 1 + 2.0/rays_num + 500/rays_num**2 # multiplier (for non-uniform theta)
    i = 0
    theta_series = beta - S*(1-q)/(1-q**n)*q**i
    #
    while V_theta < 0.0:
        # Keep a copy for linear interpolation at the end.
        z_old = z.copy(); theta_old = theta
        # Do the update using a low-order method (Euler) for the moment.
        dzdtheta = taylor_maccoll_odes(z, theta, g)
        z += dtheta * dzdtheta; theta += dtheta
        rho, V_r, V_theta, h, p = z
        if False: print("DEBUG theta=", theta, "V_r=", V_r, "V_theta=", V_theta)
        #
        if theta < theta_series: # save flow properties for desired thetas
            i += 1
            theta_series -= S*(1-q)/(1-q**n)*q**i
            V = sqrt(V_r**2 + V_theta**2)
            T = h / C_p
            a = sqrt(g*R*T)
            M.append(V/a)
            flow_dir.append(theta + atan(V_theta/V_r))
            theta_vec.append(theta)
            mu.append(asin(1/M[-1]))
    #
    # At this point, V_theta should have crossed zero so
    # we can linearly-interpolate the cone-surface conditions.
    V_theta_old = z_old[2]
    frac = (0.0 - V_theta_old)/(V_theta - V_theta_old)
    z_c = z_old*(1.0-frac) + z*frac
    theta_c = theta_old*(1.0-frac) + theta*frac
    # At the cone surface...
    rho, V_r, V_theta, h, p = z_c
    V = sqrt(V_r**2 + V_theta**2)
    T = h / C_p
    a = sqrt(g*R*T)
    M.append(V/a)
    flow_dir.append(theta + atan(V_theta/V_r))
    theta_vec.append(theta_c)
    mu.append(asin(1/M[-1]))
    #
    assert abs(V_theta) < 1.0e-6
    #
    return M, flow_dir, theta_vec, mu

# ------------------------------- end ----------------------------------
