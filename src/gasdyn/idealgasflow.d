/** idealgasflow.d
 * Quasi-one-dimensional steady-flow of an ideal gas.
 *
 * Contents:
 *
 * One-dimensional flows:
 *
 * - Isentropic flow relations.
 *   State zero (0) refers to the stagnation condition.
 *   State star is the sonic (throat) condition.
 * - 1D (Normal) Shock Relations
 *   State 1 is before the shock and state 2 after the shock.
 *   Velocities are in a shock-stationary frame.
 * - 1-D flow with heat addition (Rayleigh-line)
 *   State star is the (hypothetical) sonic condition.
 *
 * Two-dimensional flows:
 *
 * - Prandtl-Meyer functions
 * - Oblique-shock relations
 * - Taylor-Maccoll conical flow
 *
 * Authors: Peter Jacobs, Rowan Gollan, Momar Hughes
 *          Centre for Hypersonics, School of Engineering
 *          The University of Queensland
 *
 * Version:
 * 30-Sep-1994: Xplore version
 * 16-May-2004: Python equivalent adapted from the Xplore version.
 * Things happened... then.
 * 16-Feb-2015: D port of isentropic relations from the C code
 * 09-Sep-2015: Merge with the cfpylib/gasdyn/ideal_gas_flow.py module.
 * 2016-Nov-19: Port the conical flow functions
 */

module idealgasflow;

import std.conv;
import std.math;
import std.string;
import std.stdio;
import nm.bbla;
import nm.bracketing;
import nm.ridder;
import nm.linesearch;
import gasflowexception;

/// Isentropic flow

/**
 * Area ratio A/Astar for an isentropic, quasi-one-dimensional flow.
 * Input:
 *   M: Mach number at area A
 *   g: ratio of specific heats
 */
double A_Astar(double M, double g=1.4)
{
    double t1 = (g + 1.0) / (g - 1.0);
    double m2 = M^^2;
    double t2 = 1.0 / m2 * (2.0 / (g + 1.0) * (1.0 + (g - 1.0) * 0.5 * m2))^^t1;
    return sqrt(t2);
}

/**
 * Compute ratio of total (stagnation) temp to static temp
 * Input:
 *   M : Mach number
 *   g : ratio of specific heats
 */
double T0_T(double M, double g=1.4)
{
    return 1.0 + 0.5 * (g - 1.0) * M*M;
}

/**
 * Compute ratio of total (stagnation) pressure to static pressure.
 * Input :
 *   M : Mach number
 *   g : ratio of specific heats
 */
double p0_p(double M, double g=1.4)
{
    return pow(T0_T(M,g), g/(g-1.0));
}

/**
 * Compute ratio of stagnation density to free-stream density.
 * Input :
 *   M : Mach number
 *   g : ratio of specific heats
 */
double r0_r(double M, double g=1.4)
{
    return (T0_T(M, g))^^(1.0 / (g - 1.0));
}

unittest {
    double M = 2.4;
    double g = 1.4;
    assert(isClose(T0_T(M,g), 2.152), "Total temperature fail");
    assert(isClose(p0_p(M,g), 14.620), "Total pressure fail");
    assert(isClose(r0_r(M,g), 6.7937), "Total density fail");
}

//------------------------------------------------------------------

/// 1-D normal shock relations.

/**
 * Mach number M2 after a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 */
double m2_shock(double M1, double g=1.4)
{
    if (M1 < 1.0) {
        throw new GasFlowException(text("r2_r1: subsonic Mach number: ", M1));
    }
    double numer = 1.0 + (g - 1.0) * 0.5 * M1^^2;
    double denom = g * M1^^2 - (g - 1.0) * 0.5;
    return sqrt(numer / denom);
}

/**
 * Density ratio r2/r1 across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 */
double r2_r1(double M1, double g=1.4)
{
    if (M1 < 1.0) {
        throw new GasFlowException(text("r2_r1: subsonic Mach number: ", M1));
    }
    double numer = (g + 1.0) * M1^^2;
    double denom = 2.0 + (g - 1.0) * M1^^2;
    return numer / denom;
}

/**
 * Velocity ratio u2/u1 across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 */
double u2_u1(double M1, double g=1.4)
{
    if (M1 < 1.0) {
        throw new GasFlowException(text("u2_u1: subsonic Mach number: ", M1));
    }
    return 1.0 / r2_r1(M1, g);
}

/**
 * Static pressure ratio p2/p1 across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 */
double p2_p1(double M1, double g=1.4)
{
    if (M1 < 1.0) {
        throw new GasFlowException(text("p2_p1: subsonic Mach number: ", M1));
    }
    return 1.0 + 2.0 * g / (g + 1.0) * (M1^^2 - 1.0);
}

/**
 * Static temperature ratio T2/T1 across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 */
double T2_T1(double M1, double g=1.4)
{
    return p2_p1(M1, g) / r2_r1(M1, g);
}

/**
 * Stagnation pressure ratio p02/p01 across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 */
double p02_p01(double M1, double g=1.4)
{
    double t1 = (g + 1.0) / (2.0 * g * M1^^2 - (g - 1.0));
    double t2 = (g + 1.0) * M1^^2 / (2.0 + (g - 1.0) * M1^^2);
    return t1^^(1.0/(g-1.0)) * t2^^(g/(g-1.0));
}

/**
 * Nondimensional entropy change ds across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 * Returns: ds/Cv
 */
double ds_Cv(double M1, double g=1.4)
{
    double t1 = p2_p1(M1, g);
    double t2 = r2_r1(M1, g);
    return log(t1 * t2^^g);
}

/**
 * Pitot pressure for a specified Mach number free-stream flow.
 * Will shock the gas if required.
 * Input:
 *   p1: static pressure of incoming flow
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 * Returns: Pitot pressure (absolute)
 */
double pitot_p(double p1, double M1, double g=1.4)
{
    if (M1 > 1.0) {
        double p2 = p2_p1(M1,g)*p1;
        double M2 = m2_shock(M1, g);
        return p0_p(M2, g)*p2;
    } else {
        return p0_p(M1, g)*p1;
    } // end if
} // end pitot_p()

unittest {
    double M = 2.0;
    double g = 1.4;
    assert(isClose(m2_shock(M,g), 0.5774), "Mach number after shock fail");
    assert(isClose(p2_p1(M,g), 4.50), "Pressure ratio across shock fail");
    assert(isClose(T2_T1(M,g), 1.687), "Temperature ratio across shock fail");
    assert(isClose(r2_r1(M,g), 2.667), "Density ratio across shock fail");
}

//------------------------------------------------------------------

/// 1-D flow with heat addition (Rayleigh-line)

/**
 * Total temperature ratio for flow with heat addition.
 * Input:
 *   M: initial Mach number
 *   g: ratio of specific heats
 * Returns: T0/T0star where T0 is the total temperature of the initial flow
 *   and T0star is the total temperature that would be achieved
 *   if enough heat is added to get to sonic conditions.
 */
double T0_T0star(double M, double g=1.4)
{
    double term1 = (g + 1.0) * M^^2;
    double term2 = (1.0 + g * M^^2)^^2;
    double term3 = 2.0 + (g - 1.0) * M^^2;
    return term1 / term2 * term3;
}

/**
 * Computes M from Total Temperature ratio for Rayleigh-line flow.
 * Input:
 *   T0T0star: total temperature ratio (star indicating sonic conditions)
 *   g: ratio of specific heats
 * Returns: initial Mach number of flow
 *
 * Note that supersonic flow is assumed for the initial guess.
 */
double M_Rayleigh(double T0T0star, double g=1.4)
{
    auto f_to_solve = delegate(double m){return T0_T0star(m, g) - T0T0star;};
    double M1 = 1.0;
    double M2 = 2.0;
    int result_flag = bracket!f_to_solve(M1, M2);
    // [TODO] should test result_flag.
    return solve!f_to_solve(M1,M2);
}

/**
 * Static temperature ratio T/Tstar for Rayleigh-line flow.
 * Input:
 *   M: initial Mach number
 *   g: ratio of specific heats
 * Returns: T/Tstar where T is the static temperature of the initial flow
 *   and Tstar is the static temperature that would be achieved
 *   if enough heat is added to get to sonic conditions.
 */
double T_Tstar(double M, double g=1.4)
{
    return M^^2 * ( (1.0 + g) / (1.0 + g * M^^2) )^^2;
}

/**
 * Static pressure ratio p/pstar for Rayleigh-line flow.
 * Input:
 *   M: initial Mach number
 *   g: ratio of specific heats
 * Returns: p/pstar where p is the static pressure of the initial flow
 *   and pstar is the static pressure that would be achieved
 *   if enough heat is added to get to sonic conditions.
 */
double p_pstar(double M, double g=1.4)
{
    return (1.0 + g) / (1.0 + g * M^^2);
}

/**
 * Density ratio r/rstar for Rayleigh-line flow.
 * Input:
 *   M: initial Mach number
 *   g: ratio of specific heats
 * Returns: r/rstar where r is the density of the initial flow
 *   and rstar is the density that would be achieved
 *   if enough heat is added to get to sonic conditions.
 */
double r_rstar(double M, double g=1.4)
{
    return 1.0 / M^^2 / (1.0 + g) * (1.0 + g * M^^2);
}

/**
 * Stagnation pressure ratio p0/p0star for Rayleigh-line flow.
 * Input:
 *   M: initial Mach number
 *   g: ratio of specific heats
 * Returns: p0/p0star where p0 is the total pressure of the initial flow
 *   and p0star is the total pressure that would be achieved
 *   if enough heat is added to get to sonic conditions.
 */
double p0_p0star(double M, double g=1.4)
{
    double term1 = (2.0 + (g - 1.0) * M^^2) / (g + 1.0);
    double term2 = g / (g - 1.0);
    return (1.0 + g) / (1.0 + g * M^^2) * term1^^term2;
}

unittest {
    double M = 2.0;
    double g = 1.4;
    assert(isClose(T0_T0star(M,g), 0.7934), "Rayleigh-line total T0_T0star fail");
    assert(isClose(T_Tstar(M,g), 0.5289), "Rayleigh-line static T_Tstar fail");
    assert(isClose(p_pstar(M,g), 0.3636), "Rayleigh-line static p_pstar fail");
    assert(isClose(r_rstar(M,g), 0.6875), "Rayleigh-line static p_pstar fail");
    assert(isClose(M_Rayleigh(T0_T0star(M,g),g), M), "Rayleigh-line inverse fail");
}

//------------------------------------------------------------------

/// Isentropic flow turning (Prandtl-Meyer)

/**
 * Calculates the Prandtl-Meyer function (nu) given a Mach number,
 * using eqn 4.40 (p134) in Anderson's 1990 text
 * Modern Compressible Flow With Historic Perspective, 2nd edition
 *
 * Note, nu = 0.0 when M = 1.0
 * Inputs :
 * M : Mach number
 * g : ratio of specific heats
 */
double PM1(double M,double g=1.4)
{
    double temp1 = sqrt((g + 1.0) / (g - 1.0));
    double temp2 = M * M - 1.0;
    double nu = 0.0;
    if (temp2 < 0.0) {
        throw new GasFlowException(text("PM1 received a subsonic Mach number: ", M ));
    } else {
        temp2 = sqrt(temp2);
        nu = temp1 * atan(temp2 / temp1) - atan(temp2);
    } // end if
    return nu;
} // end PM1()

/**
 * Calculate the Mach number given a Prandtl-Meyer function value.
 * Note, nu = 0.0 when M = 1.0
 * Input:
 *   nu : Prandtl-Meyer function (radians)
 *   g  : ratio of specific heats
 *   tol: tolerance to satisfy, for terminating iteration.
 */
double PM2(double nu, double g=1.4, double tol=1.0e-6)
{
    if (nu < 0.0) {
        throw new GasFlowException("Given negative value for Prandtl-Meyer function.");
    } // end if
    //
    // Generate an initial guess and, it it is good, return it.
    double M_0  = MFromNu_approximate(nu, g);
    double nu_0 = PM1(M_0, g);
    double f_0  = nu - nu_0;
    if (fabs(f_0) < tol) { return M_0; }
    //
    // Make some improvements using the secant method.
    double M_1  = M_0 * 1.001;
    double nu_1 = PM1(M_1, g);
    double f_1  = nu - nu_1;
    int count = 0;
    do {
        ++count;
        // Improve the guess to Mach number.
        double slope = (f_1 - f_0) / (M_1 - M_0);
        double M_2 = M_1 - f_1 / slope;
        double nu_2 = PM1(M_2, g);
        double f_2 = nu - nu_2;
        // Prepare for next iteration.
        M_0 = M_1; nu_0 = nu_1; f_0 = f_1;
        M_1 = M_2; nu_1 = nu_2; f_1 = f_2;
    } while (fabs(f_1) > tol && count < 30);
    //
    if (fabs(f_1) > tol) {
        throw new GasFlowException(text("PM2: iteration did not converge."));
    }
    return M_1;
} // end PM2()

/**
 * Calculate the Mach number given a Prandtl-Meyer function using
 * the polynomial approximation from S.M. Fraser (1975)
 * Calculation of Mach number from given turning angle in supersonic flow,
 * The Aeronautical Journal, February 1975.
 * For large Mach numbers, use an asymptotic expansion.
 * Input:
 *   Nu : Prandtl-Meyer function (radians)
 *   g  : ratio of specific heats
 */
double MFromNu_approximate(double nu,double g=1.4)
{
    double M = 0.0;
    double nu_d = nu * 180.0 / PI;
    if (nu_d < 0.0) {
        M = 0.0;
    } else if (nu_d < 5.0) {
        M = 1.0 + 7.932e-2 * pow(nu_d, 2.0/3.0) *
            (1.0 + nu_d * (3.681e-2 + nu_d * (-5.99e-3 + nu_d * 5.719e-4)) );
    } else if ( nu_d < 65.0 ) {
        M = 1.071 + nu_d * (3.968e-2 + nu_d * (-4.615e-4 + nu_d *
                                               (1.513e-5 + nu_d * (-1.840e-7 + nu_d * 1.186e-9))));
    } else {
        // Use an asymptotic expansion for large M.
        double bigG = sqrt((g + 1.0) / (g - 1.0));
        double nu_max = PI_2 * (bigG - 1.0);
        M = (1.0 - bigG * bigG) / (nu - nu_max);
    } // end if
    return M;
} // end MFromNu_approx()

/**
 * Compute Mach angle from Mach number.
 * Input:
 *   M : Mach number (M >= 1.0)
 * Returns Mach angle in radians.
 */
double MachAngle(double M)
{
    if (M < 1.0) {
        throw new GasFlowException(text("MachAngle: subsonic Mach number: ", M));
    } else {
        return asin(1.0/M);
    } // end if
} // end MachAngle()

unittest {
    double M = 2.4;
    double g = 1.4;
    assert(isClose(PM1(M,g), 0.6413), "Prandtl-Meyer fail");
    try {
        PM1(0.8,g);
    } catch (GasFlowException e) {
        auto found_subsonic = (indexOf(e.toString(), "subsonic") != -1);
        assert(found_subsonic, "Prandtl-Meyer failed to catch subsonic M");
    }
    double nu = 0.6413479572;
    assert(isClose(PM2(nu,g), 2.4), "Inverse Prandtl-Meyer fail");
    try {
        PM2(-0.5,g);
    } catch (GasFlowException e) {
        auto found_negative = (indexOf(e.toString(), "negative") != -1);
        assert(found_negative, "Prandtl-Meyer failed to catch negative nu");
    }
    assert(isClose(MachAngle(M), 0.430), "Mach angle fail");
    try {
        MachAngle(0.8);
    } catch (GasFlowException e) {
        auto found_subsonic = (indexOf(e.toString(), "subsonic") != -1);
        assert(found_subsonic, "Mach angle failed to catch subsonic M");
    }
} // end unittest

// -----------------------------------------------------------------

/// Oblique shock relations

// beta is shock angle wrt on-coming stream direction (in radians)
// theta is flow deflection wrt on-coming stream (in radians)

/**
 * Oblique shock wave angle.
 * Input:
 *   M1: upstream Mach number
 *   theta: flow deflection angle (radians)
 * Returns: shock angle with respect to initial flow direction (radians)
 */
double beta_obl(double M1, double theta, double g=1.4, double tol=1.0e-6)
{
    if (M1 < 1.0) {
        throw new GasFlowException(text("beta_obl: subsonic Mach number: ", M1));
    }
    if (theta == 0.0) {
        throw new GasFlowException("beta_obl: zero deflection angle, so there can be no shock");
    }
    int sign_beta = (theta < 0.0) ? -1 : 1;
    theta = fabs(theta);
    auto f_to_solve = delegate(double beta){return theta_obl(M1, beta, g) - theta;};
    //
    double b1 = asin(1.0/M1); // the weakest shock is at the Mach angle
    if (theta < tol) return sign_beta*b1;
    double f1 = f_to_solve(b1);
    if (fabs(f1) < tol) return sign_beta*b1;
    double b2 = b1 * 1.05; // second guess slightly stronger
    double f2 = f_to_solve(b2);
    if (fabs(f2) < tol) return sign_beta*b2;
    int result_flag = bracket!f_to_solve(b1, b2);
    if (result_flag < 0) {
        string msg = text("beta_obl: could not bracket beta for M1=", M1, " theta=", theta);
        throw new GasFlowException(msg);
    }
    return sign_beta*solve!f_to_solve(b1,b2);
} // end beta_obl()

/**
 * Oblique shock wave angle.
 * Input:
 *   p2_p1: Static pressure ratio p2/p1 across an oblique shock
 *   theta: flow deflection angle (radians)
 * Returns: shock angle with respect to initial flow direction (radians)
 */
double beta_obl2(double M1, double p2_p1, double g=1.4)
{
    if (M1 < 1.0) {
        throw new GasFlowException(text("beta_obl2: subsonic Mach number: ",M1));
    } // end if
    if (p2_p1 < 1.0) {
        throw new GasFlowException(text("beta_obl2: invalid p2_p1: ", p2_p1));
    } // end if
    double dum1 = sqrt(((g+1.)*p2_p1+g-1.)/2./g);
    return asin(dum1/M1);
} // end beta_obl2()


/**
 * Compute the deflection angle given the shock wave angle.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: theta, flow deflection angle (radians)
 */
double theta_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (fabs(M1n - 1) < 1.0e-9) return 0.0;
    if (M1n < 1.0) {
        string msg = text("theta_obl: subsonic normal Mach number: ", M1n,
                          " for M1=", M1, " beta=", beta);
        throw new GasFlowException(msg);
    }
    int sign_beta = (beta < 0.0) ? -1 : 1;
    beta = fabs(beta);
    if (dtan_theta(M1,beta,g) < 0.0){
        string msg = text("theta_obl: assume shock is detached",
                          " for M1=", M1, " beta=", beta);
        throw new GasFlowException(msg);
    } // end if
    double t1 = 2.0 / tan(beta) * (M1n^^2 - 1.0);
    double t2 = M1^^2 * (g + cos(2.0 * beta)) + 2.0;
    return sign_beta*atan(t1/t2);
} // end theta_obl()

/**
 * Computes derivative of tan of flow deflection angle [tan(theta)]
 * wrt shock angle [beta], for oblique shock
 * derivative is negative if the shock is detached
 * Input:
 *   M: pre-shock Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 */
double dtan_theta(double M1, double beta, double g=1.4)
{
    beta = fabs(beta);
    double dum1 = M1^^2*(cos(2.0*beta)+g)+2.0;
    double dum2 = M1^^2*sin(beta)^^2-1.0;
    return 4.0*M1^^2*cos(beta)^^2/dum1
        + 4.0*M1^^2*sin(2.0*beta)*dum2/(dum1^^2)/tan(beta)
        - 2.0*dum2/dum1/(sin(beta)^^2);
} // end dtan_theta()

/**
 * Mach number after an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: Mach number in flow after the shock
 */
double M2_obl(double M1, double beta, double theta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("M2_obl: subsonic normal Mach number: ", M1n));
    }
    double numer = 1.0 + (g - 1.0) * 0.5 * M1n^^2;
    double denom = g * M1n^^2 - (g - 1.0) * 0.5;
    return sqrt(numer / denom / (sin(beta - theta))^^2 );
} // end M2_obl()

/**
 * Density ratio r2/r1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: r2/r1
 */
double r2_r1_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("r2_r1_obl: subsonic normal Mach number: ", M1n));
    }
    return r2_r1(M1n,g);
} // end r2_r1_obl()

/**
 * normal velocity ratio Vn2/Vn1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: Vn2/Vn1
 */
double Vn2_Vn1_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("Vn2_Vn1_obl: subsonic normal Mach number: ", M1n));
    }
    return u2_u1(M1n,g);
}

/**
 * Speed ratio V2/V1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: V2/V1
 */
double V2_V1_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("V2_V1_obl: subsonic normal Mach number: ", M1n));
    }
    return sqrt((sin(beta) / r2_r1_obl(M1, beta, g))^^2 + (cos(beta))^^2);
}

/**
 * Static pressure ratio p2/p1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: p2/p1
 */
double p2_p1_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("p2_p1_obl: subsonic normal Mach number: ", M1n));
    }
    return p2_p1(M1n,g);
}

/**
 * Static temperature ratio T2/T1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: T2/T1
 */
double T2_T1_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("T2_T1_obl: subsonic normal Mach number: ", M1n));
    }
    return T2_T1(M1n,g);
}

/**
 * Ratio of stagnation pressures p02/p01 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: p02/p01
 */
double p02_p01_obl(double M1, double beta, double g=1.4)
{
    double M1n = M1 * fabs(sin(beta));
    if (M1n < 1.0) {
        throw new GasFlowException(text("p02_p01_obl: subsonic normal Mach number: ", M1n));
    }
    return p02_p01(M1n,g);
}

unittest {
    double M = 2.0;
    double g = 1.4;
    double beta = 44.0 * PI / 180.0;
    double theta = 14.0 * PI / 180.0;
    assert(isClose(beta_obl(M, theta, g), beta), "Oblique shock, beta from theta fail");
    assert(isClose(beta_obl2(M, 2.088, g), beta), "Oblique shock, beta from p2_p1 fail");
    assert(isClose(theta_obl(M, beta, g), theta), "Oblique shock, theta from beta fail");
    assert(isClose(M2_obl(M, beta, theta, g), 1.482), "Oblique shock, M2 after shock fail");
    assert(isClose(T2_T1_obl(M, beta, g), 1.249), "Oblique shock, temperature ratio fail");
    assert(isClose(p2_p1_obl(M, beta, g), 2.088), "Oblique shock, pressure ratio fail");
    assert(isClose(r2_r1_obl(M, beta, g), 1.673), "Oblique shock, density ratio fail");
    assert(isClose(p02_p01_obl(M, beta, g), 0.9608), "Oblique shock, total-pressure fail");
    assert(isClose(Vn2_Vn1_obl(M, beta, g), 0.598), "Oblique shock, normal velocity ratio fail");
    assert(isClose(V2_V1_obl(M, beta, g),0.828), "Oblique shock, absolute velocity ratio fail");
    try {
        beta_obl(M,40.*PI/180.,g);
    } catch (Exception e) {
        // Since we are no longer sure what the Exception message will contain,
        // we will just assume that we have successfully caught the detached shock
        // calculation failure for good reason.
        assert(true, "beta_obl failed to catch detached shock");
    }
    try {
        T2_T1_obl(0.8,beta,g);
    } catch (Exception e) {
        assert(true, "Oblique shock relations failed to catch subsonic Mach");
    }
}

//------------------------------------------------------------------------
/// Taylor-Maccoll cone flow.

double[] taylor_maccoll_odes(double[] z, double theta, double g=1.4)
{
    /**
    The ODEs from the Taylor-Maccoll formulation.

    See PJ's workbook for Feb 2012 for details.
    We've packaged them formally so that we might one day use
    a more sophisticated ODE integrator requiring fewer steps.
    **/
    double rho=z[0]; double V_r=z[1]; double V_theta=z[2];
    double h=z[3]; double p=z[4];
    // Assemble linear system for determining the derivatives wrt theta.
    auto A = zeros!double(5,6); // Augmented matrix with rhs in last column.
    A[0,0] = V_theta; A[0,2] = rho; A[0,5] = -2.0*rho*V_r - rho*V_theta/tan(theta);
    A[1,1] = 1.0; A[1,5] = V_theta;
    A[2,1] = rho*V_r; A[2,2] = rho*V_theta; A[2,4] = 1.0;
    A[3,1] = V_r; A[3,2] = V_theta; A[3,3] = 1.0;
    A[4,0] = h*(g-1)/g; A[4,3] = rho*(g-1)/g; A[4,4] = -1.0;
    gaussJordanElimination!double(A);
    double[] dzdtheta =  A.getColumn(5);
    return dzdtheta;
}

double[] theta_cone(double V1, double p1, double T1, double beta,
                    double R=287.1, double g=1.4)
{
    /**
    Compute the cone-surface angle and conditions given the shock wave angle.

    Input:
    V1: speed of gas into shock
    p1: free-stream pressure
    T1: free-stream static temperature
    beta: shock wave angle wrt stream direction (in radians)
    R: gas constant
    g: ratio of specific heats

    Returns:
    an array of [theta_c, V_c, p_c, T_c]
    where
    theta_c is stream deflection angle in radians
    V_c is the cone-surface speed of gas in m/s
    p_c is the cone-surface pressure
    T_c is the cone-surface static temperature

    The computation starts with the oblique-shock jump and then integrates
    across theta until V_theta goes through zero.
    The cone surface corresponds to V_theta == 0.

    .. Versions:
       08-Mar-2012 : PJ ideal-gas version adapted from the cea2_gas_flow.py.
       24-Jun-2012 : RJG added checks to catch the limiting case when beta < mu
                   : and a linear interpolation when beta is only slightly larger
                   : than mu (1% larger)
       2016-11-19 ported to D by PJ
    **/
    // When beta is only this fraction larger than mu,
    // we'll apply a linear interpolation
    enum LINEAR_INTERP_SWITCH = 1.01;
    // Free-stream properties and gas model.
    double a1 = sqrt(g*R*T1);
    double M1 = V1 / a1;
    double C_p = R * g / (g-1);
    double h1 = C_p * T1;
    double rho1 = p1 / (R * T1);
    // Test beta in relation to the Mach angle, mu
    double mu = asin(1.0/M1);
    double beta2 = LINEAR_INTERP_SWITCH*mu;
    // Test for an infinitely-weak shock angle
    if (beta <= mu) { return [0.0, V1, p1, T1]; }
    //
    if (beta < beta2) {
        // It is difficult to integrate between the shock and cone body
        // when the shock angle is only slightly larger than the Mach
        // angle. In this instance, find the value at LINEAR_INTER_SWITCH*mu
        // and linearly interpolate to find the value at beta
        auto results = theta_cone(V1, p1, T1, beta2, R, g);
        double theta2=results[0]; double V2=results[1];
        double p2=results[2]; double T2=results[3];
        double frac = (beta - mu)/(beta2 - mu);
        double theta_c = frac*theta2;
        double V = (1.0 - frac)*V1 + frac*V2;
        double p = (1.0 - frac)*p1 + frac*p2;
        double T = (1.0 - frac)*T1 + frac*T2;
        return [theta_c, V, p, T];
    }
    //
    // Start at the point just downstream the oblique shock.
    double theta_s = theta_obl(M1, beta, g);
    double M2 = M2_obl(M1, beta, theta_s, g);
    if (M2 < 1.0) {
        throw new GasFlowException("gone subsonic at shock");
    }
    double rho2 = rho1 * r2_r1_obl(M1, beta, g);
    double V2 = V1 * V2_V1_obl(M1, beta, g);
    double p2 = p1 * p2_p1_obl(M1, beta, g);
    double T2 = T1 * T2_T1_obl(M1, beta, g);
    double h2 = T2 * C_p;
    //
    // Initial conditions for Taylor-Maccoll integration.
    double dtheta = -0.05 * PI / 180.0;  // fraction-of-a-degree steps
    double theta = beta;
    double V_r = V2 * cos(beta - theta_s);
    double V_theta = -V2 * sin(beta - theta_s);
    double rho = rho2; double p = p2; double h = h2;
    // For integrating across the shock layer, the state vector is:
    double[5] z = [rho2, V_r, V_theta, h2, p2];
    double[5] z_old; double theta_old= theta;
    while (V_theta < 0.0) {
        // Keep a copy for linear interpolation at the end.
        z_old[] = z[]; theta_old = theta;
        // Do the update using a low-order method (Euler) for the moment.
        double[5] dzdtheta = taylor_maccoll_odes(z, theta, g);
        z[] += dtheta * dzdtheta[]; theta += dtheta;
        rho=z[0]; V_r=z[1]; V_theta=z[2]; h=z[3]; p=z[4];
    }
    // At this point, V_theta should have crossed zero so
    // we can linearly-interpolate the cone-surface conditions.
    double V_theta_old = z_old[2];
    double frac = (0.0 - V_theta_old)/(V_theta - V_theta_old);
    double[5] z_c; z_c[] = z_old[]*(1.0-frac) + z[]*frac;
    double theta_c = theta_old*(1.0-frac) + theta*frac;
    // At the cone surface...
    rho=z_c[0]; V_r=z_c[1]; V_theta=z_c[2]; h=z_c[3]; p=z_c[4];
    double T = h / C_p;
    if (abs(V_theta) > 1.0e-6) {
        throw new GasFlowException("oops, did not correctly find cone surface");
    }
    return [theta_c, V_r, p, T];
} // end theta_cone()

double beta_cone(double V1, double p1, double T1, double theta,
                 double R=287.1, double g=1.4)
{
    /**
    Compute the conical shock wave angle given the cone-surface deflection angle.

    Input:
    V1: speed of gas into shock
    p1: free-stream pressure
    T1: free-stream static temperature
    theta: stream deflection angle (in radians)
    R: gas constant
    g: ratio of specific heats

    Returns:
    shock wave angle wrt incoming stream direction (in radians)

    This ideal-gas version adapted from the cea2_gas_flow version, 08-Mar-2012.
    and then ported to D, November 2016.
    **/
    // Free-stream properties and gas model.
    double a1 = sqrt(g*R*T1);
    double M1 = V1 / a1;
    double C_p = R * g / (g-1);
    double h1 = C_p * T1;
    double rho1 = p1 / (R * T1);
    // We guess values of beta until this error measure is (close to) zero.
    auto error_in_theta = delegate (double beta_guess)
    {
        double[] results = theta_cone(V1, p1, T1, beta_guess, R, g);
        double theta_guess = results[0]; // here, we only care about this value
        return theta_guess - theta;
    };
    // Initial guess at bracket for shock wave angle.
    double b1 = asin(1.0/M1) * 1.01; // to be stronger than a Mach wave
    double b2 = 1.05 * b1;
    bracket!error_in_theta(b1, b2, asin(1.0/M1), PI/2);
    return solve!error_in_theta(b1, b2, 1.0e-4);
} // end beta_cone()

double beta_cone2(double M1, double theta, double R=287.1, double g=1.4)
{
    /**
    Compute the conical shock wave angle,
    given the cone-surface deflection angle and free stream Mach number.

    :param M1: free stream Mach number
    :param theta: stream deflection angle (in radians)
    :param R: gas constant
    :param g: ratio of specific heats
    :returns: shock wave angle wrt incoming stream direction (in radians)

    .. This version basically delegates work to beta_cone().
    **/
    // Compute free stream velocity assuming unit value temperature
    double T1 = 1.0;
    double a1 = sqrt(g*R*T1);
    double V1 = M1*a1;
    // Set free stream pressure to unit value
    double p1 = 1.0;
    // Now ready to call beta_cone()
    return beta_cone(V1, p1, T1, theta, R, g);
} // end beta_cone2()

unittest {
    double M1 = 1.5; double p1 = 100.0e3; double T1 = 300.0;
    double R = 287.1; double g = 1.4; double rho1 = p1/(R*T1);
    double a1 = sqrt(g*R*T1);
    double V1 = M1 * a1;
    double beta = 49.0 * PI/180; // conical shock angle
    double[] results = theta_cone(V1, p1, T1, beta);
    double theta_c=results[0];  double V_c=results[1];
    double p_c=results[2]; double T_c=results[3];
    assert(isClose(theta_c*180.0/PI, 19.96), "cone flow deflection angle fail");
    assert(isClose((p_c - p1)/(0.5*rho1*V1*V1), 0.386), "cone pressure coefficient fail");
    assert(isClose(beta_cone(V1, p1, T1, 20.0*PI/180)*180/PI, 49.0),
           "cone shock angle from deflection, V, p and T fail");
    assert(isClose(beta_cone2(M1, 20.0*PI/180)*180/PI, 49.0),
           "cone shock angle from deflection and M fail");
}
