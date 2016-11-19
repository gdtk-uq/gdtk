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
 * - [TODO] Taylor-Maccoll conical flow
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
 */

module idealgasflow;

import std.conv;
import std.math;
import std.string;
import nm.bracketing;
import nm.ridder;
import nm.linesearch;

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
    assert(approxEqual(T0_T(M,g), 2.152), "Total temperature fail");
    assert(approxEqual(p0_p(M,g), 14.620), "Total pressure fail");
    assert(approxEqual(r0_r(M,g), 6.7937), "Total density fail");
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
	throw new Error(text("r2_r1: subsonic Mach number: ", M1));
    } // end if
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
	throw new Error(text("r2_r1: subsonic Mach number: ", M1));
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
	throw new Error(text("u2_u1: subsonic Mach number: ", M1));
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
	throw new Error(text("p2_p1: subsonic Mach number: ", M1));
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
 * Nodimensional entropy change ds across a normal shock.
 * Input:
 *   M1: Mach number of incoming flow
 *   g: ratio of specific heats
 * Returns: ds/Cv
 */
double DS_Cv(double M1, double g=1.4)
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
    assert(approxEqual(m2_shock(M,g), 0.5774), "Mach number after shock fail");
    assert(approxEqual(p2_p1(M,g), 4.50), "Pressure ratio across shock fail");
    assert(approxEqual(T2_T1(M,g), 1.687), "Temperature ratio across shock fail");
    assert(approxEqual(r2_r1(M,g), 2.667), "Density ratio across shock fail");
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
    assert(approxEqual(T0_T0star(M,g), 0.7934), "Rayleigh-line total T0_T0star fail");
    assert(approxEqual(T_Tstar(M,g), 0.5289), "Rayleigh-line static T_Tstar fail");
    assert(approxEqual(p_pstar(M,g), 0.3636), "Rayleigh-line static p_pstar fail");
    assert(approxEqual(r_rstar(M,g), 0.6875), "Rayleigh-line static p_pstar fail");
    assert(approxEqual(M_Rayleigh(T0_T0star(M,g),g), M), "Rayleigh-line inverse fail");
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
	throw new Error(text("PM1 received a subsonic Mach number: ", M ));
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
	throw new Error("Given negative value for Prandtl-Meyer function.");
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
	throw new Error(text("PM2: iteration did not converge."));
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
	throw new Error(text("MachAngle: subsonic Mach number: ", M));
    } else {
	return asin(1.0/M);
    } // end if
} // end MachAngle()

unittest {
    double M = 2.4;
    double g = 1.4;
    assert(approxEqual(PM1(M,g), 0.6413), "Prandtl-Meyer fail");
    try {
	PM1(0.8,g);
    }
    catch (Error e) {
	auto found_subsonic = (indexOf(e.toString(), "subsonic") != -1);
	assert(found_subsonic, "Prandtl-Meyer failed to catch subsonic M");
    }
    double nu = 0.6413479572;
    assert(approxEqual(PM2(nu,g), 2.4), "Inverse Prandtl-Meyer fail");
    try {
	PM2(-0.5,g);
    } catch (Error e) {
	auto found_negative = (indexOf(e.toString(), "negative") != -1);
	assert(found_negative, "Prandtl-Meyer failed to catch negative nu");
    }
    assert(approxEqual(MachAngle(M), 0.430), "Mach angle fail");
    try {
	MachAngle(0.8);
    } catch (Error e) {
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
double beta_obl(double M1, double theta, double g=1.4,double tol=1.0e-6)
{
    if (M1 < 1.0) {
	throw new Error(text("beta_obl: subsonic Mach number: ", M1));
    } // end if
    int sign = 1; if(theta<0.0){sign=-1;}
    theta = fabs(theta);
    auto f_to_solve = delegate(double beta){return theta_obl(M1, beta, g) - theta;};
    //    
    double b1 = asin(1.0/M1); 
    double b2 = b1 * 1.05;
    int result_flag = bracket!f_to_solve(b1, b2);
    // [TODO] should test result_flag.
    return sign*solve!f_to_solve(b1,b2);
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
	throw new Error(text("beta_obl2: subsonic Mach number: ",M1));
    } // end if
    if (p2_p1 < 1.0) {
	throw new Error(text("beta_obl2: invalid p2_p1: ", p2_p1));
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
    if (M1 < 1.0) {
	throw new Error(text("theta_obl: subsonic Mach number: ", M1));
    } // end if
    int sign = 1; if(beta<0.0){sign=-1;}
    beta = fabs(beta);
    if (dtan_theta(M1,beta,g) < 0.0){
	throw new Error(text("theta_obl: shock is detached."));
    } // end if
    double M1n = M1 * sin(beta);
    double t1 = 2.0 / tan(beta) * (M1n^^2 - 1.0); 
    double t2 = M1^^2 * (g + cos(2.0 * beta)) + 2.0;
    return sign*atan(t1/t2);
} // end theta_obl()

/**
 * Computes derivative of tan of flow deflection angle [tan(theta)]
 * wrt shock angle [beta], for oblique shock
 * derivative is negative if the shock is detached
 * Input:
 *   M: pre-shock Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 */
double dtan_theta(double M1,double beta,double g=1.4)
{
    beta = fabs(beta);
    double dum1 = M1^^2*(cos(2.0*beta)+g)+2.0;
    double dum2 = M1^^2*sin(beta)^^2-1.0;
    return 4.0*M1^^2*cos(beta)^^2/dum1 + 4.0*M1^^2*sin(2.0*beta)*dum2/(dum1^^2)/tan(beta) - 2.0*dum2/dum1/(sin(beta)^^2);
} // end dtan_theta()

/**
 * Mach number after an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: Mach number in flow after the shock
 */
double M2_obl(double M1,double beta,double theta,double g=1.4)
{
    if (M1 < 1.0) {
	throw new Error(text("M2_obl: subsonic Mach number: ", M1));
    } // end if
    double M1n = M1 * fabs(sin(beta));
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
    if (M1 < 1.0) {
	throw new Error(text("MachAngle: subsonic Mach number: ", M1));
    }
    double M1n = M1 * fabs(sin(beta));
    return r2_r1(M1n,g);
} // end r2_r1_obl()

/**
 * normal velocity ratio u2/u1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 * Returns: u2/u1
 */
double u2_u1_obl(double M1, double beta,double g=1.4)
{
    if (M1 < 1.0) {
	throw new Error(text("u2_u1_obl: subsonic Mach number: ", M1));	
    }
    double M1n = M1 * fabs(sin(beta));
    return u2_u1(M1n,g);	
}

/**
 * absolute velocity ratio V2/V1 across an oblique shock.
 * Input:
 *   M1: upstream Mach number
 *   beta: shock angle with respect to initial flow direction (radians)
 *   theta: flow deflection angle (radians)
 * Returns: V2/V1
 */
double V2_V1_obl(double M1, double beta,double theta, double g=1.4)
{
    if (M1 < 1.0) {
	throw new Error(text("v2_V1_obl: subsonic Mach number: ", M1));	
    }
    double M2 = M2_obl(M1,beta,theta,g);
    double T2_T1 = T2_T1_obl(M1,beta,g);
    return M2/M1*sqrt(T2_T1);
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
    if (M1 < 1.0) {
	throw new Error(text("p2_p1_obl: subsonic Mach number: ", M1));
    }
    double M1n = M1 * fabs(sin(beta));
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
    if (M1 < 1.0) {
	throw new Error(text("T2_T1_obl: subsonic Mach number: ", M1));
    }
    double M1n = M1 * fabs(sin(beta));
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
    if (M1 < 1.0) {
	throw new Error(text("p02_p01_obl: subsonic Mach number: ", M1));
    }
    double M1n = M1 * fabs(sin(beta));
    return p02_p01(M1n,g);
}

unittest {
    double M = 2.0;
    double g = 1.4;
    double beta = 44.0 * PI / 180.0;
    double theta = 14.0 * PI / 180.0;
    assert(approxEqual(beta_obl(M, theta, g), beta), "Oblique shock, beta from theta fail");
    assert(approxEqual(beta_obl2(M, 2.088, g), beta), "Oblique shock, beta from p2_p1 fail");
    assert(approxEqual(theta_obl(M, beta, g), theta), "Oblique shock, theta from beta fail");
    assert(approxEqual(M2_obl(M, beta, theta, g), 1.482), "Oblique shock, M2 after shock fail");
    assert(approxEqual(T2_T1_obl(M, beta, g), 1.249), "Oblique shock, temperature ratio fail");
    assert(approxEqual(p2_p1_obl(M, beta, g), 2.088), "Oblique shock, pressure ratio fail");
    assert(approxEqual(r2_r1_obl(M, beta, g), 1.673), "Oblique shock, density ratio fail");
    assert(approxEqual(p02_p01_obl(M, beta, g), 0.9608), "Oblique shock, total-pressure fail");
    assert(approxEqual(u2_u1_obl(M, beta, g), 0.598), "Oblique shock, normal velocity ratio fail");
    assert(approxEqual(V2_V1_obl(M, beta,theta, g),0.828), "Oblique shock, absolute velocity ratio fail");
    try {
	beta_obl(M,40.*PI/180.,g);
    } catch (Error e) {
	auto found_detached = (indexOf(e.toString(), "detached") != -1);
	assert(found_detached, "beta_obl failed to catch detached shock");
    }
    try {
	T2_T1_obl(0.8,beta,g);
    } catch (Error e) {
	auto found_subsonic = (indexOf(e.toString(), "subsonic") != -1);
	assert(found_subsonic, "Oblique shock relations failed to catch subsonic Mach");
    }
}

//------------------------------------------------------------------------

/// [TODO] Taylor-Maccoll cone flow.

