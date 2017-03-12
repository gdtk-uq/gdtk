/** gasflow.d
 * Quasi-one-dimensional flow calculations using GasState and GasModel objects.
 *
 * These 1-D flow functions were first abstracted from estcj.py
 * and then made a little more generic so that they could be used
 * to support other gas flow programs such as nenzfr.
 *
 * Contents:
 *   Normal shocks
 *   Isentropic Expansions
 *   Oblique shocks
 *   Conical shocks
 *
 * Authors: PA Jacobs, RJ Gollan
 * 
 * Version: (Python code)
 * 26-Feb-2012 : functions moved out of estcj.py to this module.
 * 02-May-2013 : added more expansion_to_throat_calculation function from
 *               Matt McGilvray's gun tunnel version of nenzfr. -Chris James
 * 2017-Mar-09 : Ported to D by PJ
 *
 */

module gasflow;

import std.conv;
import std.math;
import std.string;
import std.stdio;
import std.algorithm;
import nm.bbla;
import nm.bracketing;
import nm.ridder;
import nm.linesearch;
import gas.gas_model;


double[] shock_ideal(const(GasState) state1, double Vs, GasState state2, GasModel gm)
/**
 * Computes post-shock conditions in the shock frame, assuming ideal gas.
 *
 * Input:
 *   state1: reference to pre-shock Gas state (given)
 *   Vs: speed of gas coming into shock (given)
 *   state2: reference to post-shock Gas state (to be estimated)
 *   gm: the gas model in use
 *
 * Returns: the dynamic array [V2, Vg], containing the post-shock gas speeds,
 *   V2 in the shock-reference frame,
 *   Vg in the lab frame.
 */
{
    double M1 = Vs / state1.a;
    double V1 = Vs;
    double gam = gm.gamma(state1);
    double R = gm.R(state1);
    //
    state2.rho = state1.rho * (gam+1.0)*M1*M1 / (2.0+(gam-1.0)*M1*M1);
    state2.p = state1.p * (2.0*gam*M1*M1 - (gam-1.0)) / (gam+1.0);
    state2.Ttr = state2.p / (R*state2.rho);
    gm.update_thermo_from_pT(state2);
    gm.update_sound_speed(state2);
    //
    double V2 = state1.rho / state2.rho * V1;
    double Vg = V1 - V2;
    return [V2, Vg];
} // end shock_ideal()

double my_limiter(double delta, double orig, double frac=0.5)
// Limit the magnitude of delta to no more than a fraction of the original.
// It occasionally happens that the Newton iterations go badly.
// It is worth trying to take smaller steps in these situations,
// assuming that the computed direction is still a fair guess.
{
    return copysign(fmin(abs(delta),frac*abs(orig)), delta);
}

double[] normal_shock(const(GasState) state1, double Vs, GasState state2, GasModel gm)
/**
 * Computes post-shock conditions in a shock-stationary frame.
 *
 * Input:
 *   state1: reference to pre-shock Gas state (given)
 *   Vs: speed of gas coming into shock (given)
 *   state2: reference to post-shock Gas state (to be estimated)
 *   gm: the gas model in use
 *
 * Returns: the dynamic array [V2, Vg], containing the post-shock gas speeds,
 *   V2 in the shock-reference frame,
 *   Vg in the lab frame.
 *
 * Note that, at the moment, we have implemented only a simple ideal-gas starting
 * guess for the Newton iteration.  We may need to reintroduce Chris James' 
 * guessed R and gamma options for some very strong shocks and difficult gases. 
 * An alternative may be do simply reduce the temperature by some factor and 
 * raise the density by the same factor, keeping the pressure about the same.
 * I vaguely recall Malcolm McIntosh had a correlation with shock velocity or 
 * Mach number for his adjustments to the ideal-gas initial guess. 
 */
{
    // Initial guess via ideal gas relations.
    double[] velocities = shock_ideal(state1, Vs, state2, gm);
    double V2 = velocities[0]; double Vg = velocities[1];
    // We assume that state2 now contains a fair initial guess
    // and set up the target values for the Rankine-Hugoniot relations.
    double V1 = Vs;
    double momentum = state1.p + state1.rho * V1 * V1;
    double total_enthalpy = gm.enthalpy(state1) + 0.5 * V1 * V1;
    //
    auto Fvector = delegate(double rho2, double T2)
    {
	// Constraint equations for state2 from the normal shock relations.
	// The correct post-shock values allow this vector to evaluate to zeros.
	state2.rho = rho2; state2.Ttr = T2;
	gm.update_thermo_from_rhoT(state2);
	V2 = V1 * state1.rho / rho2; // mass conservation
        double f1 = momentum - state2.p - state2.rho*V2*V2;
        double f2 = total_enthalpy - gm.enthalpy(state2) - 0.5*V2*V2;
	return [f1, f2];
    };
    //
    // Augmented matrix for the linear equation coefficients.
    Matrix Ab = new Matrix(2, 3);
    //
    double rho_delta = 1.0;
    double T_delta = 1.0;
    double rho_tol = 1.0e-3; // tolerance in kg/m^3
    double T_tol = 0.25;  // tolerance in degrees K
    //
    // Update the estimates using the Newton-Raphson method.
    //
    foreach (count; 0..20) {
        double rho_save = state2.rho;
        double T_save = state2.Ttr;
        double[] f_save = Fvector(rho_save, T_save);
        // Use finite differences to compute the Jacobian.
        double d_rho = rho_save * 0.01;
        double d_T = T_save * 0.01;
        double[] f_values = Fvector(rho_save + d_rho, T_save);
        double df0drho = (f_values[0] - f_save[0]) / d_rho;
        double df1drho = (f_values[1] - f_save[1]) / d_rho;
        f_values = Fvector(rho_save, T_save + d_T);
        double df0dT = (f_values[0] - f_save[0]) / d_T;
        double df1dT = (f_values[1] - f_save[1]) / d_T;
	Ab[0,0] = df0drho; Ab[0,1] = df0dT; Ab[0,2] = -f_save[0];
	Ab[1,0] = df1drho; Ab[1,1] = df1dT; Ab[1,2] = -f_save[1];
	gaussJordanElimination(Ab);
        rho_delta = Ab[0,2]; T_delta = Ab[1,2];
        // Possibly limit the increments so that the Newton iteration is
        // less inclined to go crazy.
        rho_delta = my_limiter(rho_delta, rho_save);
        T_delta = my_limiter(T_delta, T_save);
        double rho_new = rho_save + rho_delta;
        double T_new   = T_save + T_delta;
	state2.rho=rho_new; state2.Ttr = T_new;
        gm.update_thermo_from_rhoT(state2);
        // Check convergence.
        if (abs(rho_delta) < rho_tol && abs(T_delta) < T_tol) { break; }
    } // end foreach count
    // Back-out velocities via continuity.
    V2 = V1 * state1.rho / state2.rho;
    Vg = V1 - V2;
    return [V2, Vg];
} // end normal_shock()


double[] normal_shock_p2p1(const(GasState) state1, double p2p1,
			   GasState state2, GasModel gm)
/**
 * Computes post-shock conditions, using high-temperature gas properties
 * and a shock-stationary frame.
 *
 * Input:
 *   state1: pre-shock gas state (given)
 *   p2p1: ration of pressure across the shock (given)
 *   state2: reference to the post-shock state (to be computed)
 *   gmodel: reference to the gas model
 *
 * Returns: an array containing
 *   the incident shock speed, V1
 *   the post-shock gas speed, V2 in the shock-reference frame
 *                             Vg in the lab frame.
 */
{
    state2.copy_values_from(state1);
    // Initial guess via ideal gas relations.
    double g = gm.gamma(state1);
    double Ms = sqrt(1+(g+1)/2/g*(p2p1-1.0));
    double V1ideal = Ms * state1.a;
    double[] velocities;
    // Set up error function that will be zero when we have the correct shock speed.
    auto error_in_p2p1 = delegate(double Vs) {
        velocities = normal_shock(state1, Vs, state2, gm);
        return (state2.p/state1.p - p2p1)/p2p1;
    };
    double vguess1 = V1ideal;
    double vguess2 = 1.1 * V1ideal;
    if (bracket!error_in_p2p1(vguess1, vguess2) < 0) {
	throw new Exception("normal_shock_p2p1 could not bracket the shock velocity.");
    }
    double V1 = solve!error_in_p2p1(vguess1, vguess2, 1.0e-3);
    velocities = normal_shock(state1, V1, state2, gm);
    double V2 = velocities[0]; double Vg = velocities[1];
    return [V1, V2, Vg];
} // end normal_shock_p2p1()


double reflected_shock(const(GasState) state2, double Vg,
		       GasState state5, GasModel gm)
/**
 * Computes state5 which has brought the gas to rest at the end of the shock tube.
 *
 * Input
 * state2: the post-incident-shock gas state
 * Vg: the lab-frame velocity of the gas in state 2
 * state5: reference to the stagnation state (to be computed)
 *
 * Returns: Vr, the reflected shock speed in the lab frame.
 */
{
    // As an initial guess, 
    // assume that we have a very strong shock in an ideal gas.
    double gam = gm.gamma(state2);
    double density_ratio = (gam+1.0)/(gam-1.0);
    double Vr_a = Vg / density_ratio;
    double[] velocities = normal_shock(state2, Vr_a+Vg, state5, gm);
    double V5 = velocities[0]; 
    // The objective function is the difference in speeds,
    // units are m/s.  A value of zero for this function means
    // that, as the shock propagates upstream with speed ur,
    // the processed test gas is left in the end of the tube
    // with a velocity of zero in the laboratory frame.
    double f_a = V5 - Vr_a;
    //
    // Now, update this guess using the secant method.
    //
    double Vr_b = 1.1 * Vr_a;
    velocities = normal_shock(state2, Vr_b+Vg, state5, gm); 
    V5 = velocities[0]; 
    double f_b = V5 - Vr_b;
    if (abs(f_a) < abs(f_b)) {
	swap(f_a, f_b);
	swap(Vr_a, Vr_b);
    }
    int count = 0;
    while (abs(f_b) > 0.5 && count < 20) {
        double slope = (f_b - f_a) / (Vr_b - Vr_a);
        double Vr_c = Vr_b - f_b / slope;
        velocities = normal_shock(state2, Vr_c+Vg, state5, gm);
	V5 = velocities[0];
        double f_c = V5 - Vr_c;
        if (abs(f_c) < abs(f_b)) {
            Vr_b = Vr_c; f_b = f_c;
        } else {
            Vr_a = Vr_c; f_a = f_c;
	}
        count += 1;
    }
    if (count >= 20) {
        throw new Exception("Reflected shock iteration did not converge.");
    }
    // At this point, Vr_b should be our best guess.
    // Update the gas state data and return the best-guess value.
    velocities = normal_shock(state2, Vr_b+Vg, state5, gm);
    return Vr_b;
} // end reflected_shock()


double expand_from_stagnation(double p_over_p0, const(GasState) state0,
			      GasState state1, GasModel gm)
/**
 * Given a stagnation condition state0, expand to a new pressure.
 *
 * Input:
 *   p_over_p0: pressure ratio
 *   state0: GasState object specifying stagnation conditions
 *   state1: GasState object for the expanded conditions (to be computed)
 *
 * Returns: the corresponding velocity (in m/s) of the expanded stream.
 */
{
    state1.copy_values_from(state0);
    state1.p = state0.p * p_over_p0;
    double s0 = gm.entropy(state0);
    gm.update_thermo_from_ps(state1, s0);
    // Matt McGilvray had a note about CEA giving bad entropy values
    // so we'll assert things are OK before proceeding.
    assert (abs(gm.entropy(state1) - s0)/abs(s0) < 0.001, "Bad entropy value.");
    double static_enthalpy = gm.enthalpy(state1);
    double total_enthalpy = gm.enthalpy(state0);
    double V = sqrt(2.0*(total_enthalpy - static_enthalpy));
    return V;
} // end expand_from_stagnation()


double expand_to_mach(double mach, const(GasState) state0,
		      GasState state1, GasModel gm)
/**
 * Given a stagnation condition state0, expand to a given Mach number.
 *
 * Input:
 *   mach: target mach number
 *   state0: GasState object specifying stagnation conditions
 *   state1: GasState object for the expanded conditions (to be computed)
 *
 * Returns: the corresponding velocity (in m/s) of the expanded stream.
 */
{
    double total_enthalpy = gm.enthalpy(state0);
    double s0 = gm.entropy(state0);
    state1.copy_values_from(state0);
    double p_over_p0_guess1 = 1.0;
    double p_over_p0_guess2 = 0.90;
    // [TODO] Could probably do better with ideal gas guess.
    auto error_in_mach = delegate(double p_over_p0) {
	double V = expand_from_stagnation(p_over_p0, state0, state1, gm);
	double a = state1.a;
	return mach - V/a;
    };
    if (bracket!error_in_mach(p_over_p0_guess1, p_over_p0_guess2) < 0) {
	throw new Exception("expand_to_mach() could not bracket the pressure ratio.");
    }
    double p_over_p0 = solve!error_in_mach(p_over_p0_guess1, p_over_p0_guess2, 1.0e-6);
    state1.p = state0.p * p_over_p0;
    gm.update_thermo_from_ps(state1, s0);
    double static_enthalpy = gm.enthalpy(state1);
    double V = sqrt(2.0*(total_enthalpy - static_enthalpy));
    return V;
} // end expand_to_mach()


void total_condition(const(GasState) state1, double V1, 
		     GasState state0, GasModel gm)
/**
 * Given a free-stream condition and velocity,
 * compute the corresponding stagnant condition
 * at which the gas is brought to rest isentropically.
 *
 * Input
 *   state1: Gas object specifying free-stream condition
 *   V1: free-stream velocity, m/s
 *   state0: reference to the stagnation gas state (to be computed)
 *   gm: gas model
 */
{
    double H1 = gm.enthalpy(state1) + 0.5*V1*V1;
    double s1 = gm.entropy(state1);
    state0.copy_values_from(state1);
    auto error_in_total_enthalpy = delegate(double x) {
        // The enthalpy at the stagnation condition should match
        // the total enthalpy of the stream.
	state0.p = x * state1.p;
        gm.update_thermo_from_ps(state0, s1);
	return (H1 - gm.enthalpy(state0))/abs(H1);
    };
    double x1 = 1.0; double x2 = 1.01;
    // [TODO] could probably do better with an ideal gas guess
    if (bracket!error_in_total_enthalpy(x1, x2) < 0) {
	throw new Exception("total_condition() could not bracket the pressure ratio.");
    }
    double x_total = solve!error_in_total_enthalpy(x1, x2, 1.0e-4);
    state0.p = x_total * state1.p;
    gm.update_thermo_from_ps(state0, s1);
} // end total_condition()


void pitot_condition(const(GasState) state1, double V1, 
		     GasState state2pitot, GasModel gm)
/**
 * Given a free-stream condition, compute the corresponding Pitot condition
 * at which the gas is brought to rest, possibly through a shock.
 *
 * Input
 *   state1: Gas object specifying free-stream condition
 *   V1: free-stream velocity, m/s
 *   state2: reference to the final gas state (to be computed)
 *   gm: gas model
 */
{
    if (V1 > state1.a) {
        // Supersonic free-stream; process through a shock first.
	GasState state2 = new GasState(state1);
        double[] velocities = normal_shock(state1, V1, state2, gm);
	double V2 = velocities[0];
        total_condition(state2, V2, state2pitot, gm);
    } else {
        // Subsonic free-stream
        total_condition(state1, V1, state2pitot, gm);
    }
} // end pitot_condition()


double steady_flow_with_area_change(const(GasState)state1, double V1, double A2_over_A1,
				    GasState state2, GasModel gm, double tol = 1.0e-4)
/**
 * Given station 1 condition, velocity and area-ratio A2/A1,
 * compute the steady, isentropic condition at station 2.
 *
 * Input:
 *   state1: Gas object specifying condition at station 1 (given)
 *   V1: velocity at station 1, m/s (given)
 *   A2_over_A1: area ratio between stations A2/A1 (given)
 *   state2: reference to GasState object for condition 2 (to be computed)
 *   tol: tolerance for function solver to find nozzle outlet condition
 *
 * Returns: velocity at station 2, m/s
 */
{
    state2.copy_values_from(state1);
    gm.update_sound_speed(state2);
    double M1 = abs(V1)/state2.a;
    double V2 = V1;
    if (abs(A2_over_A1 - 1.0) < 1.0e-6) { return V2; } // essentially no change
    //
    GasState total_cond = new GasState(state1);
    total_condition(state1, V1, total_cond, gm);
    double p2p1_max = total_cond.p/state1.p;
    double p2p1_min = 0.001;
    // Establish a suitable bracket for the pressure ratio.
    // [TODO] When setting up the initial guess for pressure ratio,
    // we could probably do better with the ideal relation between M and A/Astar.
    double p2p1_guess1;
    double p2p1_guess2;
    if (M1 > 1.0) {
        if (A2_over_A1 > 1.0) {
            // For a supersonic expansion, we will see a drop in presure.
	    p2p1_guess1 = 0.9;
            p2p1_guess2 = 1.0;
	} else {
            // For a supersonic compression, we will see a rise in pressure.
            p2p1_guess1 = min(1.1, 1.0+0.1*(p2p1_max-1));
            p2p1_guess2 = min(1.2, 1.0+0.9*(p2p1_max-1));
	}
    } else { // subsonic
        if (A2_over_A1 < 1.0) {
            // Subsonic nozzle will accelerate to lower pressures.
            p2p1_guess1 = 0.5;
            p2p1_guess2 = 1.0;
        } else {
            // Subsonic diffuser will decelerate to higher pressure.
            p2p1_guess1 = min(1.1, 1.0+0.1*(p2p1_max-1));
	    p2p1_guess2 = min(1.2, 1.0+0.9*(p2p1_max-1));
	}
    }
    // Set up constraint data and the error-function to be given to the solver.
    double H1 = gm.enthalpy(state1) + 0.5*V1*V1;
    double mdot1 = state1.rho * V1; // assuming unit area at station 1
    double s1 = gm.entropy(state1);
    auto error_in_mass_flux = delegate(double p2p1) {
        // The mass flux should be the same at each station.
        state2.copy_values_from(state1);
	state2.p *= p2p1;
	gm.update_thermo_from_ps(state2, s1);
        V2 = sqrt(2*(H1 - gm.enthalpy(state2)));
	double mdot2 = state2.rho * V2 * A2_over_A1;
        double mdot_error = (mdot2 - mdot1)/abs(mdot1);
        return mdot_error;
    };
    if (bracket!error_in_mass_flux(p2p1_guess1, p2p1_guess2, p2p1_min, p2p1_max) < 0) {
	throw new Exception("steady_flow_with_area_change() could not bracket" ~
			    " the pressure ratio.");
    }
    double p2p1 = solve!error_in_mass_flux(p2p1_guess1, p2p1_guess2, tol);
    state2.copy_values_from(state1);
    state2.p *= p2p1;
    gm.update_thermo_from_ps(state2, s1);
    V2 = sqrt(2*(H1 - gm.enthalpy(state2)));
    return V2;
} // end steady_flow_with_area_change()


//------------------------------------------------------------------------
// Finite-strength waves along characteristic lines.

double finite_wave_dp(string characteristic, double V1, const(GasState) state1,
		      double p2, GasState state2, GasModel gm, int steps=100)
/**
 * Process the gas isentropically, following a characteristic line.
 *
 * See Section 7.6 Finite Nonlinear Waves in JD Anderson's text
 * Modern Compressible Flow.
 *
 * Input:
 *   characteristic: is either 'cplus' or 'cminus'
 *   V1: initial gas velocity, in m/s
 *   state1: initial gas state
 *   p2: new pressure after processing, in Pa
 *   state2: reference to the final gas state (to be computed)
 *   gm: reference to the current gas model
 *   steps: number of small steps to take through the process
 *
 * Returns: flow velocity after processing.
 */
{
    double V2 = V1;
    double p1 = state1.p;
    double s1 = gm.entropy(state1);
    //
    double dV;
    double dp = (p2 - state1.p)/steps;
    // Use more than the requested number of steps if p2 < dp, 
    // to prevent an overshoot into -ve pressure. (Chris James)
    while (p2 < abs(dp)) {
        steps = to!int(steps * 1.1);
        dp = (p2 - state1.p)/steps;
    }
    state2.p = p1+0.5*dp; // effectively mid-point of next step
    gm.update_thermo_from_ps(state2, s1);        
    gm.update_sound_speed(state2);        
    foreach (i; 0 .. steps) {
        double rhoa = state2.rho * state2.a;
        if (characteristic == "cminus") {
            dV = dp / rhoa;
        } else {
            dV = -dp / rhoa;
	}
        V2 += dV;
        state2.p += dp;  // prepare for next step
        gm.update_thermo_from_ps(state2, s1);
	gm.update_sound_speed(state2);
    }
    // back up to the correct end-point
    state2.p -= 0.5 * dp;
    gm.update_thermo_from_ps(state2, s1);
    gm.update_sound_speed(state2);
    return V2;
} // end finite_wave_dp()


double finite_wave_dv(string characteristic, double V1, const(GasState) state1,
		      double V2_target, GasState state2, GasModel gm,
		      int steps=100, double Tmin=200.0)
/**
 * Process the gas isentropically, following a characteristic line.
 *
 * See Section 7.6 Finite Nonlinear Waves in JD Anderson's text
 * Modern Compressible Flow.
 *
 * Input:
 *   characteristic: is either 'cplus' or 'cminus'
 *   V1: initial gas velocity, in m/s
 *   state1: initial gas state
 *   V2_target: desired velocity after processing, in m/s
 *     Note that we may not reach the requested velocity before pressure 
 *     and temperature become too small.
 *   state2: reference to the final gas state (to be computed)
 *   gm: reference to the current gas model
 *   steps: number of small steps to take through the process
 *   Tmin: temperature (in Kelvin) below which we terminate the process.
 *     We have this minimum to avoid problems with the thermodynamic
 *     polynomials of CEA2 program.  If you really want to work with very low
 *     temperatures, it's probably best to use an ideal gas model.
 *
 * Returns: flow velocity after processing.
 */
{
    double V2 = V1;
    double dV = (V2_target - V1)/steps;
    double s1 = gm.entropy(state1);
    state2.copy_values_from(state1);
    gm.update_sound_speed(state2);
    foreach (i; 0 .. steps) {
        double rhoa = state2.rho * state2.a;
	double dp;
        if (characteristic == "cminus") {
            dp = dV * rhoa;
        } else {
            dp = -dV * rhoa;
	}
        V2 += dV;
        state2.p += dp;
        gm.update_thermo_from_ps(state2, s1);
	gm.update_sound_speed(state2);
        if (state2.Ttr < Tmin) { break; }
    }
    return V2;
} // end finite_wave_dv()


//------------------------------------------------------------------------
// Oblique shock relations

double[] theta_oblique(const(GasState) state1, double V1, double beta,
		       GasState state2, GasModel gm)
/**
 * Compute the deflection angle and post-shock conditions given the shock wave angle.
 *
 * Input:
 *   state1: upstream gas condition
 *   V1: speed of gas into shock
 *   beta: shock wave angle with respect to stream direction (in radians)
 *   state2: reference to downstream gas condition (to be computed)
 *   gm: reference to current gas model
 *
 * Returns: array of theta, V2.
 *   theta is stream deflection angle in radians
 *   V2 is post-shock speed of gas in m/s
 */
{
    double V1_n = V1 * sin(beta);
    double V_t = V1 * cos(beta);
    state2.copy_values_from(state1);
    gm.update_sound_speed(state2);
    double M1_n = V1 / state2.a; // normal Mach number coming into shock
    if (M1_n < 1.0) {
        throw new Exception(format("theta_oblique(): subsonic inflow M1_n=%e", M1_n));
    }
    double[] velocities = normal_shock(state1, V1_n, state2, gm);
    double V2_n = velocities[0]; double Vg_n = velocities[1]; 
    double V2 = sqrt(V2_n * V2_n + V_t * V_t);
    double theta = beta - atan2(V2_n, V_t);
    return [theta, V2];
} // end theta_oblique()


double beta_oblique(const(GasState) state1, double V1, double theta,
		    GasModel gm)
/**
 * Compute the oblique shock wave angle given the deflection angle.
 *
 * Input:
 *   state1: upstream gas condition
 *   V1: speed of gas into shock
 *   theta: stream deflection angle (in radians)
 *   gm: reference to the current gas model
 *
 * Returns: shock wave angle wrt incoming stream direction (in radians)
 */
{
    GasState state2 = new GasState(state1);
    gm.update_sound_speed(state2);
    double M1 = V1 / state2.a;
    double b1 = max(asin(1.0/M1), 1.1*theta);
    double b2 = b1 * 1.05;
    auto error_in_theta = delegate(double beta_guess) {
        double[] shock_results = theta_oblique(state1, V1, beta_guess, state2, gm);
        double theta_guess = shock_results[0]; double V2 = shock_results[1]; 
        double error_value = theta_guess - theta;
        return error_value;
    };
    if (bracket!error_in_theta(b1, b2, asin(1.0/M1), PI/2) < 0) {
	throw new Exception("beta_oblique(): failed to converge on a shock-wave angle.");
    }
    double beta_result = solve!error_in_theta(b1, b2, 1.0e-4);
    return beta_result;
}

//------------------------------------------------------------------------
// Taylor-Maccoll cone flow.

double[] EOS_derivatives(const(GasState) state_0, GasModel gm)
/**
 * Compute equation-of-state derivatives at the specified state.
 *
 * Input:
 *   state: a complete state (with valid data)
 *   gm: reference to current gas model
 *
 * Returns: array of approximations [drho/dp, drho/dh]
 */
{
    double rho_0 = state_0.rho;
    double h_0 = gm.enthalpy(state_0);
    double T_0 = state_0.Ttr;
    // Choose relatively-small increments in enthalpy (J/kg) and pressure (Pa).
    double dh = abs(h_0) * 0.01 + 1000.0;
    double dp = state_0.p * 0.01 + 1000.0;
    // We're actually going to work in changes of p and T.
    double Cp = gm.dhdT_const_p(state_0);
    double dT = dh/Cp;
    // Use finite-differences to get the partial derivative.
    GasState state_new = new GasState(state_0);
    state_new.p = state_0.p + dp;
    gm.update_thermo_from_pT(state_new);
    double drhodp = (state_new.rho - rho_0) / dp;
    // and again, for the change in h, holding p constant.
    state_new.copy_values_from(state_0);
    state_new.p = state_0.p;
    state_new.Ttr = state_0.Ttr + dT;
    gm.update_thermo_from_pT(state_new);
    double drhodh = (state_new.rho - rho_0) / dh;
    // Assume that these first-order differences will suffice.
    return [drhodp, drhodh];
} // end EOS_derivatives()


double[] taylor_maccoll_odes(double[] z, double theta,
			     const(GasState) gas_state, GasModel gm)
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
    auto A = zeros(5,6); // Augmented matrix with rhs in last column.
    double[] derivs = EOS_derivatives(gas_state, gm);
    double dfdp = derivs[0]; double dfdh = derivs[1];
    debug { writeln("DEBUG dfdp=", dfdp, " dfdh=", dfdh); }
    A[0,0] = V_theta; A[0,2] = rho; A[0,5] = -2.0*rho*V_r - rho*V_theta/tan(theta);
    A[1,1] = 1.0; A[1,5] = V_theta;
    A[2,1] = rho*V_r; A[2,2] = rho*V_theta; A[2,4] = 1.0;
    A[3,1] = V_r; A[3,2] = V_theta; A[3,3] = 1.0;
    A[4,0] = 1.0; A[4,3] = -dfdh; A[4,4] = -dfdp;
    gaussJordanElimination(A);
    double[] dzdtheta =  A.getColumn(5);
    return dzdtheta;
}
/+

def theta_cone(state1, V1, beta):
    """
    Compute the cone-surface angle and conditions given the shock wave angle.

    :param state1: upstream gas condition
    :param V1: speed of gas into shock
    :param beta: shock wave angle wrt stream direction (in radians)
    :returns: tuple of theta_c, V_c and state_c:
        theta_c is stream deflection angle in radians
        V_c is cone-surface speed of gas in m/s
        state_c is cone-surface gas state

    The computation starts with the oblique-shock jump and then integrates
    across theta until V_theta goes through zero.
    The cone surface corresponds to V_theta == 0.
    """
    # Start at the point just downstream the oblique shock.
    theta_s, V2, state2 = theta_oblique(state1, V1, beta)
    #
    # Initial conditions.
    dtheta = -0.5 * math.pi / 180.0  # fraction-of-a-degree steps
    theta = beta
    V_r = V2 * math.cos(beta - theta_s)
    V_theta = -V2 * math.sin(beta - theta_s)
    rho = state2.rho; h = state2.h; p = state2.p
    gas_state = state2.clone()
    # For integrating across the shock layer, the state vector is:
    z = numpy.array([rho, V_r, V_theta, h, p])
    while V_theta < 0.0:
        # Keep a copy for linear interpolation at the end.
        z_old = z.copy(); theta_old = theta
        # Do the update using a low-order method (Euler) for the moment.
        dzdtheta = taylor_maccoll_odes(z, theta, gas_state)
        z += dtheta * dzdtheta; theta += dtheta
        rho, V_r, V_theta, h, p = z
        gas_state.set_ph(p, h)
        if DEBUG_GAS_FLOW: print "DEBUG theta=", theta, "V_r=", V_r, "V_theta=", V_theta
    # At this point, V_theta should have crossed zero so
    # we can linearly-interpolate the cone-surface conditions.
    V_theta_old = z_old[2]
    frac = (0.0 - V_theta_old)/(V_theta - V_theta_old)
    z_c = z_old*(1.0-frac) + z*frac
    theta_c = theta_old*(1.0-frac) + theta*frac
    # At the cone surface...
    rho, V_r, V_theta, h, p = z_c
    gas_state.set_ph(p, h)
    assert abs(V_theta) < 1.0e-6
    #
    return theta_c, V_r, gas_state


def beta_cone(state1, V1, theta):
    """
    Compute the conical shock wave angle given the cone-surface deflection angle.

    :param state1: upstream gas condition
    :param V1: speed of gas into shock
    :param theta: stream deflection angle (in radians)
    :returns: shock wave angle wrt incoming stream direction (in radians)
    """
    M1 = V1 / state1.a
    b1 = max(math.asin(1.0/M1), theta) * 1.01 # to be stronger than a Mach wave
    b2 = b1 * 1.05
    def error_in_theta(beta_guess):
        theta_guess, V_c, state_c = theta_cone(state1, V1, beta_guess)
        return theta_guess - theta
    beta_result = secant(error_in_theta, b1, b2, tol=1.0e-4)
    if beta_result == 'FAIL':
        raise RuntimeError('beta_cone(): failed to converge on a shock-wave angle.')
    return beta_result

+/
//------------------------------------------------------------------------

unittest {
}    
