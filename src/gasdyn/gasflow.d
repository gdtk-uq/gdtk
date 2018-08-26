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
 * [TODO] change to using Brent's function solver rather than Ridder's.
 */

module gasflow;

import std.conv;
import std.math;
import std.string;
import std.stdio;
import std.algorithm;
import nm.complex;
import nm.number;
import nm.bbla;
import nm.bracketing;
import nm.ridder;
import nm.linesearch;
import gas;
import gas.cea_gas;
import gasflowexception;


number[] shock_ideal(const(GasState) state1, number Vs, GasState state2, GasModel gm)
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
    number M1 = Vs / state1.a;
    // We seem to be hitting a round-off problem with the following test.
    // The small decrement is to be sure that we really have been given
    // a truly subsonic inflow.
    if ((M1-1.0) < -1.0e-6) {
        string msg = text("M1 subsonic: Vs=", Vs, " a=", state1.a, " M1=", M1);
        throw new GasFlowException(msg);
    }
    number V1 = Vs;
    number gam = gm.gamma(state1);
    number R = gm.R(state1);
    //
    state2.rho = state1.rho * (gam+1.0)*M1*M1 / (2.0+(gam-1.0)*M1*M1);
    state2.p = state1.p * (2.0*gam*M1*M1 - (gam-1.0)) / (gam+1.0);
    state2.T = state2.p / (R*state2.rho);
    gm.update_thermo_from_pT(state2);
    gm.update_sound_speed(state2);
    //
    number V2 = state1.rho / state2.rho * V1;
    number Vg = V1 - V2;
    return [V2, Vg];
} // end shock_ideal()

number[] normal_shock(const(GasState) state1, number Vs, GasState state2,
                      GasModel gm, double rho_tol=1.0e-6, double T_tol = 0.1)
/**
 * Computes post-shock conditions in a shock-stationary frame.
 *
 * Input:
 *   state1: reference to pre-shock Gas state (given)
 *   Vs: speed of gas coming into shock (given)
 *   state2: reference to post-shock Gas state (to be estimated)
 *   gm: the gas model in use
 *   rho_tol: tolerance in kg/m^3
 *   T_tol: temperature tolerance in degrees K
 *
 * Returns: the dynamic array [V2, Vg], containing the post-shock gas speeds,
 *   V2 in the shock-reference frame,
 *   Vg in the lab frame.
 *
 */
{
    // Initial guess via ideal gas relations.
    number[] velocities = shock_ideal(state1, Vs, state2, gm);
    number V2 = velocities[0]; number Vg = velocities[1];
    if (cast(CEAGas) gm !is null) {
        // The ideal gas estimate may have an unreasonably high value for temperature.
        // I vaguely recall Malcolm McIntosh had a correlation with shock velocity or 
        // Mach number for his adjustments to the ideal-gas initial guess, however,
        // this is Chris James' simple adjustment that seems to allow the calculation
        // to proceed.
        state2.T = fmin(state2.T, 20000.0);
        gm.update_thermo_from_pT(state2);
    }
    // We assume that state2 now contains a fair initial guess
    // and set up the target values for the Rankine-Hugoniot relations.
    number V1 = Vs;
    number momentum = state1.p + state1.rho * V1 * V1;
    number total_enthalpy = gm.enthalpy(state1) + 0.5 * V1 * V1;
    //
    auto Fvector = delegate(number rho2, number T2)
    {
        // Constraint equations for state2 from the normal shock relations.
        // The correct post-shock values allow this vector to evaluate to zeros.
        state2.rho = rho2; state2.T = T2;
        gm.update_thermo_from_rhoT(state2);
        V2 = V1 * state1.rho / rho2; // mass conservation
        number f1 = momentum - state2.p - state2.rho*V2*V2;
        number f2 = total_enthalpy - gm.enthalpy(state2) - 0.5*V2*V2;
        return [f1, f2];
    };
    // Augmented matrix for the linear equation coefficients.
    Matrix!number Ab = new Matrix!number(2, 3);
    number rho_delta = 1.0; number T_delta = 1.0;
    //
    // Update the estimates for rho,T using the Newton-Raphson method.
    //
    // We put an upper limit on the number of iterations that should be plenty
    // for well behaved gas models.  The CEA2 code, however, can have inconsistent
    // thermo data at the joins of the polynomial pieces and the Newton iterations
    // may go into a limit cycle that is close to the solution but never within
    // the requested tolerances.  If this happens, we just accept whatever final
    // iteration is computed.
    foreach (count; 0 .. 20) {
        number rho_save = state2.rho;
        number T_save = state2.T;
        number[] f_save = Fvector(rho_save, T_save);
        // Use finite differences to compute the Jacobian.
        number d_rho = rho_save * 0.01;
        number d_T = T_save * 0.01;
        number[] f_values = Fvector(rho_save + d_rho, T_save);
        number df0drho = (f_values[0] - f_save[0]) / d_rho;
        number df1drho = (f_values[1] - f_save[1]) / d_rho;
        f_values = Fvector(rho_save, T_save + d_T);
        number df0dT = (f_values[0] - f_save[0]) / d_T;
        number df1dT = (f_values[1] - f_save[1]) / d_T;
        Ab[0,0] = df0drho; Ab[0,1] = df0dT; Ab[0,2] = -f_save[0];
        Ab[1,0] = df1drho; Ab[1,1] = df1dT; Ab[1,2] = -f_save[1];
        gaussJordanElimination(Ab);
        // These are the full increments according to the Newton method.
        rho_delta = Ab[0,2]; T_delta = Ab[1,2];
        // Limit the increments so that the iteration is more stable for the difficult cases.
        rho_delta = copysign(fmin(abs(rho_delta),0.2*abs(rho_save)), rho_delta);
        T_delta = copysign(fmin(abs(T_delta),0.2*abs(T_save)), T_delta);
        state2.rho = rho_save + rho_delta;
        state2.T = T_save + T_delta;
        gm.update_thermo_from_rhoT(state2);
        // Check convergence.
        if (abs(rho_delta) < rho_tol && abs(T_delta) < T_tol) { break; }
    } // end foreach count
    // Back-out velocities via continuity.
    V2 = V1 * state1.rho / state2.rho;
    Vg = V1 - V2;
    gm.update_sound_speed(state2); // be sure that state2 is complete
    return [V2, Vg];
} // end normal_shock()


number[] normal_shock_p2p1(const(GasState) state1, number p2p1,
                           GasState state2, GasModel gm)
/**
 * Computes post-shock conditions, using high-temperature gas properties
 * and a shock-stationary frame.
 *
 * Input:
 *   state1: pre-shock gas state (given)
 *   p2p1: ratio of pressure across the shock (given)
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
    number g = gm.gamma(state1);
    number Ms = sqrt(1+(g+1)/2/g*(p2p1-1.0));
    number V1ideal = Ms * state1.a;
    number[] velocities;
    // Set up error function that will be zero when we have the correct shock speed.
    auto error_in_p2p1 = delegate(number Vs) {
        velocities = normal_shock(state1, Vs, state2, gm);
        return (state2.p/state1.p - p2p1)/p2p1;
    };
    number vguess1 = V1ideal;
    number vguess2 = 1.1 * V1ideal;
    if (bracket!error_in_p2p1(vguess1, vguess2) < 0) {
        throw new GasFlowException("normal_shock_p2p1 could not bracket the shock velocity.");
    }
    number V1 = solve!error_in_p2p1(vguess1, vguess2, 1.0e-6);
    velocities = normal_shock(state1, V1, state2, gm);
    number V2 = velocities[0]; number Vg = velocities[1];
    return [V1, V2, Vg];
} // end normal_shock_p2p1()


number reflected_shock(const(GasState) state2, number Vg,
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
    number gam = gm.gamma(state2);
    number density_ratio = (gam+1.0)/(gam-1.0);
    number Vr_a = Vg / density_ratio;
    number[] velocities = normal_shock(state2, Vr_a+Vg, state5, gm);
    number V5 = velocities[0]; 
    // The objective function is the difference in speeds,
    // units are m/s.  A value of zero for this function means
    // that, as the shock propagates upstream with speed ur,
    // the processed test gas is left in the end of the tube
    // with a velocity of zero in the laboratory frame.
    number f_a = V5 - Vr_a;
    //
    // Now, update this guess using the secant method.
    //
    number Vr_b = 1.1 * Vr_a;
    velocities = normal_shock(state2, Vr_b+Vg, state5, gm); 
    V5 = velocities[0]; 
    number f_b = V5 - Vr_b;
    if (abs(f_a) < abs(f_b)) {
        swap(f_a, f_b);
        swap(Vr_a, Vr_b);
    }
    int count = 0;
    while (abs(f_b) > 0.5 && count < 20) {
        number slope = (f_b - f_a) / (Vr_b - Vr_a);
        number Vr_c = Vr_b - f_b / slope;
        velocities = normal_shock(state2, Vr_c+Vg, state5, gm);
        V5 = velocities[0];
        number f_c = V5 - Vr_c;
        if (abs(f_c) < abs(f_b)) {
            Vr_b = Vr_c; f_b = f_c;
        } else {
            Vr_a = Vr_c; f_a = f_c;
        }
        count += 1;
    }
    if (count >= 20) {
        throw new GasFlowException("Reflected shock iteration did not converge.");
    }
    // At this point, Vr_b should be our best guess.
    // Update the gas state data and return the best-guess value.
    velocities = normal_shock(state2, Vr_b+Vg, state5, gm);
    return Vr_b;
} // end reflected_shock()


number expand_from_stagnation(const(GasState) state0, number p_over_p0,
                              GasState state1, GasModel gm)
/**
 * Given a stagnation condition state0, expand to a new pressure.
 *
 * Input:
 *   state0: GasState object specifying stagnation conditions
 *   p_over_p0: pressure ratio
 *   state1: GasState object for the expanded conditions (to be computed)
 *
 * Returns: the corresponding velocity (in m/s) of the expanded stream.
 */
{
    state1.copy_values_from(state0);
    state1.p = state0.p * p_over_p0;
    number s0 = gm.entropy(state0);
    gm.update_thermo_from_ps(state1, s0);
    // Matt McGilvray had a note about CEA giving bad entropy values
    // so we'll assert things are OK before proceeding.
    assert (abs(gm.entropy(state1) - s0)/abs(s0) < 0.001, "Bad entropy value.");
    number static_enthalpy = gm.enthalpy(state1);
    number total_enthalpy = gm.enthalpy(state0);
    // In case we are given a pressure-ratio of 1.0 and the noise in the
    // CEA calculations lead to a situation where the new static enthalpy
    // is just slightly larger than the total enthalpy value, we fall back
    // to a zero value for velocity.
    number V = 0.0;
    if (total_enthalpy > static_enthalpy) {
        V = sqrt(2.0*(total_enthalpy - static_enthalpy));
    }
    return V;
} // end expand_from_stagnation()


number expand_to_mach(const(GasState) state0, number mach,
                      GasState state1, GasModel gm)
/**
 * Given a stagnation condition state0, expand to a given Mach number.
 *
 * Input:
 *   state0: GasState object specifying stagnation conditions
 *   mach: target mach number for the expanded gas
 *   state1: GasState object for the expanded conditions (to be computed)
 *
 * Returns: the corresponding velocity (in m/s) of the expanded stream.
 */
{
    number total_enthalpy = gm.enthalpy(state0);
    number s0 = gm.entropy(state0);
    state1.copy_values_from(state0);
    number p_over_p0_guess1 = 1.0;
    number p_over_p0_guess2 = 0.90;
    // [TODO] Could probably do better with ideal gas guess.
    auto error_in_mach = delegate(number p_over_p0) {
        number V = expand_from_stagnation(state0, p_over_p0, state1, gm);
        number a = state1.a;
        return mach - V/a;
    };
    if (bracket!(error_in_mach, number)(p_over_p0_guess1, p_over_p0_guess2,
                                        to!number(1.0e-6), to!number(1.0)) < 0) {
        throw new GasFlowException("expand_to_mach() could not bracket the pressure ratio.");
    }
    number p_over_p0 = solve!(error_in_mach, number)(p_over_p0_guess1, p_over_p0_guess2, 1.0e-6);
    state1.p = state0.p * p_over_p0;
    gm.update_thermo_from_ps(state1, s0);
    number static_enthalpy = gm.enthalpy(state1);
    number V = 0.0;
    if (total_enthalpy > static_enthalpy) {
        V = sqrt(2.0*(total_enthalpy - static_enthalpy));
    }
    return V;
} // end expand_to_mach()


void total_condition(const(GasState) state1, number V1, 
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
    number H1 = gm.enthalpy(state1) + 0.5*V1*V1;
    number s1 = gm.entropy(state1);
    state0.copy_values_from(state1);
    auto error_in_total_enthalpy = delegate(number x) {
        // The enthalpy at the stagnation condition should match
        // the total enthalpy of the stream.
        state0.p = x * state1.p;
        gm.update_thermo_from_ps(state0, s1);
        return (H1 - gm.enthalpy(state0))/abs(H1);
    };
    number x1 = 1.0; number x2 = 1.01;
    // [TODO] could probably do better with an ideal gas guess
    if (bracket!error_in_total_enthalpy(x1, x2) < 0) {
        throw new GasFlowException("total_condition() could not bracket the pressure ratio.");
    }
    number x_total = solve!error_in_total_enthalpy(x1, x2, 1.0e-4);
    state0.p = x_total * state1.p;
    gm.update_thermo_from_ps(state0, s1);
} // end total_condition()


void pitot_condition(const(GasState) state1, number V1, 
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
        number[] velocities = normal_shock(state1, V1, state2, gm);
        number V2 = velocities[0];
        total_condition(state2, V2, state2pitot, gm);
    } else {
        // Subsonic free-stream
        total_condition(state1, V1, state2pitot, gm);
    }
} // end pitot_condition()


number steady_flow_with_area_change(const(GasState)state1, number V1, number A2_over_A1,
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
 *   tol: tolerance for function solver to find condition 2
 *
 * Returns: velocity at station 2, m/s
 */
{
    state2.copy_values_from(state1);
    gm.update_sound_speed(state2);
    number M1 = abs(V1)/state2.a;
    number V2 = V1;
    if (abs(A2_over_A1 - 1.0) < 1.0e-6) { return V2; } // essentially no change
    //
    GasState total_cond = new GasState(state1);
    total_condition(state1, V1, total_cond, gm);
    number p2p1_max = total_cond.p/state1.p;
    number p2p1_min = 0.001;
    // Establish a suitable bracket for the pressure ratio.
    // [TODO] When setting up the initial guess for pressure ratio,
    // we could probably do better with the ideal relation between M and A/Astar.
    number p2p1_guess1;
    number p2p1_guess2;
    if (M1 > 1.0) {
        if (A2_over_A1 > 1.0) {
            // For a supersonic expansion, we will see a drop in presure.
            p2p1_guess1 = 0.9;
            p2p1_guess2 = 1.0;
        } else {
            // For a supersonic compression, we will see a rise in pressure.
            p2p1_guess1 = fmin(to!number(1.1), 1.0+0.1*(p2p1_max-1));
            p2p1_guess2 = fmin(to!number(1.2), 1.0+0.9*(p2p1_max-1));
        }
    } else { // subsonic
        if (A2_over_A1 < 1.0) {
            // Subsonic nozzle will accelerate to lower pressures.
            p2p1_guess1 = 0.5;
            p2p1_guess2 = 1.0;
        } else {
            // Subsonic diffuser will decelerate to higher pressure.
            p2p1_guess1 = fmin(to!number(1.1), 1.0+0.1*(p2p1_max-1));
            p2p1_guess2 = fmin(to!number(1.2), 1.0+0.9*(p2p1_max-1));
        }
    }
    // Set up constraint data and the error-function to be given to the solver.
    number H1 = gm.enthalpy(state1) + 0.5*V1*V1;
    number mdot1 = state1.rho * V1; // assuming unit area at station 1
    number s1 = gm.entropy(state1);
    auto error_in_mass_flux = delegate(number p2p1) {
        // The mass flux should be the same at each station.
        state2.copy_values_from(state1);
        state2.p *= p2p1;
        gm.update_thermo_from_ps(state2, s1);
        V2 = sqrt(2*(H1 - gm.enthalpy(state2)));
        number mdot2 = state2.rho * V2 * A2_over_A1;
        number mdot_error = (mdot2 - mdot1)/abs(mdot1);
        return mdot_error;
    };
    if (bracket!(error_in_mass_flux,number)(p2p1_guess1, p2p1_guess2, p2p1_min, p2p1_max) < 0) {
        throw new GasFlowException("steady_flow_with_area_change() could not bracket" ~
                                   " the pressure ratio.");
    }
    number p2p1 = solve!(error_in_mass_flux,number)(p2p1_guess1, p2p1_guess2, tol);
    state2.copy_values_from(state1);
    state2.p *= p2p1;
    gm.update_thermo_from_ps(state2, s1);
    V2 = sqrt(2*(H1 - gm.enthalpy(state2)));
    return V2;
} // end steady_flow_with_area_change()


//------------------------------------------------------------------------
// Finite-strength waves along characteristic lines.

number finite_wave_dp(const(GasState) state1, number V1,
                      string characteristic, number p2,
                      GasState state2, GasModel gm, int steps=100)
/**
 * Process the gas isentropically, following a characteristic line.
 *
 * See Section 7.6 Finite Nonlinear Waves in JD Anderson's text
 * Modern Compressible Flow.
 *
 * Input:
 *   state1: initial gas state
 *   V1: initial gas velocity, in m/s
 *   characteristic: is either 'cplus' or 'cminus'
 *   p2: target pressure after processing, in Pa
 *   state2: reference to the final gas state (to be computed)
 *   gm: reference to the current gas model
 *   steps: number of small steps to take through the process
 *
 * Returns: flow velocity after processing.
 */
{
    number V2 = V1;
    number p1 = state1.p;
    number s1 = gm.entropy(state1);
    //
    number dV;
    number dp = (p2 - state1.p)/steps;
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
        number rhoa = state2.rho * state2.a;
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


number finite_wave_dv(const(GasState) state1, number V1,
                      string characteristic, number V2_target,
                      GasState state2, GasModel gm,
                      int steps=100, double Tmin=200.0)
/**
 * Process the gas isentropically, following a characteristic line.
 *
 * See Section 7.6 Finite Nonlinear Waves in JD Anderson's text
 * Modern Compressible Flow.
 *
 * Input:
 *   state1: initial gas state
 *   V1: initial gas velocity, in m/s
 *   characteristic: is either 'cplus' or 'cminus'
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
    number V2 = V1;
    number dV = (V2_target - V1)/steps;
    number s1 = gm.entropy(state1);
    state2.copy_values_from(state1);
    gm.update_sound_speed(state2);
    foreach (i; 0 .. steps) {
        number rhoa = state2.rho * state2.a;
        number dp;
        if (characteristic == "cminus") {
            dp = dV * rhoa;
        } else {
            dp = -dV * rhoa;
        }
        V2 += dV;
        state2.p += dp;
        gm.update_thermo_from_ps(state2, s1);
        gm.update_sound_speed(state2);
        if (state2.T < Tmin) { break; }
    }
    return V2;
} // end finite_wave_dv()


//------------------------------------------------------------------------
// Oblique shock relations

number[] theta_oblique(const(GasState) state1, number V1, number beta,
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
    number V1_n = V1 * sin(beta);
    number V_t = V1 * cos(beta);
    state2.copy_values_from(state1);
    gm.update_sound_speed(state2);
    number M1_n = V1 / state2.a; // normal Mach number coming into shock
    if ((M1_n-1.0) < -1.0e-6) {
        throw new GasFlowException(format("theta_oblique(): subsonic inflow M1_n=%e", M1_n));
    }
    number[] velocities = normal_shock(state1, V1_n, state2, gm);
    number V2_n = velocities[0]; number Vg_n = velocities[1]; 
    number V2 = sqrt(V2_n * V2_n + V_t * V_t);
    number theta = beta - atan2(V2_n, V_t);
    return [theta, V2];
} // end theta_oblique()


number beta_oblique(const(GasState) state1, number V1, number theta,
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
    number M1 = V1 / state2.a;
    number b1 = max(asin(1.0/M1), 1.1*theta);
    number b2 = b1 * 1.05;
    auto error_in_theta = delegate(number beta_guess) {
        number[] shock_results = theta_oblique(state1, V1, beta_guess, state2, gm);
        number theta_guess = shock_results[0]; number V2 = shock_results[1]; 
        number error_value = theta_guess - theta;
        return error_value;
    };
    if (bracket!(error_in_theta,number)(b1, b2, to!number(asin(1.0/M1)), to!number(PI/2)) < 0) {
        throw new GasFlowException("beta_oblique(): failed to bracket a shock-wave angle.");
    }
    number beta_result = solve!(error_in_theta,number)(b1, b2, 1.0e-6);
    return beta_result;
}

//------------------------------------------------------------------------
// Taylor-Maccoll cone flow.

number[2] EOS_derivatives(const(GasState) state_0, GasModel gm)
/**
 * Compute equation-of-state derivatives at the specified state.
 *
 * Input:
 *   state: a complete state (with valid data)
 *   gm: reference to current gas model
 *
 * Returns: array of approximate derivatives [drho/dp, drho/du]
 */
{
    number rho_0 = state_0.rho;
    // Choose relatively-small increments in energy (J/kg) and pressure (Pa).
    number du = abs(state_0.u) * 0.01 + 1000.0;
    number dp = state_0.p * 0.01 + 1000.0;
    // We're actually going to work in changes of p and T.
    // The function name for Cv evaluation does look a bit odd,
    // in that the increment in internal energy with p constant
    // will not be constant volume.  However, what we really want
    // is the change in temperature for a change in internal energy.
    number Cv = gm.dudT_const_v(state_0);
    number dT = du/Cv;
    // Use finite-differences to get the partial derivative.
    GasState state_new = new GasState(state_0);
    state_new.p = state_0.p + dp;
    gm.update_thermo_from_pT(state_new);
    number drhodp = (state_new.rho - rho_0) / dp;
    // and again, for the change in u, holding p constant.
    state_new.copy_values_from(state_0);
    state_new.p = state_0.p;
    state_new.T = state_0.T + dT;
    gm.update_thermo_from_pT(state_new);
    number drhodu = (state_new.rho - rho_0) / du;
    // Assume that these first-order differences will suffice.
    return [drhodp, drhodu];
} // end EOS_derivatives()


number[5] taylor_maccoll_odes(number[5] z, number theta,
                              const(GasState) gas_state, GasModel gm)
{
    /**
    The ODEs from the Taylor-Maccoll formulation.

    See PJ's workbook for March 2017 for details.
    We've packaged them formally so that we might one day use
    a more sophisticated ODE integrator requiring fewer steps.
    **/
    number rho=z[0]; number V_r=z[1]; number V_theta=z[2];
    number u=z[3]; number p=z[4];
    // Assume gas_state is current.
    // Assemble linear system for determining the derivatives wrt theta.
    auto A = zeros!number(5,6); // Augmented matrix with rhs in last column.
    number[2] derivs = EOS_derivatives(gas_state, gm);
    number drhodp = derivs[0]; number drhodu = derivs[1];
    number one = 1.0;
    A[0,0] = V_theta; A[0,2] = rho; A[0,5] = -2.0*rho*V_r - rho*V_theta/tan(theta);
    A[1,1] = one; A[1,5] = V_theta;
    A[2,1] = rho*V_r; A[2,2] = rho*V_theta; A[2,4] = one;
    A[3,0] = -p/(rho^^2); A[3,1] = V_r; A[3,2] = V_theta; A[3,3] = one; A[3,4] = one/rho;
    A[4,0] = one; A[4,3] = -drhodu; A[4,4] = -drhodp;
    gaussJordanElimination!number(A);
    number[5] dzdtheta =  A.getColumn(5);
    return dzdtheta;
}

number[2] theta_cone(const(GasState) state1, number V1, number beta,
                     GasState state_c, GasModel gm)
/**
 * Compute the cone-surface angle and conditions given the shock wave angle.
 *
 * Input:
 *   state1: upstream gas condition
 *   V1: speed of gas into shock
 *   beta: shock wave angle wrt stream direction (in radians)
 *   state_c: reference to the gas state at the cone surface (to be computed)
 *   gm: reference to the current gas model
 *
 * Returns: array of [theta_c, V_c]:
 *   theta_c is stream deflection angle in radians
 *   V_c is cone-surface speed of gas in m/s
 *
 * The computation starts with the oblique-shock jump and then integrates
 * across theta until V_theta goes through zero.
 * The cone surface corresponds to V_theta == 0.
 */
{
    // [TODO] Implement Rowan's linear interpolation for weak shocks.
    //
    // Start at the point just downstream the oblique shock.
    GasState state2 = new GasState(state1);
    number[] shock_results = theta_oblique(state1, V1, beta, state2, gm);
    number theta_s = shock_results[0]; number V2 = shock_results[1];
    //
    // Initial conditions.
    number dtheta = -0.5 * PI/180.0;  // fraction-of-a-degree steps
    number theta = beta;
    number V_r = V2 * cos(beta - theta_s);
    number V_theta = -V2 * sin(beta - theta_s);
    number rho = state2.rho; number u = state2.u; number p = state2.p;
    //
    GasState gas_state = new GasState(state2);
    gas_state.rho = rho; gas_state.u = u;
    gm.update_thermo_from_rhou(gas_state);
    // For integrating across the shock layer, the state vector is:
    number[5] z = [rho, V_r, V_theta, u, p];
    number[5] z_old; number theta_old; number V_theta_old;
    while (V_theta < 0.0) {
        // Keep a copy for linear interpolation at the end.
        z_old[] = z[]; theta_old = theta; V_theta_old = V_theta;
        // Do the update using a low-order method (Euler) for the moment.
        number[5] dzdtheta = taylor_maccoll_odes(z, theta, gas_state, gm);
        foreach (i; 0 .. z.length) { z[i] += dtheta * dzdtheta[i]; }
        theta += dtheta;
        rho=z[0]; V_r=z[1]; V_theta=z[2]; u=z[3]; p=z[4];
        gas_state.rho = rho; gas_state.u = u;
        gm.update_thermo_from_rhou(gas_state);
        assert(abs(gas_state.p - p)/p < 0.001, "pressure diverging");
    }
    // At this point, V_theta should have crossed zero so
    // we can linearly-interpolate the cone-surface conditions.
    number frac = (0.0 - V_theta_old)/(V_theta - V_theta_old);
    number[5] z_c;
    foreach (i; 0 .. z_c.length) { z_c[i] = z_old[i]*(1.0-frac) + z[i]*frac; }
    number theta_c = theta_old*(1.0-frac) + theta*frac;
    // At the cone surface...
    rho=z_c[0]; V_r=z_c[1]; V_theta=z_c[2]; u=z_c[3]; p=z_c[4];
    state_c.rho = rho; state_c.u = u;
    gm.update_thermo_from_rhou(state_c);
    gm.update_sound_speed(state_c);
    assert(abs(V_theta) < 1.0e-6, "V_theta should be very small");
    return [theta_c, V_r];
} // end theta_cone()


number beta_cone(const(GasState) state1, number V1, number theta, GasModel gm)
/**
 * Compute the conical shock wave angle given the cone-surface deflection angle.
 *
 * Input:
 * state1: upstream gas condition
 * V1: speed of gas into shock
 * theta: stream deflection angle (in radians)
 * gm: reference to current gas model
 *
 * Returns: shock wave angle wrt incoming stream direction (in radians)
 */
{
    GasState state2 = new GasState(state1);
    gm.update_sound_speed(state2);
    number M1 = V1 / state2.a;
    number b1 = max(asin(1.0/M1), 1.1*theta); // to be stronger than a Mach wave
    number b2 = b1 * 1.05;
    auto error_in_theta = delegate(number beta_guess) {
        number[2] shock_results = theta_cone(state1, V1, beta_guess, state2, gm);
        number theta_guess = shock_results[0]; number V2 = shock_results[1]; 
        number error_value = theta_guess - theta;
        return error_value;
    };
    if (bracket!(error_in_theta,number)(b1, b2, to!number(asin(1.0/M1)), to!number(PI/2)) < 0) {
        throw new GasFlowException("beta_cone(): failed to converge on a shock-wave angle.");
    }
    number beta_result = solve!(error_in_theta,number)(b1, b2, 1.0e-6);
    return beta_result;
} // end beta_cone()

//--------------------------------------------------------------------------------

version (test_gasflow) {
    // [TODO]
}
