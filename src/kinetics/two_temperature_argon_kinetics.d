/**
 * two_temperature_argon_kinetics.d
 *
 * Two-temperature reacting argon based on the model of:
 * Martin I. Hoffert and Hwachii Lien
 * "Quasi-One-Dimensional, Nonequilibrium Gas Dynamics of Partially Ionised Two-Temperature Argon"
 * The Physics of Fluids 10, 1769 (1967); doi 10.1063/1.1762356
 *
 * Authors: Daniel Smith, Rory Kelly and Peter J.
 * Version: 19-July-2017: initial cut.
 *          11-Feb-2019: rebuild to remove use of Matrix class
 *          23-Jan-2020: adaptive time step
 */

module kinetics.two_temperature_argon_kinetics;

import std.stdio : writeln, writefln;
import std.format;
import std.math;
import std.conv;
import ntypes.complex;
import nm.number;
import nm.bbla;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

final class UpdateArgonFrac : ThermochemicalReactor {

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        // We need to pick a number of pieces out of the gas-model file, again.
        // Although some already exist in the GasModel object, they are private.
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "TwoTemperatureReactingArgon");
        _ion_tol = getDoubleWithDefault(L, -1, "ion_tol", 0.0);
        _integration_method = getStringWithDefault(L, -1, "integration_method", "Backward_Euler");
        _T_min_for_reaction = getDoubleWithDefault(L, -1, "T_min_for_reaction", 3000.0);
        _Te_default = getDoubleWithDefault(L, -1, "Te_default", 10000.0);
        _n_step_suggest = getIntWithDefault(L, -1, "n_step_suggest", 10);
        _Newton_Raphson_tol = getDoubleWithDefault(L, -1, "Newton_Raphson_tol", 1.0e-6); // relative error
        _max_iter_newton = getIntWithDefault(L, -1, "max_iter_newton", 30);
        _argon_argon_rate_multiplier = getIntWithDefault(L, -1, "argon_argon_rate_multiplier", 1);
        lua_close(L);
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        if (dtSuggest < 0.0) {
            // A negative value indicated that the flow solver did not have
            // a suggested time-step size.
            dtSuggest = tInterval/_n_step_suggest;
        }
        // There are changes only if the gas is hot enough.
        if (Q.T > _T_min_for_reaction) {
            //
            // Determine the current number densities.
            number[3] numden;
            _gmodel.massf2numden(Q, numden);
            number n_e = numden[2]; // number density of electrons
            number n_Ar = numden[0]; // number density of Ar
            //
            // This is a model of an isolated reactor so the total number
            // of heavy particles and the energy in the reactor remain constant.
            // Note, also, that n_e == n_Ar_plus.
            _n_total = n_e + n_Ar;
            number alpha = n_e/_n_total;
            _u_total = Q.u + Q.u_modes[0];
            //
            // Given the amount of energy, we have a limit to the number of ions
            // that can be formed before driving the energy of the heavy particles
            // too low.
            // Let's put a lower limit of 200K for the temperature of the heavy particles.
            number u_available = _u_total - 3.0/2.0*_Rgas*200.0;
            number alpha_max = fmin(1.0, u_available/(_Rgas*(Q.T_modes[0]+_theta_ion)));
            _n_e_max = alpha_max * _n_total;
            _u_min_heavy_particles = 3.0/2.0*_Rgas*200.0;
            //
            // Retain a copy of the initial state, in case we restart
            // the integration with more, smaller time time steps.
            number initial_n_e = n_e;
            number initial_translational_energy = Q.u;
            //
            // Start with the suggested time step size.
            int NumberSteps = cast(int) fmax(floor(tInterval/dtSuggest), 1.0);
            number[2] y;
            int integration_attempt = 0;
            bool finished_integration = false;
            while (!finished_integration) {
                ++integration_attempt;
                try {
                    // Integrate across the interval with fixed time steps.
                    _chem_dt = tInterval/NumberSteps;
                    // Pack the state vector, ready for the integrator.
                    y[0] = initial_n_e;
                    y[1] = initial_translational_energy;
                    //
                    switch (_integration_method) {
                    case "Forward_Euler":
                        foreach (n; 0 .. NumberSteps) {
                            number[2] myF = F(y, Q);
                            foreach (i; 0 .. 2) { y[i] += _chem_dt * myF[i]; }
                            if (!state_vector_is_within_limits(y)) {
                                string msg = "State vector not within limits";
                                debug { msg ~= format("\n    y=[%g,%g]", y[0], y[1]); }
                                throw new ThermochemicalReactorUpdateException(msg);
                            }
                        } // end foreach n
                        break;
                    case "Backward_Euler":
                        // Perturbation sizes for the finite-difference Jacobian.
                        number[2] h = [1.0e-6*_n_total, 1.0e-6*_u_total];
                        number[2] y_prev; y_prev[0] = y[0]; y_prev[1] = y[1];
                        foreach (n; 0 .. NumberSteps) {
                            double norm_error = 1.0e50; // something large that will be replaced
                            int iter = 0;
                            do {
                                number[2][2] J = Jacobian(y, y_prev, h, Q);
                                number[2] rhs = BackEuler_F(y, y_prev, Q);
                                number[2] dy; solve2(J, rhs, dy);
                                foreach (i; 0 .. 2) { y[i] -= dy[i]; }
                                if (!state_vector_is_within_limits(y)) {
                                    string msg = "State vector not within limits";
                                    debug { msg ~= format("\n    y=[%g,%g] iter=%d", y[0], y[1], iter); }
                                    throw new ThermochemicalReactorUpdateException(msg);
                                }
                                number[2] error = BackEuler_F(y, y_prev, Q);
                                norm_error = fabs(error[0].re/(_n_total.re+1.0)) + fabs(error[1].re/(_u_total.re+1.0));
                                ++iter;
                            } while (norm_error > _Newton_Raphson_tol && iter < _max_iter_newton);
                            if (norm_error > _Newton_Raphson_tol) {
                                string msg = "Backward-Euler Newton iteration did not converge.";
                                debug { msg ~= format("\n    norm_error=%g iter=%d", norm_error, iter); }
                                throw new ThermochemicalReactorUpdateException(msg);
                            }
                            foreach (i; 0 .. 2) { y_prev[i] = y[i]; } // needed for next time step
                        } // end foreach n
                        break;
                    case "RK4":
                        number[2] k1, k2_in, k2, k3_in, k3, k4_in, k4;
                        foreach (n; 0 .. NumberSteps) {
                            number[2] myF = F(y, Q);
                            foreach (i; 0 .. 2) { k1[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 2) { k2_in[i] = y[i] + k1[i]/2.0; }
                            myF = F(k2_in, Q);
                            foreach (i; 0 .. 2) { k2[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 2) { k3_in[i] = y[i] + k2[i]/2.0; }
                            myF = F(k3_in, Q);
                            foreach (i; 0 .. 2) { k3[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 2) { k4_in[i] = y[i] + k3[i]; }
                            myF = F(k4_in, Q);
                            foreach (i; 0 .. 2) { k4[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 2) { y[i] += 1.0/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]); }
                            if (!state_vector_is_within_limits(y)) {
                                string msg = "State vector not within limits";
                                debug { msg ~= format("\n    y=[%g,%g]", y[0], y[1]); }
                                throw new ThermochemicalReactorUpdateException(msg);
                            }
                        } // end foreach n
                        break;
                    default:
                        throw new Exception("Invalid ODE update selection.");
                    } // end switch _integration_method
                    finished_integration = true;
                    if (integration_attempt == 1) {
                        // Success on the first attempt.
                        // Increase the time step a little, for next call.
                        _chem_dt *= 1.1;
                    }
                } catch(Exception e) {
                    if (integration_attempt < 10) {
                        // We will allow a retry with a smaller time step.
                        NumberSteps *= 3;
                    } else {
                        // We have taken too many attempts and have continued to fail.
                        string msg = "Thermochemical update fails after several attempts.";
                        debug { msg ~= text("\n Previous exception message: ", e.msg); }
                        throw new ThermochemicalReactorUpdateException(msg);
                    }
                }
            } // end while !finished_integration

            // Unpack the results, clipping round-off errors.
            n_e = fmax(y[0], 0.0);
            Q.u = fmin(y[1], _u_total);
            // Reconstruct the other parts of the flow state.
            if (n_e > 0.0) {
                n_Ar = _n_total - n_e;
                Q.u_modes[0] = fmax(_u_total-Q.u, 0.0);
            } else {
                // Do not let the number of electrons go negative.
                // Force the numbers back to something physically realizable.
                Q.u = _u_total;
                Q.u_modes[0] = 0.0;
                //mass
                n_Ar = _n_total;
                n_e = 0.0;
                alpha = n_e/(n_e + n_Ar);
                //energy
                //temperature
                Q.T = 2.0/3.0*Q.u/_Rgas;
                Q.T_modes[0] = _Te_default;
            }

            // Must update the mass fractions before updating the gas state from rho and u...
            // Number density of Argon+ is the same as electron number density.
            numden[0] = n_Ar; numden[1] = n_e; numden[2] = n_e;
            _gmodel.numden2massf(numden, Q);
            _gmodel.update_thermo_from_rhou(Q);
            _gmodel.update_sound_speed(Q);
            //
            // Remember the size of the successful time step for the next call.
            dtSuggest = _chem_dt;
        } // end if Q.T > _T_min_for_reaction
    } // end opCall()

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source)
    {
        number[2] myF = rates(Q);

        number electron_formation_rate_in_numden    = myF[0];
        number energy_to_translation_per_kg_mixture = myF[1];

        if (Q.T > _T_min_for_reaction) {
            source[0] = -1.0*electron_formation_rate_in_numden*_m_Ar;
            source[1] =  1.0*electron_formation_rate_in_numden*_m_Arp;
            source[2] =  1.0*electron_formation_rate_in_numden*_m_e;
            source[3] = -1.0*energy_to_translation_per_kg_mixture*Q.rho;
        } else {
            source[0] = 0.0;
            source[1] = 0.0;
            source[2] = 0.0;
            source[3] = 0.0;
        }
    }

    private:
    @nogc
    number[2] F(ref const(number[2]) y, ref GasState Q)
    {
        update_Q_from_state_vector(y, Q);
        return rates(Q);
    }

    @nogc
    void update_Q_from_state_vector(ref const(number[2]) y, ref GasState Q)
    {
        // Definition of our state vector:
        // [0] number density of electrons and
        // [1] translational energy of heavy particles
        number[2] y_dash; y_dash[0] = 0.0; y_dash[1] = 0.0;

        // Update the GasState from state vector.
        //
        // We do this so that we always work the rate calculations from
        // a physically realizable state.
        number n_e = y[0];
        if (n_e < 0.0) {
            // Do not let the number of electrons go negative.
            string msg = "Electron number density tried to go negative.";
            throw new Exception(msg);
        }
        number n_Ar = _n_total - n_e;
        // There might be round-off error.
        Q.u = fmin(y[1], _u_total);
        Q.u_modes[0] = fmax(_u_total-Q.u, 0.0);
        //
        number alpha = n_e/(_n_total);
        number Te = (alpha > _ion_tol) ? Q.T_modes[0] : to!number(_Te_default);
        number T = Q.T;
        // Must update the mass fractions before updating the gas state from rho and u...
        // Number density of Argon+ is the same as electron number density.
        number[3] nden; nden[0] = n_Ar; nden[1] = n_e; nden[2] = n_e;
        _gmodel.numden2massf(nden, Q);
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

    @nogc
    number[2] rates(in GasState Q)
    // Compute the rate of change of the state vector for the ionisation reactions.
    // It no longer has the side effect of keeping the GasState up to date with y vector.
    {
        // Definition of our state vector:
        // [0] number density of electrons and
        // [1] translational energy of heavy particles
        number[2] y_dash; y_dash[0] = 0.0; y_dash[1] = 0.0;

        number[3] numden;
        _gmodel.massf2numden(Q, numden);
        number n_e = numden[2]; // number density of electrons
        number n_Arp= numden[1]; // number density of Ar
        number n_Ar = numden[0]; // number density of Ar
        number n_total = n_e + n_Ar;
        number alpha = n_e/(n_total);
        number Te = (alpha > _ion_tol) ? Q.T_modes[0] : to!number(_Te_default);
        number T = Q.T;

        //=====================================================================
        // Rate constants for ionisation reactions.
        // PJ, 2020-01-22, Having Te drop low is a problem for kre,
        // so we might over-write it locally.
        Te = fmax(Te, 3000.0);
        number kfA = 1.68e-26*T*sqrt(T)*(_theta_A1star/T+2)*exp(-_theta_A1star/T)*_argon_argon_rate_multiplier;
        number kfe = 3.75e-22*Te*sqrt(Te)*(_theta_A1star/Te+2)*exp(-_theta_A1star/Te);
        number krA = 5.8e-49*(_theta_A1star/T+2)*exp((_theta_ion-_theta_A1star)/T)*_argon_argon_rate_multiplier;
        number kre = 1.29e-44*(_theta_A1star/Te+2)*exp((_theta_ion-_theta_A1star)/Te);

        // Determine the current rate of ionisation.
        // Production rate due to argon collisions.
        number n_dot_A = kfA*pow(n_Ar,2) - krA*n_Ar*pow(n_e,2);
        // Production rate due to electron collisions.
        number n_dot_e = kfe*n_Ar*n_e - kre*pow(n_e,3);
        number n_dot = n_dot_A + n_dot_e;

        y_dash[0] = n_dot;

        //=====================================================================
        // Energy moving to electronic mode due to reactions.
        number alpha_dot = n_dot/n_total;
        number u_dot_reac = 3.0/2.0*_Rgas*alpha_dot*Te
            + alpha_dot*_Rgas*_theta_ion - n_dot_e*_Kb*_theta_ion/Q.rho;
        //
        number Q_ea, Q_ei, v_ea, v_ei;
        if (alpha > _ion_tol) {
            // Calculate Qea
            if (Te < 10.0e3) {
                Q_ea = 0.39 + Te*(-0.551e-4 + 0.595e-8*Te);
            } else if (Te < 500.0e3) {
                Q_ea = -0.35 + 0.775e-4*Te;
            } else {
                Q_ea = -0.35 + 0.775e-4*500.0e3;
            }
            Q_ea *= 1.0e-20;
            //
            // Calculate Q_ei
            Q_ei = 1.95e-10/Te/Te*log(1.53e8*pow(Te,3)/(n_e/1.0e6));
            if (Q_ei < 0.0) { Q_ei = 0.0; }
            //
            // Find the collision frequencies
            // electron-Ar collisions
            v_ea = (1.0-alpha)*Q.rho/_m_Ar*sqrt(8*_Kb*Te/to!number(PI)/_m_e)*Q_ea;
            // electron-Ar+ collisions
            v_ei = alpha*Q.rho/_m_Ar*sqrt(8*_Kb*Te/to!number(PI)/_m_e)*Q_ei;
            //
            alpha_dot = n_dot/(n_e+n_Ar);
            // energy transferred to electron mode through collisions
            number u_dot = 3*n_e*_m_e/_m_Ar*(v_ea+v_ei)*_Kb*(T-Te)/Q.rho;
            y_dash[1] = -u_dot - u_dot_reac;
        } else {
            y_dash[1] = 0.0 - u_dot_reac;
        }
        //
        return y_dash;
    } // end F()

    @nogc
    number[2] BackEuler_F(ref const(number[2]) y, ref const(number[2]) y_prev, ref GasState Q)
    {
        // When the value of y satisfies the backward Euler-update formula,
        // the output vector will have zero for each element.
        number[2] myF = F(y, Q);
        number[2] output;
        foreach (i; 0 .. 2) { output[i] = y_prev[i] - y[i] + _chem_dt*myF[i]; }
        return output;
    }

    @nogc
    number[2][2] Jacobian(ref const(number[2]) y, ref const(number[2]) y_prev,
                          ref const(number[2]) h, ref GasState Q)
    {
        // Compute dG/dy, where G is the backward Euler update function (above).
        number[2] myF0 = BackEuler_F(y, y_prev, Q);
        number[2] yph0; yph0[0] = y[0]+h[0]; yph0[1] = y[1]; // y plus increment in index 0
        number[2] myF1 = BackEuler_F(yph0, y_prev, Q);
        number[2] col1; foreach (i; 0 .. 2) { col1[i] = (myF1[i] - myF0[i])/h[0]; }
        number[2] yph1; yph1[0] = y[0]; yph1[1] = y[1]+h[1]; // y plus increment in index 1
        number[2] myF2 = BackEuler_F(yph1, y_prev, Q);
        number[2] col2; foreach (i; 0 .. 2) { col2[i] = (myF2[i] - myF0[i])/h[1]; }
        number[2][2] J;
        J[0][0] = col1[0]; J[0][1] = col2[0];
        J[1][0] = col1[1]; J[1][1] = col2[1];
        number det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (fabs(det) < 1.0e-20) {
            string msg = "Jacobian appears singular.";
            debug { msg ~= text(" J=", J, " y=", y, " y_prev=", y_prev, "current Q=", Q); }
            throw new ThermochemicalReactorUpdateException(msg);
        }
        return J;
    }

    @nogc
    void solve2(ref const(number[2][2]) a, ref const(number[2]) b, ref number[2] x)
    {
        // Solves the matrix equation Ax=b.
        number det = a[0][0]*a[1][1] - a[0][1]*a[1][0];
        x[0] = -(a[0][1]*b[1] - a[1][1]*b[0])/det;
        x[1] = (a[0][0]*b[1] - a[1][0]*b[0])/det;
    }

    @nogc
    bool state_vector_is_within_limits(ref number[2] y)
    {
        bool result = true;
        if (y[0] < 0.0) { result = false; }
        if (y[0] > _n_e_max) { result = false; }
        if (y[1] < _u_min_heavy_particles) { result = false; }
        // Allow for finite-difference perturbations when evaluating Jacobian.
        if (y[1] > _u_total*1.01) { result = false; }
        return result;
    }

private:
    immutable double _m_Ar = 6.6335209e-26; // mass of argon (kg)
    immutable double _m_e = 9.10938e-31; // mass of electron (kg)
    immutable double _m_Arp = _m_Ar - _m_e; // mass of argon ion (kg)
    immutable double _Kb = Boltzmann_constant;
    immutable double _Rgas = 208.0;
    immutable double _theta_ion = 183100.0;
    immutable double _theta_A1star = 135300.0;

    // Adjustable parameters, will be set in the constructor.
    string _integration_method;
    double _T_min_for_reaction; // degrees K
    double _Te_default;
    double _ion_tol;
    int _n_step_suggest; // number of substeps to take if dtSuggest<0.0
    double _Newton_Raphson_tol;
    int _max_iter_newton; // limit the number of Newton iterations for backard Euler
    double _argon_argon_rate_multiplier; // increases the argon-argon rate constant by this factor

    // The following items are constants for the duration
    // of the update for our abstract isolated reactor.
    number _n_total; // number density of atoms and ions combined
    number _u_total; // energy within the reactor
    // Since the energy in the reactor is limited,
    // the number of ionized particles also limited.
    number _n_e_max;
    // We don't want the translational energy of the heavy particles
    // dropping too low, so keep limit for their internal energy.
    number _u_min_heavy_particles;

    double _chem_dt;
} // end class UpdateArgonFrac


version(two_temperature_argon_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.two_temperature_reacting_argon;
    void main() {
        string modelFileName = "sample-input/two-temperature-reacting-argon-model.lua";
        lua_State* L = init_lua_State();
        doLuaFile(L, modelFileName);
        auto gm = new TwoTemperatureReactingArgon(L);
        lua_close(L);
        auto gd = GasState(3, 1);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.T_modes[0] = 300.0;
        gd.massf[0] = 1.0; gd.massf[1] = 0.0; gd.massf[2] = 0.0;

        // Need molar masses to determine alpha
        double[3] _mol_masses;
        _mol_masses[0] = 39.948e-3; // Units are kg/mol
        _mol_masses[2] = 5.485799e-7; // Units are kg/mol
        _mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
        double theta_ion = 183100.0;
        double alpha;
        double new_vel;

        double Ru = R_universal;
        double m_Ar = 6.6335209e-26; //mass of argon (kg)
        double M_Ar = Avogadro_number*m_Ar;

        // This is the benchmark case presented in the Hoffert & Lien paper,
        // Section IV. Normal-shock relaxation zone: the local steady-state approximation
        //
        // Pre-shock conditions
        //==============================
        double rho1 = 0.0213657;
        double T1 = 300.0;
        double p1 = 1333.22;
        double u1 = 6.1e3;
        double M1 = 18.91;
        //
        // Inital gas properties, immediate post-shock
        auto gs = GasState(3, 1);
        gs.p = p1 * 446.735124;
        gs.T = T1 * 112.620756;
        gs.T_modes[0] = T1; // presume that the electron temperature has not jumped
        gs.massf[0] = 1.0; gs.massf[1] = 0.0; gs.massf[2] = 0.0;
        gm.update_thermo_from_pT(gs);
        gm.update_sound_speed(gs); // (not necessary)

        // Some space- and time-stepping information
        double vel = 1548.98; // m/s
        double x = 0.0;
        double maxtime = 4.0e-6; // max time for the simulation
        double dt = 1.0e-9; // time interval for chemistry update
        int maxsteps = to!int(maxtime/dt + 1); // number of steps in which to iterate through
        int printfreq = maxsteps/10;
        int writefreq = maxsteps/1000;

        // initialise the storage arrays
        double[] t_list;
        double[] x_list;
        double[] T_list;
        double[] T_modes_list;
        double[] P_list;
        double[] u_list;
        double[] alpha_list;

        // Initial values for storage arrays
        t_list ~= 0.0;
        x_list ~= x;
        T_list ~= gs.T;
        T_modes_list ~= gs.T_modes[0];
        P_list ~= gs.p;
        u_list ~= vel;
        alpha_list ~= 0.0;

        // Now, step along, allowing the reactions to proceed.
        //
        auto argonChemUpdate = new UpdateArgonFrac(modelFileName, gm);
        double[maxParams] params; // ignore
        double dtSuggest = dt; // To give 100 steps per integration interval.
        //
        foreach (i; 1 .. maxsteps) {
            // Perform the chemistry update
            argonChemUpdate(gs, dt, dtSuggest, params);
            // New ionisation fraction
            alpha = (gs.massf[2]/_mol_masses[2]) /
                ((gs.massf[2]/_mol_masses[2])+(gs.massf[0]/_mol_masses[0]));
            // Update x position
            new_vel = u1/8.*(5 + 3.0/pow(M1,2) - sqrt(9*pow(1-1.0/pow(M1,2),2) +
                                                      96*alpha/5./pow(M1,2)*theta_ion/T1));
            x += (vel + new_vel)/2*dt;
            vel = new_vel;
            // Update rho, P and u.
            gs.rho = rho1*u1/vel;
            gs.p = rho1*((Ru/M_Ar)*T1 + pow(u1,2)) - gs.rho*pow(vel,2);
            // New internal energy
            double e_new = 5*(Ru/M_Ar)*T1/2. + pow(u1,2)/2 - pow(vel,2)/2 - gs.p/gs.rho;
            // Give all of the new energy to the heavy particles.
            gs.u = gs.u + (e_new - (gs.u + gs.u_modes[0]));
            // Update the temperature of the heavy paritcles based on this.
            gm.update_thermo_from_rhou(gs);
            //
            if (writefreq == 0) { writefreq = 1; }
            if ((i % writefreq) == 0) { // only save 1000 points
                t_list ~= i*dt;
                x_list ~= x;
                T_list ~= gs.T;
                T_modes_list ~= gs.T_modes[0];
                P_list ~= gs.p;
                u_list ~= vel;
                alpha_list ~= alpha;
            }
        } // end foreach i

        // writeln("writing to Data file... please wait");
        double[][] collateddata = [t_list,x_list,T_list,T_modes_list,P_list,u_list,alpha_list];
        //
        File file = File("two_temperature_argon_kinetics_test_results.data","w");
        foreach (i; 0..t_list.length) {
            file.writeln(collateddata[0][i], " ", collateddata[1][i], " ",
                         collateddata[2][i], " ", collateddata[3][i], " ",
                         collateddata[4][i], " ", collateddata[5][i], " ",
                         collateddata[6][i]);
        }
    assert(isClose(gs.T, 14730, 1.0e-3), failedUnitTest());
    assert(isClose(gs.p, 705738, 1.0e-3), failedUnitTest());
    assert(isClose(vel, 700.095, 1.0e-3), failedUnitTest());
    } // end main()
} // end two_temperature_argon_kinetics_test



