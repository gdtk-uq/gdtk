// shock_tube_nozzle.d: The principal analysis code is here.
//
// Authors: Peter J., Nick Gibbons, Rowan Gollan.
//
module shock_tube_nozzle;

import std.stdio;
import std.array;
import std.string;
import std.conv;
import std.getopt;
import std.file;
import std.math;
import dyaml;
import nm.secant;
import nm.schedule;
import gas;
import gas.cea_gas;
import gas.therm_perf_gas;
import gas.equilibrium_gas;
import gas.therm_perf_gas_equil;
import kinetics;
import gasdyn.gasflow;
import configuration;

// Storage for the calculated gas state at the end of the facility nozzle.
struct Result{
    double x;
    double area_ratio;
    double velocity;
    double Mach_number;
    double p_pitot;
    double rayleigh_pitot;
    double pressure, density, temperature, viscosity;
    double[] T_modes;
    double[] massf;
    double massflux_rel_err, enthalpy_rel_err, pitot_rel_err;
}

class Nenzf1dAnalysisException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

Result analyse(int verbosityLevel, Config config)
{
    // The core of the nenzf1d program is a sequence of state-to-state calculations
    // for the test gas that initially fills the shock tube.
    //
    // Start by unpacking the config structure.
    // We may be in the stand-alone nenzf1d program or
    // we may have been called by a Python program,
    // hence the wrapping of the config data.
    string gm1_filename = config.gm1_filename;
    auto gm1 = init_gas_model(gm1_filename);
    string gm2_filename = config.gm2_filename;
    string reactions_filename = config.reactions_filename;
    string reactions_filename2 = config.reactions_filename2;
    string[] species = config.species;
    double[] molef = config.molef;
    double T1 = config.T1;
    double p1 = config.p1;
    double Vs = config.Vs;
    double pe = config.pe;
    double Te = config.Te;
    double meq_throat = config.meq_throat;
    double ar = config.ar;
    double pp_ps = config.pp_ps;
    double C = config.C;
    double[] xi = config.xi;
    double[] di = config.di;
    bool check_supersonic_expansion = true;
    if (verbosityLevel >= 2) {
        writeln("Input data.");
        writeln("  gas-model-1= ", gm1_filename);
        writeln("  gas-model-2= ", gm2_filename);
        writeln("  reactions-file= ", reactions_filename);
        writeln("  reactions_file2= ", reactions_filename2);
        writeln("  species= ", species);
        writeln("  molef= ", molef);
        writeln("  T1= ", T1);
        writeln("  p1= ", p1);
        writeln("  Vs= ", Vs);
        writeln("  pe= ", pe);
        writeln("  Te= ", Te);
        writeln("  meq_throat= ", meq_throat);
        writeln("  ar= ", ar);
        writeln("  pp_ps= ", pp_ps);
        writeln("  C= ", C);
        writeln("  xi= ", xi);
        writeln("  di= ", di);
	writeln("  check_supersonic_expansion= ", check_supersonic_expansion);
    }
    //
    // --------------------------------------------
    // ANALYSIS PART A, SHOCK TUBE to NOZZLE-SUPPLY
    // --------------------------------------------
    //
    // Set up equilibrium-gas or frozen-kinetics flow analysis of shock tube.
    // Usually, this will involve a cea2 gas model that allows high-temperature effects.
    //
    void write_cea_state(ref const(GasState) gs, string fill="  ")
    {
        writefln("%spressure    %g kPa", fill, gs.p/1000.0);
        writefln("%sdensity     %g kg/m^3", fill, gs.rho);
        writefln("%stemperature %g K", fill, gs.T);
    }
    auto state5s = GasState(gm1);
    if (Te > 0.0 && pe > 0.0) {
        // We skip straight to setting the nozzle-supply condition.
        if (verbosityLevel >= 1) { writeln("Directly set nozzle-supply (stagnation) state 5s."); }
        state5s.p = pe;
        state5s.T = Te;
        gm1.update_thermo_from_pT(state5s);
        gm1.update_sound_speed(state5s);
        double H5s = gm1.internal_energy(state5s) + state5s.p/state5s.rho; // stagnation enthalpy
        if (verbosityLevel >= 1) {
            writefln("  H5s         %g MJ/kg", H5s/1.0e6);
            write_cea_state(state5s);
            if (verbosityLevel >= 3) { writeln("  state5s= ", state5s); }
        }
    } else {
        // Process the gas through initial and reflected shocks
        // to get the nozzle-supply condition.
        if (verbosityLevel >= 1) { writeln("Initial gas in shock tube (state 1)."); }
        auto state1 = GasState(gm1);
        state1.p = p1; state1.T = T1; state1.massf = [1.0,];
        gm1.update_thermo_from_pT(state1);
        gm1.update_sound_speed(state1);
        double H1 = gm1.internal_energy(state1) + state1.p/state1.rho;
        if (verbosityLevel >= 1) {
            write_cea_state(state1);
            writefln("  H1          %g MJ/kg", H1/1.0e6);
            if (verbosityLevel >= 3) { writeln("  state1: ", state1); }
        }
        //
        if (verbosityLevel >= 1) { writeln("Incident-shock process to state 2."); }
        auto state2 = GasState(gm1);
        double[2] velocities = normal_shock(state1, Vs, state2, gm1);
        double V2 = velocities[0];
        double Vg = velocities[1];
        if (verbosityLevel >= 1) {
            writefln("  V2          %g m/s", V2);
            writefln("  Vg          %g m/s", Vg);
            write_cea_state(state2);
            if (verbosityLevel >= 3) { writeln("  state2: ", state2); }
        }
        //
        if (verbosityLevel >= 1) { writeln("Reflected-shock process to state 5."); }
        GasState state5 = GasState(gm1);
        double Vr = reflected_shock(state2, Vg, state5, gm1);
        //
        if (verbosityLevel >= 1) { writeln("Isentropic relaxation to state 5s."); }
        state5s.copy_values_from(state5);
        // Entropy is set, then pressure is relaxed via an isentropic process.
        state5s.p = (pe > 0.0) ? pe : state5.p;
        double entropy5 = gm1.entropy(state5);
        gm1.update_thermo_from_ps(state5s, entropy5);
        double H5s = gm1.internal_energy(state5s) + state5s.p/state5s.rho; // stagnation enthalpy
        if (verbosityLevel >= 1) {
            writefln("  entropy     %g J/kg/K", entropy5);
            writefln("  H5s         %g MJ/kg", H5s/1.0e6);
            writefln("  H5s-H1      %g MJ/kg", (H5s-H1)/1.0e6);
            write_cea_state(state5s);
            if (verbosityLevel >= 3) { writeln("  state5s= ", state5s); }
        }
    } // end shock processing analysis
    //
    // At this point, we have the stagnation condition that drives the flow
    // through the nozzle expansion.
    //
    if (verbosityLevel >= 1) {
        writeln("Isentropic flow to throat to state 6 (Mach 1).");
    }
    double error_at_throat(double x)
    {
        // Returns Mach number error as pressure is changed.
        GasState state = GasState(gm1);
        double V = expand_from_stagnation(state5s, x, state, gm1);
        gm1.update_sound_speed(state);
        double err = (V/state.a) - 1.0;
        return err;
    }
    double x6 = 1.0;
    try {
        // Note that with the CEA2 calculations, we really cannot get the
        // following iteration to converge better then 1.0e-4.
        x6 = nm.secant.solve!(error_at_throat, double)(0.95, 0.90, 1.0e-4);
    } catch (Exception e) {
        writeln("Failed to find throat conditions iteratively.");
        throw e;
    }
    auto state6 = GasState(gm1);
    double V6 = expand_from_stagnation(state5s, x6, state6, gm1);
    double mflux6 = state6.rho * V6;  // mass flux per unit area, at throat
    if (verbosityLevel >= 1) {
        writefln("  V6          %g m/s", V6);
        writefln("  mflux6      %g", mflux6);
        write_cea_state(state6);
        if (verbosityLevel >= 3) { writeln("  state6= ", state6); }
    }
    //
    if (verbosityLevel >= 1) {
        writefln("Isentropic flow to slightly-expanded state 6e (Mach %g).", meq_throat);
    }
    double error_at_small_expansion(double x)
    {
        // Returns Mach number error as pressure is changed.
        GasState state = GasState(gm1);
        double V = expand_from_stagnation(state5s, x, state, gm1);
        gm1.update_sound_speed(state);
        double err = (V/state.a) - meq_throat;
        return err;
    }
    double x6e = x6;
    try {
        // Note that with the CEA2 calculations, we really cannot get the
        // following iteration to converge better then 1.0e-4.
        x6e = nm.secant.solve!(error_at_small_expansion, double)(0.95*x6, 0.90*x6, 1.0e-4);
    } catch (Exception e) {
        writeln("Failed to find slightly-expanded conditions iteratively.");
        throw e;
    }
    auto state6e = GasState(gm1);
    double V6e = expand_from_stagnation(state5s, x6e, state6e, gm1);
    double mflux6e = state6e.rho * V6e;  // mass flux per unit area, at slightly-expanded state
    // Compute the area-ratio at this slightly-expanded state.
    double ar6e = mflux6 / mflux6e;
    if (verbosityLevel >= 1) {
        writefln("  V6e         %g m/s", V6e);
        writefln("  mflux6e     %g", mflux6e);
        writefln("  ar6e        %g", ar6e);
        write_cea_state(state6e);
        if (verbosityLevel >= 3) { writeln("  state6e= ", state6e); }
    }
    //
    if (verbosityLevel >= 1) {
        writeln("Isentropic expansion to nozzle exit of given area (state 7).");
    }
    // The mass flux going through the nozzle exit has to be the same
    // as that going through the nozzle throat.
    double error_at_exit(double x)
    {
        // Returns mass_flux error as for a given exit pressure."
        auto state = GasState(gm1);
        double V = expand_from_stagnation(state5s, x, state, gm1);
        double mflux = state.rho * V * ar;
        return (mflux-mflux6)/mflux6;
    }
    // It appears that we need a pretty good starting guess for the pressure ratio.
    // Maybe a low value is OK.
    double x7 = x6;
    try {
        x7 = nm.secant.solve!(error_at_exit, double)(0.001*x6, 0.00005*x6, 1.0e-4, 1.0/state5s.p, 1.0);
    } catch (Exception e) {
        writeln("Failed to find exit conditions iteratively. Proceeding with expansion...");
        // Note that we have the throat conditions and may proceed
        // with a nonequilibrium expansion calculation.
        x7 = x6;
    }
    auto state7 = GasState(gm1);
    double V7 = expand_from_stagnation(state5s, x7, state7, gm1);
    double mflux7 = state7.rho * V7 * ar;
    auto state7_pitot = GasState(gm1);
    pitot_condition(state7, V7, state7_pitot, gm1);
    if (verbosityLevel >= 1) {
        writefln("  area_ratio  %g", ar);
        writefln("  V7          %g m/s", V7);
        write_cea_state(state7);
        writefln("  mflux7      %g", mflux7);
        writefln("  pitot7      %g kPa", state7_pitot.p/1000.0);
        if (verbosityLevel >= 3) { writeln("  state7= ", state7); }
        writeln("End of part A: shock-tube and frozen/eq nozzle analysis.");
    }
    //
    // -------------------------------------
    // ANALYSIS PART B, SUPERSONIC EXPANSION
    // -------------------------------------
    //
    auto gm2 = init_gas_model(gm2_filename);
    auto reactor = init_thermochemical_reactor(gm2, reactions_filename, reactions_filename2);
    double[10] reactor_params; // An array that passes extra parameters to the reactor.
    //
    void write_tp_state(ref const(GasState) gs, string fill="  ")
    {
        writefln("%spressure    %g kPa", fill, gs.p/1000.0);
        writefln("%sdensity     %g kg/m^3", fill, gs.rho);
        writefln("%stemperature %g K", fill, gs.T);
        foreach (i; 0 .. gs.T_modes.length) {
            string label = format("T_modes[%d]", i);
            writefln("%s%-12s%g K", fill, label, gs.T_modes[i]);
        }
        foreach (name; species) {
            string label = format("massf[%s]", name);
            writefln("%s%-12s%g", fill, label, gs.massf[gm2.species_index(name)]);
        }
    }
    //
    GasState init_tp_state_from_cea(ref const(GasState) state6e, ref const(GasState) state6, GasModel gm2){
        if (state6e.ceaSavedData is null) {
            throw new Exception("Failed to find ceaSavedData in state6e");
        }
        if (verbosityLevel >= 1) {
            writeln("Initializing gas state using CEA saved data");
        }
        if (verbosityLevel >= 2) {
            writeln("Start part B state mass fractions from CEA.");
            writeln("massf=", state6e.ceaSavedData.massf);
        }
        GasState gs = GasState(gm2);
        gs.p = state6e.p; gs.T = state6e.T;
        foreach (ref Tmode; gs.T_modes) { Tmode = state6e.T; }
        foreach (name; species) {
            gs.massf[gm2.species_index(name)] = state6.ceaSavedData.massf[name];
        }
        // CEA2 should be good to 0.01 percent.
        // If it is not, we probably have something significant wrong.
        scale_mass_fractions(gs.massf, 0.0, 0.0001);
        gm2.update_thermo_from_pT(gs);
        gm2.update_sound_speed(gs);
        return gs;
    }
    //
    GasState init_tp_state_from_eqgas(ref GasState state6e, EquilibriumGas gm_eq, GasModel gm2){
        // Although EquilibriumGas has officially one species, it keeps an internal
        // thermally perfect gas model for doing calculations. Here, do an apparently
        // pointless update_thermo_from_pT on state 6e to set the internal gas state
        // mass fractions, and then copy them into the new GasState, based on gm2.
        // @author: Nick Gibbons
        gm_eq.update_thermo_from_pT(state6e);
        GasState tpgs = GasState(gm_eq.savedGasState);
        ThermallyPerfectGasEquilibrium tpgm = gm_eq.savedGasModel;
        if (verbosityLevel >= 1) {
            writeln("Initializing gas state using equilibrium gas saved data");
        }
        if (verbosityLevel >= 2) {
            writeln("Start part B state mass fractions from ceq.");
            writeln("massf=", tpgs.massf);
        }
        GasState gs = GasState(gm2);
        gs.p = tpgs.p; gs.T = tpgs.T;
        foreach (ref Tmode; gs.T_modes) { Tmode = tpgs.T; }
        // We're assuming that the species order is the same in both models,
        // which could be a problem if the user doesn't keep them in order.
        foreach (name; species) {
            gs.massf[gm2.species_index(name)] = tpgs.massf[tpgm.species_index(name)];
        }
        gm2.update_thermo_from_pT(gs);
        gm2.update_sound_speed(gs);
        return gs;
    }
    //
    // We will continue with a non-equilibrium chemistry expansion
    // only if we have the correct gas models in play.
    auto gm_cea = cast(CEAGas) gm1;
    auto gm_eq = cast(EquilibriumGas) gm1;
    GasState gas0 = GasState(gm2);
    if (gm_cea !is null){
        try {
            gas0 = init_tp_state_from_cea(state6e, state6, gm2);
        } catch (Exception e) {
            writeln("Error in initialising tp state from cea");
            throw e;
        }
    } else if (gm_eq !is null) {
        gas0 = init_tp_state_from_eqgas(state6e, gm_eq, gm2);
    } else {
        writeln("Cannot continue with nonequilibrium expansion.");
        throw new Exception("  Gas model 1 is not of class CEAGas or EquilibriumGas.");
    }
    if (verbosityLevel >= 1) {
        writeln("Begin part B: continue supersonic expansion with finite-rate chemistry.");
    }
    //
    // Geometry of nozzle expansion.
    //
    auto diameter_schedule = new Schedule!double(xi, di);
    double x = xi[0];
    double d = diameter_schedule.interpolate_value(x);
    double area_at_throat = 0.25*d*d*std.math.PI; // for later normalizing the exit area
    //
    // Since we start slightly supersonic,
    // we need to determin where we are along the nozzle profile.
    double error_in_area_ratio(double x)
    {
        double d = diameter_schedule.interpolate_value(x);
        double area = 0.25*d*d*std.math.PI;
        return area/area_at_throat - ar6e;
    }
    // It appears that we need a pretty good starting guess for the pressure ratio.
    // Maybe a low value is OK.
    double xstart = xi[0];
    try {
        xstart = nm.secant.solve!(error_in_area_ratio, double)(xi[0], xi[1], 1.0e-9, xi[0], xi[$-1]);
    } catch (Exception e) {
        writeln("Failed to find starting position iteratively. Defaulting to xi[0]...");
        xstart = xi[0];
    }
    d = diameter_schedule.interpolate_value(xstart);
    double area_at_start = 0.25*d*d*std.math.PI;
    if (verbosityLevel >= 1) {
        writeln("Start position:");
        writefln("  x           %g m", xstart);
        writefln("  area-ratio  %g", area_at_start/area_at_throat);
    }
    //
    //
    // Make sure that we start the supersonic expansion with a velocity
    // that is slightly higher than the speed of sound for the thermally-perfect gas.
    // This sound speed for the finite-rate chemistry seems slightly higher than
    // the sound speed for the equilibrium gas.
    double v = fmax(V6e, 1.001*gas0.a);
    if (verbosityLevel >= 1) {
        writeln("Start condition:");
        writefln("  velocity    %g m/s", v);
        writefln("  sound-speed %g m/s", gas0.a);
        writefln("  (v-V6e)/V6e %g", (v-V6e)/V6e);
        write_tp_state(gas0);
        if (verbosityLevel >= 3) { writeln("gas0=", gas0); }
    }
    //
    string sample_header = "x(m) A(m**2) rho(kg/m**3) p(Pa) T(degK)";
    foreach (i; 0 .. gm2.n_modes) { sample_header ~= format(" T_modes[%d](degK)", i); }
    sample_header ~= " u(J/kg)";
    foreach (i; 0 .. gm2.n_modes) { sample_header ~= format(" u_modes[%d](J/kg)", i); }
    foreach (name; species) { sample_header ~= format(" massf_%s", name); }
    sample_header ~= " v(m/s) dt_suggest(s) mdot(kg/s)";
    //
    string sample_data(double x, double area, double v, ref const(GasState) gas, double dt_suggest)
    {
        string txt = format("%g %g %g %g %g", x, area, gas.rho, gas.p, gas.T);
        foreach (i; 0 .. gas.T_modes.length) { txt ~= format(" %g", gas.T_modes[i]); }
        txt ~= format(" %g", gas.u);
        foreach (i; 0 .. gas.u_modes.length) { txt ~= format(" %g", gas.u_modes[i]); }
        foreach (name; species) { txt ~= format(" %g", gas.massf[gm2.species_index(name)]); }
        txt ~= format(" %g %g %g", v, dt_suggest, gas.rho*v*area);
        return txt;
    }
    //
    double[2] eos_derivatives(ref const(GasState) gas0, double tol=0.0001)
    {
        // Finite difference evaluation, assuming that gas0 is valid state already.
        // Only the trans-rotational internal energy will be perturbed with any
        // other internal energy modes being left unperturbed.
        double p0 = gas0.p; double rho0 = gas0.rho; double u0 = gas0.u;
        //
        auto gas1 = GasState(gm2);
        gas1.copy_values_from(gas0);
        double drho = rho0 * tol; gas1.rho = rho0 + drho;
        gm2.update_thermo_from_rhou(gas1);
        double dpdrho = (gas1.p - p0)/drho;
        //
        gas1.copy_values_from(gas0);
        double du = u0 * tol; gas1.u = u0 + du;
        gm2.update_thermo_from_rhou(gas1);
        double dpdu = (gas1.p - p0)/du;
        //
        return [dpdrho, dpdu];
    }
    //
    // Supersonic expansion process for gas.
    //
    x = xstart;
    double area = area_at_start;
    double massflux_at_start = area*gas0.rho*v;
    double H_at_start = gm2.enthalpy(gas0) + 0.5*v*v;
    // We may use the Pitot pressure as a stopping criteria.
    double p_pitot = C * gas0.rho*v*v;
    //
    double dt_suggest = 1.0e-12;  // suggested starting time-step for chemistry
    if (verbosityLevel >= 2) {
        writeln(sample_header);
        writeln(sample_data(xi[0], area, v, gas0, dt_suggest));
        writeln("Start reactions...");
    }
    double t = 0;  // time is in seconds
    double x_end = config.x_end;  // Nozzle-exit position for terminating expansion.
    double t_final = config.t_final;
    double t_inc = config.t_inc;
    double t_inc_factor = config.t_inc_factor;
    double t_inc_max = config.t_inc_max;
    if (verbosityLevel >= 2) {
        writeln("Stepping parameters:");
        writefln("  x_end= %g", x_end);
        writefln("  t_final= %g", t_final);
        writefln("  t_inc= %g", t_inc);
        writefln("  t_inc_factor= %g", t_inc_factor);
        writefln("  t_inc_max= %g", t_inc_max);
    }
    //
    while ((x < x_end) &&
           (area < ar) &&
           (t < t_final) &&
           (p_pitot > pp_ps*state5s.p)) {
        // At the start of the step...
        double rho = gas0.rho; double T = gas0.T; double p = gas0.p;
        double u = gm2.internal_energy(gas0);
        //
        // Do the chemical increment.
        auto gas1 = GasState(gm2); // we need an update state
        gas1.copy_values_from(gas0);
        reactor(gas1, t_inc, dt_suggest, reactor_params);
        gm2.update_thermo_from_rhou(gas1);
        //
        double du_chem = gm2.internal_energy(gas1) - u;
        double dp_chem = gas1.p - p;
        if (verbosityLevel >= 3) {
            writeln("du_chem=", du_chem, " dp_chem=", dp_chem);
        }
        //
        // Update the independent variables for the end point of this step.
        // Simple, Euler update for the spatial stepper.
        double t1 = t + t_inc;
        double x1 = x + v*t_inc;
        double d1 = diameter_schedule.interpolate_value(x1);
        double area1 = 0.25*d1*d1*std.math.PI;
        //
        // Do the gas-dynamic accommodation after the chemical change.
        double darea = area1 - area;
        double etot = u + 0.5*v*v;
        double[2] df = eos_derivatives(gas1);
        double dfdr = df[0]; double dfdu = df[1];
        // double A = area+0.5*darea; // Slight improvement of mass-flux error.
        double A = area; // Original formulation of the linear constraint equations uses this.
        if (verbosityLevel >= 3) {
            writeln("dfdr=", dfdr, " dfdu=", dfdu);
            writeln("x=", x, " v=", v, " diam=", d, " A=", A, " dA=", darea);
        }
        // Linear solve to get the accommodation increments.
        //   [v*A,      rho*A,          0.0, 0.0    ]   [drho  ]   [-rho*v*dA                          ]
        //   [0.0,      rho*v,          1.0, 0.0    ] * [dv    ] = [-dp_chem                           ]
        //   [v*etot*A, (rho*etot+p)*A, 0.0, rho*v*A]   [dp_gda]   [-rho*v*A*du_chem -(rho*etot+p)*v*dA]
        //   [dfdr,     0.0,           -1.0, dfdu   ]   [du_gda]   [0.0                                ]
        //
        // Compute the accommodation increments using expressions from Maxima.
        double denom = A*(rho*rho*v*v - dfdr*rho*rho - dfdu*p);
        double drho = (A*(dp_chem - du_chem*dfdu)*rho*rho - darea*rho^^3*v*v) / denom;
        double dv = (A*(du_chem*dfdu - dp_chem)*rho +
                     darea*(dfdr*rho*rho + dfdu*p))*v / denom;
        double dp_gda = -((darea*dfdr*rho^^3 + A*dfdu*du_chem*rho*rho + darea*dfdu*p*rho)*v*v
                          - A*dfdr*dp_chem*rho*rho - A*dfdu*dp_chem*p) / denom;
        double du_gda = -(A*(du_chem*rho*rho*v*v - du_chem*dfdr*rho*rho - dp_chem*p)
                          + darea*p*rho*v*v) / denom;
        if (verbosityLevel >= 3) {
            writeln("drho=", drho, " dv=", dv, " dp_gda=", dp_gda, " du_gda=", du_gda);
            writefln("residuals= %g %g %g %g",
                     v*area*drho + rho*area*dv + rho*v*darea,
                     rho*v*dv + (dp_gda + dp_chem),
                     v*etot*drho*area + (rho*etot+p)*area*dv + rho*v*area*(du_gda + du_chem) + v*(rho*etot+p)*darea,
                     dfdr*drho - dp_gda + dfdu*du_gda);
        }
        // Add the gas-dynamic accommodation increments.
        gas1.rho += drho;
        double v1 = v + dv;
        double p1_check = gas1.p + dp_gda;
        gas1.u += du_gda;
        gm2.update_thermo_from_rhou(gas1);
        gm2.update_sound_speed(gas1);
        if (verbosityLevel >= 3) {
            writeln("At new point x1=", x1, " v1=", v1,
                    ": gas1.p=", gas1.p, " p1_check=", p1_check,
                    " rel_error=", fabs(gas1.p-p1_check)/p1_check);
        }
        // Have now finished the chemical and gas-dynamic update.
        if (verbosityLevel >= 2) {
            writeln(sample_data(x1, area1, v1, gas1, dt_suggest));
        }
	// Check that we are going down the supersonic-branch of the expansion.
	if (check_supersonic_expansion && (dv < 0.0 || du_gda > 0.0)) {
            string msg = format("We seem to be going down the subsonic branch dv=%g du_gda=%g", dv, du_gda);
	    msg ~= "\nMaybe the equilibrium-flow analysis did not continue far enough past M=1 at the throat.";
	    msg ~= "\nTry increasing the value of meq_throat.";
            msg ~= format(" Present value: meq_throat=%g", meq_throat);
	    throw new Nenzf1dAnalysisException(msg);
	}
        // House-keeping for the next step.
        v = v1; t = t1; x = x1; d =  d1; area = area1;
        gas0.copy_values_from(gas1);
        p_pitot = C * gas0.rho*v*v;
        t_inc = fmin(t_inc*t_inc_factor, t_inc_max);
    } // end while
    //
    // Build a structure of the computed data for printing upstairs in "main"
    Result result;
    result.x = x;
    result.area_ratio = area/area_at_throat;
    result.velocity = v;
    result.Mach_number = v/gas0.a;
    result.p_pitot = p_pitot;
    //
    double rayleigh_pitot = 0.0;
    try {
        // Note that we need the gas model to be able to compute entropy
        // in order to apply the Rayleigh-Pitot formula.
        // Not all gas models provide and implementation, so test by calling it first.
        GasState gs_pitot = GasState(gas0);
        auto s = gm2.entropy(gs_pitot);
        gm2.update_thermo_from_ps(gs_pitot, s);
        pitot_condition(gas0, v, gs_pitot, gm2);
        rayleigh_pitot = gs_pitot.p;
    } catch (Exception e) {
        debug {
            writeln("The calculation of the Rayleigh Pitot pressure has failed.");
            writeln("The exception message is: ", e.msg);
        }
        // Leave zero value for Rayleigh Pitot pressure and continue anyway.
    }
    result.rayleigh_pitot = rayleigh_pitot;
    result.pressure = gas0.p;
    result.density = gas0.rho;
    result.temperature = gas0.T;
    foreach (Tmode; gas0.T_modes) result.T_modes ~= Tmode;
    foreach (name; species) result.massf ~= gas0.massf[gm2.species_index(name)];
    gm2.update_trans_coeffs(gas0);
    result.viscosity = gas0.mu;
    //
    double massflux = area * gas0.rho * v;
    result.massflux_rel_err = fabs(massflux - massflux_at_start)/massflux_at_start;
    double H = gm2.enthalpy(gas0) + 0.5*v*v;
    result.enthalpy_rel_err = fabs(H - H_at_start)/H_at_start;
    if (result.rayleigh_pitot > 0.0) {
        result.pitot_rel_err = fabs(p_pitot - result.rayleigh_pitot)/p_pitot;
    }
    //
    if (verbosityLevel >= 1) { writeln("End."); }
    return result;
} // end analyse()

