// main.d for nenzf1d
// Reflected shock tube analysis followed by supersonic nozzle expansion.
//
// Peter J.
// School of Mechanical and Mining Engineering
// The University of Queensland
//
// 2020-09-26 Initial code built from Python prototype.
// 2020-10-16 Ready for use with 1T thermally-perfect gas model, I believe.
// 2020-12-09 Allow 2T gas models.

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
import kinetics;
import gasflow;


int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.
    // Be careful with the usageMsg string; it has embedded newline characters.
    string usageMsg = "Usage: nenzf1d <input-file>
Options:
   --verbosity=<int>   0 == very terse output
                       1 == key results printed
                       2 == echo input as well as printing more detailed results
                       3 == debug printing as well
   --help              Print this help message.";
    if (args.length < 2) {
        writeln("Too few arguments. You need to specify the input file.");
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    int verbosityLevel = 1; // default to commenting on major steps
    bool helpWanted = false;
    try {
        getopt(args,
               "verbosity", &verbosityLevel,
               "help", &helpWanted
               );
    } catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed: ");
        args = args[1 .. $]; // Dispose of program name in first argument.
        foreach (myarg; args) writeln("    arg: ", myarg);
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    if (verbosityLevel >= 1) {
        writeln("NENZF1D 2020-12-10: shock-tunnel with nonequilibrium nozzle flow.");
    }
    if (verbosityLevel >= 2) {
        writeln("Revision: PUT_REVISION_STRING_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
    }
    if (helpWanted) {
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 0;
        return exitFlag;
    }
    // Read the input file.
    string inputFile = args[1].strip();
    if (!exists(inputFile)) {
        writeln("Did not find input-file: ", inputFile);
        exitFlag = 2;
        return exitFlag;
    }
    // Extract our job parameters from the YAML input file.
    auto config = dyaml.Loader.fromFile(inputFile).load();
    if (verbosityLevel >= 1) {
        writeln("  "~config["title"].as!string);
    }
    //
    // The first gas model is for the shock-tube analysis,
    // assuming frozen reactions or full equilibrium.
    string gm1_filename = config["gas-model-1"].as!string;
    auto gm1 = init_gas_model(gm1_filename);
    //
    // The second gas model is for the expansion process that has finite-rate 1T chemistry.
    string gm2_filename = config["gas-model-2"].as!string;
    string reactions_filename = config["reactions"].as!string;
    string reactions_filename2 = "";
    try {
        reactions_filename2 = config["reactions-file2"].as!string;
    } catch (YAMLException e) {
        // We cannot set the second reactions file so assume 1T chemistry.
    }
    //
    string[] species;
    foreach(string name; config["species"]) { species ~= name; }
    double[] molef;
    foreach(name; species) {
        double mf = 0.0;
        try { mf = to!double(config["molef"][name].as!string); } catch (YAMLException e) {}
        molef ~= mf;
    }
    //
    // Initial gas state in shock tube.
    double T1 = to!double(config["T1"].as!string);
    double p1 = to!double(config["p1"].as!string);
    double Vs = to!double(config["Vs"].as!string);
    // Observed relaxation pressure for reflected-shock, nozzle-supply region.
    // A value of 0.0 indicates that we should use the ideal shock-reflection pressure.
    double pe = 0.0;
    try { pe = to!double(config["pe"].as!string); } catch (YAMLException e) {}
    // Mach number (in equilibrium gas) at nozzle throat.
    double meq_throat = 1.0;
    // Nominally, it would be 1.0 for sonic flow.
    // It may be good to expand a little more so that,
    // on changing to the frozen-gas sound speed in the nonequilibrium gas model,
    // the flow remains slightly supersonic.
    try { meq_throat = to!double(config["meq_throat"].as!string); } catch (YAMLException e) {}
    // Nozzle-exit area ratio for terminating expansion.
    double ar = 1.0;
    try { ar = to!double(config["ar"].as!string); } catch (YAMLException e) {}
    // Alternatively, we might stop on pPitot/pSupply becoming less than pp_ps.
    double pp_ps = 0.0;
    try { pp_ps = to!double(config["pp_ps"].as!string); } catch (YAMLException e) {}
    // pPitot = C * rho*V^^2
    // A value of C=1.0 is a good default.
    // In a number of sphere simulations, for flows representative of T4 flow conditions,
    // values of pPitot/(rho*v^^2) appeared to be in the range 0.96 to 1.0.
    double C = 1.0;
    try { C = to!double(config["C"].as!string); } catch (YAMLException e) {}
    // Nozzle x,diameter schedule.
    double[] xi; foreach(string val; config["xi"]) { xi ~= to!double(val); }
    double[] di; foreach(string val; config["di"]) { di ~= to!double(val); }
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
        writeln("  meq_throat= ", meq_throat);
        writeln("  ar= ", ar);
        writeln("  pp_ps= ", pp_ps);
        writeln("  C= ", C);
        writeln("  xi= ", xi);
        writeln("  di= ", di);
    }
    //
    // ---------------------------
    // ANALYSIS PART A, SHOCK TUBE
    // ---------------------------
    //
    void write_cea_state(GasState gs, string fill="  ")
    {
        writefln("%spressure    %g kPa", fill, gs.p/1000.0);
        writefln("%sdensity     %g kg/m^3", fill, gs.rho);
        writefln("%stemperature %g K", fill, gs.T);
    }
    // Set up equilibrium-gas flow analysis of shock tube.
    // Let's assume a cea2 gas model.
    if (verbosityLevel >= 1) { writeln("Initial gas in shock tube (state 1)."); }
    auto state1 = new GasState(gm1);
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
    auto state2 = new GasState(gm1);
    double[2] velocities = normal_shock(state1, Vs, state2, gm1);
    double V2 = velocities[0];
    double Vg = velocities[1];
    if (verbosityLevel >= 1) {
        writefln("  V2          %g km/s", V2/1000.0);
        writefln("  Vg          %g km/s", Vg/1000.0);
        write_cea_state(state2);
        if (verbosityLevel >= 3) { writeln("  state2: ", state2); }
    }
    //
    if (verbosityLevel >= 1) { writeln("Reflected-shock process to state 5."); }
    GasState state5 = new GasState(gm1);
    double Vr = reflected_shock(state2, Vg, state5, gm1);
    //
    if (verbosityLevel >= 1) { writeln("Isentropic relaxation to state 5s."); }
    auto state5s = new GasState(gm1);
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
    //
    if (verbosityLevel >= 1) {
        writefln("Isentropic flow to throat to state 6 (Mach %g).", meq_throat);
    }
    double error_at_throat(double x)
    {
        // Returns Mach number error as pressure is changed.
        GasState state = new GasState(gm1);
        double V = expand_from_stagnation(state5s, x, state, gm1);
        gm1.update_sound_speed(state);
        return (V/state.a) - meq_throat;
    }
    double x6 = 1.0;
    try {
        x6 = nm.secant.solve!(error_at_throat, double)(0.95, 0.90, 1.0e-4);
    } catch (Exception e) {
        writeln("Failed to find throat conditions iteratively.");
        writeln("  Exception message: %s", e.msg);
        exitFlag = 2;
        return exitFlag;
    }
    auto state6 = new GasState(gm1);
    double V6 = expand_from_stagnation(state5s, x6, state6, gm1);
    double mflux6 = state6.rho * V6;  // mass flux per unit area, at throat
    if (verbosityLevel >= 1) {
        writefln("  V6          %g km/s", V6);
        writefln("  mflux6      %g", mflux6);
        write_cea_state(state6);
        if (verbosityLevel >= 3) { writeln("  state6= ", state6); }
    }
    //
    if (verbosityLevel >= 1) { writeln("Isentropic expansion to nozzle exit of given area (state 7)."); }
    // The mass flux going through the nozzle exit has to be the same
    // as that going through the nozzle throat.
    double error_at_exit(double x)
    {
        // Returns mass_flux error as for a given exit pressure."
        auto state = new GasState(gm1);
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
        writeln("Failed to find exit conditions iteratively.");
        // Note that we have the throat conditions and may proceed
        // with a nonequilibrium expansion calculation.
        x7 = x6;
    }
    auto state7 = new GasState(gm1);
    double V7 = expand_from_stagnation(state5s, x7, state7, gm1);
    double mflux7 = state7.rho * V7 * ar;
    auto state7_pitot = new GasState(gm1);
    pitot_condition(state7, V7, state7_pitot, gm1);
    if (verbosityLevel >= 1) {
        writefln("  area_ratio  %g", ar);
        writefln("  V7          %g km/s", V7/1000.0);
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
    // We will continue with a non-equilibrium chemistry expansion
    // only if we have the correct gas models in play.
    auto gm_cea = cast(CEAGas) gm1;
    if (gm_cea is null) {
        writeln("Cannot continue with nonequilibrium expansion.");
        writeln("  Gas model 1 is not of class CEAGas.");
        exitFlag = 2;
        return exitFlag;
    }
    auto gm2 = init_gas_model(gm2_filename);
    auto reactor = init_thermochemical_reactor(gm2, reactions_filename, reactions_filename2);
    double[10] reactor_params; // An array that passes extra parameters to the reactor.
    //
    void write_tp_state(GasState gs, string fill="  ")
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
    if (state6.ceaSavedData is null) {
        exitFlag = 3;
        return exitFlag;
    }
    if (verbosityLevel >= 2) {
        writeln("Throat state mass fractions from CEA.");
        writeln("massf=", state6.ceaSavedData.massf);
    }
    GasState gas0 = new GasState(gm2);
    gas0.p = state6.p; gas0.T = state6.T;
    foreach (ref Tmode; gas0.T_modes) { Tmode = state6.T; }
    try {
        foreach (name; species) {
            gas0.massf[gm2.species_index(name)] = state6.ceaSavedData.massf[name];
        }
        // CEA2 should be good to 0.01 percent.
        // If it is not, we probably have something significant wrong.
        scale_mass_fractions(gas0.massf, 0.0, 0.0001);
    } catch (Exception e) {
        writeln("Failed to transfer mass fractions to finite-rate-chemistry gas.");
        writeln(e.msg);
        exitFlag = 4;
        return exitFlag;
    }
    gm2.update_thermo_from_pT(gas0);
    gm2.update_sound_speed(gas0);
    if (verbosityLevel >= 1) {
        writeln("Begin part B: supersonic expansion with finite-rate chemistry.");
    }
    //
    // Make sure that we start the supersonic expansion with a velocity
    // that is slightly higher than the speed of sound for the thermally-perfect gas.
    // This speed seems slightly higher than for the equilibrium gas.
    double v = 1.001 * gas0.a;
    if (verbosityLevel >= 1) {
        writeln("Throat condition:");
        writefln("  velocity    %g km/s", v/1000.0);
        writefln("  sound-speed %g km/s", gas0.a/1000.0);
        writefln("  (v-V6)/V6   %g", (v-V6)/V6);
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
    string sample_data(double x, double area, double v, GasState gas, double dt_suggest)
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
    double[2] eos_derivatives(ref GasState gas0, double tol=0.0001)
    {
        // Finite difference evaluation, assuming that gas0 is valid state already.
        // Only the trans-rotational internal energy will be perturbed with any
        // other internal energy modes being left unperturbed.
        double p0 = gas0.p; double rho0 = gas0.rho; double u0 = gas0.u;
        //
        auto gas1 = new GasState(gm2);
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
    // Supersonic expansion.
    auto diameter_schedule = new Schedule(xi, di);
    double x = xi[0];
    double d = diameter_schedule.interpolate_value(x);
    double area = 0.25*d*d*std.math.PI;
    double area_at_throat = area; // for later normalizing the exit area
    double massflux_at_throat = area*gas0.rho*v;
    double H_at_throat = gm2.enthalpy(gas0) + 0.5*v*v;
    // We may use the Pitot pressure as a stopping criteria.
    double p_pitot = C * gas0.rho*v*v;
    //
    double dt_suggest = 1.0e-12;  // suggested starting time-step for chemistry
    double dt_therm = dt_suggest;
    if (verbosityLevel >= 2) {
        writeln(sample_header);
        writeln(sample_data(xi[0], area, v, gas0, dt_suggest));
        writeln("Start reactions...");
    }
    double t = 0;  // time is in seconds
    double x_end = xi[$-1];  // Nozzle-exit position for terminating expansion.
    try { x_end = to!double(config["x_end"].as!string); } catch (YAMLException e) {}
    double t_final = 2.0 * x_end / v;  // long enough to convect past exit
    try { t_final = to!double(config["t_final"].as!string); } catch (YAMLException e) {}
    double t_inc = 1.0e-10; // start small
    try { t_inc = to!double(config["t_inc"].as!string); } catch (YAMLException e) {}
    double t_inc_factor = 1.0001;
    try { t_inc_factor = to!double(config["t_inc_factor"].as!string); } catch (YAMLException e) {}
    double t_inc_max = 1.0e-7;
    try { t_inc_max = to!double(config["t_inc_max"].as!string); } catch (YAMLException e) {}
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
        auto gas1 = new GasState(gm2); // we need an update state
        gas1.copy_values_from(gas0);
        reactor(gas1, t_inc, dt_suggest, dt_therm, reactor_params);
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
        // House-keeping for the next step.
        v = v1; t = t1; x = x1; d =  d1; area = area1;
        gas0.copy_values_from(gas1);
        p_pitot = C * gas0.rho*v*v;
        t_inc = fmin(t_inc*t_inc_factor, t_inc_max);
    } // end while
    //
    writeln("Exit condition:");
    writefln("  x           %g m", x);
    writefln("  area-ratio  %g", area/area_at_throat);
    writefln("  velocity    %g km/s", v/1000.0);
    writefln("  Mach        %g", v/gas0.a);
    writefln("  p_pitot     %g kPa (C.rho.V^2)", p_pitot/1000);
    double rayleigh_pitot = 0.0;
    try {
        // We need the gas model to be able to compute entropy in order
        // to apply the Rayleigh-Pitot formula.
        GasState gs_pitot = new GasState(gas0);
        pitot_condition(gas0, v, gs_pitot, gm2);
        rayleigh_pitot = gs_pitot.p;
        writefln("  p_pitot     %g kPa (Rayleigh-Pitot, frozen)", rayleigh_pitot/1000);
    } catch (Exception e) { /* Do nothing. */ }
    write_tp_state(gas0);
    gm2.update_trans_coeffs(gas0);
    writefln("  viscosity   %g Pa.s", gas0.mu);
    //
    writeln("Expansion error-indicators:");
    double massflux = area * gas0.rho * v;
    writefln("  relerr-mass %g", fabs(massflux - massflux_at_throat)/massflux_at_throat);
    double H = gm2.enthalpy(gas0) + 0.5*v*v;
    writefln("  relerr-H    %g", fabs(H - H_at_throat)/H_at_throat);
    if (rayleigh_pitot > 0.0) {
        writefln("  relerr-pitot %g", fabs(p_pitot - rayleigh_pitot)/p_pitot);
    }
    //
    if (verbosityLevel >= 1) { writeln("End."); }
    return exitFlag;
} // end main
