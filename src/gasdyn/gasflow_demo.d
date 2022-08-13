/** gasflow_demo.d
 * Example of computing a reflected-shock tunnel flow, followed by other cases.
 *
 * PJ, 2017-03-09
 */

import std.stdio, std.math;
import gas;
import gasflow;
import idealgasflow; // for the oblique shock relations

void main(){
    writeln("Begin gasflow demo (reflected-shock tunnel calculation)");
    GasModel gm = init_gas_model("../gas/sample-data/cea-lut-air-version-test.lua");
    // GasModel gm = init_gas_model("../gas/sample-data/cea-air5species-gas-model.lua");
    GasState s1 = GasState(gm, 125.0e3, 300.0);
    writeln("    s1: ", s1);
    //
    writeln("Incident shock (ideal gas)");
    GasState s2 = GasState(s1);
    double[] velocities = shock_ideal(s1, 2414.0, s2, gm);
    double V2 = velocities[0]; double Vg = velocities[1];
    writeln("    V2=", V2, " Vg=", Vg);
    writeln("    s2: ", s2);
    //
    writeln("Incident shock (full gas model) with rho,T iteration");
    velocities = normal_shock(s1, 2414.0, s2, gm);
    V2 = velocities[0]; Vg = velocities[1];
    writeln("    V2=", V2, " Vg=", Vg);
    writeln("    s2: ", s2);
    //
    writeln("Incident shock (full gas model) with p,T iteration");
    velocities = normal_shock_1(s1, 2414.0, s2, gm);
    V2 = velocities[0]; Vg = velocities[1];
    writeln("    V2=", V2, " Vg=", Vg);
    writeln("    s2: ", s2);
    //
    writeln("Incident shock computed from pressure ratio (full gas model)");
    velocities = normal_shock_p2p1(s1, 7.314e6/125.0e3, s2, gm);
    double V1 = velocities[0]; V2 = velocities[1]; Vg = velocities[2];
    writeln("    V1=", V1);
    writeln("    V2=", V2, " Vg=", Vg);
    writeln("    s2: ", s2);
    //
    writeln("Reflected shock");
    GasState s5 = GasState(s1);
    double Vr_b = reflected_shock(s2, Vg, s5, gm);
    writeln("    Vr_b=", Vr_b);
    writeln("    s5:", s5);
    //
    writeln("Expand from stagnation (with ratio of pressure to match observation)");
    GasState s5s = GasState(s5);
    double V5s = expand_from_stagnation(s5, 34.37/59.47, s5s, gm);
    writeln("    V5s=", V5s, " Mach=", V5s/s5s.a);
    writeln("    s5s:", s5s);
    writeln("    (h5s-h1)=", gm.enthalpy(s5s) - gm.enthalpy(s1));
    //
    writeln("Expand to throat conditions (Mach 1.0001)");
    GasState s6 = GasState(s5s);
    double V6 = expand_to_mach(s5s, 1.0001, s6, gm);
    writeln("    V6=", V6, " Mach=", V6/s6.a);
    writeln("    s6:", s6);
    //
    writeln("Something like a Mach 4 nozzle.");
    GasState s7 = GasState(s6);
    double V7 = steady_flow_with_area_change(s6, V6, 27.0, s7, gm);
    writeln("    V7=", V7, " Mach=", V7/s7.a);
    writeln("    s7:", s7);
    //
    writeln("Total condition");
    GasState s8 = GasState(s7);
    total_condition(s7, V7, s8, gm);
    writeln("    s8:", s8);
    //
    writeln("Pitot condition from state 7");
    GasState s9 = GasState(s7);
    pitot_condition(s7, V7, s9, gm);
    writeln("    pitot-p/total-p=", s9.p/s8.p);
    writeln("    s9:", s9);
    //
    writeln("\nSteady, isentropic flow with area change. (more checks)");
    GasState s10a = GasState(gm, 1.0e5, 320.0); // ideal air, not high T
    gm.update_sound_speed(s10a);
    double V10a = 1.001 * s10a.a;
    GasState s10b = GasState(s10a);
    writeln("something like M4 nozzle with ideal air");
    double V10b = steady_flow_with_area_change(s10a, V10a, 10.72, s10b, gm);
    writeln("    M=", V10b/s10b.a, " expected 4");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.006586/0.5283);
    //
    writeln("slightly supersonic start");
    V10b = steady_flow_with_area_change(s10a, V10a, 1.030, s10b, gm);
    writeln("    M=", V10b/s10b.a, " expected 1.2");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.4124/0.5283);
    //
    writeln("sonic to M=0.2");
    V10a = 0.999 * s10a.a;
    V10b = steady_flow_with_area_change(s10a, V10a, 2.9635, s10b, gm);
    writeln("    M=", V10b/s10b.a, " expected 0.2");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.9725/0.5283);
    //
    writeln("M=0.2 to M=0.5");
    V10a = 0.2 * s10a.a;
    V10b = steady_flow_with_area_change(s10a, V10a, 1.3398/2.9635, s10b, gm);
    writeln("    M=", V10b/s10b.a, " expected 0.5");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.8430/0.9725);
    //
    writeln("\nFinite wave process along a cplus characteristic, stepping in pressure.");
    V1 = 0.0; s1.p = 1.0e5; s1.T = 320.0; // ideal air, not high T
    gm.update_sound_speed(s1);
    double Jplus = V1 + 2*s1.a/(1.4-1);
    V2 = finite_wave_dp(s1, V1, "cplus", 60.0e3, s2, gm);
    writeln("    V2=", V2);
    writeln("    s2:", s2);
    writeln("    ideal V2=", Jplus - 2*s2.a/(1.4-1));
    //
    writeln("\nFinite wave process along a cplus characteristic, stepping in velocity.");
    V1 = 0.0; s1.p = 1.0e5; s1.T = 320.0; // ideal air, not high T
    gm.update_sound_speed(s1);
    Jplus = V1 + 2*s1.a/(1.4-1);
    V2 = finite_wave_dv(s1, V1, "cplus", 125.0, s2, gm);
    writeln("    V2=", V2);
    writeln("    s2:", s2);
    writeln("    ideal Jplus=", Jplus, " actual Jplus=", V2 + 2*s2.a/(1.4-1));
    //
    double M1 = 1.5;
    writefln("\nOblique-shock demo for M1=%g.", M1);
    s1.p = 1.0e5; s1.T = 300.0; // ideal air, not high T
    gm.update_thermo_from_pT(s1);
    gm.update_sound_speed(s1);
    double beta = 45.0 * PI/180.0;
    writeln("    given beta(degrees)=", beta*180/PI);
    V1 = 1.5 * s1.a;
    writeln("    s1:", s1);
    double[] shock_results = theta_oblique(s1, V1, beta, s2, gm);
    double theta = shock_results[0]; V2 = shock_results[1];
    writeln("    theta=", theta);
    writeln("    V2=", V2);
    writeln("    s2:", s2);
    writeln("    c.f. ideal gas angle=", theta_obl(M1, beta));
    //
    writeln("Oblique shock angle from deflection.");
    double beta2 = beta_oblique(s1, V1, theta, gm);
    writeln("    beta2(degrees)=", beta2*180/PI);
    //
    M1 = 1.5;
    s1.p = 1.0e5; s1.T = 300.0; // ideal air, not high T
    gm.update_thermo_from_pT(s1);
    gm.update_sound_speed(s1);
    writefln("\nTaylor-Maccoll cone flow demo with M1=%g", M1);
    writeln("for M1=1.5, beta=49deg, expect theta=20deg from NACA1135.");
    V1 = M1 * s1.a;
    beta = 49.0 * PI/180.0;
    GasState state_c = GasState(s1);
    double[2] cone_results = theta_cone(s1, V1, beta, state_c, gm);
    double theta_c = cone_results[0]; double V_c = cone_results[1];
    writeln("    theta_c(deg)=", theta_c*180.0/PI);
    writeln("    expected 20deg, surface speed V_c=", V_c);
    writeln("    surface pressure coefficient=",
            (state_c.p - s1.p)/(0.5*s1.rho*V1*V1), " expected 0.385");
    writeln("    state_c:", state_c);
    //
    M1 = 1.5;
    s1.p = 1.0e5; s1.T = 300.0; // ideal air, not high T
    gm.update_thermo_from_pT(s1);
    gm.update_sound_speed(s1);
    writefln("\nTaylor-Maccoll cone flow demo with M1=%g", M1);
    writeln("for M1=1.5, beta=49.0404423512deg, expect theta=20deg from NACA1135.");
    V1 = M1 * s1.a;
    beta = 49.0404423512 * PI/180.0;
    cone_results = theta_cone(s1, V1, beta, state_c, gm);
    theta_c = cone_results[0]; V_c = cone_results[1];
    writeln("    theta_c(deg)=", theta_c*180.0/PI);
    writeln("    expected 20deg, surface speed V_c=", V_c);
    writeln("    surface pressure coefficient=",
            (state_c.p - s1.p)/(0.5*s1.rho*V1*V1), " expected 0.385");
    writeln("    state_c:", state_c);
    //
    M1 = 1.8;
    s1.p = 1.0e5; s1.T = 300.0; // ideal air, not high T
    gm.update_thermo_from_pT(s1);
    gm.update_sound_speed(s1);
    writefln("\nTaylor-Maccoll cone flow demo with M1=%g", M1);
    writeln("for M1=1.8, beta=45deg, expect theta=24deg from NACA1135.");
    V1 = M1 * s1.a;
    beta = 45.0 * PI/180.0;
    cone_results = theta_cone(s1, V1, beta, state_c, gm);
    theta_c = cone_results[0]; V_c = cone_results[1];
    writeln("    theta_c(deg)=", theta_c*180.0/PI);
    writeln("    expected 24deg, surface speed V_c=", V_c);
    writeln("    surface pressure coefficient=",
            (state_c.p - s1.p)/(0.5*s1.rho*V1*V1), " expected 0.466");
    writeln("    state_c:", state_c);
    //
    M1 = 1.5;
    writeln("\nConical shock from cone with half-angle 20deg in M1=", M1);
    V1 = M1 * s1.a;
    beta = beta_cone(s1, V1, 20.0*PI/180, gm);
    writeln("sigma(deg)=", beta*180/PI, " expected 49deg");
    //
    writeln("\nFinished gasflow demo.");
} // end main()
