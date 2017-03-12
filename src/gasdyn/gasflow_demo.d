/** gasflow_demo.d
 * Example of computing a reflected-shock tunnel flow.
 *
 * PJ, 2017-03-09
 */

import std.stdio, std.math;
import gas.gas_model;
import gasflow;

void main(){
    writeln("Begin gasflow demo (reflected-shock tunnel calculation)");
    GasModel gm = init_gas_model("../gas/sample-data/cea-air13species-gas-model.lua");
    GasState s1 = new GasState(gm, 125.0e3, 300.0);
    writeln("    s1: ", s1);
    //
    writeln("Incident shock (ideal gas)");
    GasState s2 = new GasState(s1);
    double[] vel_results = shock_ideal(s1, 2414.0, s2, gm);
    double vel2 = vel_results[0]; double velg = vel_results[1];
    writeln("    vel2=", vel2, " velg=", velg);
    writeln("    s2: ", s2);
    //
    writeln("Incident shock (full gas model)");
    vel_results = normal_shock(s1, 2414.0, s2, gm);
    vel2 = vel_results[0]; velg = vel_results[1];
    writeln("    vel2=", vel2, " velg=", velg);
    writeln("    s2: ", s2);
    //
    writeln("Incident shock computed from pressure ratio (full gas model)");
    vel_results = normal_shock_p2p1(s1, 7.314e6/125.0e3, s2, gm);
    double vel1 = vel_results[0]; vel2 = vel_results[1]; velg = vel_results[2];
    writeln("    vel1=", vel1);
    writeln("    vel2=", vel2, " velg=", velg);
    writeln("    s2: ", s2);
    //
    writeln("Reflected shock");
    GasState s5 = new GasState(s1);
    double velr_b = reflected_shock(s2, velg, s5, gm);
    writeln("    velr_b=", velr_b);
    writeln("    s5:", s5);
    //
    writeln("Expand from stagnation (with ratio of pressure to match observation)");
    GasState s5s = new GasState(s5);
    double vel5s = expand_from_stagnation(34.37/59.47, s5, s5s, gm);
    writeln("    vel5s=", vel5s, " Mach=", vel5s/s5s.a);
    writeln("    s5s:", s5s);
    writeln("    (h5s-h1)=", gm.enthalpy(s5s) - gm.enthalpy(s1)); 
    //
    writeln("Expand to throat conditions (Mach 1)");
    GasState s6 = new GasState(s5s);
    double vel6 = expand_to_mach(1.0, s5s, s6, gm);
    writeln("    vel6=", vel6, " Mach=", vel6/s6.a);
    writeln("    s6:", s6);
    //
    writeln("Something like a Mach 4 nozzle.");
    GasState s7 = new GasState(s6);
    double vel7 = steady_flow_with_area_change(s6, vel6, 27.0, s7, gm);
    writeln("    vel7=", vel7);
    writeln("    s7:", s7);
    //
    writeln("Total condition");
    GasState s8 = new GasState(s7);
    total_condition(s7, vel7, s8, gm);
    writeln("    s8:", s8);
    //
    writeln("Pitot condition from state 7");
    GasState s9 = new GasState(s7);
    pitot_condition(s7, vel7, s9, gm);
    writeln("    pitot-p/total-p=", s9.p/s8.p);
    writeln("    s9:", s9);
    //
    writeln("\nSteady, isentropic flow with area change. (more checks)");
    GasState s10a = new GasState(gm, 1.0e5, 320.0); // ideal air, not high T
    gm.update_sound_speed(s10a);
    double vel10a = 1.001 * s10a.a;
    GasState s10b = new GasState(s10a);
    writeln("something like M4 nozzle with ideal air");
    double vel10b = steady_flow_with_area_change(s10a, vel10a, 10.72, s10b, gm);
    writeln("    M=", vel10b/s10b.a, " expected 4");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.006586/0.5283);
    //
    writeln("slightly supersonic start");
    vel10b = steady_flow_with_area_change(s10a, vel10a, 1.030, s10b, gm); 
    writeln("    M=", vel10b/s10b.a, " expected 1.2");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.4124/0.5283);
    //
    writeln("sonic to M=0.2");
    vel10a = 0.999 * s10a.a;
    vel10b = steady_flow_with_area_change(s10a, vel10a, 2.9635, s10b, gm);
    writeln("    M=", vel10b/s10b.a, " expected 0.2");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.9725/0.5283);
    //
    writeln("M=0.2 to M=0.5");
    vel10a = 0.2 * s10a.a;
    vel10b = steady_flow_with_area_change(s10a, vel10a, 1.3398/2.9635, s10b, gm);
    writeln("    M=", vel10b/s10b.a, " expected 0.5");
    writeln("    p2/p1=", s10b.p/s10a.p, " expected ", 0.8430/0.9725);
/+
    #
    print "\nFinite wave process along a cplus characteristic, stepping in p."
    V1 = 0.0
    s9 = Gas({'Air':1.0})
    s9.set_pT(1.0e5, 320.0)
    Jplus = V1 + 2*s9.a/(1.4-1)
    V2, s10 = finite_wave_dp('cplus', V1, s9, 60.0e3)
    print "V2=", V2, "s10:"
    s10.write_state(sys.stdout)
    print "ideal V2=", Jplus - 2*s10.a/(1.4-1)
    #
    print "\nFinite wave process along a cplus characteristic, stepping in V."
    V1 = 0.0
    s9.set_pT(1.0e5, 320.0)
    Jplus = V1 + 2*s9.a/(1.4-1)
    V2, s10 = finite_wave_dv('cplus', V1, s9, 125.0)
    print "V2=", V2, "s10:"
    s10.write_state(sys.stdout)
    print "ideal Jplus=", Jplus, " actual Jplus=", V2 + 2*s10.a/(1.4-1)
    #
    M1 = 1.5
    print "\nOblique-shock demo for M1=%g." % M1
    from ideal_gas_flow import theta_obl
    s1.set_pT(100.0e3, 300.0)
    beta = 45.0 * math.pi/180
    V1 = 1.5 * s1.a
    print "s1:"
    s1.write_state(sys.stdout)
    theta, V2, s2 = theta_oblique(s1, V1, beta)
    print "theta=", theta, "V2=", V2, "s2:"
    s2.write_state(sys.stdout)
    print "c.f. ideal gas angle=", theta_obl(M1, beta)
    #
    print "Oblique shock angle from deflection."
    beta2 = beta_oblique(s1, V1, theta)
    print "beta2(degrees)=", beta2*180/math.pi
    #
    M1 = 1.5
    print "\nTaylor-Maccoll cone flow demo with M1=%g" % M1
    print "for M1=1.5, beta=49deg, expect theta=20deg from NACA1135."
    V1 = M1 * s1.a
    beta = 49.0 * math.pi/180
    theta_c, V_c, s_c = theta_cone(s1, V1, beta)
    print "theta_c(deg)=", theta_c*180.0/math.pi, "expected 20deg, surface speed V_c=", V_c
    print "surface pressure coefficient=", (s_c.p - s1.p)/(0.5*s1.rho*V1*V1), "expected 0.385"
    print "s_c:"
    s_c.write_state(sys.stdout)
    #
    M1 = 1.5
    print "\nTaylor-Maccoll cone flow demo with M1=%g" % M1
    print "for M1=1.5, beta=49.0404423512deg, expect theta=20deg from NACA1135."
    V1 = M1 * s1.a
    beta = 49.0404423512 * math.pi/180
    theta_c, V_c, s_c = theta_cone(s1, V1, beta)
    print "theta_c(deg)=", theta_c*180.0/math.pi, "expected 20deg, surface speed V_c=", V_c
    print "surface pressure coefficient=", (s_c.p - s1.p)/(0.5*s1.rho*V1*V1), "expected 0.385"
    print "s_c:"
    s_c.write_state(sys.stdout)
    #
    M1 = 1.8
    print "\nTaylor-Maccoll cone flow demo with M1=%g" % M1
    print "for M1=1.8, beta=45deg, theta=24deg from NACA1135."
    V1 = M1 * s1.a
    beta = 45.0 * math.pi/180
    theta_c, V_c, s_c = theta_cone(s1, V1, beta)
    print "theta_c(deg)=", theta_c*180.0/math.pi, "expected 24deg, surface speed V_c=", V_c
    print "surface pressure coefficient=", (s_c.p - s1.p)/(0.5*s1.rho*V1*V1), "expected 0.466"
    print "s_c:"
    s_c.write_state(sys.stdout)
    #
    M1 = 1.5
    print "\nConical shock from cone with half-angle 20deg in M1=", M1
    V1 = M1 * s1.a
    beta = beta_cone(s1, V1, 20.0*math.pi/180)
    print "sigma(deg)=", beta*180/math.pi, "expected 49deg"
    #
    print "Done."
+/
	
    writeln("Finished gasflow demo.");
} // end main()
