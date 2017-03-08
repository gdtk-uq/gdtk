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
    GasModel gm = init_gas_model("../gas/sample-data/cea-air5species-gas-model.lua");
    GasState s1 = new GasState(gm, 100.0e3, 300.0);
    writeln("s1:", s1);
    writeln("Incident shock");
    GasState s2 = new GasState(s1);
    double[] vel_results = shock_ideal(s1, 3000.0, s2, gm);
    double vel2 = vel_results[0];
    double velg = vel_results[1];
    writeln("vel2=", vel2, "velg=", velg);
    writeln("s2:", s2);
/+
    #
    print "Reflected shock"
    s5 = s1.clone()
    Vr_b = reflected_shock(s2, Vg, s5)
    print "Vr_b=", Vr_b
    print "s5:"
    s5.write_state(sys.stdout)
    #
    print "Expand from stagnation"
    s6, V = expand_from_stagnation(0.0025, s5)
    print "V=", V, "Mach=", V/s6.a, "s6:"
    s6.write_state(sys.stdout)
    #
    print "Total condition"
    s7 = total_condition(s6, V)
    print "s7:"
    s7.write_state(sys.stdout)
    print "Pitot condition from state 6"
    s8 = pitot_condition(s6, V)
    print "pitot-p/total-p=", s8.p/s5.p, "s8:"
    s8.write_state(sys.stdout)
    #
    print "\nSteady, isentropic flow with area change."
    s8a = Gas({'Air':1.0})
    s8a.set_pT(1.0e5, 320.0)
    V8a = 1.001 * s8a.a
    V8b, s8b = steady_flow_with_area_change(s8a, V8a, 10.72) # something like M4 nozzle
    print "M=", V8b/s8b.a, "expected 4,  p2/p1=", s8b.p/s8a.p, "expected", 0.006586/0.5283
    V8b, s8b = steady_flow_with_area_change(s8a, V8a, 1.030) # slightly supersonic
    print "M=", V8b/s8b.a, "expected 1.2,  p2/p1=", s8b.p/s8a.p, "expected", 0.4124/0.5283
    V8a = 0.999 * s8a.a
    V8b, s8b = steady_flow_with_area_change(s8a, V8a, 2.9635) # sonic to M=0.2
    print "M=", V8b/s8b.a, "expected 0.2,  p2/p1=", s8b.p/s8a.p, "expected", 0.9725/0.5283
    V8a = 0.2 * s8a.a
    V8b, s8b = steady_flow_with_area_change(s8a, V8a, 1.3398/2.9635) # M=0.2 to M=0.5
    print "M=", V8b/s8b.a, "expected 0.5,  p2/p1=", s8b.p/s8a.p, "expected", 0.8430/0.9725
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
