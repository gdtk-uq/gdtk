/*
 * file: idealgasflow_demo.d
 * author: Momar Hughes and Peter J.
 * date: 5 Apr 2015 initial code
 *       2016-11-19 added Taylor-Maccoll demo
 */

import std.stdio, std.math;
import gasdyn.idealgasflow;

void main(){
    writeln("Begin gasdynamic demo...");
    double M = 2.4, g = 1.4;
    writefln("Use M=%g and gamma=%g",M,g);
        
    writeln("Test: T0_T...");
    writeln("...returns T0/T = ",T0_T(M,g));
    writeln("Literature value = 2.152");
        
    writeln("Test: p0_p...");
    writeln("...returns p0/p = ",p0_p(M,g));
    writeln("Literature value = 14.620");
        
    writeln("Test: NuFromM...");
    writeln("...returns nu = ",PM1(M,g));
    writeln("Literature value = 0.6413");
    writeln("...provided M=0.8:");
    try {
        PM1(0.8,g);
    } catch (Exception e) {
        writeln("Caught subsonic value.");
    }
        
    double nu = 0.6413479572;
    writeln("Test: MFromNu...");
    writeln("...returns M = ",PM2(nu,g));
    writeln("Literature value = 2.4");
    writeln("...provided nu=-0.5:");
    try {
        PM2(-0.5,g);
    } catch (Exception e) {
        writeln("Caught subsonic value.");
    }
        
    writeln("Test: MachAngle...");
    writeln("...returns mu = ",MachAngle(M));
    writeln("Literature value = 0.430");
    writeln("...provided M=0.8:");
    try {
        MachAngle(0.8);
    } catch (Exception e) {
        writeln("Caught subsonic value.");
    }

    M = 2.0;
    writeln("Oblique shock relations may not quite match (data is from chart)...");
    double beta = 44.0*PI/180.0; double theta = 14.0*PI/180.0; // from chart, M=2
    writefln("Computed: M1=%g, theta(beta=%g)=%g, beta(theta=%g)=%g",
             M, beta, theta_obl(M, beta), theta, beta_obl(M, theta));
    writeln("Conditions behind shock:");
    writefln("M2=%g, expected 1.482 (from chart, 14 degree deflection)",
             M2_obl(M, beta, theta));
    writefln("Computed: T2/T1=%g, p2/p1=%g, r2/r1=%g",
             T2_T1_obl(M, beta), p2_p1_obl(M, beta), r2_r1_obl(M, beta));
    writeln("Expected: T2/T1=1.249, p2/p1=2.088, r2/r1=1.673",
            " (approx. normal-shock table M=1.390)");
    writefln("V2/V1=%g, p02/p01=%g", V2_V1_obl(M, beta), p02_p01_obl(M, beta));
    writeln("Expected: V2/V1=0.8304=sin(B)/sin(B-d)*r1/r2");
    writeln("");

    double M1 = 1.5; double p1 = 100.0e3; double T1 = 300.0;
    double R = 287.1; g = 1.4; double rho1 = p1/(R*T1);
    writefln("Taylor-Maccoll cone flow demo with M1=%g", M1);
    writeln("for M1=1.5, beta=49deg, expect theta=20deg from NACA1135.");
    double a1 = sqrt(1.4*287*T1);
    double V1 = M1 * a1;
    beta = 49.0 * PI/180;
    double[] results = theta_cone(V1, p1, T1, beta);
    double theta_c=results[0];  double V_c=results[1];
    double p_c=results[2]; double T_c=results[3];
    writeln("theta_c(deg)=", theta_c*180.0/PI, " expected 20deg, surface speed V_c=", V_c);
    writeln("surface pressure coefficient=", (p_c - p1)/(0.5*rho1*V1*V1), " expected 0.385");
    writefln("p_c: %g, T_c: %g", p_c, T_c);
    writeln("");
    writeln("Conical shock from cone with half-angle 20deg in M1=", M1);
    beta = beta_cone(V1, p1, T1, 20.0*PI/180);
    writeln("sigma(deg)=", beta*180/PI, " expected 49deg");
    writeln("Repeat above test, but call beta_cone2()");
    beta = beta_cone2(M1, 20.0*PI/180);
    writeln("sigma(deg)=", beta*180/PI, " expected 49deg");
        
    writeln("Finished gasdynamic demo.");
} // end main()
