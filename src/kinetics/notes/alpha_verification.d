/**
  * Author: Kyle Damm
  * Date: 26-Oct-2015
  * Title: Verification of the Alpha-QSS ODE update
  * The Alpha-Qss method is 2nd order accurate.
  * 
  * EXECUTION:    make -f make.alpha > /dev/null
  *               ./alpha_verify 
  */
import std.stdio;
import std.conv;
import std.math;
import nm.ridder;
import std.string;
import gas;
import gas.therm_perf_gas;
import kinetics.reaction_mechanism;
import kinetics.chemistry_update;
import std.algorithm;

void main()
{
    //============================== Analytical ==================================
    // Let's solve the analytical function for the steady state value we expect
    // the numerical solver to approximate
    double k_f = 3.1080121901939430e-05;
    double k_b = 5.4707058051099997e-07;
    double tMax = 60000.0;
    double c_0 = 4.54;

    auto zeroFun = delegate (double c_HI) {
        double tmp_a = sqrt(k_f * k_b);
        double tmp_b = 1.0 / (4.0* c_0 * tMax);
        double tmp_c = 2.0*c_0 + c_HI * (2.0*sqrt(k_b/k_f) - 1);
        double tmp_d = 2.0*c_0 - c_HI * (2.0*sqrt(k_b/k_f) + 1);
        return tmp_a - tmp_b * to!double(log(tmp_c/tmp_d));
    };

    double C0 = 7.0;
    double C1 = 7.15;
    double C_HI = solve!(zeroFun,double)(C0, C1);
    writefln("Analytical Solution = %16.16f", C_HI);

    //============================== Numerical ==================================
    // Now we will numerically solve the problem for a number of dt sizes and
    // check if we get the expected solution and convergence. Note that dt is
    // initially set to 4 seconds as oppose to 4000 seconds, this is due to
    // the alpha-Qss method returning NAN for dt values over 4 seconds.
    auto gm = new ThermallyPerfectGas("../sample-input/H2-I2-HI.lua");
    auto rmech = createReactionMechanism("../sample-input/H2-I2-inp.lua", gm);
    auto alphaStep = new AlphaQssStep(gm, rmech);

    auto gd = GasState(3, 1);
    gd.T[0] = 700.0;
    double c0 = 4.54;
    gd.p = 2.0*c0*R_universal*gd.T[0];
    double[] molef = [0.5, 0.5, 0.0];
    gm.molef2massf(molef, gd);
    gm.update_thermo_from_pT(gd);
    double[] conc0 = [c0, c0, 0.0];
    rmech.eval_rate_constants(gm, gd);

    tMax = 60000.0;
    double dtInit = 4.0;
    double[] dtVals = [dtInit];
    // To get an error ratio reduction of a factor
    // of 4, we would reduce the timestep by: 0.5;
    foreach ( i; 0..11 ) dtVals ~= 0.5*dtVals[$-1];
    double analyticalVal = 7.1420197868416215;
    double[] numVals;
    double[] err;

    foreach ( dt; dtVals ) {
        numVals ~= numericalEstimate(dt, tMax, conc0, alphaStep);
        err ~= analyticalVal - numVals[$-1];
    }

    writeln("|    dt    |  numerical value  |         error         |       ratio       |");
    writeln("|----------+-------------------+-----------------------+-------------------+");
    writefln("| %8.3f | %16.14f  | % 16.14e |                   |", dtVals[0], numVals[0], err[0]);
    foreach ( i; 1..dtVals.length ) {
        writefln("| %8.3f | %16.14f  | % 16.14e | % 16.14f |", dtVals[i], numVals[i], err[i], err[i-1]/err[i]);
    }

    auto f = File("rkf-verification-results.dat", "w");
    f.writeln("# dt   value    error    error-ratio");
    foreach ( i; 1..dtVals.length ) {
        f.writefln("%12.6f  %20.16e  %20.16e %20.16e", dtVals[i], numVals[i], err[i], err[i-1]/err[i]);
    }
    f.close();
}

  

double numericalEstimate(double dt, double tMax, double[] conc0, AlphaQssStep step)
{
    // Function for performing the ODE integration stepping, separated to allow for
    // multiple dt values to be evaluated successively.
    double t = 0.0;
    double[] conc1;
    conc1.length = conc0.length;
    double dtDummy;
    while ( (tMax - t) > 1.0e-9 ) {
        dt = min(dt, tMax - t);
        step(conc0, dt, conc1, dtDummy);
        t += dt;
        conc0 = conc1.dup;
    }
    return conc1[2];
}
