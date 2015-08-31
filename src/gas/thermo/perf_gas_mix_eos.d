/**
 * perf_gas_mix_eos.d
 * Implements a mixture of perfect gases equation of state.
 * This module provides simple functions for the
 * the p-v-T behaviour of a mixture perfect gases.
 * 
 * Author: Rowan G. and Peter J.
 * Version: 2014-09-07 -- first cut
 */

module gas.thermo.perf_gas_mix_eos;

import gas.gas_model;
import gas.thermo.pvt_eos;

/++
 PerfectGasMixEOS is a thermal equation of state.
 
 The perfect gas mixture model assumes point masses and
 perfectly elastic collisions.
+/
class PerfectGasMixEOS : PVT_EOS {
public:
    this(in double[] R) {
	_R = R.dup;
    }

    /++
      Compute the pressure assuming density and temperature
      are up-to-date in GasState Q.
    +/
    @nogc override void update_pressure(ref GasState Q) const {
	double Rmix = mass_average(Q, _R);
	Q.p = Q.rho*Rmix*Q.T[0];
    }

    /++
      Compute the density assuming pressure and temperature
      are up-to-date in GasState Q.
    +/
    @nogc override void update_density(ref GasState Q) const {
	double Rmix = mass_average(Q, _R);
	Q.rho = Q.p/(Rmix*Q.T[0]);
    }

    /++
      Compute the temperature assuming density and pressure
      are up-to-date in GasState Q.
    +/
    @nogc override void update_temperature(ref GasState Q) const {
	double Rmix = mass_average(Q, _R);
	Q.T[0] = Q.p/(Rmix*Q.rho);
    }

private:
    double[] _R; /// specific gas constants in J/(kg.K)
}	  

unittest {
    import std.math;
    import std.stdio;
    import util.msg_service;
    double[] R = [297.0, 260.0]; // N2, O2
    auto pg = new PerfectGasMixEOS(R);
    auto gd = GasState(2, 1);
    gd.T[0] = 300.0;
    gd.rho = 1.2;
    gd.massf[0] = 0.78;
    gd.massf[1] = 0.22;
    pg.update_pressure(gd);
    assert(approxEqual(gd.p, 103989.6), failedUnitTest());
    gd.p = 103989.6;
    gd.rho = 0.0;
    pg.update_density(gd);
    assert(approxEqual(gd.rho, 1.2), failedUnitTest());
    gd.rho = 1.2;
    gd.T[0] = 0.0;
    pg.update_temperature(gd);
    assert(approxEqual(gd.T[0], 300.0), failedUnitTest());
}



	   


