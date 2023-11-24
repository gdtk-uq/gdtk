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

import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
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

    this(in double[] R, bool separateElectronTemperature, int electronSpeciesIdx, int electronModeIdx) {
        _R = R.dup;
        _separateElecTemp = separateElectronTemperature;
        _eSpIdx = electronSpeciesIdx;
        _eMIdx = electronModeIdx;
    }

    /++
      Compute the pressure assuming density and temperature
      are up-to-date in GasState Q.
    +/
    @nogc override void update_pressure(ref GasState Q) const {
        number Rmix = heavyParticleGasConstant(Q);
        Q.p = Q.rho*Rmix*Q.T;
        if (_separateElecTemp) {
            Q.p += Q.rho*Q.massf[_eSpIdx]*_R[_eSpIdx]*Q.T_modes[_eMIdx];
        }
    }

    /++
      Compute the density assuming pressure and temperature
      are up-to-date in GasState Q.
    +/
    @nogc override void update_density(ref GasState Q) const {
        number Rmix = heavyParticleGasConstant(Q);
        number denom = Rmix*Q.T;
        if (_separateElecTemp) {
            denom += Q.massf[_eSpIdx]*_R[_eSpIdx]*Q.T_modes[_eMIdx];
        }
        Q.rho = Q.p/denom;
    }

    /++
      Compute the temperature assuming density and pressure
      are up-to-date in GasState Q.
      We can update the transrotaional temperature only if we
      assume the electron temperature is already known.
    +/
    @nogc override void update_temperature(ref GasState Q) const {
        number Rmix = heavyParticleGasConstant(Q);
        number p = Q.p;
        if (_separateElecTemp) {
            p -= Q.rho*Q.massf[_eSpIdx]*_R[_eSpIdx]*Q.T_modes[_eMIdx];
        }
        Q.T = p/(Rmix*Q.rho);
    }

private:
    double[] _R; /// specific gas constants in J/(kg.K)
    bool _separateElecTemp = false;
    int _eSpIdx = -1;
    int _eMIdx = -1;

    @nogc
    number heavyParticleGasConstant(ref GasState Q) const
    {
        number Rmix = 0.0;
        foreach (isp; 0 .. Q.massf.length) {
            if (isp == _eSpIdx && _separateElecTemp) continue;
            Rmix += Q.massf[isp]*_R[isp];
        }
        return Rmix;
    }
}

version(perf_gas_mix_eos_test) {
    import std.math;
    import std.stdio;
    import util.msg_service;
    int main() {
        double[] R = [297.0, 260.0]; // N2, O2
        auto pg = new PerfectGasMixEOS(R, false, -1, -1);
        auto gd = GasState(2, 1);
        gd.T = 300.0;
        gd.rho = 1.2;
        gd.massf[0] = 0.78;
        gd.massf[1] = 0.22;
        pg.update_pressure(gd);
        assert(isClose(gd.p, 103989.6, 1.0e-6), failedUnitTest());
        gd.p = 103989.6;
        gd.rho = 0.0;
        pg.update_density(gd);
        assert(isClose(gd.rho, 1.2, 1.0e-6), failedUnitTest());
        gd.rho = 1.2;
        gd.T = 0.0;
        pg.update_temperature(gd);
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());

        return 0;
    }
}
