/**
 * Interface and implementation of energy exchange system.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.energy_exchange_system;

import std.stdio;
import std.math;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;

import kinetics.energy_exchange_mechanism;
import kinetics.relaxation_time;

interface EnergyExchangeSystem {
    @nogc void evalRelaxationTimes(in GasState gs);
    @nogc void evalRates(in GasState gs, ref number[] rates);
}

class TwoTemperatureEnergyExchange : EnergyExchangeSystem {
public:

    this(GasModel gmodel)
    {
        mGmodel = gmodel;
        mGsEq = new GasState(gmodel);
        mVibRelaxers.length = 1;
        mVibRelaxers[0] = 0;
        mHeavyParticles.length = 2;
        mHeavyParticles[0] = 0; mHeavyParticles[1] = 1;
        mMolef.length = 2;
        mVT.length = 1;
        mVT[0].length = 2;
        double M_N2 = 2.80134000e-02 * 1000.0;
        double M_N =  1.40067000e-02 * 1000.0;
        double mu0 = (M_N2*M_N2)/(M_N2 + M_N2);
        double mu1 = (M_N2*M_N)/(M_N2 + M_N);
        double theta_v = 3354.0;
        // N2-N2
        double a = 1.16e-3*sqrt(mu0)*pow(theta_v, 4./3.);
        double b = 0.015*pow(mu0, 1./4.);
        writeln("a= ", a, " b= ", b);
        auto rt = new MillikanWhiteVT(0, a, b, mu0);
        mVT[0][0] = new LandauTellerVT(0, 0, rt, gmodel);
        a = 1.16e-3*sqrt(mu1)*pow(theta_v, 4./3.);
        b = 0.015*pow(mu1, 1./4.);
        writeln("a= ", a, " b= ", b);
        rt = new MillikanWhiteVT(1, a, b, mu1);
        mVT[0][1] = new LandauTellerVT(0, 0, rt, gmodel);
    }
    
    @nogc
    void evalRelaxationTimes(in GasState gs)
    {
        mGmodel.massf2molef(gs, mMolef);
        foreach (p; mVibRelaxers) {
            foreach (q; mHeavyParticles) {
                mVT[p][q].evalRelaxationTime(gs, mMolef);
            }
        }
    }
    
    @nogc
    void evalRates(in GasState gs, ref number[] rates)
    {
        // Compute a star state at transrotational equilibrium.
        mGsEq.copy_values_from(gs);
        mGsEq.T_modes[0] = gs.T;
        mGmodel.update_thermo_from_pT(mGsEq);
        number rate = 0.0;
        // Compute rate change from VT exchange
        foreach (p; mVibRelaxers) {
            foreach (q; mHeavyParticles) {
                rate += mVT[p][q].rate(gs, mGsEq, mMolef);
            }
        }
        rates[0] = rate;
        // TODO: Compute rate change from ET exchange
        // TODO: Compute rate change from chemistry coupling.
    }

private:
    int[] mVibRelaxers;
    int[] mHeavyParticles;
    number[] mMolef;
    GasModel mGmodel;
    GasState mGsEq;
    EnergyExchangeMechanism[][] mVT;
    // EnergyExchangeMechanism[] mET;
    // EnergyExchangeMechanism[] mChemCoupling;
    
}

