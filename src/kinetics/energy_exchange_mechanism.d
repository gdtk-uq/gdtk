/**
 * Interface and implementation for energy exchange mechanisms.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.energy_exchange_mechanism;

import std.math;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;


interface EnergyExchangeMechanism {
    @nogc number rate(in GasState gs, in GasState gsEq, number[] molef);
}

class LandauTellerVT {
public:

    @nogc
    override number rate(in GasState gs, in GasState gsEq, number[] molef)
    {
        number evStar = mGmodel.energyPerSpeciesInMode(gsEq, m_p, m_mode);
        number ev = mGmodel.energyPerSpeciesInMode(gs, m_p, m_mode);
        number tau = mRT.eval(gs, molef);
        // NOTE 1. tau has already been weighted by colliding mole fractions.
        //         This is taken care of as part of by using bath pressure in
        //         calculation of relaxation time.
        // NOTE 2. massf scaling is applied here to convert J/s/kg-of-species-ip to J/s/kg-of-mixture
        return gs.massf[m_p] * (evStar - ev)/tau;
    }

private:
    GasModel mGmodel;
    int m_p;
    int m_mode;
    RelaxationTime mRT;
}
