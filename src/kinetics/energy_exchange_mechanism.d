/**
 * Interface and implementation for energy exchange mechanisms.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.energy_exchange_mechanism;

import std.math;
import std.stdio;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;

import kinetics.relaxation_time;


class EnergyExchangeMechanism {
    @property @nogc number tau() const { return m_tau; }
    @nogc void evalRelaxationTime(in GasState gs, number[] molef)
    {
        m_tau = mRT.eval(gs, molef);
    }
    @nogc abstract number rate(in GasState gs, in GasState gsEq, number[] molef);

private:
    number m_tau;
    RelaxationTime mRT;
    GasModel mGmodel;
}

class LandauTellerVT : EnergyExchangeMechanism {
public:

    this(int p, int mode, RelaxationTime RT, GasModel gmodel)
    {
        m_p = p;
        m_mode = mode;
        mRT = RT.dup();
        mGmodel = gmodel;
    }
    
    @nogc
    override number rate(in GasState gs, in GasState gsEq, number[] molef)
    {
        number evStar = mGmodel.energyPerSpeciesInMode(gsEq, m_p, m_mode);
        number ev = mGmodel.energyPerSpeciesInMode(gs, m_p, m_mode);
        // NOTE 1. tau has already been weighted by colliding mole fractions.
        //         This is taken care of as part of by using bath pressure in
        //         calculation of relaxation time.
        // NOTE 2. massf scaling is applied here to convert J/s/kg-of-species-ip to J/s/kg-of-mixture
        return gs.massf[m_p] * (evStar - ev)/m_tau;
    }

private:
    int m_p;
    int m_mode;
}
