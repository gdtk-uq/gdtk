/**
 * Interface and implementation for energy exchange mechanisms.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.energy_exchange_mechanism;

import std.format;
import std.math;
import std.stdio;
import std.conv;

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

    this(lua_State *L, int p, int q, int mode, GasModel gmodel)
    {
	m_p = p;
        m_q = q;
	m_mode = mode;
	mGmodel = gmodel;
	lua_getfield(L, -1, "relaxation_time");
	mRT = createRelaxationTime(L, q);
	lua_pop(L, 1);
    }
    
    this(int p, int q, int mode, RelaxationTime RT, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_mode = mode;
        mRT = RT.dup();
        mGmodel = gmodel;
    }
    
    @nogc
    override number rate(in GasState gs, in GasState gsEq, number[] molef)
    {
        // If bath pressure is very small, then practically no relaxation
        // occurs from collisions with particle q.
        if (m_tau < 0.0) // signals very small mole fraction of q
            return to!number(0.0);
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
    int m_q;
    int m_mode;
}

EnergyExchangeMechanism createVTMechanism(lua_State *L, int p, int q, int mode, GasModel gmodel)
{
    auto rateModel = getString(L, -1, "rate");
    switch (rateModel) {
    case "Landau-Teller":
	return new LandauTellerVT(L, p, q, mode, gmodel);
    default:
	string msg = format("The VT mechanism rate model: %s is not known.", rateModel);
	throw new Error(msg);
    }
    
    

    
}
