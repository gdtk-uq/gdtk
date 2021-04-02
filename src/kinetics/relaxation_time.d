/**
 * Interface and implementation for relaxation time calculations for energy exchange mechanisms.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.relaxation_time;

import std.string;
import std.math;
import std.conv;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;


interface RelaxationTime {
    RelaxationTime dup();
    @nogc number eval(in GasState gs, number[] molef);
}

class MillikanWhiteVT : RelaxationTime {
public:
    this(int q, double a, double b, double mu)
    {
        m_q = q;
        m_a = a;
        m_b = b;
        m_mu = mu;
    }

    this(lua_State *L, int q)
    {
        m_q = q;
        m_a = getDouble(L, -1, "a");
        m_b = getDouble(L, -1, "b");
        m_mu = getDouble(L, -1, "mu");
    }

    MillikanWhiteVT dup()
    {
        return new MillikanWhiteVT(m_q, m_a, m_b, m_mu);
    }
    
    @nogc
    override number eval(in GasState gs, number[] molef)
    {
        // If bath pressure is very small, then practically no relaxation
        // occurs from collisions with particle q.
        if (molef[m_q] <= SMALL_MOLE_FRACTION)
            return to!number(-1.0);
        // Set the bath pressure at that of the partial pressure of the 'q' colliders
        // and compute in atm for use in Millikan-White expression
        number pBath = molef[m_q]*gs.p/P_atm;
        number tau = (1.0/pBath) * exp(m_a * (pow(gs.T, -1./3.) - m_b) - 18.42);
        return tau;
    }

private:
    int m_q;
    double m_a;
    double m_b;
    double m_mu;
}

RelaxationTime createRelaxationTime(lua_State *L, int q)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "Millikan-White":
	return new MillikanWhiteVT(L, q);
    default:
	string msg = format("The relaxation time model: %s is not known.", model);
	throw new Error(msg);
    }

}

