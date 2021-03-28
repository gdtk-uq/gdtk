/**
 * Interface and implementation for relaxation time calculations for energy exchange mechanisms.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.relaxation_time;

import std.math;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;


interface RelaxationTime {
    @nogc numer eval(in GasState gs, number[] molef);
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

    this(lua_State *L)
    {
        m_q = getInt(L, -1, "q");
        m_a = getDouble(L, -1, "a");
        m_b = getDouble(L, -1, "b");
        m_mu = getDouble(L, -1, "mu");
    }

    @nogc
    number eval(GasState gs, number[] molef)
    {
        // If bath pressure is very small, then practically no relaxation
        // occurs from collisions with particle q.
        if (molef[m_q] <= SMALL_MOLE_FRACTION)
            return -1.0;
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
