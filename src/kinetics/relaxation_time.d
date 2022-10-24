/**
 * Interface and implementation for relaxation time calculations for energy exchange mechanisms.
 *
 *  - "Conservation Equations and Physical Models for Hypersonic Air Flows in Thermal and Chemical Nonequilibrium"
 *     Peter A. Gnoffo and Roop N. Gupta and Judy L. Shinn, NASA Technical Paper 2867, 1989
 *
 *  - "Review of Chemical-Kinetic Problems for Future NASA Missions, I: Earth Entries"
 *     Chul Park, Journal of Thermophysics and Heat Transfer, Volume 7, Number 3, July-Sept 1993
 *
 *  - "Modification of chemical-kinetic parameters for 11-air species in re-entry flows"
 *     Jae Gang Kim and Sung Min Jo, International Journal of Heat and Mass Transfer, Volume 169, 2021
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.relaxation_time;

import std.string;
import std.math;
import std.conv;
import std.stdio;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;


interface RelaxationTime {
    RelaxationTime dup();
    @nogc number eval(in GasState gs, number[] molef, number[] numden);
}

class MillikanWhiteVT : RelaxationTime {
public:
    this(int q, double a, double b)
    {
        m_q = q;
        m_a = a;
        m_b = b;
    }

    this(lua_State *L, int q)
    {
        m_q = q;
        m_a = getDouble(L, -1, "a");
        m_b = getDouble(L, -1, "b");
    }

    MillikanWhiteVT dup()
    {
        return new MillikanWhiteVT(m_q, m_a, m_b);
    }

    @nogc
    override number eval(in GasState gs, number[] molef, number[] numden)
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
}

class ParkHTCVT : RelaxationTime {
/*
    High temperature correction taken from Gnoffo, 1989 (equations 56-58). Originally
    from Park, 1985.

    @author: Nick Gibbons
*/
    this(int p, int q, double sigma, double mu, RelaxationTime rt, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_sigma = sigma;
        m_mu = mu;
        m_rt = rt.dup();
        m_gmodel = gmodel;
    }

    this(lua_State *L, int p, int q, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_sigma = getDouble(L, -1, "sigma");

        double mq = gmodel.mol_masses[m_q]/Avogadro_number;
        double mp = gmodel.mol_masses[m_p]/Avogadro_number;
        m_mu = mq*mp/(mq+mp); // Reduced mass of the collision pair

        lua_getfield(L, -1, "submodel");
        m_rt = createRelaxationTime(L, p, q, gmodel);
        lua_pop(L, 1);
        m_gmodel = gmodel;
    }

    ParkHTCVT dup()
    {
        return new ParkHTCVT(m_p, m_q, m_sigma, m_mu, m_rt, m_gmodel);
    }

    @nogc
    override number eval(in GasState gs, number[] molef, number[] numden)
    {
        /*
        TODO: Park's original derivation is for the relaxation of a single
        species gas and it's not entirely clear what to do with the expression
        in the multi-species case.

        I've chosen to use the hard shell collision frequency expression, with
        the number density of the colliding particle as the relevant parameter.
        This is because we want the rate *per unit mass* of the relaxing
        species, although we should probably look into a more rigorous
        derivation of this so I can be sure that this is right.
        */
        if (molef[m_q] <= SMALL_MOLE_FRACTION)
            return to!number(-1.0);

        number tau_submodel = m_rt.eval(gs, molef, numden);

        number sigma = compute_cross_section(gs);
        number cs = sqrt(8.0*Boltzmann_constant*gs.T/pi/m_mu); // Mean particle velocity
        number tau_park = 1.0/(sigma*cs*numden[m_q]);        // Notice q is the colliding particle
        return tau_submodel + tau_park;
    }

protected:
    immutable double pi = to!double(PI);
    int m_p;
    int m_q;
    double m_sigma;
    double m_mu;
    RelaxationTime m_rt;
    GasModel m_gmodel;

    @nogc
    number compute_cross_section(ref const(GasState) gs) {
        return to!number(m_sigma);
    }
}

class ParkHTC2VT : ParkHTCVT {
/*
    High temperature correction taken Park 1993 (page 387).
    @author: Nick Gibbons
*/
    this(int p, int q, double sigma, double mu, RelaxationTime rt, GasModel gmodel)
    {
        super(p, q, sigma, mu, rt, gmodel);
    }

    this(lua_State *L, int p, int q, GasModel gmodel)
    {
        super(L, p, q, gmodel);
    }

    override ParkHTC2VT dup()
    {
        return new ParkHTC2VT(m_p, m_q, m_sigma, m_mu, m_rt, m_gmodel);
    }

protected:
    @nogc
    override number compute_cross_section(ref const(GasState) gs) {
        return m_sigma*(50e3/gs.T)*(50e3/gs.T);
    }
}

class KimHTCVT : ParkHTCVT {
/*
    High temperature correction taken from Kim 2021 (Table I).
    @author: Nick Gibbons
*/
    this(int p, int q, double sigma, double exponent, double mu, RelaxationTime rt, GasModel gmodel)
    {
        m_exponent = exponent;
        super(p, q, sigma, mu, rt, gmodel);
    }

    this(lua_State *L, int p, int q, GasModel gmodel)
    {
        m_exponent = getDouble(L, -1, "exponent");
        super(L, p, q, gmodel);
    }

    override KimHTCVT dup()
    {
        return new KimHTCVT(m_p, m_q, m_sigma, m_exponent, m_mu, m_rt, m_gmodel);
    }

protected:
    double m_exponent;

    @nogc
    override number compute_cross_section(ref const(GasState) gs) {
        return m_sigma*pow(gs.T, m_exponent);
    }
}


RelaxationTime createRelaxationTime(lua_State *L, int p, int q, GasModel gmodel)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "Millikan-White":
	return new MillikanWhiteVT(L, q);
    case "ParkHTC":
	return new ParkHTCVT(L, p, q, gmodel);
    case "ParkHTC2":
	return new ParkHTC2VT(L, p, q, gmodel);
    case "KimHTC":
	return new KimHTCVT(L, p, q, gmodel);
    default:
	string msg = format("The relaxation time model: %s is not known.", model);
	throw new Error(msg);
    }
}

