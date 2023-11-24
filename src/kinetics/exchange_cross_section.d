/*
Definitions for electron-heavy particle energy exchange models.

References:
 - "Conservation Equations and Physical Models for Hypersonic Air Flows in Thermal and Chemical Nonequilibrium"
   Peter A. Gnoffo and Roop N. Gupta and Judy L. Shinn, NASA Technical Paper 2867, 1989

 - "Effects of hydrogen impurities on shock structure and stability in ionizing monatomic gases. Part 1. Argon"
   I. I. Glass and W. S. Lio, Journal of Fluid Mechanics vol. 85, part 1, pp 55-77, 1978
 
@author: Nick Gibbons (21/04/20)
*/

module kinetics.exchange_cross_section;

import std.string;
import std.math;
import std.conv;
import std.stdio;

import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;


interface ExchangeCrossSection {
    ExchangeCrossSection dup();
    @nogc number opCall(in GasState gs, number[] numden);
}

class GnoffoNeutralECS : ExchangeCrossSection {
    /*
        For collision cross sections specified with a simple three-parameter fit,
        as in Gnoffo, 1989, equation 66

        Notes: A limit of 30,000 K is adopted for the curve-fit, as that is the 
        value mentioned in the report that the numbers are valid to. The exchange
        rate will still scale above this temperature however, because of the 
        Te dependance in the rate equation.
    */
    this(int e, int q, double a, double b, double c, int mode) {
        this.e = e;
        this.q = q;
        this.a = a;
        this.b = b;
        this.c = c;
        this.mode = mode;
    }

    this(lua_State *L, int e, int q, int mode) {
        this.e = e;
        this.q = q;
        this.a = getDouble(L, -1, "a");
        this.b = getDouble(L, -1, "b");
        this.c = getDouble(L, -1, "c");
        this.mode = mode;
    }

    GnoffoNeutralECS dup() {
        return new GnoffoNeutralECS(e, q, a, b, c, mode);
    }
    
    @nogc
    override number opCall(in GasState gs, number[] numden) {
        number Te = gs.T_modes[mode];
        number Tfit = fmin(30e3, Te);
        return a + b*Tfit + c*Tfit*Tfit;
    }

private:
    int e,q,mode;
    double a,b,c;
}

class CoulombIonECS : ExchangeCrossSection {
    /*
        For an electron colliding with a positively charged ion, an
        analytical expression exists for the cross section.

        See Glass and Liu 1978, Equation (16)

        Notes: For some reason the expression for this quantity in Glass and
        Liu has slightly different numerical constants in it compared to the
        one in Gnoffo.  This results in relaxation rates that are ~10%
        different. This is probably just due different approximations in the
        derivation of the formulas, but it would be nice to know which one is
        less approximate... (See notes from 21/04/28)
    */
    this(int e, int q, int mode) {
        this.e = e;
        this.q = q;
        this.mode = mode;
    }

    this(lua_State *L, int e, int q, int mode) {
        this.e = e;
        this.q = q;
        this.mode = mode;
    }

    CoulombIonECS dup() {
        return new CoulombIonECS(e, q, mode);
    }
    
    @nogc
    override number opCall(in GasState gs, number[] numden) {
        number ne = fmax(1.0, numden[e]);
        number ir0 = 4.0*pi*Boltzmann_constant*gs.T_modes[mode]*vacuum_permittivity/elementary_charge/elementary_charge;
        number Lambda = log(9.0/4.0*ir0*ir0*ir0/pi/ne);
        return 2.0*pi/9.0/ir0/ir0*Lambda;
    }

private:
    int e,q,mode;
    immutable double pi = to!double(PI);
}

ExchangeCrossSection createExchangeCrossSection(lua_State *L, int e, int q, int mode)
{
    auto type = getString(L, -1, "type");
    switch (type) {
    case "GnoffoNeutral":
        return new GnoffoNeutralECS(L, e, q, mode);
    case "Coulomb":
        return new CoulombIonECS(L, e, q, mode);
    default:
        string msg = format("The exchange cross section type: %s is not known.", type);
        throw new Error(msg);
    }

}

