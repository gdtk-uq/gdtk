/**
 * Models for gas conductivity. Perhaps this should be a gas model
 *
 * Author: Nick Gibbons
 * Version: 2021-05-24
 */

module efieldconductivity;

import std.stdio;
import std.format;
import std.math;
import std.conv;

import geom;
import nm.number;
import mass_diffusion;
import gas.gas_state;
import gas.gas_model;
import gas.physical_constants;

interface ConductivityModel{
    @nogc number opCall(ref const(GasState) gs, const Vector3 pos, GasModel gm);
}


class TestConductivity : ConductivityModel{
    this() {}
    final const number opCall(ref const(GasState) gs, const Vector3 pos, GasModel gm){
        double sigma = -1.0*exp(pos.x.re)*cos(pos.y.re);
        return to!number(sigma);
    }
}

class ConstantConductivity : ConductivityModel{
/*
    Test with a just constant conductivity.
*/
    this() {}
    @nogc final number opCall(ref const(GasState) gs, const Vector3 pos, GasModel gm){
        return to!number(1.0);
    }
}

class RaizerConductivity : ConductivityModel{
/*
    Test with the formula from: Y. P. Razier, Gas Discharge Physics (Springer-Verlag, 1991)
     - Valid for air, nitrogen and argon when weakly ionised.
*/
    this() {}
    @nogc final number opCall(ref const(GasState) gs, const Vector3 pos, GasModel gm){
        version(multi_T_gas) {
            double Tref = gs.T_modes[0].re; // Hmmm. This will crash in single temp
        } else {
            double Tref = gs.T.re;
        }
        number sigma = 8300.0*exp(-36000.0/Tref);
        //debug{writefln(" gs: %s sigma: %e ", gs, sigma);}
        return sigma;
    }
}

class DiffusionConductivity : ConductivityModel{
/*
    Compute the electrical conductivity using an expression derived by NNG
    based on chapter 6 of Seshadri, 1925

     "Fundamentals of Plasma Physics"
     S. R. Seshadri
     American Elsivier Publishing Company, Inc. NU 10017

    @author: Nick Gibbons (09/12/21)
*/
    this(GasModel gm) {
        if (!gm.is_plasma) throw new Error("DiffusionConductivity model requires a GasModel with is_plasma=true");

        nsp = gm.n_species;
        number_density.length = nsp;
        Davg.length = nsp;
        // We lie to BinaryDiffusion that our is_plasma is false to prevent it from enforcing ambipolar diffusion
        bd = new BinaryDiffusion(nsp, false, gm.charge);
    }

    @nogc final number opCall(ref const(GasState) gs, const Vector3 pos, GasModel gm){
    /*
        We assume that the Einstein relations are valid for a multicomponent plasma.
          - Also assume no applied magnetic field and weak species gradients.
    */
        gm.massf2numden(gs, number_density);
        bd.computeAvgDiffCoeffs(gs, gm, Davg);

        number sigma = 0.0;
        foreach(i; 0 .. nsp){
            double Z = gm.charge[i];
            number n = number_density[i];
            number T = gs.T; // TODO: Consider using the electron temperature for i==[e-]
            number D = Davg[i];
            sigma +=  Z*Z*D*n*elementary_charge*elementary_charge/Boltzmann_constant/T;
        }
        //debug{writefln(" gs: %s sigma: %e ", gs, sigma);}
        return sigma;
    }
private:
    size_t nsp;
    number[] number_density;
    number[] Davg;
    BinaryDiffusion bd;
}

ConductivityModel create_conductivity_model(string name, GasModel gm){
    ConductivityModel conductivity_model;
    switch (name) {
    case "test":
        conductivity_model = new TestConductivity();
        break;
    case "constant":
        conductivity_model = new ConstantConductivity();
        break;
    case "raizer":
        conductivity_model = new RaizerConductivity();
        break;
    case "diffusion":
        conductivity_model = new DiffusionConductivity(gm);
        break;
    case "none":
        throw new Error("User has asked for solve_electric_field but failed to specify a conductivity model.");
    default:
        string errMsg = format("The conductivity model '%s' is not available.", name);
        throw new Error(errMsg);
    }
    return conductivity_model;
}
