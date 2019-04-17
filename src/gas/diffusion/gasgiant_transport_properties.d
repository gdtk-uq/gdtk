/**
 * Authors: Igor Segrovets and Rowan G.
 *
 */

module gas.diffusion.gasgiant_transport_properties;

import std.math;
import std.conv : to;

import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;

enum Species {He=0, H2, H }	


// Keep order consistent with above.


class GasGiantViscosity : Viscosity {
    this () {
        _a_22.length = 3;
	foreach (ref e; _a_22) {
	    e.length = 3;
	    foreach (ref ee; e) {
               ee.length = 7;
	    }
	}
	foreach (i; 0 .. _a_22.length) {
            string key = speciesNames[i] ~ ":" ~ speciesNames[i];
            foreach (j, ref e; _a_22[i][i]) {
                e = a_22[key][j];
            }
        }


    }
    override GasGiantViscosity dup() const {
        return new GasGiantViscosity();
    }

    @nogc
    number eval(in GasState Q) {
        return selfViscosity(Q.T, 0, 0);
    }

private:
    number[][][] _a_22;

    @nogc
    number selfViscosity(number T, int isp, int jsp)
    {
        // FIX ME: Igor, check on this when isp != jsp.
        number sig2Omega_22 = collisionIntegral_22(T, isp, jsp);
        number mu = 2.6693e-6*PI*sqrt(M[isp]*T)/(PI*sig2Omega_22);
        
        return mu;
    }

    @nogc
    number collisionIntegral_22(number T, int isp, int jsp)
    {
        auto x = log(T);
        auto a1 = _a_22[isp][jsp][0];
        auto a2 = _a_22[isp][jsp][1];
        auto a3 = _a_22[isp][jsp][2];
        auto a4 = _a_22[isp][jsp][3];
        auto a5 = _a_22[isp][jsp][4];
        auto a6 = _a_22[isp][jsp][5];
        auto a7 = _a_22[isp][jsp][6];
        auto tmpA = exp((x-a3)/a4);
        auto tmpB = exp((a3-x)/a4);
        auto tmpC = exp((x-a6)/a7);
        auto tmpD = exp((a6-x)/a7);

        auto sig2Omega_22 = (a1+a2*x)*tmpA / (tmpA + tmpB);
        sig2Omega_22 += a5*tmpC / (tmpC + tmpD);
        return sig2Omega_22;
    }

}

/*
class GasGiantThermalConductivity : ThermalConductivity {
    GasGiantThermalConductivity dup() const {}
    
    @nogc
    number eval(in GasState Q) {}
}
*/

version(gasgiant_transport_properties_test) {
    import std.stdio;
    int main()
    {
        GasGiantViscosity ggv = new GasGiantViscosity();
        auto Q = new GasState(1, 1);
        auto f = File("gasgiant-visc-test.dat", "w");
        f.writeln("# T(K)   mu-He  mu-H2  mu-H");
        double T = 300.0;
        double dT = 100.0;
        double Tmax = 10000.0;
        while (T <= Tmax) {
            Q.T = T;
            auto mu0 = ggv.selfViscosity(Q.T, 0, 0);
            auto mu1 = ggv.selfViscosity(Q.T, 1, 1);
            auto mu2 = ggv.selfViscosity(Q.T, 2, 2);
            f.writefln("%f %.8e %.8e %.8e", T, mu0, mu1, mu2);
            T += dT;
        }
        f.close();
        return 0;
    }
}

static double[3] M;
static double[3] m;
static double[7][string] a_22;
static string[] speciesNames;
static this()
{
    speciesNames ~= "He";
    speciesNames ~= "H2";
    speciesNames ~= "H";

    M[Species.He] = 4.002602; 
    m[Species.He] = (M[Species.He]/Avogadro_number)*10e-3;

    a_22["He:He"] = [15.9887272, -4.3159308, 8.0, -15.52249873, 14.35288253, 4.0, 4.80040319];
    
    M[Species.H2] = 2.015882;
    m[Species.H2] = (M[Species.H2]/Avogadro_number)*10e-3;

    a_22["H2:H2"] = [27.54387526,-1.98253166,3.88885724,-12.91940775,0.3470796,8.72131306,-0.88296275];

    M[Species.H] = 2.015882/2;
    m[Species.H] = (M[Species.H2]/Avogadro_number)*10e-3;

    a_22["H:H"] = [22.08948804,-1.85066626,8.50932055,-7.66943974,0.77454531,9.69545318,-0.62104466];
}



