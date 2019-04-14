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

enum Species {He=0}

class GasGiantViscosity : Viscosity {
    this () {
        _a_22.length = 1;
        _a_22[0].length = 1;
        _a_22[0][0].length = 7;
        foreach (i, ref e; _a_22[0][0]) {
            e = a_22["He:He"][i];
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
        auto ggv = new GasGiantViscosity();
        auto Q = new GasState(1, 1);
        auto f = File("He-visc-test.dat", "w");
        f.writeln("# T(K)   mu");
        double T = 300.0;
        double dT = 100.0;
        double Tmax = 10000.0;
        while (T <= Tmax) {
            Q.T = T;
            auto mu = ggv.eval(Q);
            f.writefln("%f %.8e", T, mu);
            T += dT;
        }
        f.close();
        return 0;
    }
}

static double[1] M;
static double[1] m;
static double[7][string] a_22;

static this()
{
    M[Species.He] = 4.002602; 
    m[Species.He] = M[Species.He]/Avogadro_number;

    a_22["He:He"] = [15.9887272, -4.3159308, 8.0, -15.52249873, 14.35288253, 4.0, 4.80040319];
    
}



