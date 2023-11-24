/**
 * wilke_mixing_therm_cond.d
 * Implements Wilke's mixing rule to compute the
 * thermal conductivity of a mixture of gases.
 * The notation follows that used by White (2006).
 *
 * References:
 * Wilke, C.R. (1950)
 * A Viscosity Equation for Gas Mixtures.
 * Journal of Chemical Physics, 18:pp. 517--519
 *
 * White, F.M. (2006)
 * Viscous Fluid Flow, Third Edition
 * NcGraw Hill, New York
 * (see page 34)
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-09-08 -- initial cut
 * Version: 2023-06-28 -- Optimisations by NNG
 */

module gas.diffusion.wilke_mixing_therm_cond;

import std.math;
import ntypes.complex;
import nm.number;
import util.msg_service;

import gas.gas_model;
import gas.gas_state;
import gas.diffusion.therm_cond;

class WilkeMixingThermCond : ThermalConductivityMixtureModel {
public:
    this(in ThermalConductivity[] tcms, in double[] mol_masses)
    in {
        assert(tcms.length == mol_masses.length,
               brokenPreCondition("tcms.length and mol_masses.length", __LINE__, __FILE__));
    }
    do {
        foreach (tcm; tcms) {
            _tcms ~= tcm.dup;
        }
        nsp = mol_masses.length;
        _mol_masses = mol_masses.dup;
        _x.length = _mol_masses.length;
        _k.length = _mol_masses.length;
        _numer.length = _mol_masses.length;
        _denom.length = _mol_masses.length;
        foreach (ref n; _numer) n.length = _mol_masses.length;
        foreach (ref d; _denom) d.length = _mol_masses.length;

        foreach(i; 0 .. nsp){
            foreach(j; 0 .. nsp){
                _numer[i][j] = pow(_mol_masses[j]/_mol_masses[i], 0.25);
                _denom[i][j] = sqrt(8.0 + 8.0*_mol_masses[i]/_mol_masses[j]);
            }
        }
    }
    this(in WilkeMixingThermCond src) {
        this(src._tcms, src._mol_masses);
    }
    override WilkeMixingThermCond dup() const {
        return new WilkeMixingThermCond(this);
    }

    override number eval(ref const(GasState) Q, int imode) {
        // 1. Evaluate the mole fractions
        massf2molef(Q.massf, _mol_masses, _x);
        // 2. Calculate the component thermoconductivities
        number T = Q.T;
        number logT = log(T);
        foreach(isp; 0 .. nsp) {
            _k[isp] = _tcms[isp].eval(T, logT);
        }
        // 3. Interaction potentials are now precalculated
        // 4. Apply mixing formula
        number sum;
        number k = 0.0;
        foreach (i; 0 .. nsp) {
            if ( _x[i] < SMALL_MOLE_FRACTION ) continue;
            sum = 0.0;
            foreach (j; 0 .. nsp) {
                if ( _x[j] < SMALL_MOLE_FRACTION ) continue;
                number phiij = 1.0 + sqrt(_k[i]/_k[j])*_numer[i][j];
                phiij *= phiij;
                phiij /= _denom[i][j];
                sum += _x[j]*phiij;
            }
            k += _k[i]*_x[i]/sum;
        }
        return k;
    }

private:
    size_t nsp;
    ThermalConductivity[] _tcms; // component viscosity models
    double[] _mol_masses; // component molecular weights
    number[][] _denom, _numer; // precomputed interation potential components
    // Working array space
    number[] _x;
    number[] _k;
}


version(wilke_mixing_therm_cond_test) {
    int main()
    {
        import std.stdio;
        import gas.diffusion.sutherland_therm_cond;
        // Placeholder test. Redo with CEA curves.
        number T = 300.0;
        auto tcm_N2 = new SutherlandThermCond(273.0, 0.0242, 150.0);
        auto tcm_O2 = new SutherlandThermCond(273.0, 0.0244, 240.0);
        auto tcm = new WilkeMixingThermCond([tcm_N2, tcm_O2], [28.0e-3, 32.0e-3]);

        auto gd = GasState(2, 0);
        gd.T = T;
        gd.massf[0] = 0.8;
        gd.massf[1] = 0.2;
        tcm.update_thermal_conductivity(gd);
        assert(isClose(0.0263063, gd.k, 1.0e-3), failedUnitTest());

        return 0;
    }
}

