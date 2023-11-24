/**
 * Ideal helium model.
 *
 * Although we could treat helium as part of the general ideal gas model,
 * this model allows us to implement the viscosity model to match that used
 * in Peter's C&F paper on expansion tube simulation.
 *
 * Jacobs, P.A. (1994)
 * Numerical Simulation of Transient Hypervelocity Flow in an Expansion Tube
 * Computers and Fluids, 23:1, pp. 77--101
 *
 * Author: Rowan G. and Peter J.
 */

module gas.ideal_helium;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

class IdealHelium: GasModel {
public:

    this()
    {
        type_str = "IdealHelium";
        _n_species = 1;
        _n_modes = 0;
        _species_names.length = 1;
        _species_names[0] = "helium";
        create_species_reverse_lookup();

        _mol_masses.length = 1;
        _mol_masses[0] = 4.002602e-3;
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "IdealHelium =()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q) const
    {
        if (Q.T <= 0.0 || Q.p <= 0.0) {
            string msg = "Temperature and/or pressure was negative for update_thermo_from_pT.";
            throw new GasModelException(msg);
        }
        Q.rho = Q.p/(Q.T*_Rgas);
        Q.u = _Cv*Q.T;
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        if (Q.u <= 0.0 || Q.rho <= 0.0) {
            string msg = "Internal energy and/or density was negative for update_thermo_from_rhou.";
            throw new GasModelException(msg);
        }
        Q.T = Q.u/_Cv;
        Q.p = Q.rho*_Rgas*Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        if (Q.T <= 0.0 || Q.rho <= 0.0) {
            string msg = "Temperature and/or density was negative for update_thermo_from_rhoT.";
            throw new GasModelException(msg);
        }
        Q.p = Q.rho*_Rgas*Q.T;
        Q.u = _Cv*Q.T;
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        if (Q.p <= 0.0 || Q.rho <= 0.0) {
            string msg = "Pressure and/or density was negative for update_thermo_from_rhop.";
            throw new GasModelException(msg);
        }
        Q.T = Q.p/(Q.rho*_Rgas);
        Q.u = _Cv*Q.T;
    }

    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        Q.T = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
        update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        Q.T = h / _Cp;
        Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.T/_T1)));
        update_thermo_from_pT(Q);
    }
    override void update_sound_speed(ref GasState Q) const
    {
        if (Q.T <= 0.0) {
            string msg = "Temperature was negative for update_sound_speed.";
            throw new GasModelException(msg);
        }
        Q.a = sqrt(_gamma*_Rgas*Q.T);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        Q.mu = 5.023e-7 * pow(Q.T, 0.647);
        Q.k = _Cp*Q.mu/_Prandtl;
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
        throw new Exception("not implemented");
    }
    */
    override number dudT_const_v(in GasState Q) const
    {
        return to!number(_Cv);
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return to!number(_Cp);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        number R = gas_constant(Q);
        return R*Q.T;
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(_Rgas);
    }
    override @nogc number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q) const
    {
        return Q.u + Q.p/Q.rho;
    }
    override number entropy(in GasState Q) const
    {
        return _s1 + _Cp * log(Q.T/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    // Thermodynamic constants
    double _Rgas = 2077; // J/kg/K
    double _gamma = 5./3.;   // ratio of specific heats
    double _Cv = 3144; // J/kg/K
    double _Cp = 5191; // J/kg/K
    // Reference values for entropy
    double _s1 = 31517;  // J/kg/K , from NIST
    double _T1 = 298.15;  // K
    double _p1 = 101.325e3;  // Pa
    // Transport-related properties
    double _Prandtl = 0.67;

} // end class Ideal_helium

version(ideal_helium_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        auto gm = new IdealHelium();
        auto gd = GasState(1, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 1.0;
        assert(approxEqualNumbers(gm.R(gd), to!number(2077), 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(approxEqualNumbers(gd.p, to!number(1.0e5), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[0], to!number(1.0), 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        assert(approxEqualNumbers(gd.rho, to!number(0.160488), 1.0e-4), failedUnitTest());
        assert(approxEqualNumbers(gd.u, to!number(943200), 1.0e-4), failedUnitTest());
        assert(approxEqualNumbers(gd.a, to!number(1019.06820), 1.0e-4), failedUnitTest());
        gm.update_trans_coeffs(gd);
        assert(approxEqualNumbers(gd.mu, to!number(2.012151e-5), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.k, to!number(0.155897), 1.0e-6), failedUnitTest());

        return 0;
    }
}
