/**
 * ideal_gas_ab.d
 *
 * Two-component reacting gas as described in.
 * JM Powers and TD Aslam (2006)
 * Exact solution for multidimensional compressible reactive flow
 * for verifying numerical algorithms.
 * AIAA Journal Vol. 44 No. 2 pages 337-344
 *
 * This gas model is useful as a demonstration of building a custom
 * reacting gas model and as a flow-solver verification tool.
 *
 * This binary reacting gas also appears in other works:
 *
 * Yee, H.C., Kotov, D.V., Wang, W. and Shu, C.-W. (2013)
 * Spurious behavior of shock-capturing methods by the fractional
 * step approach: Problems containing stiff source terms and discontinuities.
 * Journal of Computational Physics, 241: pp. 266--291
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2017-01-07: initial cut.
 *          2020-05-08: renamed IdealGasAB
 */

module gas.ideal_gas_ab;

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

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;

// The basic gas model.

class IdealGasAB: GasModel {
public:
    string modelType = "Powers-Aslam";
    this(lua_State *L) {
        type_str = "IdealGasAB";
        // Some parameters are fixed and some come from the gas model file.
        _n_species = 2;
        _n_modes = 0;
        _species_names.length = 2;
        _species_names[0] = "A";
        _species_names[1] = "B";
        // Bring table to TOS
        lua_getglobal(L, "IdealGasAB");
        // [TODO] test that we actually have the table as item -1
        lua_getfield(L, -1, "type");
        if (!lua_isnil(L, -1)) {
            modelType = to!string(luaL_checkstring(L, -1));
        }
        lua_pop(L, 1);
        // Now, pull out the remaining numeric value parameters.
        _Rgas = getDouble(L, -1, "R");
        _mol_masses.length = 2;
        _mol_masses[0] = R_universal / _Rgas;
        _mol_masses[1] = _mol_masses[0];
        _gamma = getDouble(L, -1, "gamma");
        // Heat of reaction
        _q = getDouble(L, -1, "q");
        lua_pop(L, 1); // dispose of the table
        // Entropy reference, same as for IdealAir
        _s1 = 0.0;
        _T1 = 298.15;
        _p1 = 101.325e3;
        // Compute derived parameters
        _Cv = _Rgas / (_gamma - 1.0);
        _Cp = _Rgas*_gamma/(_gamma - 1.0);
        create_species_reverse_lookup();
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= "IdealGasAB =(";
        repr ~= "species=[\"A\", \"B\"]";
        repr ~= ", Mmass=[" ~ to!string(_mol_masses[0]);
        repr ~= "," ~ to!string(_mol_masses[1]) ~ "]";
        repr ~= ", gamma=" ~ to!string(_gamma);
        repr ~= ", q=" ~ to!string(_q);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q) const
    {
        Q.rho = Q.p/(Q.T*_Rgas);
        Q.u = _Cv*Q.T - Q.massf[1]*_q;
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        Q.T = (Q.u + Q.massf[1]*_q)/_Cv;
        Q.p = Q.rho*_Rgas*Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        Q.p = Q.rho*_Rgas*Q.T;
        Q.u = _Cv*Q.T - Q.massf[1]*_q;
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        Q.T = Q.p/(Q.rho*_Rgas);
        Q.u = _Cv*Q.T - Q.massf[1]*_q;
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
        Q.a = sqrt(_gamma*_Rgas*Q.T);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        // The gas is inviscid for the test cases described in the AIAA paper.
        Q.mu = 0.0;
        Q.k = 0.0;
    }
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
    override number internal_energy(in GasState Q) const
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
    double _Rgas; // J/kg/K
    double _gamma;   // ratio of specific heats
    double _Cv; // J/kg/K
    double _Cp; // J/kg/K
    // Reference values for entropy
    double _s1;  // J/kg/K
    double _T1;  // K
    double _p1;  // Pa
    // Molecular transport coefficents are zero.
    // Heat of reaction.
    double _q; // J/kg
} // end class IdealGasAB

// Unit test of the basic gas model...

version(ideal_gas_ab_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/ideal-gas-ab-model.lua");
        auto gm = new IdealGasAB(L);
        lua_close(L);
        auto gd = GasState(2, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 0.75; gd.massf[1] = 0.25;
        assert(isClose(gm.R(gd), 287.0, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 2, failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[0], 0.75, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[1], 0.25, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        number my_rho = 1.0e5 / (287.0 * 300.0);
        assert(isClose(gd.rho, my_rho, 1.0e-4), failedUnitTest());
        number my_Cv = gm.dudT_const_v(gd);
        number my_u = my_Cv*300.0 - 0.25*300000.0;
        assert(isClose(gd.u, my_u, 1.0e-3), failedUnitTest());
        number my_Cp = gm.dhdT_const_p(gd);
        number my_a = sqrt(my_Cp/my_Cv*287.0*300.0);
        assert(isClose(gd.a, my_a, 1.0e-3), failedUnitTest());
        gm.update_trans_coeffs(gd);
        assert(isClose(gd.mu, 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.k, 0.0, 1.0e-6), failedUnitTest());

        return 0;
    }
}
