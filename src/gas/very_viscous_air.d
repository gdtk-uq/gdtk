/**
 * very_viscous_air.d
 * A special ideal air model used for the Method of Manufactured Solutions
 * viscous case.
 *
 * Author: Rowan G.
 * Version: 2015-05-05
 */

module gas.very_viscous_air;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import ntypes.complex;
import nm.number;
import util.msg_service;
import util.lua;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;

class VeryViscousAir: GasModel {
public:
    this() {
        type_str = "VeryViscousAir";
        // Default model is mostly initialized in the private data below.
        _n_species = 1;
        _n_modes = 0;
        _species_names ~= "very viscous";
        _Rgas = 287.0;
        _mol_masses ~= R_universal/_Rgas;
        _gamma = 1.4;
        _Cv = _Rgas / (_gamma - 1.0);
        _Cp = _Rgas*_gamma/(_gamma - 1.0);
        _mu = 10.0;
        double Pr = 1.0;
        _k = _mu * _Cp / Pr;
        create_species_reverse_lookup();
    }

    this(lua_State* L)
    {
        this();
        lua_getglobal(L, "VeryViscousAir");
        // Possibly override k,  mu and number of species
        lua_getfield(L, -1, "mu");
        if ( !lua_isnil(L, -1) ) {
            _mu = to!double(lua_tonumber(L, -1));
        }
        lua_pop(L, 1);
        lua_getfield(L, -1, "k");
        if ( !lua_isnil(L, -1) ) {
            _k = to!double(lua_tonumber(L, -1));
        }
        lua_pop(L, 1);
        lua_getfield(L, -1, "number_species");
        if ( !lua_isnil(L, -1) ) {
            _n_species = luaL_checkint(L, -1);
        }
        lua_pop(L, 1);

        if (_n_species > 1) {
            // Let's rename the species.
            _species_names.length = _n_species;
            foreach (isp; 0 .. _n_species) {
                _species_names[isp] = format("air%d", isp);
            }
            // Let's fix the _mol_masses array.
            _mol_masses.length = _n_species;
            foreach (isp; 0 .. _n_species) _mol_masses[isp] = R_universal/_Rgas;
        }
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "VeryViscousAir =()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q) const
    {
        Q.rho = Q.p/(Q.T*_Rgas);
        Q.u = _Cv*Q.T;
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        Q.T = Q.u/_Cv;
        Q.p = Q.rho*_Rgas*Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        Q.p = Q.rho*_Rgas*Q.T;
        Q.u = _Cv*Q.T;
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        Q.T = Q.p/(Q.rho*_Rgas);
        Q.u = _Cv*Q.T;
    }

    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        throw new Exception("Not implemented.");
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        throw new Exception("Not implemented.");
    }
    override void update_sound_speed(ref GasState Q) const
    {
        Q.a = sqrt(_gamma*_Rgas*Q.T);
    }
    override void update_trans_coeffs(ref GasState Q) const
    {
        Q.mu = _mu;
        Q.k = _k;
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
    double _Rgas = 287.0; // J/kg/K
    double _gamma = 1.4;   // ratio of specific heats
    double _Cv = R_universal/0.02896 / 0.4; // J/kg/K
    double _Cp = R_universal/0.02896 * 1.4/0.4; // J/kg/K
    // Reference values for entropy
    double _s1 = 0.0; // J/kg/K
    double _T1 = 298.15; // K
    double _p1 = 101.325e3; // Pa
    // Molecular transport coefficent constants.
    double _mu;
    double _k;

} // end class Very_viscous_air

version(very_viscous_air_test) {
    int main() {
        auto gm = new VeryViscousAir();
        auto gs = GasState(gm, 100.0e3, 300.0);
        assert(isClose(gm.R(gs), 287.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gs);
        gm.update_sound_speed(gs);
        assert(isClose(gs.rho, 1.16144, 1.0e-6), failedUnitTest());
        assert(isClose(gs.u, 215250, 1.0e-6), failedUnitTest());
        assert(isClose(gs.a, 347.189, 1.0e-6), failedUnitTest());
        gm.update_trans_coeffs(gs);
        assert(isClose(gs.mu, 10.0, 1.0e-6), failedUnitTest());
        assert(isClose(gs.k, 10045, 1.0e-6), failedUnitTest());

        return 0;
    }
}
