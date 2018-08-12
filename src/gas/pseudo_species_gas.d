/**
 * pseudo_species_gas.d
 * Authors: Pierre Mariotto, Rowan G. and Peter J. 
 *
 * This model describes a gas in which the populations of microstates (in terms
 * of their internal energy state) are tracked. We might group together a few
 * microstates that are closely distributed in terms of energy. We will end
 * up with a number of groups. Some with a few microstates per group, some might
 * have only a single microstate in the groups. These groups we call pseudo-species.
 *
 */

module gas.pseudo_species_gas;

import std.algorithm.iteration : sum;
import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.pseudo_species;

class PseudoSpeciesGas : GasModel {
public:
    this(lua_State *L)
    {
        auto nPseudoSpecies = getInt(L, LUA_GLOBALSINDEX, "number_pseudo_species");
        _n_species = cast(uint) nPseudoSpecies;
        _n_modes = 0;
        _species_names.length = _n_species;

        lua_getfield(L, LUA_GLOBALSINDEX, "pseudo_species");
        foreach (isp; 0 .. _n_species) {
            lua_rawgeti(L, -1, isp);
            if (lua_isnil(L, -1)) {
                string msg = format("There was an error when attempting to information about pseudo-species %d.\n", isp);
                throw new Error(msg);
            }
            _pseudoSpecies ~= createPseudoSpecies(L);
            lua_pop(L, 1);
            _species_names[isp] = _pseudoSpecies[$-1].name;
            _mol_masses ~= _pseudoSpecies[$-1].mol_mass;
            _R ~= R_universal/_pseudoSpecies[$-1].mol_mass;
            _dof ~= _pseudoSpecies[$-1].DOF;
        }
        lua_pop(L, 1);
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "PseudoSpeciesGas=()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.rho = Q.p/(R_mix*Q.T);

        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T + uNoneq;
    }

    override void update_thermo_from_rhou(GasState Q)
    {
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.T = (Q.u - uNoneq)/Cv_mix; 
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
    }

    override void update_thermo_from_rhoT(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T + uNoneq;
    }

    override void update_thermo_from_rhop(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.T = Q.p/(Q.rho*R_mix);
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T + uNoneq;
    }

    override void update_thermo_from_ps(GasState Q, number s)
    {
        throw new Error("ERROR: PseudoSpeciesGas:update_thermo_from_ps NOT IMPLEMENTED.");
    }

    override void update_thermo_from_hs(GasState Q, number h, number s)
    {
        throw new Error("ERROR: PseudoSpeciesGas:update_thermo_from_hs NOT IMPLEMENTED.");
    }

    override void update_sound_speed(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        auto Cv_mix = Cv(Q);
        auto gamma_mix = (R_mix/Cv_mix) + 1;
        Q.a = sqrt(gamma_mix*R_mix*Q.T);
    }

    override void update_trans_coeffs(GasState Q)
    {
        // THINK ABOUT: PM/RJG
        // We will need to implement this when we've made a decision.
        // For the present, only inviscid simulations are valid.
        Q.mu = 0.0;
        Q.k = 0.0;
    }

    override number dudT_const_v(in GasState Q) const
    {
        number Cv = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv += Q.massf[isp] * (_dof[isp]/2.0) * _R[isp];
        }
        return Cv;
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return gas_constant(Q) + dudT_const_v(Q);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        throw new Error("ERROR: PseudoSpeciesGas:dpdrho_const_T NOT IMPLEMENTED.");
    }
    override number gas_constant(in GasState Q) const
    {
        number R_mix = 0.0;
        foreach (isp; 0 .. _n_species) {
            R_mix += Q.massf[isp]*_R[isp];
        }
        return R_mix;
    }
    override number internal_energy(in GasState Q) const
    {
        return Q.u + sum(Q.u_modes);
    }
    override number enthalpy(in GasState Q) const
    {
        return Q.u + sum(Q.u_modes) + Q.p/Q.rho;
    }
    override number entropy(in GasState Q) const
    {
        throw new Error("ERROR: PseudoSpeciesGas:entropy NOT IMPLEMENTED.");
    }

private:
    PseudoSpecies[] _pseudoSpecies;
    double[] _R;
    int[] _dof;

    number energyInNoneq(GasState Q) const {
        number uNoneq = 0.0;
        foreach (isp; 0 .. _n_species) {
            uNoneq += (Avogadro_number/_mol_masses[isp]) * Q.massf[isp] * _pseudoSpecies[isp].energy(Q);
        }
        return uNoneq;
    }
}

version(pseudo_species_gas_test) {
    int main()
    {
        import util.msg_service;
        import std.math : FloatingPointControl; 
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        doLuaFile(L, "sample-data/pseudo-species-3-components.lua");
        auto gm = new PseudoSpeciesGas(L);
        auto gd = new GasState(3, 3);

        return 0;
    }
    
}
