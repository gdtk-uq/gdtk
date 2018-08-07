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
        _n_modes = _n_species;
        _species_names.length = _n_species;
        _energy_mode_names.length = _n_species;

        lua_getfield(L, LUA_GLOBALSINDEX, "pseudo_species");
        foreach (isp; 0 .. _n_species) {
            lua_rawgeti(L, -1, isp);
            if (lua_isnil(L, -1)) {
                string msg = format("There was an error when attempting to information about pseudo-species %d.\n", isp);
                throw new Error(msg);
            }
            _pseudoSpecies ~= createPseudoSpecies(L);
            lua_pop(L, 1);
            _species_names[isp] = _pseudoSpecies[$-1].name();
        }
        lua_pop(L, 1);
        _energy_mode_names[] = _species_names[];
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "PseudoSpeciesGas=()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const
    {
        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        Q.rho = Q.p/(_R_N2*Q.T);

        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        foreach (imode; 0 .. _n_modes) {
            Q.u_modes[imode] = (Avogadro_number/_M_N2) * Q.massf[imode] * _pseudoSpecies[imode].energy(Q);
        }
        
        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        Q.u = 2.5*_R_N2*Q.T;
    }

    override void update_thermo_from_rhou(GasState Q) const
    {
        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        Q.T = Q.u/(2.5*_R_N2);
        Q.p = Q.rho*_R_N2*Q.T;

        // THINK ABOUT: PM/RJG
        // Presently, I have not filled in T_modes. It should not be needed anywhere else in the
        // simulation and it seems expensive to do this calculation if it isn't needed.
        // So we just assign Q.T to all modes.
        Q.T_modes[] = Q.T;
    }

    override void update_thermo_from_rhoT(GasState Q) const
    {
        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        Q.p = Q.rho*_R_N2*Q.T;
        Q.u = 2.5*_R_N2*Q.T;
        
        // THINK ABOUT: PM/RJG
        // Here we assume we know something about the temperatures in each mode
        foreach (imode; 0 .. _n_modes) {
            Q.u_modes[imode] = (Avogadro_number/_M_N2) * Q.massf[imode] * _pseudoSpecies[imode].energy(Q);
        }
    }

    override void update_thermo_from_rhop(GasState Q) const
    {
        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        Q.T = Q.p/(Q.rho*_R_N2);
        Q.u = 2.5*_R_N2*Q.T;

        /// THINK ABOUT: PM/RJG
        // Here we assume we know something about the temperatures in each mode
        foreach (imode; 0 .. _n_modes) {
            Q.u_modes[imode] = (Avogadro_number/_M_N2) * Q.massf[imode] * _pseudoSpecies[imode].energy(Q);
        }
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
        // THINK ABOUT: PM/RJG
        // We will need to generalise this expression when we have more than one regular species.
        // For the present, assume we are working with N2.
        Q.a = sqrt(_gamma*_R_N2*Q.T);
    }

    override void update_trans_coeffs(GasState Q)
    {
        // THINK ABOUT: PM/RJG
        // We will need to implement this when we've made a decision.
        // For the present, only inviscid simulations are valid.
        Q.mu = 0.0;
        Q.k = 0.0;
        Q.k_modes[] = to!number(0.0);
    }

        override number dudT_const_v(in GasState Q) const
    {
        return to!number(_R_N2/(_gamma - 1.0));
    }
    override number dhdT_const_p(in GasState Q) const
    {
        throw new Error("ERROR: PseudoSpeciesGas:dhdT_const_p NOT IMPLEMENTED.");
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        throw new Error("ERROR: PseudoSpeciesGas:dpdrho_const_T NOT IMPLEMENTED.");
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(_R_N2);
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
    double _R_N2 = 296.805; // gas constant for N2
    double _M_N2 = 0.0280134;
    double _gamma = 7./5.; // ratio of specific heats.

    PseudoSpecies[] _pseudoSpecies;
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
