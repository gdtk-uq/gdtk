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
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.pseudo_species;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.cea_viscosity;
import gas.diffusion.cea_therm_cond;
import gas.diffusion.sutherland_viscosity;
import gas.diffusion.sutherland_therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;

class PseudoSpeciesGas : GasModel {

// This class must use _n_species inherited from GasModel

public:
    this(lua_State *L)
    {
        type_str = "PseudoSpeciesGas";
        auto nPseudoSpecies = getInt(L, "number_pseudo_species");
        _n_species = cast(uint) nPseudoSpecies;
        _n_modes = 0;
        _species_names.length = _n_species;
        auto nParents = getInt(L, "number_parents_species");
        _n_parents = cast(uint) nParents;
        _parents_names.length = _n_parents;

        // Create Parents Species
        lua_getglobal(L, "db");
        _parentSpecies = new ParentsSpecies(L);
        _parentsQ = GasState(_parentSpecies);
        lua_pop(L, 1);

        // Create Pseudo Species
        lua_getglobal(L, "pseudo_species");
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
            _parentIdx ~= _pseudoSpecies[$-1].parentIdx;
        }
        lua_pop(L, 1);
        create_species_reverse_lookup();

    }

    override string toString() const
    {
        char[] repr;
        repr ~= "PseudoSpeciesGas=()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.rho = Q.p/(R_mix*Q.T);

        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T + uNoneq;
    }

    override void update_thermo_from_rhou(ref GasState Q)
    {
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.T = (Q.u - uNoneq)/Cv_mix;
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
    }

    override void update_thermo_from_rhoT(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T + uNoneq;
    }

    override void update_thermo_from_rhop(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.T = Q.p/(Q.rho*R_mix);
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T + uNoneq;
    }

    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        throw new Error("ERROR: PseudoSpeciesGas:update_thermo_from_ps NOT IMPLEMENTED.");
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        throw new Error("ERROR: PseudoSpeciesGas:update_thermo_from_hs NOT IMPLEMENTED.");
    }

    override void update_sound_speed(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        auto Cv_mix = Cv(Q);
        auto gamma_mix = (R_mix/Cv_mix) + 1;
        Q.a = sqrt(gamma_mix*R_mix*Q.T);
    }

    override void update_trans_coeffs(ref GasState Q)
    {
        // Q is the GasState for the set of pseudo-species
        // Determine the parentsGasState from Q
        updateParentsGastate(Q, _parentsQ);
        _parentSpecies._viscModel.update_viscosity(_parentsQ);
        _parentSpecies._thermCondModel.update_thermal_conductivity(_parentsQ);
        // Assign to each pseudo species viscosity and conductivity of the parent
        // species
        Q.mu = _parentsQ.mu;
        Q.k = _parentsQ.k;
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

    @nogc void updateParentsGastate(ref const(GasState) Q, ref GasState parentsQ)
    {
      // Update translational temperature
      parentsQ.T = Q.T;
      // Update parents mass fractions by suming pseudo species massf
      parentsQ.massf[] = to!number(0.0);
      foreach (int psp ; 0 .. _n_species) {
        int prtIdx = _pseudoSpecies[psp].parentIdx;
        parentsQ.massf[prtIdx] += Q.massf[psp];
      }
    }

private:

    int _n_parents;
    string[] _parents_names;
    PseudoSpecies[] _pseudoSpecies;
    double[] _R;
    int[] _dof;
    int[] _parentIdx;
    ParentsSpecies _parentSpecies;
    GasState _parentsQ;

    @nogc number energyInNoneq(ref GasState Q) const {
        number uNoneq = 0.0;
        foreach (isp; 0 .. _n_species) {
            uNoneq += Q.massf[isp] * _pseudoSpecies[isp].energy(Q);
        }
        return uNoneq;
    }
}


// Class to describe the parent of a pseudo species in charge of the transport
// properties
// It is inspired by the class ThermallyPerfectGas

class ParentsSpecies: GasModel {
public:
    this(lua_State* L)
    // Construct the model from parameters that are contained in a Lua interpreter.
    {
        getArrayOfStrings(L, "parents", _species_names);
        _n_species = cast(uint) _species_names.length;

        // 0. Initialise private work arrays

        // 1. Initialise gas constants from molecular mass
        _R.length = _n_species;
        _mol_masses.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _mol_masses[isp] = getDouble(L, -1, "M");
            _R[isp] = R_universal/_mol_masses[isp];
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        // 1b. Gather L-J parameters
        _LJ_sigmas.length = _n_species;
        _LJ_epsilons.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _LJ_sigmas[isp] = getDouble(L, -1, "sigma");
            _LJ_epsilons[isp] = getDouble(L, -1, "epsilon");
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        // 4. Set the viscosity model
        Viscosity[] vms;
        foreach ( isp; 0.._n_species ) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            lua_getfield(L, -1, "viscosity");
            lua_getfield(L, -1, "model");
            string model = to!string(luaL_checkstring(L, -1));
            lua_pop(L, 1);
            switch (model) {
            case "CEA":
                vms ~= createCEAViscosity(L);
                break;
            case "Sutherland":
                vms ~= createSutherlandViscosity(L);
                break;
            default:
                string errMsg = format("The viscosity model '%s' is not available.\n", model);
                errMsg ~= format("This error occurred for species %d when constructing "
                                 ~"the thermally perfect gas mix.\n", isp);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        _viscModel = new WilkeMixingViscosity(vms, _mol_masses);
        // 5. Set the thermal conductivity model
        ThermalConductivity[] tcms;
        foreach (isp; 0.._n_species ) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            lua_getfield(L, -1, "thermal_conductivity");
            lua_getfield(L, -1, "model");
            string model = to!string(luaL_checkstring(L, -1));
            lua_pop(L, 1);
            switch (model) {
            case "CEA":
                tcms ~= createCEAThermalConductivity(L);
                break;
            case "Sutherland":
                tcms ~=  createSutherlandThermalConductivity(L);
                break;
            default:
                string errMsg = format("The thermal conductivity model '%s' "
                                       ~"is not available.\n", model);
                errMsg ~= format("This error occurred for species %d when constructing "
                                 ~"the thermally perfect gas mix.\n", isp);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        _thermCondModel = new WilkeMixingThermCond(tcms, _mol_masses);
        create_species_reverse_lookup();
        //
        version(complex_numbers) {
            throw new Error("Do not use with complex numbers.");
        }
    } // end constructor using Lua interpreter

    override void update_trans_coeffs(ref GasState Q)
    {
        _viscModel.update_viscosity(Q);
        _thermCondModel.update_thermal_conductivity(Q);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
    }
    override void update_thermo_from_rhou(ref GasState Q)
    {
    }
    override void update_thermo_from_rhoT(ref GasState Q)
    {
    }
    override void update_thermo_from_rhop(ref GasState Q)
    {
    }
    override void update_thermo_from_ps(ref GasState Q, number s)
    {
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
    }
    override void update_sound_speed(ref GasState Q)
    {
    }
    override number dudT_const_v(in GasState Q)
    {
        return to!number(0.);
    }
    override number dhdT_const_p(in GasState Q)
    {
        return to!number(0.);
    }
    override number dpdrho_const_T(in GasState Q)
    {
        return to!number(0.);
    }
    override number gas_constant(in GasState Q)
    {
        return to!number(0.);
    }
    override number internal_energy(in GasState Q)
    {
        return to!number(0.);
    }
    override number enthalpy(in GasState Q)
    {
        return to!number(0.);
    }
    override number enthalpy(in GasState Q, int isp)
    {
        return to!number(0.);
    }
    override number entropy(in GasState Q)
    {
        return to!number(0.);
    }
    override number entropy(in GasState Q, int isp)
    {
        return to!number(0.);
    }

private:
    double[] _R;
    WilkeMixingViscosity _viscModel;
    WilkeMixingThermCond _thermCondModel;

} // end class PseudoSpeciesGas


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
        doLuaFile(L, "sample-data/pseudo-species-49-components.lua");
        auto gm = new PseudoSpeciesGas(L);
        // writeln("gm.n_species = ", gm.n_species);
        assert(gm.n_species == 49, failedUnitTest());
        // writeln("gm._n_parents = ", gm._n_parents);
        assert(gm._n_parents == 2, failedUnitTest());

        auto gd = GasState(gm);
        // writeln("gd.massf.length = ", gd.massf.length);
        assert(gd.massf.length == 49, failedUnitTest());
        gd.massf[] = 0.0;
        gd.massf[1] = 0.20908519640652;
        gd.massf[2] = 0.12949211438053;
        gd.massf[0] = 1.0 - (gd.massf[1] + gd.massf[2]);
        gd.p = 101325.0;
        gd.T = 7000.0;
        // writeln(gd.massf);
        // writeln(gm.Cv(gd));
        assert(isClose(gm.Cv(gd), 840.445, 1.0e-3), failedUnitTest());
        gm.update_thermo_from_pT(gd);
        // writeln("gd.u = ", gd.u);
        assert(isClose(gd.u, -5.36274e+06, 1.0e-3), failedUnitTest());

        gm.update_trans_coeffs(gd);
        // writeln("gd.mu = ", gd.mu);
        // writeln("gd.k = ", gd.k);
        assert(isClose(gd.mu, 0.000187585, 1.0e-3), failedUnitTest());
        assert(isClose(gd.k, 0.403616, 1.0e-3), failedUnitTest());

        return 0;
    } // end main()
}
