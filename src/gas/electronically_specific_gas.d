/**
 * Author: Brad Semple
 *
 */

module gas.electronically_specific_gas;

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
import gas.electronic_species;

class ElectronicallySpecificGas : GasModel {
public:
    this(lua_State *L) 
    {
        auto nElecSpecies = getInt(L, LUA_GLOBALSINDEX, "number_electronic_species");
        _n_species = cast(uint) nElecSpecies; //go to 15 or 16?
        _n_modes = 0; //dont know what to do with this?
        _species_names.length = _n_species; 

        lua_getfield(L, LUA_GLOBALSINDEX, "electronic_species");
        foreach (isp; 0 .. _n_species) {
            lua_rawgeti(L, -1, isp);
            if (lua_isnil(L, -1)) {
                string msg = format("There was an error when attempting to information about pseudo-species %d.\n", isp);
                throw new Error(msg);
            }
            _electronicSpecies ~= createElectronicSpecies(L);
            lua_pop(L, 1);
            _species_names[isp] = _electronicSpecies[$-1].name;
            _mol_masses ~= _electronicSpecies[$-1].mol_mass;
            _R ~= R_universal/_electronicSpecies[$-1].mol_mass;
            _level ~= _electronicSpecies[$-1].level;
            _group_degeneracy ~= _electronicSpecies[$-1].group_degeneracy;
        }
        lua_pop(L, 1);
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "ElectronicSpeciesGas=()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.rho = Q.p/(R_mix*Q.T);

        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T - uNoneq;
    }

    override void update_thermo_from_rhou(GasState Q)
    {
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.T = (Q.u + uNoneq)/Cv_mix; 
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
    }

    override void update_thermo_from_rhoT(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T - uNoneq;
    }

    override void update_thermo_from_rhop(GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.T = Q.p/(Q.rho*R_mix);
        auto uNoneq = energyInNoneq(Q);
        auto Cv_mix = Cv(Q);
        Q.u = Cv_mix*Q.T - uNoneq;
    }

    override void update_thermo_from_ps(GasState Q, number s)
    {
        throw new Error("ERROR: ElectronicallySpecificGas:update_thermo_from_ps NOT IMPLEMENTED.");
    }

    override void update_thermo_from_hs(GasState Q, number h, number s)
    {
        throw new Error("ERROR: ElectronicallySpecificGas:update_thermo_from_hs NOT IMPLEMENTED.");
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
        Q.mu = 0.0;
        Q.k = 0.0;
    }

    override number dudT_const_v(in GasState Q) const
    {
        throw new Error("ERROR: ElectronicallySpecificGas:dudT_const_v NOT IMPLEMENTED.");
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return gas_constant(Q) + dudT_const_v(Q);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        throw new Error("ERROR: ElectronicallySpecificGas:dpdrho_const_T NOT IMPLEMENTED.");
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
        throw new Error("ERROR: ElecontricallySpecificGas:entropy NOT IMPLEMENTED.");
    }

private:
    ElectronicSpecies[] _electronicSpecies;
    double[] _R;
    int[] _level;
    int[] _group_degeneracy;

    number energyInNoneq(GasState Q) const {
        number uNoneq = 0.0;
        foreach (isp; 0 .. _n_species) {
            uNoneq += Q.massf[isp] * _electronicSpecies[isp].energy(Q);
        }
        return uNoneq;
    }
}


version(electronically_specific_gas_test) {
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
        doLuaFile(L, "sample-data/brad_elec_composition.lua");
        auto gm = new ElectronicSpeciesGas(L);
        auto gd = new GasState(19,19);

        gd.massf[] = 0.0;
        gd.massf[0] = 0.81403036047055; //initialises massf of NI
        gd.massf[9] = 0.185968105968037; //initialises massf of OI
        gd.massf[18] = 1.0 - (gd.massf[0] + gd.massf[9]); //tiny massf for free electron

        gd.p = 101325.0;
        gd.T = 7000.0;
        gd.T_modes[0]=10000.0;

        gm.update_thermo_from_pT(gd);

        return 0;
    }
    
}
