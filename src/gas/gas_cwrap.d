// gas_cwrap.d
//
// This particular module will be compiled into a loadable library and
// provides a C-level API for the GasModel, GasState and ThermochemicalReactor
// classes that may be called from any language with a C-foreign-function-interface.
//
// PJ 2019-07-24: just enough to try integration with the gas makefile.
//    2023-06-03: add function to get all thermo scalars at once.
//

import core.runtime;
import core.stdc.string;
import std.string;
import std.stdio;
import std.format;
import std.conv;
import std.algorithm;
import std.math;

import gas.gas_model;
import gas.gas_state;
import gas.init_gas_model;

import kinetics.thermochemical_reactor;
import kinetics.init_thermochemical_reactor;
import kinetics.reaction_mechanism;

import gasdyn.gasflow;


// We will accumulate GasModel and GasState objects in these arrays
// and use the indices as handles in the scripting language.
GasModel[] gas_models;
GasState*[] gas_states;
ThermochemicalReactor[] thermochemical_reactors;
ReactionMechanism[] reaction_mechanisms;

extern (C) int cwrap_gas_init()
{
    Runtime.initialize();
    version(enable_fp_exceptions) {
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
    }
    gas_models.length = 0;
    gas_states.length = 0;
    thermochemical_reactors.length = 0;
    reaction_mechanisms.length = 0;
    return 0;
}

shared static this()
{
    // writeln("libgasmodule.so shared static this");
}

shared static ~this()
{
    // writeln("libgasmodule.so shared static ~this");
}

//---------------------------------------------------------------------------

extern (C) int gas_model_new(const char* file_name)
{
    try {
        GasModel gm = init_gas_model(to!string(file_name));
        gas_models ~= gm;
        return to!int(gas_models.length - 1);
    } catch (Exception e) {
        stderr.writeln("Failed to construct a new GasModel.");
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_type_str(int gm_i, char* dest_str, int n)
{
    // The gas model's type string will be copied into char* array.
    // It is presumed that sufficient space (n chars, including \0) was allocated previously.
    try {
        const char* src_str = cast(const char*) gas_models[gm_i].type_str.toStringz;
        strncpy(dest_str, src_str, n);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_n_species(int gm_i)
{
    try {
        int n = to!int(gas_models[gm_i].n_species);
        return n;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_n_modes(int gm_i)
{
    try {
        int n = to!int(gas_models[gm_i].n_modes);
        return n;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_species_name_and_length(int gm_i, int isp, char* dest_name, int* n)
{
    // The isp-th species name will be copied into char* array.
    // It is presumed that sufficient space (n chars, including \0) was allocated previously.
    try {
        string sname = gas_models[gm_i].species_name(isp);
        const char* src_name = cast(const char*) sname.toStringz;
        strncpy(dest_name, src_name, sname.length);
        *n = cast(int) sname.length;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_species_name(int gm_i, int isp, char* dest_name, int n)
{
    // The isp-th species name will be copied into char* array.
    // It is presumed that sufficient space (n chars, including \0) was allocated previously.
    try {
        char* src_name = cast(char*) gas_models[gm_i].species_name(isp).toStringz;
        strncpy(dest_name, src_name, n);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_mol_masses(int gm_i, double* mm)
{
    // Return the species molecular mass as an array.
    try {
        foreach(i, mass; gas_models[gm_i].mol_masses) { mm[i] = mass; }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern (C) int gas_state_new(int gm_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = new GasState(gm);
        // For a single-species gas, the user often ignores the mass-fraction array.
        // Let's set the single mass fraction to 1,
        // so that we do not have to NaN being printed.
        if (gm.n_species == 1) { gs.massf[0] = 1.0; }
        gas_states ~= gs;
        return to!int(gas_states.length - 1);
    } catch (Exception e) {
        stderr.writeln("Failed to construct a new GasState.");
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_set_scalar_field(int gs_i, const char* field_name, double value)
{
    try {
        GasState* gs = gas_states[gs_i];
        string name = to!string(field_name);
        switch (name) {
        case "rho":
            gs.rho = value;
            break;
        case "p":
            gs.p = value;
            break;
        case "p_e":
            gs.p_e = value;
            break;
        case "T":
            gs.T = value;
            break;
        case "u":
            gs.u = value;
            break;
        default:
            string msg = format("Cannot set field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_scalar_field(int gs_i, const char* field_name, double* value)
{
    try {
        GasState* gs = gas_states[gs_i];
        string name = to!string(field_name);
        *value = 0.0;
        switch (name) {
        case "rho":
            *value = gs.rho;
            break;
        case "p":
            *value = gs.p;
            break;
        case "p_e":
            *value = gs.p_e;
            break;
        case "T":
            *value = gs.T;
            break;
        case "u":
            *value = gs.u;
            break;
        case "a":
            *value = gs.a;
            break;
        case "k":
            *value = gs.k;
            break;
        case "mu":
            *value = gs.mu;
            break;
        default:
            string msg = format("Unavailable field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_thermo_scalars(int gs_i, double* values)
{
    try {
        GasState* gs = gas_states[gs_i];
        values[0] = gs.rho;
        values[1] = gs.p;
        values[2] = gs.T;
        values[3] = gs.u;
        values[4] = gs.a;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_set_array_field(int gs_i, const char* field_name, double* values, int n)
{
    try {
        GasState* gs = gas_states[gs_i];
        string name = to!string(field_name);
        switch (name) {
        case "massf":
            foreach (i; 0 .. n) { gs.massf[i] = values[i]; }
            break;
        case "u_modes":
            foreach (i; 0 .. n) { gs.u_modes[i] = values[i]; }
            break;
        case "T_modes":
            foreach (i; 0 .. n) { gs.T_modes[i] = values[i]; }
            break;
        default:
            string msg = format("Cannot set field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_array_field(int gs_i, const char* field_name, double* values, int n)
{
    try {
        GasState* gs = gas_states[gs_i];
        string name = to!string(field_name);
        switch (name) {
        case "massf":
            foreach (i; 0 .. n) { values[i] = gs.massf[i]; }
            break;
        case "u_modes":
            foreach (i; 0 .. n) { values[i] = gs.u_modes[i]; }
            break;
        case "T_modes":
            foreach (i; 0 .. n) { values[i] = gs.T_modes[i]; }
            break;
        case "k_modes":
            foreach (i; 0 .. n) { values[i] = gs.k_modes[i]; }
            break;
        default:
            string msg = format("Unavailable field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_ceaSavedData_field(int gs_i, const char* field_name, double* value)
{
    try {
        GasState* gs = gas_states[gs_i];
        if (gs.ceaSavedData is null) {
            throw new Exception("No available ceaSavedData.");
        }
        string name = to!string(field_name);
        *value = 0.0;
        switch (name) {
        case "rho":
            *value = gs.ceaSavedData.rho;
            break;
        case "p":
            *value = gs.ceaSavedData.p;
            break;
        case "T":
            *value = gs.ceaSavedData.T;
            break;
        case "u":
            *value = gs.ceaSavedData.u;
            break;
        case "h":
            *value = gs.ceaSavedData.h;
            break;
        case "Mmass":
            *value = gs.ceaSavedData.Mmass;
            break;
        case "Rgas":
            *value = gs.ceaSavedData.Rgas;
            break;
        case "gamma":
            *value = gs.ceaSavedData.gamma;
            break;
        case "a":
            *value = gs.ceaSavedData.a;
            break;
        case "Cp":
            *value = gs.ceaSavedData.Cp;
            break;
        case "s":
            *value = gs.ceaSavedData.s;
            break;
        case "k":
            *value = gs.ceaSavedData.k;
            break;
        case "mu":
            *value = gs.ceaSavedData.mu;
            break;
        default:
            string msg = format("Unavailable field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_ceaSavedData_massf(int gs_i, const char* species_name, double* value)
{
    try {
        GasState* gs = gas_states[gs_i];
        string name = to!string(species_name);
        *value = gs.ceaSavedData.massf[name];
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_ceaSavedData_species_names(int gs_i, char* dest_str, int n)
{
    // All of the species names will be copied into char* array.
    // It is presumed that sufficient space (n chars, including \0) was allocated previously.
    try {
        string src_str = gas_states[gs_i].ceaSavedData.massf.keys().join("\t");
        strncpy(dest_str, src_str.toStringz, n);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_copy_values(int gs_to_i, int gs_from_i)
{
    try {
        gas_states[gs_to_i].copy_values_from(*gas_states[gs_from_i]);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern (C) int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(gs.p)) { valid = false; }
        if (!isFinite(gs.T)) { valid = false; }
        foreach (i; 0 .. gm.n_modes) {
            if (!isFinite(gs.T_modes[i])) { valid = false; }
        }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for pT update."); }
        gm.update_thermo_from_pT(*gs);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(gs.rho)) { valid = false; }
        if (!isFinite(gs.u)) { valid = false; }
        foreach (i; 0 .. gm.n_modes) {
            if (!isFinite(gs.u_modes[i])) { valid = false; }
        }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for rhou update."); }
        gm.update_thermo_from_rhou(*gs);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(gs.rho)) { valid = false; }
        if (!isFinite(gs.T)) { valid = false; }
        foreach (i; 0 .. gm.n_modes) {
            if (!isFinite(gs.T_modes[i])) { valid = false; }
        }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for rhoT update."); }
        gm.update_thermo_from_rhoT(*gs);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(gs.rho)) { valid = false; }
        if (!isFinite(gs.p)) { valid = false; }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for rhop update."); }
        gm.update_thermo_from_rhop(*gs);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_ps(int gm_i, int gs_i, double s)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(s)) { valid = false; }
        if (!isFinite(gs.p)) { valid = false; }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for ps update."); }
        gm.update_thermo_from_ps(*gs, s);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_hs(int gm_i, int gs_i, double h, double s)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(s)) { valid = false; }
        if (!isFinite(h)) { valid = false; }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for hs update."); }
        gm.update_thermo_from_hs(*gs, h, s);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_sound_speed(int gm_i, int gs_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(gs.rho)) { valid = false; }
        if (!isFinite(gs.p)) { valid = false; }
        if (!isFinite(gs.T)) { valid = false; }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for sound_speed update."); }
        gm.update_sound_speed(*gs);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_trans_coeffs(int gm_i, int gs_i)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        bool valid = true;
        if (!isFinite(gs.rho)) { valid = false; }
        if (!isFinite(gs.p)) { valid = false; }
        if (!isFinite(gs.T)) { valid = false; }
        if (gm.n_species > 1) {
            foreach (i; 0 .. gm.n_species) {
                if (!isFinite(gs.massf[i])) { valid = false; }
            }
        }
        if (!valid) { throw new Exception("Gas state is not valid for trans_coeffs update."); }
        gm.update_trans_coeffs(*gs);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_Cv(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].Cv(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_Cp(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].Cp(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_dpdrho_const_T(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].dpdrho_const_T(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_R(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].R(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_gamma(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].gamma(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_Prandtl(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].Prandtl(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_internal_energy(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].internal_energy(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_enthalpy(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].enthalpy(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_entropy(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].entropy(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_molecular_mass(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].molecular_mass(*(gas_states[gs_i]));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_binary_diffusion_coefficients(int gm_i, int gs_i, double* dij)
{
    try {
        GasModel gm = gas_models[gm_i];
        GasState* gs = gas_states[gs_i];
        size_t nsp = gm.n_species;
        double[][] D;
        D.length = nsp;
        foreach(ref Di; D) Di.length=nsp;

        gm.binary_diffusion_coefficients(*gs, D);
        foreach (i; 0 .. nsp) {
            foreach(j; 0 .. nsp) {
                dij[i*nsp + j] = D[i][j];
            }
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}


extern (C) int gas_model_gas_state_enthalpy_isp(int gm_i, int gs_i, int isp,
                                                double* result)
{
    try {
        *result = gas_models[gm_i].enthalpy(*(gas_states[gs_i]), isp);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_entropy_isp(int gm_i, int gs_i, int isp,
                                               double* result)
{
    try {
        *result = gas_models[gm_i].entropy(*gas_states[gs_i], isp);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_gibbs_free_energy_isp(int gm_i, int gs_i, int isp,
                                                         double* result)
{
    try {
        *result = gas_models[gm_i].gibbs_free_energy(*(gas_states[gs_i]), isp);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_massf2molef(int gm_i, double* massf, double* molef)
{
    try {
        double[] my_massf, my_molef;
        auto nsp = gas_models[gm_i].n_species;
        my_massf.length = nsp;
        my_molef.length = nsp;
        foreach (i; 0 .. nsp) { my_massf[i] = massf[i]; }
        massf2molef(my_massf, gas_models[gm_i].mol_masses, my_molef);
        foreach (i; 0 .. nsp) { molef[i] = my_molef[i]; }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_molef2massf(int gm_i, double* molef, double* massf)
{
    try {
        double[] my_massf, my_molef;
        auto nsp = gas_models[gm_i].n_species;
        my_massf.length = nsp;
        my_molef.length = nsp;
        foreach (i; 0 .. nsp) { my_molef[i] = molef[i]; }
        molef2massf(my_molef, gas_models[gm_i].mol_masses, my_massf);
        foreach (i; 0 .. nsp) { massf[i] = my_massf[i]; }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_get_molef(int gm_i, int gs_i, double* molef)
{
    try {
        double[] my_molef;
        auto nsp = gas_models[gm_i].n_species;
        my_molef.length = nsp;
        massf2molef(gas_states[gs_i].massf, gas_models[gm_i].mol_masses, my_molef);
        foreach (i; 0 .. nsp) { molef[i] = my_molef[i]; }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_get_conc(int gm_i, int gs_i, double* conc)
{
    try {
        double[] my_conc;
        auto nsp = gas_models[gm_i].n_species;
        my_conc.length = nsp;
        gas_models[gm_i].massf2conc(*(gas_states[gs_i]), my_conc);
        foreach (i; 0 .. nsp) { conc[i] = my_conc[i]; }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern (C) int thermochemical_reactor_new(int gm_i, const char* filename1, const char* filename2)
{
    try {
        auto cr = init_thermochemical_reactor(gas_models[gm_i], to!string(filename1),
                                              to!string(filename2));
        thermochemical_reactors ~= cr;
        return to!int(thermochemical_reactors.length - 1);
    } catch (Exception e) {
        stderr.writeln("Failed to construct a new ThermochemicalReactor.");
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int thermochemical_reactor_gas_state_update(int cr_i, int gs_i,
                                                       double t_interval,
                                                       double* dt_suggest)
{
    try {
        // Extra parameters are not considered, presently. PJ 2017-04-22, 2019-11-26
        double[maxParams] params;
        double my_dt_suggest = *dt_suggest;
        thermochemical_reactors[cr_i](*(gas_states[gs_i]), t_interval, my_dt_suggest, params);
        *dt_suggest = my_dt_suggest;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int thermochemical_reactor_eval_source_terms(int cr_i, int gm_i, int gs_i,
                                                        int nsp, int nmodes, double* source)
{
    try {
        auto reactor = thermochemical_reactors[cr_i];

        double[] sourceTerms;
        sourceTerms.length = nsp + nmodes;
        reactor.eval_source_terms(gas_models[gm_i], *(gas_states[gs_i]), sourceTerms);
        foreach (i,v; sourceTerms) {
            source[i] = v;
        }
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}


//---------------------------------------------------------------------------
extern (C) int reaction_mechanism_new(int gm_i, const char* filename)
{
    try {
        auto rm = createReactionMechanism(to!string(filename), gas_models[gm_i], 300.0, 30000.0);
        reaction_mechanisms ~= rm;
        return to!int(reaction_mechanisms.length - 1);
    } catch (Exception e) {
        stderr.writeln("Failed to construct a new ReactionMechanism.");
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int reaction_mechanism_n_reactions(int rm_i)
{
    try {
        int n = to!int(reaction_mechanisms[rm_i].n_reactions);
        return n;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int reaction_mechanism_tickrates(int rm_i, int gm_i, int gs_i, double* forwardrates, double* backwardrates)
{
    // Return the individual reaction rates as an array
    try {
        auto rm = reaction_mechanisms[rm_i];
        auto gm = gas_models[gm_i];
        auto gs = gas_states[gs_i];
        
        double[] forward_tickrates, backward_tickrates, concentrations;
        forward_tickrates.length = rm.n_reactions;
        backward_tickrates.length = rm.n_reactions;
        concentrations.length = gm.n_species;

        rm.eval_rate_constants(gm, *gs);
        gm.massf2conc(*gs, concentrations);
        rm.eval_tickrates(concentrations, forward_tickrates, backward_tickrates);
        foreach(i; 0 .. rm.n_reactions) forwardrates[i] = forward_tickrates[i];
        foreach(i; 0 .. rm.n_reactions) backwardrates[i] = backward_tickrates[i];

        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}


//---------------------------------------------------------------------------
extern(C)
int gasflow_shock_ideal(int state1_id, double vs, int state2_id, int gm_id,
                        double* results)
{
    try {
        double[] my_results = shock_ideal(*(gas_states[state1_id]), vs,
                                          *(gas_states[state2_id]), gas_models[gm_id]);
        results[0] = my_results[0]; // v2
        results[1] = my_results[1]; // vg
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_normal_shock(int state1_id, double vs, int state2_id, int gm_id,
                         double* results, double rho_tol, double T_tol)
{
    try {
        double[] my_results = normal_shock(*(gas_states[state1_id]), vs,
                                           *(gas_states[state2_id]), gas_models[gm_id],
                                           rho_tol, T_tol);
        results[0] = my_results[0]; // v2
        results[1] = my_results[1]; // vg
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_normal_shock_1(int state1_id, double vs, int state2_id, int gm_id,
                           double* results, double p_tol, double T_tol)
{
    try {
        double[] my_results = normal_shock_1(*(gas_states[state1_id]), vs,
                                             *(gas_states[state2_id]), gas_models[gm_id],
                                             p_tol, T_tol);
        results[0] = my_results[0]; // v2
        results[1] = my_results[1]; // vg
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_normal_shock_p2p1(int state1_id, double p2p1, int state2_id, int gm_id,
                              double* results)
{
    try {
        double[] my_results = normal_shock_p2p1(*(gas_states[state1_id]), p2p1,
                                                *(gas_states[state2_id]), gas_models[gm_id]);
        results[0] = my_results[0]; // v1
        results[1] = my_results[1]; // v2
        results[2] = my_results[2]; // vg
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_reflected_shock(int state2_id, double vg, int state5_id, int gm_id,
                            double* results)
{
    try {
        double vr_b = reflected_shock(*(gas_states[state2_id]), vg,
                                      *(gas_states[state5_id]), gas_models[gm_id]);
        results[0] = vr_b;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_expand_from_stagnation(int state0_id, double p_over_p0, int state1_id,
                                   int gm_id, double* results)
{
    try {
        double v = expand_from_stagnation(*(gas_states[state0_id]), p_over_p0,
                                          *(gas_states[state1_id]), gas_models[gm_id]);
        results[0] = v;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_expand_to_mach(int state0_id, double mach, int state1_id,
                           int gm_id, double* results)
{
    try {
        double v = expand_to_mach(*(gas_states[state0_id]), mach,
                                  *(gas_states[state1_id]), gas_models[gm_id]);
        results[0] = v;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_total_condition(int state1_id, double v1, int state0_id, int gm_id)
{
    try {
        total_condition(*(gas_states[state1_id]), v1,
                        *(gas_states[state0_id]), gas_models[gm_id]);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_pitot_condition(int state1_id, double v1, int state2pitot_id, int gm_id)
{
    try {
        pitot_condition(*(gas_states[state1_id]), v1,
                        *(gas_states[state2pitot_id]), gas_models[gm_id]);
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_steady_flow_with_area_change(int state1_id, double v1, double a2_over_a1,
                                         int state2_id, int gm_id, double tol, double p2p1_min,
                                         double* results)
{
    try {
        double v2 = steady_flow_with_area_change(*(gas_states[state1_id]), v1, a2_over_a1,
                                                 *(gas_states[state2_id]), gas_models[gm_id],
                                                 tol, p2p1_min);
        results[0] = v2;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_finite_wave_dp(int state1_id, double v1, const char* characteristic, double p2,
                           int state2_id, int gm_id, int steps, double* results)
{
    try {
        double v2 = finite_wave_dp(*(gas_states[state1_id]), v1, to!string(characteristic), p2,
                                   *(gas_states[state2_id]), gas_models[gm_id], steps);
        results[0] = v2;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_finite_wave_dv(int state1_id, double v1, const char* characteristic, double v2_target,
                           int state2_id, int gm_id, int steps, double Tmin, double* results)
{
    try {
        double v2 = finite_wave_dv(*(gas_states[state1_id]), v1, to!string(characteristic), v2_target,
                                   *(gas_states[state2_id]), gas_models[gm_id], steps, Tmin);
        results[0] = v2;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_osher_riemann(int stateL_id, int stateR_id, double velL, double velR,
                          int stateLstar_id, int stateRstar_id,
                          int stateX0_id, int gm_id, double* results)
{
    try {
        double[5] my_results = osher_riemann(*(gas_states[stateL_id]), *(gas_states[stateR_id]), velL, velR,
                                             *(gas_states[stateLstar_id]), *(gas_states[stateRstar_id]),
                                             *(gas_states[stateX0_id]), gas_models[gm_id]);
        results[0] = my_results[0]; // pstar
        results[1] = my_results[1]; // wstar
        results[2] = my_results[2]; // wL
        results[3] = my_results[3]; // wR
        results[4] = my_results[4]; // velX0
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_osher_flux(int stateL_id, int stateR_id, double velL, double velR,
                       int gm_id, double* results)
{
    try {
        double[3] F = osher_flux(*(gas_states[stateL_id]), *(gas_states[stateR_id]), velL, velR, gas_models[gm_id]);
        results[0] = F[0]; // mass
        results[1] = F[1]; // x-momentum
        results[2] = F[2]; // energy
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_roe_flux(int stateL_id, int stateR_id, double velL, double velR,
                     int gm_id, double* results)
{
    try {
        double[3] F = roe_flux(*(gas_states[stateL_id]), *(gas_states[stateR_id]), velL, velR, gas_models[gm_id]);
        results[0] = F[0]; // mass
        results[1] = F[1]; // x-momentum
        results[2] = F[2]; // energy
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_lrivp(int stateL_id, int stateR_id, double velL, double velR,
                  int gmL_id, int gmR_id, double* wstar, double* pstar)
{
    try {
        double my_wstar = *wstar;
        double my_pstar = *pstar;
        lrivp(*(gas_states[stateL_id]), *(gas_states[stateR_id]), velL, velR,
              gas_models[gmL_id], gas_models[gmR_id], my_wstar, my_pstar);
        *wstar = my_wstar;
        *pstar = my_pstar;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_piston_at_left(int stateR_id, double velR, int gm_id,
                           double wstar, double* pstar)
{
    try {
        double my_pstar = *pstar;
        piston_at_left(*(gas_states[stateR_id]), velR, gas_models[gm_id], wstar, my_pstar);
        *pstar = my_pstar;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_piston_at_right(int stateL_id, double velL, int gm_id,
                            double wstar, double* pstar)
{
    try {
        double my_pstar = *pstar;
        piston_at_right(*(gas_states[stateL_id]), velL, gas_models[gm_id], wstar, my_pstar);
        *pstar = my_pstar;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_theta_oblique(int state1_id, double v1, double beta,
                          int state2_id, int gm_id, double* results)
{
    try {
        double[] my_results = theta_oblique(*(gas_states[state1_id]), v1, beta,
                                            *(gas_states[state2_id]), gas_models[gm_id]);
        results[0] = my_results[0]; // theta
        results[1] = my_results[1]; // v2
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_beta_oblique(int state1_id, double v1, double theta,
                          int gm_id, double* results)
{
    try {
        double beta = beta_oblique(*(gas_states[state1_id]), v1, theta, gas_models[gm_id]);
        results[0] = beta;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_theta_cone(int state1_id, double v1, double beta,
                       int state_c_id, int gm_id, double dtheta, double* results)
{
    try {
        double[2] my_results = theta_cone(*(gas_states[state1_id]), v1, beta,
                                          *(gas_states[state_c_id]), gas_models[gm_id], dtheta);
        results[0] = my_results[0]; // theta_c
        results[1] = my_results[1]; // v2_c
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C)
int gasflow_beta_cone(int state1_id, double v1, double theta,
                      int gm_id, double dtheta, double* results)
{
    try {
        double beta = beta_cone(*(gas_states[state1_id]), v1, theta, gas_models[gm_id], dtheta);
        results[0] = beta;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}
