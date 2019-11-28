// gas_cwrap.d
//
// This particular module will be compiled into a loadable library and
// provides a C-level API for the GasModel, GasState and ThermochemicalReactor
// classes that may be called from any language with a C-foreign-function-interface.
//
// PJ 2019-07-24: just enough to try integration with the gas makefile.
//

import core.runtime;
import core.stdc.string;
import std.string;
import std.stdio;
import std.format;
import std.conv;
import std.algorithm;

import gas.gas_model;
import gas.gas_state;
import gas.init_gas_model;

import kinetics.thermochemical_reactor;
import kinetics.init_thermochemical_reactor;

import gasflow;


// We will accumulate GasModel and GasState objects in these arrays
// and use the indices as handles in the scripting language.
GasModel[] gas_models;
GasState[] gas_states;
ThermochemicalReactor[] thermochemical_reactors;

extern (C) int cwrap_gas_init()
{
    Runtime.initialize();
    gas_models.length = 0;
    gas_states.length = 0;
    thermochemical_reactors.length = 0;
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

extern (C) int gas_model_new(char* file_name)
{
    try {
        GasModel gm = init_gas_model(to!string(file_name));
        gas_models ~= gm;
        return to!int(gas_models.length - 1);
    } catch (Exception e) {
        writeln("Failed to construct a new GasModel.");
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_n_species(int gm_i)
{
    try {
        int n = to!int(gas_models[gm_i].n_species);
        return n;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_n_modes(int gm_i)
{
    try {
        int n = to!int(gas_models[gm_i].n_modes);
        return n;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
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
        writeln("Exception message: ", e.msg);
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
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern (C) int gas_state_new(int gm_i)
{
    try {
        GasState gs = new GasState(gas_models[gm_i]);
        gas_states ~= gs;
        return to!int(gas_states.length - 1);
    } catch (Exception e) {
        writeln("Failed to construct a new GasState.");
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_set_scalar_field(int gs_i, char* field_name, double value)
{
    try {
        GasState gs = gas_states[gs_i];
        string name = to!string(field_name); 
        switch (name) {
        case "rho":
            gs.rho = value;
            break;
        case "p":
            gs.p = value;
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
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_scalar_field(int gs_i, char* field_name, double* value)
{
    try {
        GasState gs = gas_states[gs_i];
        string name = to!string(field_name);
        *value = 0.0;
        switch (name) {
        case "rho":
            *value = gs.rho;
            break;
        case "p":
            *value = gs.p;
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
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_set_array_field(int gs_i, char* field_name, double* values, int n)
{
    try {
        GasState gs = gas_states[gs_i];
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
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_state_get_array_field(int gs_i, char* field_name, double* values, int n)
{
    try {
        GasState gs = gas_states[gs_i];
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
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern (C) int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i)
{
    try {
        gas_models[gm_i].update_thermo_from_pT(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i)
{
    try {
        gas_models[gm_i].update_thermo_from_rhou(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i)
{
    try {
        gas_models[gm_i].update_thermo_from_rhoT(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i)
{
    try {
        gas_models[gm_i].update_thermo_from_rhop(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_ps(int gm_i, int gs_i, double s)
{
    try {
        gas_models[gm_i].update_thermo_from_ps(gas_states[gs_i], s);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_thermo_from_hs(int gm_i, int gs_i, double h, double s)
{
    try {
        gas_models[gm_i].update_thermo_from_hs(gas_states[gs_i], h, s);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_sound_speed(int gm_i, int gs_i)
{
    try {
        gas_models[gm_i].update_sound_speed(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_update_trans_coeffs(int gm_i, int gs_i)
{
    try {
        gas_models[gm_i].update_trans_coeffs(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_Cv(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].Cv(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_Cp(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].Cp(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_dpdrho_const_T(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].dpdrho_const_T(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_R(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].R(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_internal_energy(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].internal_energy(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_enthalpy(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].enthalpy(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_entropy(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].entropy(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_molecular_mass(int gm_i, int gs_i, double* result)
{
    try {
        *result = gas_models[gm_i].molecular_mass(gas_states[gs_i]);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_enthalpy_isp(int gm_i, int gs_i, int isp,
                                                double* result)
{
    try {
        *result = gas_models[gm_i].enthalpy(gas_states[gs_i], isp);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_entropy_isp(int gm_i, int gs_i, int isp,
                                               double* result)
{
    try {
        *result = gas_models[gm_i].entropy(gas_states[gs_i], isp);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_gibbs_free_energy_isp(int gm_i, int gs_i, int isp,
                                                         double* result)
{
    try {
        *result = gas_models[gm_i].gibbs_free_energy(gas_states[gs_i], isp);
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
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
        writeln("Exception message: ", e.msg);
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
        writeln("Exception message: ", e.msg);
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
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int gas_model_gas_state_get_conc(int gm_i, int gs_i, double* conc)
{
    try {
        double[] my_conc;
        auto nsp = gas_models[gm_i].n_species;
        my_conc.length = nsp;
        gas_models[gm_i].massf2conc(gas_states[gs_i], my_conc);
        foreach (i; 0 .. nsp) { conc[i] = my_conc[i]; }
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern (C) int thermochemical_reactor_new(int gm_i, char* filename1, char* filename2)
{
    try {
        auto cr = init_thermochemical_reactor(gas_models[gm_i], to!string(filename1),
                                              to!string(filename2));
        thermochemical_reactors ~= cr;
        return to!int(thermochemical_reactors.length - 1);
    } catch (Exception e) {
        writeln("Failed to construct a new ThermochemicalReactor.");
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int thermochemical_reactor_gas_state_update(int cr_i, int gs_i,
                                                       double t_interval,
                                                       double* dt_suggest)
{
    try {
        double dummyDouble;
        // Extra parameters are not considered, presently. PJ 2017-04-22, 2019-11-26
        double[maxParams] params;
        double my_dt_suggest = *dt_suggest;
        thermochemical_reactors[cr_i](gas_states[gs_i], t_interval, my_dt_suggest,
                                      dummyDouble, params);
        *dt_suggest = my_dt_suggest;
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

//---------------------------------------------------------------------------

extern(C) int gasflow_shock_ideal(int state1_id, double Vs, int state2_id, int gm_id,
                                  double* results)
{
    try {
        double[] my_results = shock_ideal(gas_states[state1_id], Vs,
                                          gas_states[state2_id], gas_models[gm_id]);
        results[0] = my_results[0];
        results[1] = my_results[1];
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern(C) int gasflow_normal_shock(int state1_id, double Vs, int state2_id, int gm_id,
                                   double* results, double rho_tol, double T_tol)
{
    try {
        double[] my_results = normal_shock(gas_states[state1_id], Vs,
                                           gas_states[state2_id], gas_models[gm_id],
                                           rho_tol, T_tol);
        results[0] = my_results[0];
        results[1] = my_results[1];
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}
    
