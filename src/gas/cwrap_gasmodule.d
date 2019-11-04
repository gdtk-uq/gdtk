// cwrap_gasmodule.d
//
// This particular module will be compiled into a loadable library and
// provides a C-level API for the GasModel and GasState classes
// that may be called from any language with a C-foreign-function-interface.
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

// We will accumulate GasModel and GasState objects in these arrays
// and use the indices as handles in the scripting language.
GasModel[] gas_models;
GasState[] gas_states;

extern (C) int cwrap_gas_module_init()
{
    Runtime.initialize();
    gas_models.length = 0;
    gas_states.length = 0;
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

extern (C) char* gas_model_species_name(int gm_i, int isp, char* dest_name, int n)
{
    // The isp-th species name will be copied into char* array.
    // It is presumed that sufficient space (n chars, including \0) was allocated previously.
    try {
        char* src_name = cast(char*) gas_models[gm_i].species_name(isp).toStringz;
        return strncpy(dest_name, src_name, n);
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return cast(char*) 0;
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
            string msg = format("Unknown field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) double gas_state_get_scalar_field(int gs_i, char* field_name)
{
    try {
        GasState gs = gas_states[gs_i];
        string name = to!string(field_name);
        double value = 0.0;
        switch (name) {
        case "rho":
            value = gs.rho;
            break;
        case "p":
            value = gs.p;
            break;
        case "T":
            value = gs.T;
            break;
        case "u":
            value = gs.u;
            break;
        default:
            string msg = format("Unknown field name: %s", name);
            throw new Exception(msg);
        }
        return value;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return 0.0;
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
        default:
            string msg = format("Unknown field name: %s", name);
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
        default:
            string msg = format("Unknown field name: %s", name);
            throw new Exception(msg);
        }
        return 0;
    } catch (Exception e) {
        writeln("Exception message: ", e.msg);
        return -1;
    }
}

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
