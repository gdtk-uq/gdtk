// cwrap_gasmodule.d
//
// This particular module will be compiled into a loadable library and
// provides a C-level API for the GasModel and GasState classes
// that may be called from any language with a C-foreign-function-interface.
//
// PJ 2019-07-24: just enough to try integration with the gas makefile.
//

import std.stdio;
import std.format;
import std.conv;
import core.runtime;

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
