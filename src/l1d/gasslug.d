// gasslug.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
//
// The GasSlug is the principal dynamic component of a simulation.
//
// PA Jacobs
// 2020-04-08
//
module gasslug;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.format;
import std.range;

import json_helper;
import geom;
import gas;
import gasflow;
import config;
import lcell;
import endcondition;
import simcore; // has the core data arrays
import misc;

class GasSlug {
public:
    size_t indx;
    string label;
    size_t gmodel_id;
    GasModel gmodel;
    size_t ncells;
    int viscous_effects;
    bool adiabatic;
    int ecL_id;
    int ecR_id;
    EndCondition ecL;
    EndCondition ecR;
    size_t nhcells;
    size_t[] hcells;

    LFace[] faces;
    LCell[] cells;

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct slug[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        label = getJSONstring(jsonData, "label", "");
        gmodel_id = getJSONint(jsonData, "gmodel_id", 0);
        gmodel = gmodels[gmodel_id];
        ncells = getJSONint(jsonData, "ncells", 0);
        viscous_effects = getJSONint(jsonData, "viscous_effects", 0);
        adiabatic = getJSONbool(jsonData, "adiabatic", false);
        ecL_id = getJSONint(jsonData, "ecL_id", -1);
        ecR_id = getJSONint(jsonData, "ecR_id", -1);
        nhcells = getJSONint(jsonData, "nhcells", 1);
        hcells = to!(size_t[])(getJSONintarray(jsonData, "hcells", [0,]));
        if (L1dConfig.verbosity_level >= 1) {
            writeln("GasSlug[", indx, "]:");
            writefln("  label= \"%s\"", label);
            writeln("  gmodel_id= ", gmodel_id);
            writeln("  ncells= ", ncells);
            writeln("  viscous_effects= ", viscous_effects);
            writeln("  adiabatic= ", adiabatic);
            writeln("  ecL_id= ", ecL_id);
            writeln("  ecR_id= ", ecR_id);
            writeln("  hcells= ", hcells);
        }
        //
        foreach (i; 0 .. ncells+1) { faces ~= new LFace(); }
        foreach (i; 0 .. ncells) { cells ~= new LCell(gmodel); }
    } // end constructor
    
    void read_face_data(File fp, int tindx=0)
    {
        skip_to_data_at_tindx(fp, tindx);
        foreach (i; 0 .. ncells+1) {
            string text = fp.readln().chomp();
            text.formattedRead!"%e %e"(faces[i].x, faces[i].area);
        }
    } // end read_face_data

    void write_face_data(File fp, int tindx=0)
    {
        if (tindx == 0) { fp.writeln("#   x   area"); }
        fp.writeln(format("# tindx %d", tindx));
        foreach (i; 0 .. ncells+1) {
            fp.writeln(format("%e %e", faces[i].x, faces[i].area));
        }
        fp.writeln("# end");
    } // end write_face_data()

    void read_cell_data(File fp, int tindx=0)
    {
        skip_to_data_at_tindx(fp, tindx);
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        foreach (j; 0 .. ncells) {
            LCell c = cells[j];
            string text = fp.readln().chomp();
            string[] items = text.split();
            int k = 0;
            c.xmid = to!double(items[k]); k++;
            c.volume = to!double(items[k]); k++;
            c.vel = to!double(items[k]); k++;
            c.L_bar = to!double(items[k]); k++;
            c.gas.rho = to!double(items[k]); k++;
            c.gas.p = to!double(items[k]); k++;
            c.gas.T = to!double(items[k]); k++;
            c.gas.u = to!double(items[k]); k++;
            c.gas.a = to!double(items[k]); k++;
            c.shear_stress = to!double(items[k]); k++;
            c.heat_flux = to!double(items[k]); k++;
            foreach (i; 0 .. nsp) {
                c.gas.massf[i] = to!double(items[k]); k++;
            }
            if (nsp > 1) { c.dt_chem = to!double(items[k]); k++; }
            foreach (i; 0 .. nmodes) {
                c.gas.T_modes[i] = to!double(items[k]); k++;
                c.gas.u_modes[i] = to!double(items[k]); k++;
            }
            if (nmodes > 0) { c.dt_therm = to!double(items[k]); k++; }
        }
    } // end read_cell_data()

    void write_cell_data(File fp, int tindx=0)
    {
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        if (tindx == 0) {
            fp.write("# xmid  volume  vel  L_bar  rho  p  T  u  a");
            fp.write("  shear_stress  heat_flux");
            foreach (i; 0 .. nsp) { fp.write(format("  massf[%d]", i)); }
            if (nsp > 1) { fp.write("  dt_chem"); }
            foreach (i; 0 .. nmodes) {
                fp.write(format("  T_modes[%d]  u_modes[%d]", i, i));
            }
            if (nmodes > 0) { fp.write("  dt_therm"); }
            fp.write("\n");
        }
        fp.writeln(format("# tindx %d", tindx));
        foreach (j; 0 .. ncells) {
            LCell c = cells[j];
            fp.write(format("%e %e %e %e", c.xmid, c.volume, c.vel, c.L_bar));
            fp.write(format(" %e %e %e %e", c.gas.rho, c.gas.p, c.gas.T, c.gas.u));
            fp.write(format(" %e %e %e", c.gas.a, c.shear_stress, c.heat_flux));
            foreach (i; 0 .. nsp) { fp.write(format(" %e", c.gas.massf[i])); }
            if (nsp > 1) { fp.write(format(" %e", c.dt_chem)); }
            foreach (i; 0 .. nmodes) {
                fp.write(format(" %e %e", c.gas.T_modes[i], c.gas.u_modes[i]));
            }
            if (nmodes > 0) { fp.write(format(" %e", c.dt_therm)); }
            fp.write("\n");
        }
        fp.writeln("# end");
    } // end write_cell_data()

    @nogc
    void compute_areas_and_volumes()
    {
        double[4] daKT;
        foreach (f; faces) {
            daKT = tube1.eval(f.x);
            f.area = daKT[1];
        }
        foreach (i; 0 .. cells.length) {
            double xL = faces[i].x;
            double xR = faces[i+1].x;
            LCell c = cells[i];
            c.volume = 0.5*(faces[i].area+faces[i+1].area)*(xR-xL);
            c.xmid = 0.5*(xR+xL);
            daKT = tube1.eval(c.xmid);
            c.diam = daKT[0]; c.K_over_L = daKT[2]; c.Twall = daKT[3];
        }
        return;
    } // end compute_areas_and_volumes()

    @nogc
    void encode_conserved()
    {
        foreach (c; cells) { c.encode_conserved(gmodel); }
        return;
    }

    @nogc
    void decode_conserved()
    {
        foreach (c; cells) { c.decode_conserved(gmodel); }
        return;
    }

    @nogc
    void source_terms()
    {
        // [TODO]
        return;
    }

    @nogc
    void apply_rivp()
    {
        // [TODO]
        // For the moment, use cell-centre values as the left and right states.
        // Need to consider the end conditions.
        return;
    } // end apply_rivp()

    @nogc
    void time_derivatives()
    {
        // [TODO]
        return;
    }

    @nogc
    void record_state()
    {
        // [TODO]
        return;
    }

    @nogc
    void restore_state()
    {
        // [TODO]
        return;
    }

    @nogc
    void predictor_step(double dt)
    {
        // [TODO]
        return;
    }

    @nogc
    void corrector_step(double dt)
    {
        // [TODO]
        return;
    }

    @nogc
    void chemical_increment(double dt)
    {
        // [TODO]
        return;
    }
    
    @nogc
    void check_cells()
    {
        // [TODO]
        return;
    }

    @nogc
    double compute_stable_time_step()
    {
        double dt_allowed;
        return dt_allowed;
    }

} // end class GasSlug
