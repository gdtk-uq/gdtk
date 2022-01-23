// stream.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module stream;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.format;
import std.range;
import std.math;
import std.algorithm;

import nm.schedule;
import json_helper;
import geom;
import gas;
import kinetics;
import config;
import flow;
import cell;


enum BCCode {wall=0, exchange=1}; // To decide what to do at the lower and upper boundary.


class StreamTube {
public:
    int indx;
    GasModel gmodel;
    GasState gs;
    Vector3 vel;
    int ncells;
    //
    Schedule!double y_lower, y_upper;
    Schedule!int bc_lower, bc_upper;
    //
    Vector3[] vertices_west;
    Vector3[] vertices_east;
    Face2D[] ifaces_west;
    Face2D[] ifaces_east;
    Face2D[] jfaces;
    Cell2D[] cells;
    //
    FlowState2D[] flowstates_west;

    this(int indx, JSONValue configData)
    {
        this.indx = indx;
        ncells = getJSONint(configData, format("ncells_%d", indx), 0);
        gmodel = init_gas_model(Config.gas_model_file);
        auto inflowData = configData[format("inflow_%d", indx)];
        double p = getJSONdouble(inflowData, "p", 100.0e3);
        double T = getJSONdouble(inflowData, "T", 300.0);
        double[] default_massf = [1.0, ];
        foreach (i; 1 .. gmodel.n_species) { default_massf ~= 0.0; }
        double[] massf = getJSONdoublearray(inflowData, "massf", default_massf);
        gs = new GasState(gmodel);
        gs.p = p; gs.T = T; gs.massf[] = massf[];
        gmodel.update_thermo_from_pT(gs);
        gmodel.update_sound_speed(gs);
        double velx = getJSONdouble(inflowData, "velx", 0.0);
        double vely = getJSONdouble(inflowData, "vely", 0.0);
        vel.set(velx, vely);
        //
        int n_bc;
        double dx_bc;
        double[] xs, y0s, y1s;
        int[] bc0s, bc1s;
        string fileName = format("%s/streamtube-%d.data", Config.job_name, indx);
        auto lines = File(fileName, "r").byLine();
        foreach (line; lines) {
            auto txt = line.strip();
            if (canFind(txt, "n_bc")) {
                n_bc = to!int(txt.split("=")[1]);
                continue;
            }
            if (canFind(txt, "dx_bc")) {
                dx_bc = to!double(txt.split("=")[1]);
                continue;
            }
            if (canFind(txt, "#")) continue;
            auto items = txt.split();
            if (items.length >= 5) {
                xs ~= to!double(items[0]);
                y0s ~= to!double(items[1]);
                y1s ~= to!double(items[2]);
                bc0s ~= to!int(items[3]);
                bc1s ~= to!int(items[4]);
            }
        }
        if (n_bc != xs.length) {
            writeln("WARNING: n_bc=", n_bc, " length=", xs.length);
        }
        if ((xs[1]-xs[0]-dx_bc)/dx_bc > 1.0e-9) {
            writeln("WARNING: dx=", xs[1]-xs[0], " dx_bc=", dx_bc);
        }
        y_lower = new Schedule!double(xs, y0s);
        y_upper = new Schedule!double(xs, y1s);
        bc_lower = new Schedule!int(xs, bc0s);
        bc_upper = new Schedule!int(xs, bc1s);
    } // end constructor


    override string toString()
    {
        string repr = "StreamTube(";
        repr ~= format("indx=%d", indx);
        repr ~= format(", gmodel=%s, gs=%s", gmodel, gs);
        repr ~= format(", ncells=%d", ncells);
        repr ~= format(", y_lower=%s, y_upper=%s", y_lower, y_upper);
        repr ~= format(", bc_lower=%s, bc_upper=%s", bc_lower, bc_upper);
        repr ~= ")";
        return repr;
    }


    int set_up_data_storage()
    // Set up the storage and make connections to the vertices.
    {
        foreach (j; 0 .. ncells+1) {
            vertices_west ~= Vector3();
            vertices_east ~= Vector3();
        }
        foreach (j; 0 .. ncells) {
            ifaces_west ~= new Face2D();
            ifaces_east ~= new Face2D();
        }
        foreach (j; 0 .. ncells+1) {
            jfaces ~= new Face2D();
        }
        foreach (j; 0 .. ncells) {
            cells ~= new Cell2D(gmodel);
            flowstates_west ~= new FlowState2D(gmodel);
        }
        foreach (j; 0 .. ncells) {
            ifaces_west[j].p0 = &(vertices_west[j]);
            ifaces_west[j].p1 = &(vertices_west[j+1]);
            ifaces_east[j].p0 = &(vertices_east[j]);
            ifaces_east[j].p1 = &(vertices_east[j+1]);
        }
        foreach (j; 0 .. ncells+1) {
            jfaces[j].p0 = &(vertices_east[j]);
            jfaces[j].p1 = &(vertices_west[j]);
        }
        foreach (j; 0 .. ncells) {
            auto c = cells[j];
            c.p00 = &(vertices_west[j]);
            c.p10 = &(vertices_east[j]);
            c.p11 = &(vertices_east[j+1]);
            c.p01 = &(vertices_west[j+1]);
        }
        return 0;
    } // end set_up_data_storage()


    int set_up_inflow_boundary()
    {
        // Compute the vertex positions.
        double x = 0.0;
        double y0 = y_lower.interpolate_value(x);
        double y1 = y_upper.interpolate_value(x);
        foreach (j; 0 .. ncells+1) {
            auto frac = to!double(j)/to!double(ncells);
            double y = (1.0-frac)*y0 + frac*y1;
            vertices_west[j].set(x, y);
        }
        // Geometry of faces.
        foreach (j; 0 .. ncells) {
            ifaces_west[j].compute_geometry();
        }
        // Inflow states
        foreach (j; 0 .. ncells) {
            auto fs = flowstates_west[j];
            fs.gas.copy_values_from(gs);
            fs.vel.set(vel);
        }
        return 0;
    } // end set_up_inflow_boundary()


    int set_up_slice(double x)
    {
        // Locations for the east vertices.
        double y0 = y_lower.interpolate_value(x);
        double y1 = y_upper.interpolate_value(x);
        foreach (j; 0 .. ncells+1) {
            auto frac = to!double(j)/to!double(ncells);
            double y = (1.0-frac)*y0 + frac*y1;
            vertices_east[j].set(x, y);
        }
        // Compute the face and cell geometric properties.
        foreach (j; 0 .. ncells) {
            ifaces_east[j].compute_geometry();
            cells[j].compute_geometry();
        }
        foreach (j; 0 .. ncells+1) {
            jfaces[j].compute_geometry();
        }
        // Finally, initialize the flow by propagating from the west.
        foreach (j; 0 .. ncells) {
            cells[j].fs.copy_values_from(flowstates_west[j]);
        }
        return 0;
    } // end set_up_slice()


    int shuffle_data_west()
    {
        foreach (j; 0 .. ncells+1) {
            vertices_west[j].set(vertices_east[j]);
        }
        foreach (j; 0 .. ncells) {
            ifaces_west[j].compute_geometry();
        }
        foreach (j; 0 .. ncells) {
            flowstates_west[j].copy_values_from(cells[j].fs);
        }
        return 0;
    }


    int write_flow_data(bool write_header)
    {
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        File fp;
        string fileName = format("%s/flow-%d.data", Config.job_name, indx);
        if (write_header) {
            fp = File(fileName, "w");
            fp.write("# x  y  velx vely rho  p  T  u  a");
            foreach (i; 0 .. nsp) { fp.write(format(" massf-%d", i)); }
            foreach (i; 0 .. nmodes) { fp.write(format(" T_modes-%d u_modes-%d", i, i)); }
            fp.write("\n");
        } else {
            fp = File(fileName, "a");
        }
        foreach (j; 0 .. ncells) {
            auto face = ifaces_west[j];
            auto fs = flowstates_west[j];
            auto g = fs.gas;
            fp.write(format("%e %e %e %e", face.pos.x, face.pos.y, fs.vel.x, fs.vel.y));
            fp.write(format(" %e %e %e %e %e", g.rho, g.p, g.T, g.u, g.a));
            foreach (i; 0 .. nsp) { fp.write(format(" %e", g.massf[i])); }
            foreach (i; 0 .. nmodes) { fp.write(format(" %e %e", g.T_modes[i], g.u_modes[i])); }
            fp.write("\n");
        }
        fp.close();
        return 0;
    }
} // end class StreamTube
