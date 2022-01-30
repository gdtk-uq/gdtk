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
    CQIndex cqi;
    GasState gs;
    Vector3 vel;
    int ncells;
    bool axiFlag = false;
    double cfl = 0.5;
    FluxCalcCode flux_calc = FluxCalcCode.hanel;
    int x_order = 2;
    //
    Schedule!double y_lower, y_upper;
    Schedule!int bc_lower, bc_upper;
    //
    Vector3[] vertices_west;
    Vector3[] vertices_east;
    Face2D[] ifaces_west;
    Face2D[] ifaces_east;
    Face2D[] jfaces;
    //
    Cell2D[] cells;
    Cell2D[] ghost_cells_left;
    Cell2D[] ghost_cells_right;
    //
    FlowState2D[] flowstates_west;
    //
    // Scratch space for taking a gas-dynamic step.
    FlowState2D fsL, fsR;


    this(int indx, JSONValue configData)
    {
        this.indx = indx;
        ncells = getJSONint(configData, format("ncells_%d", indx), 0);
        //
        gmodel = init_gas_model(Config.gas_model_file);
        cqi = new CQIndex(gmodel.n_species, gmodel.n_modes);
        axiFlag = Config.axisymmetric;
        cfl = Config.cfl;
        flux_calc = Config.flux_calc;
        x_order = Config.x_order;
        //
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
        //
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
        //
        // Scratch space
        fsL = new FlowState2D(gmodel);
        fsR = new FlowState2D(gmodel);
    } // end constructor

    override string toString()
    {
        string repr = "StreamTube(";
        repr ~= format("indx=%d", indx);
        repr ~= format(", gmodel=%s, gs=%s, vel=%s", gmodel, gs, vel);
        repr ~= format(", cqi=%s", cqi);
        repr ~= format(", axiFlag=%s", axiFlag);
        repr ~= format(", flux_calc=%s", to!string(flux_calc));
        repr ~= format(", x_order=%d", x_order);
        repr ~= format(", ncells=%d", ncells);
        repr ~= format(", y_lower=%s, y_upper=%s", y_lower, y_upper);
        repr ~= format(", bc_lower=%s, bc_upper=%s", bc_lower, bc_upper);
        repr ~= ")";
        return repr;
    }

    void set_up_data_storage()
    // Set up the storage and make connections to the vertices.
    {
        foreach (j; 0 .. ncells+1) {
            vertices_west ~= Vector3();
            vertices_east ~= Vector3();
        }
        foreach (j; 0 .. ncells) {
            ifaces_west ~= new Face2D(gmodel, cqi);
            ifaces_east ~= new Face2D(gmodel, cqi);
        }
        foreach (j; 0 .. ncells+1) {
            jfaces ~= new Face2D(gmodel, cqi);
        }
        foreach (j; 0 .. ncells) {
            cells ~= new Cell2D(gmodel, cqi);
            flowstates_west ~= new FlowState2D(gmodel);
        }
        foreach (j; 0 .. 2) {
            ghost_cells_left ~= new Cell2D(gmodel, cqi);
            ghost_cells_right ~= new Cell2D(gmodel, cqi);
        }
        foreach (j; 0 .. ncells) {
            ifaces_west[j].p0 = &(vertices_west[j]);
            ifaces_west[j].p1 = &(vertices_west[j+1]);
            //
            ifaces_east[j].p0 = &(vertices_east[j]);
            ifaces_east[j].p1 = &(vertices_east[j+1]);
        }
        foreach (j; 0 .. ncells+1) {
            auto f = jfaces[j];
            f.p0 = &(vertices_east[j]);
            f.p1 = &(vertices_west[j]);
            if (j == 0) {
                f.left_cells[1] = ghost_cells_left[1];
                f.left_cells[0] = ghost_cells_left[0];
            } else if (j == 1) {
                f.left_cells[1] = ghost_cells_left[0];
                f.left_cells[0] = cells[j-1];
            } else {
                f.left_cells[1] = cells[j-2];
                f.left_cells[0] = cells[j-1];
            }
            if (j == ncells-1) {
                f.right_cells[0] = cells[j];
                f.right_cells[1] = ghost_cells_right[0];
            } else if (j ==  ncells) {
                f.right_cells[0] = ghost_cells_right[0];
                f.right_cells[1] = ghost_cells_right[1];
            } else {
                f.right_cells[0] = cells[j];
                f.right_cells[1] = cells[j+1];
            }
        }
        foreach (j; 0 .. ncells) {
            auto c = cells[j];
            c.p00 = &(vertices_west[j]);
            c.p10 = &(vertices_east[j]);
            c.p11 = &(vertices_east[j+1]);
            c.p01 = &(vertices_west[j+1]);
            //
            c.faceN = jfaces[j+1];
            c.faceE = ifaces_east[j];
            c.faceS = jfaces[j];
            c.faceW = ifaces_west[j];
        }
        return;
    } // end set_up_data_storage()

    @nogc
    void set_up_inflow_boundary()
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
            ifaces_west[j].compute_geometry(axiFlag);
        }
        // Inflow states
        foreach (j; 0 .. ncells) {
            auto fs = flowstates_west[j];
            fs.gas.copy_values_from(gs);
            fs.vel.set(vel);
        }
        return;
    } // end set_up_inflow_boundary()

    @nogc
    void set_up_slice(double x)
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
            ifaces_east[j].compute_geometry(axiFlag);
            cells[j].compute_geometry(axiFlag);
        }
        foreach (j; 0 .. ncells+1) {
            jfaces[j].compute_geometry(axiFlag);
        }
        // Finally, initialize the flow by propagating from the west.
        foreach (j; 0 .. ncells) {
            cells[j].fs.copy_values_from(flowstates_west[j]);
        }
        return;
    } // end set_up_slice()

    @nogc
    void shuffle_data_west()
    // At the end of relaxation process (and reaching steady-state),
    // we copy the geometry and flow states over to the west side
    // of the slice.  We need to do this before setting up a new slice.
    {
        foreach (j; 0 .. ncells+1) {
            vertices_west[j].set(vertices_east[j]);
        }
        foreach (j; 0 .. ncells) {
            ifaces_west[j].compute_geometry(axiFlag);
        }
        foreach (j; 0 .. ncells) {
            flowstates_west[j].copy_values_from(cells[j].fs);
        }
        return;
    } // end shuffle_data_west()

    void write_flow_data(bool write_header)
    {
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        File fp;
        string fileName = format("%s/flow-%d.data", Config.job_name, indx);
        if (write_header) {
            fp = File(fileName, "w");
            fp.write("# x  y  velx vely  M  rho  p  T  u  a");
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
            double Vx = fs.vel.x;
            double Vy = fs.vel.y;
            double M = sqrt(Vx*Vx+Vy*Vy)/g.a;
            fp.write(format("%e %e %e %e %e", face.pos.x, face.pos.y, Vx, Vy, M));
            fp.write(format(" %e %e %e %e %e", g.rho, g.p, g.T, g.u, g.a));
            foreach (i; 0 .. nsp) { fp.write(format(" %e", g.massf[i])); }
            foreach (i; 0 .. nmodes) { fp.write(format(" %e %e", g.T_modes[i], g.u_modes[i])); }
            fp.write("\n");
        }
        fp.close();
        return;
    } // end  write_flow_data()

    @nogc
    void encode_conserved(size_t ftl)
    {
        foreach (c; cells) { c.encode_conserved(ftl, gmodel); }
        return;
    }

    @nogc
    double estimate_allowable_dt()
    {
        double dt = cells[0].estimate_local_dt(cfl);
        foreach (j; 1 .. ncells) { dt = fmin(dt, cells[j].estimate_local_dt(cfl)); }
        return dt;
    }

    @nogc
    void predictor_step(double dt)
    {
        foreach (c; cells) { c.estimate_local_dt(cfl); }
        foreach (j; 0 .. ncells) {
            ifaces_west[j].simple_flux(flowstates_west[j], gmodel);
            ifaces_east[j].simple_flux(cells[j].fs, gmodel);
        }
        foreach (j; 0 .. ncells+1) {
            auto f = jfaces[j];
            if (x_order == 1) {
                // Low-order reconstruction is just copy the nearest cell-centre flowstate.
                fsL.copy_values_from(f.left_cells[0].fs);
                fsR.copy_values_from(f.right_cells[0].fs);
            } else {
                // High-order reconstruction from the left_cells and right_cells stencil.
                f.interp_l2r2(fsL, fsR, gmodel, false);
            }
            f.calculate_flux(fsL, fsR, gmodel, flux_calc);
        }
        foreach (c; cells) {
            c.eval_dUdt(0, axiFlag);
            c.U[1][] = c.U[0][] + dt*c.dUdt[0][];
            c.decode_conserved(1, gmodel);
        }
        return;
    } // end predictor_step()

    @nogc
    void corrector_step(double dt)
    {
        foreach (j; 0 .. ncells) {
            ifaces_east[j].simple_flux(cells[j].fs, gmodel);
        }
        foreach (j; 0 .. ncells+1) {
            auto f = jfaces[j];
            if (x_order == 1) {
                // Low-order reconstruction is just copy the nearest cell-centre flowstate.
                fsL.copy_values_from(f.left_cells[0].fs);
                fsR.copy_values_from(f.right_cells[0].fs);
            } else {
                // High-order reconstruction from the left_cells and right_cells stencil.
                f.interp_l2r2(fsL, fsR, gmodel, false);
            }
            f.calculate_flux(fsL, fsR, gmodel, flux_calc);
        }
        foreach (c; cells) {
            c.eval_dUdt(1, axiFlag);
            c.U[2][] = c.U[0][] + 0.5*dt*(c.dUdt[0][] + c.dUdt[1][]);
            c.decode_conserved(2, gmodel);
        }
        return;
    } // end corrector_step()

    @nogc
    void transfer_conserved_quantities(size_t from, size_t dest)
    {
        foreach (c; cells) { c.U[dest][] = c.U[from][]; }
        return;
    }
} // end class StreamTube
