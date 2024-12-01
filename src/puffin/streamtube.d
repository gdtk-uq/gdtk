// streamtube.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module streamtube;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.format;
import std.range;
import std.math;
import std.algorithm;

import nm.schedule;
import util.json_helper;
import geom;
import gas;
import kinetics;
import config;
import flow;
import cell;
import face;
import flux;


enum BCCode {wall=0, exchange=1}; // To decide what to do at the lower and upper boundary.


class StreamTube {
public:
    int indx;
    GasModel gmodel;
    ThermochemicalReactor thermochemUpdate;
    CQIndex cqi;
    GasState gs;
    Vector3 vel;
    int ncells;
    bool axiFlag;
    bool uSSTFlag;
    double[] dUdt_usst;
    double cfl;
    FluxCalcCode flux_calc;
    int x_order;
    double compression_tol;
    double shear_tol;
    //
    Schedule!double y_lower, y_upper, add_rho, add_xmom, add_ymom, add_tE;
    Schedule!int bc_lower, bc_upper, active;
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
        cqi = CQIndex(gmodel.n_species, gmodel.n_modes);
        axiFlag = Config.axisymmetric;
        uSSTFlag = Config.add_user_supplied_source_terms;
        dUdt_usst.length = cqi.n;
        foreach (ref du; dUdt_usst) du = 0.0;
        cfl = Config.cfl;
        x_order = Config.x_order;
        flux_calc = Config.flux_calc;
        compression_tol = Config.compression_tol;
        shear_tol = Config.shear_tol;
        //
        auto inflowData = configData[format("inflow_%d", indx)];
        double p = getJSONdouble(inflowData, "p", 100.0e3);
        double T = getJSONdouble(inflowData, "T", 300.0);
        double[] default_massf = [1.0, ];
        foreach (i; 1 .. gmodel.n_species) { default_massf ~= 0.0; }
        double[] massf = getJSONdoublearray(inflowData, "massf", default_massf);
        gs = GasState(gmodel);
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
        double[] xs, y0s, y1s, add_r, add_xm, add_ym, add_tEE;
        int[] bc0s, bc1s, act;
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
            if (items.length >= 10) {
                xs ~= to!double(items[0]);
                y0s ~= to!double(items[1]);
                y1s ~= to!double(items[2]);
                bc0s ~= to!int(items[3]);
                bc1s ~= to!int(items[4]);
                act ~= to!int(items[5]);
                add_r ~= to!double(items[6]);
                add_xm ~= to!double(items[7]);
                add_ym ~= to!double(items[8]);
                add_tEE ~= to!double(items[9]);
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
        active = new Schedule!int(xs, act);
        add_rho = new Schedule!double(xs, add_r);
        add_xmom = new Schedule!double(xs, add_xm);
        add_ymom = new Schedule!double(xs, add_ym);
        add_tE = new Schedule!double(xs, add_tEE);
        //
        // Scratch space
        fsL = FlowState2D(gmodel);
        fsR = FlowState2D(gmodel);
    } // end constructor

    override string toString()
    {
        string repr = "StreamTube(";
        repr ~= format("indx=%d", indx);
        repr ~= format(", gmodel=%s, gs=%s, vel=%s", gmodel, gs, vel);
        repr ~= format(", cqi=%s", cqi);
        repr ~= format(", axiFlag=%s", axiFlag);
        repr ~= format(", flux_calc=%s", to!string(flux_calc));
        repr ~= format(", compression_tol=%g, shear_tol=%g", compression_tol, shear_tol);
        repr ~= format(", x_order=%d", x_order);
        repr ~= format(", ncells=%d", ncells);
        repr ~= format(", y_lower=%s, y_upper=%s", y_lower, y_upper);
        repr ~= format(", bc_lower=%s, bc_upper=%s", bc_lower, bc_upper);
        repr ~= format(", active=%s", active);
        repr ~= format(", add_rho=%s", add_rho);
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
            flowstates_west ~= FlowState2D(gmodel);
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
            auto fs = &(flowstates_west[j]);
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

    void write_flow_data(bool write_header, bool write_data)
    {
        int nsp = gmodel.n_species;
        int nmodes = gmodel.n_modes;
        File fp;
        string fileName = format("%s/flow-%d.data", Config.job_name, indx);
        if (write_header) {
            fp = File(fileName, "w");
            fp.write("# x  y  velx vely  M  rho  p  T  u  a  shock");
            foreach (i; 0 .. nsp) { fp.write(format(" massf-%d", i)); }
            foreach (i; 0 .. nmodes) { fp.write(format(" T_modes-%d u_modes-%d", i, i)); }
            fp.write("\n");
        } else {
            fp = File(fileName, "a");
        }
        if (write_data) {
            foreach (j; 0 .. ncells) {
                auto face = ifaces_west[j];
                auto fs = flowstates_west[j];
                GasState* g = &(fs.gas);
                double Vx = fs.vel.x;
                double Vy = fs.vel.y;
                double M = sqrt(Vx*Vx+Vy*Vy)/g.a;
                double shock = (cells[j].shock_flag) ? 1.0 : 0.0;
                fp.write(format("%e %e %e %e %e", face.pos.x, face.pos.y, Vx, Vy, M));
                fp.write(format(" %e %e %e %e %e %f", g.rho, g.p, g.T, g.u, g.a, shock));
                foreach (i; 0 .. nsp) { fp.write(format(" %e", g.massf[i])); }
                foreach (i; 0 .. nmodes) { fp.write(format(" %e %e", g.T_modes[i], g.u_modes[i])); }
                fp.write("\n");
            }
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
    void mark_shock_cells()
    {
        foreach (c; cells) { c.shock_flag = false; }
        foreach (c; ghost_cells_left) { c.shock_flag = false; }
        foreach (c; ghost_cells_right) { c.shock_flag = false; }
        foreach (f; jfaces) {
            if (f.is_shock(compression_tol, shear_tol)) {
                f.left_cells[0].shock_flag = true;
                f.right_cells[0].shock_flag = true;
            }
        }
        return;
    }

    @nogc
    void predictor_step(double dt, double xmid)
    {
        foreach (c; cells) { c.estimate_local_dt(cfl); }
        foreach (j; 0 .. ncells) {
            ifaces_west[j].simple_flux(flowstates_west[j], gmodel, cqi);
            ifaces_east[j].simple_flux(cells[j].fs, gmodel, cqi);
        }
        foreach (j; 0 .. ncells+1) {
            jfaces[j].calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi);
        }
        if (uSSTFlag) {
            // Ingo's user-supplied source terms.
            dUdt_usst[cqi.mass] = add_rho.get_value(xmid);
            dUdt_usst[cqi.xMom] = add_xmom.get_value(xmid);
            dUdt_usst[cqi.yMom] = add_ymom.get_value(xmid);
            dUdt_usst[cqi.totEnergy] = add_tE.get_value(xmid);
        } else {
            dUdt_usst[cqi.mass] = 0.0;
            dUdt_usst[cqi.xMom] = 0.0;
            dUdt_usst[cqi.yMom] = 0.0;
            dUdt_usst[cqi.totEnergy] = 0.0;
        }
        foreach (c; cells) {
            c.update_conserved_for_stage(1, dt, axiFlag, dUdt_usst, uSSTFlag, gmodel);
        }
        return;
    } // end predictor_step()

    @nogc
    void corrector_step(double dt, double xmid)
    {
        foreach (j; 0 .. ncells) {
            ifaces_east[j].simple_flux(cells[j].fs, gmodel, cqi);
        }
        foreach (j; 0 .. ncells+1) {
            jfaces[j].calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi);
        }
        if (uSSTFlag) {
            // Ingo's user-supplied source terms.
            dUdt_usst[cqi.mass] = add_rho.get_value(xmid);
            dUdt_usst[cqi.xMom] = add_xmom.get_value(xmid);
            dUdt_usst[cqi.yMom] = add_ymom.get_value(xmid);
            dUdt_usst[cqi.totEnergy] = add_tE.get_value(xmid);
        } else {
            dUdt_usst[cqi.mass] = 0.0;
            dUdt_usst[cqi.xMom] = 0.0;
            dUdt_usst[cqi.yMom] = 0.0;
            dUdt_usst[cqi.totEnergy] = 0.0;
        }
        foreach (c; cells) {
            c.update_conserved_for_stage(2, dt, axiFlag, dUdt_usst, uSSTFlag, gmodel);
        }
        return;
    } // end corrector_step()

    @nogc
    void thermochemical_increment(double dt)
    {
        foreach (c; cells) { c.thermochemical_increment(dt, gmodel, thermochemUpdate); }
    }

    @nogc
    void transfer_conserved_quantities(size_t from, size_t dest)
    {
        foreach (c; cells) { c.U[dest][] = c.U[from][]; }
        return;
    }
} // end class StreamTube
