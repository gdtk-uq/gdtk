// fluidblock.d -- Part of the Lorikeet transient-flow calculator.
//
// PA Jacobs
// 2022-12-12: Adapt from the Puffin and Chicken codes.
//
module fluidblock;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.format;
import std.range;
import std.math;
import std.algorithm;
import core.stdc.math: HUGE_VAL;

import nm.schedule;
import json_helper;
import geom;
import gas;
import kinetics;
import config;
import flow;
import cell;
import face;
import flux;

enum BCCode {wall_with_slip=0, exchange=1, inflow=2, outflow=3};

struct BC {
    int code = BCCode.wall_with_slip;
    FlowState2D* fs;
}

class FluidBlock {
public:
    int indx;
    int i, j;
    bool active;
    GasModel gmodel;
    CQIndex cqi;
    int nic, njc;
    bool axiFlag;
    double cfl;
    FluxCalcCode flux_calc;
    int x_order;
    double compression_tol;
    double shear_tol;
    BC bc_west, bc_east, bc_south, bc_north;
    //
    Vector3[] vertices;
    Face2D[] ifaces;
    Face2D[] jfaces;
    Cell2D[] cells;
    // Storage of the halo of ghost flow states around the active cells.
    Cell2D[] ghost_cells;
    int n_vertices, n_cells, n_ifaces, n_jfaces, n_ghost_cells;
    //
    // Scratch space for taking a gas-dynamic step.
    FlowState2D fsL, fsR;


    this(int indx, JSONValue configData)
    // Configure the fluid block from the blob of JSON data associated with it.
    {
        this.indx = indx;
        nic = getJSONint(configData, format("ncells_%d", indx), 0);
        //
        gmodel = init_gas_model(Config.gas_model_file);
        cqi = CQIndex(gmodel.n_species, gmodel.n_modes);
        axiFlag = Config.axisymmetric;
        cfl = Config.cfl;
        x_order = Config.x_order;
        flux_calc = Config.flux_calc;
        compression_tol = Config.compression_tol;
        shear_tol = Config.shear_tol;
        //
        /+
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
        +/
        // Scratch space
        fsL = FlowState2D(gmodel);
        fsR = FlowState2D(gmodel);
    } // end constructor

    override string toString()
    {
        string repr = "FluidBlock(";
        repr ~= format("i=%d, j=%d", i, j);
        repr ~= format(", gmodel=%s", gmodel);
        repr ~= format(", cqi=%s", cqi);
        repr ~= format(", axiFlag=%s", axiFlag);
        repr ~= format(", flux_calc=%s", to!string(flux_calc));
        repr ~= format(", compression_tol=%g, shear_tol=%g", compression_tol, shear_tol);
        repr ~= format(", x_order=%d", x_order);
        repr ~= format(", n_cells=%d", n_cells);
        repr ~= format(", bc_west=%s, bc_east=%s", bc_west, bc_east);
        repr ~= format(", bc_south=%s, bc_north=%s", bc_south, bc_north);
        repr ~= ")";
        return repr;
    }

    @nogc
    int cell_index(int i, int j)
    {
        return j*nic + i;
    }

    @nogc
    int vertex_index(int i, int j)
    {
        return j*(nic+1) + i;
    }

    @nogc
    int iface_index(int i, int j)
    {
        return i*njc + j;
    }

    @nogc
    int jface_index(int i, int j)
    {
        return j*nic + j;
    }

    void set_up_data_storage()
    // Set up the storage and make connections to the vertices.
    {
        /+
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
        +/
        return;
    } // end set_up_data_storage()

    void read_grid_data()
    {
        return;
    }

    @nogc
    void set_up_geometry()
    {
        foreach (f; ifaces) { f.compute_geometry(axiFlag); }
        foreach (f; jfaces) { f.compute_geometry(axiFlag); }
        foreach (c; cells) { c.compute_geometry(axiFlag); }
    }

    void read_flow_data(int tindx)
    {
        return;
    }

    void write_flow_data(int tindx)
    {
        /+
        bool write_header = false;
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
        foreach (j; 0 .. n_cells) {
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
        fp.close();
        +/
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
        double dt = HUGE_VAL;
        foreach (c; cells) { dt = fmin(dt, c.estimate_local_dt(cfl)); }
        return dt;
    }

    @nogc
    void mark_shock_cells()
    {
        foreach (c; cells) { c.shock_flag = false; }
        foreach (c; ghost_cells) { c.shock_flag = false; }
        foreach (f; jfaces) {
            if (f.is_shock(compression_tol, shear_tol)) {
                f.left_cells[0].shock_flag = true;
                f.right_cells[0].shock_flag = true;
            }
        }
        return;
    }

    @nogc
    void predictor_step(double dt)
    {
        foreach (f; ifaces) {
            f.calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi);
        }
        foreach (f; jfaces) {
            f.calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi);
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
        foreach (f; ifaces) {
            f.calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi);
        }
        foreach (f; jfaces) {
            f.calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi);
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
} // end class FluidBlock
