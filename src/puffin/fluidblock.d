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

int BC_code_from_name(string name)
{
    if (name == "wall_with_slip") return BCCode.wall_with_slip;
    if (name == "exchange") return BCCode.exchange;
    if (name == "inflow") return BCCode.inflow;
    if (name == "outflow") return BCCode.outflow;
    return BCCode.wall_with_slip;
}

struct BC {
    int code = BCCode.wall_with_slip;
    FlowState2D* fs;

    string toString() {
        return format("BC(code=%d, fs=%s)", code, ((fs) ? to!string(*fs) : "null"));
    }
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
        gmodel = init_gas_model(Config.gas_model_file);
        cqi = CQIndex(gmodel.n_species, gmodel.n_modes);
        axiFlag = Config.axisymmetric;
        cfl = Config.cfl;
        x_order = Config.x_order;
        flux_calc = Config.flux_calc;
        compression_tol = Config.compression_tol;
        shear_tol = Config.shear_tol;
        //
        // writeln("DEBUG configData.type=", configData.type, " configData=", configData);
        nic = getJSONint(configData, "nic", 0);
        njc = getJSONint(configData, "njc", 0);
        active = getJSONbool(configData, "active", true);
        JSONValue jsonBCs = configData["bcs"];
        JSONValue jsonBC = jsonBCs["iminus"];
        // writeln("DEBUG iminus jsonBC=", jsonBC);
        bc_west.code = BC_code_from_name(jsonBC["tag"].str);
        if (bc_west.code == BCCode.inflow) {
            // [TODO] Push this into a function so that we can reuse it below. PJ 2022-12-12
            auto fs = new FlowState2D(gmodel);
            JSONValue jsonFlow = jsonBC["flow_state"];
            double p = getJSONdouble(jsonFlow, "p", 100.0e3);
            double T = getJSONdouble(jsonFlow, "T", 300.0);
            double[] default_massf = [1.0, ];
            foreach (i; 1 .. gmodel.n_species) { default_massf ~= 0.0; }
            double[] massf = getJSONdoublearray(jsonFlow, "massf", default_massf);
            fs.gas.p = p; fs.gas.T = T; fs.gas.massf[] = massf[];
            gmodel.update_thermo_from_pT(fs.gas);
            gmodel.update_sound_speed(fs.gas);
            //
            double velx = getJSONdouble(jsonFlow, "velx", 0.0);
            double vely = getJSONdouble(jsonFlow, "vely", 0.0);
            fs.vel.set(velx, vely);
            bc_west.fs = fs;
        }
        jsonBC = jsonBCs["iplus"];
        bc_east.code = BC_code_from_name(jsonBC["tag"].str);
        jsonBC = jsonBCs["jminus"];
        bc_south.code = BC_code_from_name(jsonBC["tag"].str);
        jsonBC = jsonBCs["jplus"];
        bc_north.code = BC_code_from_name(jsonBC["tag"].str);
        //
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
        // Allocate the actual data space.
        n_vertices = (nic+1)*(njc+1);
        vertices.length = n_vertices;
        //
        n_ifaces = (nic+1)*njc;
        foreach (indx; 0 .. n_ifaces) { ifaces ~= new Face2D(gmodel, cqi); }
        //
        n_jfaces = nic*(njc+1);
        foreach (indx; 0 .. n_jfaces) { jfaces ~= new Face2D(gmodel, cqi); }
        //
        n_cells = nic*njc;
        foreach (indx; 0 .. n_cells) { cells ~= new Cell2D(gmodel, cqi); }
        //
        // Connect cells to faces and vertices.
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                auto c = cells[cell_index(i,j)];
                c.p00 = &(vertices[vertex_index(i,j)]);
                c.p10 = &(vertices[vertex_index(i+1,j)]);
                c.p11 = &(vertices[vertex_index(i+1,j+1)]);
                c.p01 = &(vertices[vertex_index(i,j+1)]);
                //
                c.faceW = ifaces[iface_index(i,j)];
                c.faceE = ifaces[iface_index(i+1,j)];
                c.faceS = jfaces[jface_index(i,j)];
                c.faceN = jfaces[jface_index(i,j+1)];
            }
        }
        // Now, add ghost cells and make connection faces to their neighbour cells.
        foreach (j; 0 .. njc) {
            auto ghost_cell_left_1 = new Cell2D(gmodel, cqi);
            auto ghost_cell_left_0 = new Cell2D(gmodel, cqi);
            auto ghost_cell_right_0 = new Cell2D(gmodel, cqi);
            auto ghost_cell_right_1 = new Cell2D(gmodel, cqi);
            foreach (i; 0 ..nic) {
                Face2D f = ifaces[iface_index(i,j)];
                // Want unit normal to right and pointing in i-direction.
                f.p0 = &(vertices[vertex_index(i,j)]);
                f.p1 = &(vertices[vertex_index(i,j+1)]);
                if (i == 0) {
                    // Left-most face in block.
                    f.left_cells[1] = ghost_cell_left_1;
                    f.left_cells[0] = ghost_cell_left_0;
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = cells[cell_index(i+1,j)];
                } else if (i == 1 && i == nic-1) {
                    // Middle face when there is only 2 cells across block.
                    f.left_cells[1] = ghost_cell_left_0;
                    f.left_cells[0] = cells[cell_index(i-1,j)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = ghost_cell_right_0;
                } else if (i == 1) {
                    // First face in from left in a large block.
                    f.left_cells[1] = ghost_cell_left_0;
                    f.left_cells[0] = cells[cell_index(i-1,j)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = cells[cell_index(i+1,j)];
                } else if (i == nic-1) {
                    // First face in from right in a large block.
                    f.left_cells[1] = cells[cell_index(i-2,j)];
                    f.left_cells[0] = cells[cell_index(i-1,j)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = ghost_cell_right_0;
                } else if (i == nic) {
                    // Last face on block edge.
                    f.left_cells[1] = cells[cell_index(i-2,j)];
                    f.left_cells[0] = cells[cell_index(i-1,j)];
                    f.right_cells[0] = ghost_cell_right_0;
                    f.right_cells[1] = ghost_cell_right_1;
                } else {
                    // Interior face in a large block.
                    f.left_cells[1] = cells[cell_index(i-2,j)];
                    f.left_cells[0] = cells[cell_index(i-1,j)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = cells[cell_index(i+1,j)];
                }
            } // end for i
        } // end for j
        foreach (i; 0 .. nic) {
            auto ghost_cell_left_1 = new Cell2D(gmodel, cqi);
            auto ghost_cell_left_0 = new Cell2D(gmodel, cqi);
            auto ghost_cell_right_0 = new Cell2D(gmodel, cqi);
            auto ghost_cell_right_1 = new Cell2D(gmodel, cqi);
            foreach (j; 0 ..njc+1) {
                Face2D f = jfaces[jface_index(i,j)];
                // Want unit normal to right and pointing in j-direction.
                f.p0 = &(vertices[vertex_index(i+1,j)]);
                f.p1 = &(vertices[vertex_index(i,j)]);
                if (j == 0) {
                    // Left-most face in block.
                    f.left_cells[1] = ghost_cell_left_1;
                    f.left_cells[0] = ghost_cell_left_0;
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = cells[cell_index(i,j+1)];
                } else if (j == 1 && j == njc-1) {
                    // Middle face when there is only 2 cells across block.
                    f.left_cells[1] = ghost_cell_left_0;
                    f.left_cells[0] = cells[cell_index(i,j-1)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = ghost_cell_right_0;
                } else if (j == 1) {
                    // First face in from left in a large block.
                    f.left_cells[1] = ghost_cell_left_0;
                    f.left_cells[0] = cells[cell_index(i,j-1)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = cells[cell_index(i,j+1)];
                } else if (j == njc-1) {
                    // First face in from right in a large block.
                    f.left_cells[1] = cells[cell_index(i,j-2)];
                    f.left_cells[0] = cells[cell_index(i,j-1)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = ghost_cell_right_0;
                } else if (j == njc) {
                    // Last face on block edge.
                    f.left_cells[1] = cells[cell_index(i,j-2)];
                    f.left_cells[0] = cells[cell_index(i,j-1)];
                    f.right_cells[0] = ghost_cell_right_0;
                    f.right_cells[1] = ghost_cell_right_1;
                } else {
                    // Interior face in a large block.
                    f.left_cells[1] = cells[cell_index(i,j-2)];
                    f.left_cells[0] = cells[cell_index(i,j-1)];
                    f.right_cells[0] = cells[cell_index(i,j)];
                    f.right_cells[1] = cells[cell_index(i,j+1)];
                }
            } // end for j
        } // end for i
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
