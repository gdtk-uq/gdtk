// fluidblock.d -- Part of the Lorikeet transient-flow calculator.
//
// PA Jacobs
// 2022-12-12: Adapt from the Puffin and Chicken codes.
//
module fluidblock;

import std.conv;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.format;
import std.range;
import std.math;
import std.algorithm;
import core.stdc.math: HUGE_VAL;
import gzip;

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

enum BCCode {wall_with_slip=0, exchange=1, inflow=2, outflow=3};
string[] BCCodeNames = ["wall_with_slip", "exchange", "inflow", "outflow"];

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
        return format("BC(type=%s, fs=%s)", BCCodeNames[code], ((fs) ? to!string(*fs) : "null"));
    }
}

class FluidBlock {
public:
    int indx;
    int i, j;
    bool active;
    GasModel gmodel;
    ThermochemicalReactor thermochemUpdate;
    CQIndex cqi;
    int nic, njc;
    bool axiFlag;
    bool uSSTFlag;
    double[] dUdt_usst;
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
        uSSTFlag = Config.add_user_supplied_source_terms;
        dUdt_usst.length = cqi.n;
        foreach (ref du; dUdt_usst) du = 0.0;
        cfl = Config.cfl;
        x_order = Config.x_order;
        flux_calc = Config.flux_calc;
        compression_tol = Config.compression_tol;
        shear_tol = Config.shear_tol;
        //
        i = getJSONint(configData, "i", -1);
        j = getJSONint(configData, "j", -1);
        if (i < 0 || j < 0) {
            throw new Exception("Invalid indices for block.");
        }
        nic = Config.nics[i];
        njc = Config.njcs[j];
        active = getJSONbool(configData, "active", true);
        if (active) {
            if (Config.blk_ids[i][j] != indx) {
                writefln("indx=%d i=%d j=%d blk_ids=%s", indx, i, j, Config.blk_ids);
                throw new Exception("Block seems out of place.");
            }
        }
        //
        void setBC(ref BC bc, JSONValue jsonBC) {
            bc.code = BC_code_from_name(jsonBC["tag"].str);
            if (bc.code == BCCode.inflow) {
                JSONValue jsonFlow = jsonBC["flow_state"];
                JSONValue jsonGas = jsonFlow["gas"];
                double p = getJSONdouble(jsonGas, "p", 100.0e3);
                double T = getJSONdouble(jsonGas, "T", 300.0);
                double[] default_massf = [1.0, ];
                foreach (i; 1 .. gmodel.n_species) { default_massf ~= 0.0; }
                double[] massf = getJSONdoublearray(jsonGas, "massf", default_massf);
                bc.fs = new FlowState2D(gmodel);
                bc.fs.gas.p = p; bc.fs.gas.T = T; bc.fs.gas.massf[] = massf[];
                gmodel.update_thermo_from_pT(bc.fs.gas);
                gmodel.update_sound_speed(bc.fs.gas);
                double[] vel = getJSONdoublearray(jsonFlow, "vel", [0.0,0.0]);
                bc.fs.vel.set(vel[0], vel[1]);
            }
        }
        JSONValue jsonBCs = configData["bcs"];
        setBC(bc_west, jsonBCs["iminus"]);
        setBC(bc_east, jsonBCs["iplus"]);
        setBC(bc_south, jsonBCs["jminus"]);
        setBC(bc_north, jsonBCs["jplus"]);
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
        repr ~= format(", uSSTFlag=%s", uSSTFlag);
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
        return j*nic + i;
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
            foreach (i; 0 .. nic+1) {
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
            foreach (j; 0 .. njc+1) {
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
    } // end set_up_data_storage()

    void read_grid_data()
    {
        string fileName = format("%s/grid/grid-%04d-%04d.gz", Config.job_name, i, j);
        if (!(exists(fileName) && isFile(fileName))) {
            writefln("Grid file name: %s", fileName);
            throw new Exception("Grid file cannot be found.");
        }
        auto grid = new StructuredGrid(fileName, "gziptext");
        if (!(nic+1 == grid.niv && njc+1 == grid.njv && 1 == grid.nkv)) {
            writefln("nic=%d njc=%d grid.niv=%d grid.njv=%d grid.nkv=%d",
                     nic, njc, grid.niv, grid.njv, grid.nkv);
            throw new Exception("Grid size did not match.");
        }
        foreach (j; 0 .. njc+1) {
            foreach(i; 0 .. nic+1) {
                vertices[vertex_index(i,j)].set(grid[i,j]);
            }
        }
    } // end read_grid_data()

    @nogc
    void set_up_geometry()
    {
        foreach (f; ifaces) { f.compute_geometry(axiFlag); }
        foreach (f; jfaces) { f.compute_geometry(axiFlag); }
        foreach (c; cells) { c.compute_geometry(axiFlag); }
    }

    void read_flow_data(int tindx)
    {
        string fileName = format("%s/flow/t%04d/flow-%04d-%04d.gz", Config.job_name, tindx, i, j);
        if (!(exists(fileName) && isFile(fileName))) {
            writefln("Flow file name: %s", fileName);
            throw new Exception("Flow file cannot be found.");
        }
        auto byLine = new GzipByLine(fileName);
        foreach (varName; Config.iovar_names) {
            foreach (i; 0 .. n_cells) {
                auto line = byLine.front; byLine.popFront();
                double value = to!double(line.strip());
                // Do not overwrite the cell's geometry data.
                if (varName == "vol" || varName == "pos.x" || varName == "posy") continue;
                cells[i].iovar_set(varName, value);
            }
        }
        byLine.range.f.close();
        // Since not all thermo data may have been set,
        // we update the thermo assuming that p and T are correct.
        foreach (c; cells) { gmodel.update_thermo_from_pT(c.fs.gas); }
    } // end read_flow_data()

    void write_flow_data(int tindx)
    {
        string dirName = format("%s/flow/t%04d", Config.job_name, tindx);
        string fileName = dirName ~ format("/flow-%04d-%04d.gz", i, j);
        if (!(exists(dirName) && isDir(dirName))) {
            mkdir(dirName);
        }
        auto outfile = new GzipOut(fileName);
        foreach (varName; Config.iovar_names) {
            foreach (i; 0 .. n_cells) {
                outfile.compress(format("%e\n", cells[i].iovar_get(varName)));
            }
        }
        outfile.finish();
    } // end write_flow_data()

    @nogc
    void encode_conserved(size_t ftl)
    {
        foreach (c; cells) { c.encode_conserved(ftl, gmodel); }
    }

    @nogc
    void decode_conserved(size_t ftl)
    {
        foreach (c; cells) { c.decode_conserved(ftl, gmodel); }
    }

    @nogc
    void thermochemical_increment(double dt)
    {
        foreach (c; cells) { c.thermochemical_increment(dt, gmodel, thermochemUpdate); }
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
    }

    @nogc
    void update_conserved_for_stage(int stage, double dt)
    {
        foreach (f; ifaces) { f.calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi); }
        foreach (f; jfaces) { f.calculate_flux(fsL, fsR, gmodel, flux_calc, x_order, cqi); }
        // [TODO] One day, we should implement user-defined source terms.
        foreach (c; cells) { c.update_conserved_for_stage(stage, dt, axiFlag, dUdt_usst, uSSTFlag, gmodel); }
    }

    @nogc
    void transfer_conserved_quantities(size_t from, size_t dest)
    {
        foreach (c; cells) { c.U[dest][] = c.U[from][]; }
    }
} // end class FluidBlock
