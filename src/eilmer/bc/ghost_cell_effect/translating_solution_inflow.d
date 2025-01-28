// translating_solution_inflow.d

module bc.ghost_cell_effect.translating_solution_inflow;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.json;

import geom;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvcell;
import fluidblock;
import gas;
import bc;
import flowsolution;
import json_helper;

/*
    Custom inflow boundary condition for NNG and VW's combustion work.

    The goal is to take a freeze frame from another Eilmer solution, called sim 1,
    and imagine feeding it transiently into a second simulation, called sim 2. 

     _____________________     _____
    |                     |   /     *_____
    |  Flow solution 1    |  /solut       *________
    |_____________________| /_____ ion 1->         />
                                  *______         />
                                         *_______/>

*/

immutable string fakejsoncontent="
{
\"block_0\": {
    \"omegaz\": 0.000000000000000000e+00,
},
}
";

/*
	Implementation of python-style modulo function, from the D forums.
*/
@nogc pure int mod(int n, int d){
    int r = n % d;
    return sgn(r) == -(sgn(d)) ? r + d : r;
 }

class TranslatingSolutionInflowEffect : GhostCellEffect {
public:
    this(int id, int boundary, string fileName, string jsonFileName)
    {
        super(id, boundary, "translatingSolutionInflowEffect");
        this.fileName = fileName;
        this.jsonFileName = jsonFileName;

        string content;
        try {
            content = readText(jsonFileName);
        } catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to read config file: ", jsonFileName));
        }
        JSONValue jsondata;
        try {
            jsondata = parseJSON!string(content);
        } catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to parse JSON from file: ", jsonFileName));
        }

        JSONValue blkjson = jsondata["block_" ~ to!string(blk.id)];
        jbins = getJSONintarray(blkjson, "jbins", []);
        if (jbins.length==0) throw new Error("Failed to find jbins array for block");
        kbins = getJSONintarray(blkjson, "kbins", []);
        if (jbins.length==0) throw new Error("Failed to find kbins array for block");

        writefln("mpi %d blk %d read in %d jbins %s and %d kbins", GlobalConfig.mpi_rank_for_local_task, blk.id, jbins.length, jbins, kbins.length);

        JSONValue fakeJsonData;
        try {
            fakeJsonData = parseJSON!string(fakejsoncontent);
        } catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to parse dummy JSON content "));
        }

        this.xs0  = getJSONdouble(jsondata, "xs0", -1.0);
        this.xs1  = getJSONdouble(jsondata, "xs1", -1.0);
        this.niv  = getJSONint(jsondata, "niv", 0);

        this.translation_velocity  = getJSONdouble(jsondata, "translation_velocity", -1.0);
        this.upstream_grid_rotation= getJSONdouble(jsondata, "upstream_grid_rotation", -1.0);
        this.upstream_grid_shift   = Vector3(getJSONdoublearray(jsondata, "upstream_grid_shift", []));

        double upstream_frame_velx = cos(upstream_grid_rotation)*translation_velocity;
        double upstream_frame_vely = sin(upstream_grid_rotation)*translation_velocity;

        FluidBlockLite fbl = new FluidBlockLite(fileName, 0, fakeJsonData, Grid_t.structured_grid, "rawbinary");

        // TODO: We need to make sure that the ordering of faces in eilmer matches the order we assumed
        // in the code that builds the jbins and kbins array
        fss.length = jbins.length;
        foreach(fidx; 0 .. jbins.length){
            int j = jbins[fidx];
            int k = kbins[fidx];
            foreach(i; 0 .. fbl.nic){
                fss[fidx] ~= FlowState(GlobalConfig.gmodel_master, GlobalConfig.turb_model.nturb);
                fss[fidx][i].vel.x   = fbl["vel.x",i,j,k];
                fss[fidx][i].vel.y   = fbl["vel.y",i,j,k];
                fss[fidx][i].vel.z   = fbl["vel.z",i,j,k];
                fss[fidx][i].gas.p   = fbl["p",i,j,k];
                fss[fidx][i].gas.T   = fbl["T",i,j,k];
                fss[fidx][i].gas.rho = fbl["rho",i,j,k];
                fss[fidx][i].gas.a   = fbl["a",i,j,k];
                fss[fidx][i].gas.u   = fbl["u",i,j,k];
                fss[fidx][i].gas.mu  = fbl["mu",i,j,k];
                fss[fidx][i].gas.k   = fbl["k",i,j,k];
                fss[fidx][i].S       = fbl["S",i,j,k];
                fss[fidx][i].mu_t    = fbl["mu_t",i,j,k];
                fss[fidx][i].k_t     = fbl["k_t",i,j,k];
                foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                    string sp_name = GlobalConfig.gmodel_master.species_name(sp);
                    fss[fidx][i].gas.massf[sp] = fbl[format("massf[%d]-%s", sp, sp_name),i,j,k];
                }
                // Change reference frame of velocity vector
                fss[fidx][i].vel.x += upstream_frame_velx;
                fss[fidx][i].vel.y += upstream_frame_vely;
            }
        }
    }

    override string toString() const
    {
        return format("translatingSolutionInflow(fileName=\"%s\", jsonFileName=\"%s\")",
                      fileName, jsonFileName);
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto ghost0 = (bc.outsigns[f.i_bndry] == 1) ? f.right_cell : f.left_cell;
    }

    // not @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            auto ghost0 = (bc.outsigns[i] == 1) ? f.right_cell : f.left_cell;
        }
    } // end apply_unstructured_grid()

    // not @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        double c = cos(-1.0*upstream_grid_rotation);
        double s = sin(-1.0*upstream_grid_rotation);

        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            size_t i = f.i_bndry;
            auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
            double gcx = ghost.pos[0].x.re;
            double gcy = ghost.pos[0].y.re;
            double gcz = ghost.pos[0].z.re;

            double gcx_prime = gcx*c - gcy*s - upstream_grid_shift.x - translation_velocity*t;
            double gcy_prime = gcx*s + gcy*c - upstream_grid_shift.y;
            double gcz_prime = gcz - upstream_grid_shift.z;
            int ibin = get_bin_that_ghost_cells_are_inside_of(gcx_prime, xs0, xs1, niv);
            ghost.fs.copy_values_from(fss[i][ibin]);
        }
    } // end apply_for_interface_structured_grid()

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        double c = cos(-1.0*upstream_grid_rotation);
        double s = sin(-1.0*upstream_grid_rotation);

        size_t ighost = 0;
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
                double gcx = ghost.pos[0].x.re;
                double gcy = ghost.pos[0].y.re;
                double gcz = ghost.pos[0].z.re;

                double gcx_prime = gcx*c - gcy*s - upstream_grid_shift.x - translation_velocity*t;
                double gcy_prime = gcx*s + gcy*c - upstream_grid_shift.y;
                double gcz_prime = gcz - upstream_grid_shift.z;
                int ibin = get_bin_that_ghost_cells_are_inside_of(gcx_prime, xs0, xs1, niv);
                ghost.fs.copy_values_from(fss[ighost][ibin]);
                ighost += 1;
            }
        }
    } // end apply_structured_grid()

private:
    string fileName, jsonFileName;
    int niv;
    double xs1, xs0;
    double upstream_grid_rotation, translation_velocity;
    Vector3 upstream_grid_shift;
    FlowState[][] fss;
    int[] jbins, kbins;

    int get_bin_that_ghost_cells_are_inside_of(double gcx, double llim, double ulim, int n){
    /*
        Taking a ghost cell position in the sim 1 reference frame, figure out which of
        sim 1's cells that ghost cell would be inside of.
    */ 
        double float_idx = (gcx - llim)/(ulim-llim)*(n-1);
        int int_idx = to!int(floor(float_idx));
        int wrap_idx = mod(int_idx, n-1);
        return wrap_idx;
    }
}
