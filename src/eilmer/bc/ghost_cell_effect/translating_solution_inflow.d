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
import sfluidblock;

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

        // We're going to load in FluidBlockLites, as if they were ordinary fluid blocks
        // The FluidBlockLite constructor expects a json table, so we make a nearby empty
        // one to keep it happy.
        JSONValue fakeJsonData;
        try {
            fakeJsonData = parseJSON!string(fakejsoncontent);
        } catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to parse dummy JSON content "));
        }

        // This json file has specific data needed for this boundary condition. It's created
        // by the python script that prepared the stream files for reading in.
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

        this.xs0  = getJSONdouble(jsondata, "xs0", -1.0);
        this.xs1  = getJSONdouble(jsondata, "xs1", -1.0);
        this.nic  = getJSONint(jsondata, "nic", 0);

        this.translation_velocity  = getJSONdouble(jsondata, "translation_velocity", -1.0);
        this.upstream_grid_rotation= getJSONdouble(jsondata, "upstream_grid_rotation", -1.0);
        this.upstream_grid_shift   = Vector3(getJSONdoublearray(jsondata, "upstream_grid_shift", []));

        double c = cos(-1.0*upstream_grid_rotation);
        double s = sin(-1.0*upstream_grid_rotation);

        JSONValue blkjson = jsondata["block_" ~ to!string(blk.id)];
        nfaces = getJSONint(blkjson, "nfaces", 0);
        nghost = getJSONint(blkjson, "nghost", 0);

        FluidBlockLite[2] fbl;
        fbl[0] = new FluidBlockLite(format(fileName, blk.id, 0), 0, fakeJsonData, Grid_t.structured_grid, "rawbinary");
        fbl[1] = new FluidBlockLite(format(fileName, blk.id, 1), 0, fakeJsonData, Grid_t.structured_grid, "rawbinary");

        fss.length = nfaces*nghost;

        SFluidBlock this_blk = cast(SFluidBlock) blk;
        if (!this_blk) throw new Error("FlowBlock must be a structured-grid block.");

        size_t nkc = this_blk.nkc;
        size_t njc = this_blk.njc;
        size_t ighost = 0;
        foreach (k; 0 .. nkc) {
          foreach (j; 0 .. njc) {
            foreach (n; 0 .. nghost) {
                assert (fbl[n].nic == nic);
                foreach (i; 0 .. fbl[n].nic) {
                    fss[ighost] ~= FlowState(GlobalConfig.gmodel_master, GlobalConfig.turb_model.nturb);
                    fss[ighost][i].gas.p   = fbl[n]["p",i,j,k];
                    fss[ighost][i].gas.T   = fbl[n]["T",i,j,k];
                    fss[ighost][i].gas.rho = fbl[n]["rho",i,j,k];
                    fss[ighost][i].gas.a   = fbl[n]["a",i,j,k];
                    fss[ighost][i].gas.u   = fbl[n]["u",i,j,k];
                    fss[ighost][i].gas.mu  = fbl[n]["mu",i,j,k];
                    fss[ighost][i].gas.k   = fbl[n]["k",i,j,k];
                    fss[ighost][i].S       = fbl[n]["S",i,j,k];
                    fss[ighost][i].mu_t    = fbl[n]["mu_t",i,j,k];
                    fss[ighost][i].k_t     = fbl[n]["k_t",i,j,k];
                    foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                        string sp_name = GlobalConfig.gmodel_master.species_name(sp);
                        fss[ighost][i].gas.massf[sp] = fbl[n][format("massf[%d]-%s", sp, sp_name),i,j,k];
                    }
                    // Change reference frame of velocity vector. Does the order of shift and rotation matter?
                    double velx  = fbl[n]["vel.x",i,j,k];
                    double vely  = fbl[n]["vel.y",i,j,k];
                    double velz  = fbl[n]["vel.z",i,j,k];
                    velx += translation_velocity;
                    fss[ighost][i].vel.x = velx*cos(upstream_grid_rotation) - velz*sin(upstream_grid_rotation);
                    fss[ighost][i].vel.y = vely;
                    fss[ighost][i].vel.z = velx*sin(upstream_grid_rotation) + velz*cos(upstream_grid_rotation);
                }
                ighost += 1;
            }
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

        int l,u;
        double w0;
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            size_t i = f.i_bndry;
            size_t ighost = i*blk.n_ghost_cell_layers + n;
            auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
            double gcx = ghost.pos[0].x.re;
            double gcy = ghost.pos[0].y.re;
            double gcz = ghost.pos[0].z.re;

            double gcx_prime = gcx*c - gcz*s - upstream_grid_shift.x.re - translation_velocity*t;
            //double gcy_prime = gcy - upstream_grid_shift.y.re;
            //double gcz_prime = gcx*s + gcz*c - upstream_grid_shift.z.re;

            get_interp_idxs_and_weight(gcx_prime, l, u, w0);
            ghost.fs.copy_average_values_from(fss[ighost][l], fss[ighost][u], w0);
            ighost += 1;
        }
    } // end apply_for_interface_structured_grid()

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        double c = cos(-1.0*upstream_grid_rotation);
        double s = sin(-1.0*upstream_grid_rotation);

        assert(nfaces==bc.faces.length);
        assert(nghost==blk.n_ghost_cell_layers);

        size_t ighost = 0;
        int l,u;
        double w0;
        foreach (fidx, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                auto ghost = (bc.outsigns[fidx] == 1) ? f.right_cells[n] : f.left_cells[n];
                double gcx = ghost.pos[0].x.re;
                double gcy = ghost.pos[0].y.re;
                double gcz = ghost.pos[0].z.re;

                double gcx_prime = gcx*c - gcz*s - upstream_grid_shift.x.re - translation_velocity*t;
                //double gcy_prime = gcy - upstream_grid_shift.y.re;
                //double gcz_prime = gcx*s + gcz*c - upstream_grid_shift.z.re;

                get_interp_idxs_and_weight(gcx_prime, l, u, w0);
                ghost.fs.copy_average_values_from(fss[ighost][l], fss[ighost][u], w0);
                ighost += 1;
            }
        }
    } // end apply_structured_grid()

private:
    string fileName, jsonFileName;
    int nfaces,nghost; // These are for this current block
    int nic;          // This is for the upstream sim
    double xs1, xs0;  // These are also the upstream sim
    double upstream_grid_rotation, translation_velocity;
    Vector3 upstream_grid_shift;
    FlowState[][] fss;

    void get_interp_idxs_and_weight(double gcxp, out int l, out int u, out double w0){
    /*
        We're given a position in the x direction, and need to figure out the two
        flowstates on either side of us, as well as a weight parameter to help
        interpolate between them.
    */
        double float_idx = (gcxp - xs0)/(xs1-xs0)*(nic-1);
        int int_l = to!int(floor(float_idx));
        int int_u = to!int(floor(float_idx))+1;
        l = mod(int_l, nic);
        u = mod(int_u, nic);
        double w1 = float_idx-int_l;
        w0 = 1.0-w1;
    }

}
