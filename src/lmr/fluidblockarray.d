// fluidblockarray.d
// A coordinating class for FluidBlock objects that
// is used in the shock-fitting calculations.
//
// PJ, 2021-Feb

module fluidblockarray;

import std.algorithm;
import std.conv;
import std.stdio;
import std.string;
import std.math;
import std.json;
import std.format;
import ntypes.complex;
import nm.number;
import geom;
import globalconfig;
import globaldata;
import json_helper;
version(mpi_parallel) {
    import mpi;
}


class FBArray {
    // Shock fitting is coordinated across arrays of FluidBlock objects.
    //
    // Here is the storage for that coordination data.
    int nib, njb, nkb; // Numbers of blocks in each index direction.
    int niv, njv, nkv; // Numbers of vertices (overall) in each index direction.
    int[] nics, njcs, nkcs; // Numbers of cells (per index direction) within each sub-block
    int[] blockIds; // Block ids in a single list.
    int[][][] blockArray; // Block ids arranged in an array with indices ib,jb,kb.
    bool shock_fitting; // Flag to indicate that we want to do shock fitting and move the grid.
    //
    // Shock-fitting data storage follows.
    double[][][] velocity_weights; // Fraction of shock-front velocity.
    Vector3[][] p_east; // East-most point on rail.
    Vector3[][] p_west; // West-most point on rail.
    number[][] face_ws; // Wave speeds at face centres.
    number[][] face_a;  // Sound-speed immediately after shock.
    Vector3[][] face_pos; // Positions of face centres.
    Vector3[][] vtx_vel; // Computed velocities for the vertices on the shock boundary.
    Vector3[][] vtx_dir; // Rail directions for the vertices on the shock boundary.
    Vector3[][] vtx_pos; // Positions of the vertices on the shock boundary.
    version(mpi_parallel) {
        // The tasks associated with this FBarray with have their own communicator
        // so that the can synchronize the content of the data storage arrays above.
        MPI_Comm mpicomm;
        double[] buffer;
    }

    this(int nib, int njb, int nkb, const(int[]) ids,
         int niv, int njv, int nkv,
         const(int[]) nics, const(int[]) njcs, const(int[]) nkcs,
         bool sf_flag)
    {
        // Construct with just configuration data.
        // The shock-fitting data will be added later, if relevant.
        assert(nib*njb*nkb == ids.length, "Incorrect number of block Ids");
        this.nib = nib; this.njb = njb; this.nkb = nkb;
        blockIds.length = ids.length;
        foreach (i; 0 .. ids.length) { blockIds[i] = ids[i]; }
        // There is a particular order of definition blocks
        // in the FBArray:new() function over in prep.lua.
        size_t list_index = 0;
        blockArray.length = nib;
        foreach (i; 0 .. nib) {
            blockArray[i].length = njb;
            foreach (j; 0 .. njb) {
                blockArray[i][j].length = nkb;
                foreach (k; 0 .. nkb) {
                    blockArray[i][j][k] = blockIds[list_index];
                    list_index++;
                }
            }
        }
        // Information about the underlying grid and its subdivision.
        this.niv = niv; this.njv = njv; this.nkv = nkv;
        this.nics.length = nics.length;
        foreach (i; 0 .. nics.length) { this.nics[i] = nics[i]; }
        this.njcs.length = njcs.length;
        foreach (i; 0 .. njcs.length) { this.njcs[i] = njcs[i]; }
        this.nkcs.length = nkcs.length;
        foreach (i; 0 .. nkcs.length) { this.nkcs[i] = nkcs[i]; }
        // Other bits.
        this.shock_fitting = sf_flag;
        if (sf_flag) {
            // Make space for the shock-fitting intermediate data.
            int njc = sum(njcs); int nkc = sum(nkcs);
            assert((njv==njc+1) &&
                   ((GlobalConfig.dimensions==2 && nkv==1 && nkc==1) ||
                    (GlobalConfig.dimensions==3 && nkv==nkc+1)),
                   "Mismatch in cells and vertices at shock boundary.");
            face_ws.length = njc; foreach (j; 0 .. njc) { face_ws[j].length = nkc; }
            face_a.length = njc; foreach (j; 0 .. njc) { face_a[j].length = nkc; }
            face_pos.length = njc; foreach (j; 0 .. njc) { face_pos[j].length = nkc; }
            vtx_vel.length = njv; foreach (j; 0 .. njv) { vtx_vel[j].length = nkv; }
            vtx_dir.length = njv; foreach (j; 0 .. njv) { vtx_dir[j].length = nkv; }
            vtx_pos.length = njv; foreach (j; 0 .. njv) { vtx_pos[j].length = nkv; }
            //
            version(mpi_parallel) {
                // Size buffers for the MPI exchange of the largest sub-block surface.
                // The largest array will be needed for the vertex positions or velocities.
                int n = 0;
                foreach (nk; nkcs) {
                    foreach (nj; njcs) {
                        n = max(n, (nj+1)*((GlobalConfig.dimensions == 2) ? 1 : nk+1));
                    }
                }
                this.buffer.length = n * 3;
                // Note that the 3 represents the maximum number of items per point,
                // to be sent in an MPI transfer.
            }
        }
    }
    this(const(FBArray) other)
    {
        this(other.nib, other.njb, other.nkb, other.blockIds,
             other.niv, other.njv, other.nkv,
             other.nics, other.njcs, other.nkcs,
             other.shock_fitting);
    }
    this(JSONValue json_data)
    {
        int nib = getJSONint(json_data, "nib", 0);
        int njb = getJSONint(json_data, "njb", 0);
        int nkb = getJSONint(json_data, "nkb", 0);
        int[] oops; oops.length = nib*njb*nkb; foreach(ref item; oops) { item = -1; }
        int[] ids = getJSONintarray(json_data, "idflatlist", oops);
        bool all_positive = true;
        foreach (item; ids) { if (item < 0) { all_positive = false; } }
        assert(all_positive, "One or more blocks ids are not as expected.");
        int niv = getJSONint(json_data, "niv", 0);
        int njv = getJSONint(json_data, "njv", 0);
        int nkv = getJSONint(json_data, "nkv", 1);
        int[] oops2; oops2.length = nib; foreach(ref item; oops2) { item = -1; }
        int[] nics = getJSONintarray(json_data, "nics", oops2);
        int[] oops3; oops3.length = njb; foreach(ref item; oops3) { item = -1; }
        int[] njcs = getJSONintarray(json_data, "njcs", oops3);
        int[] oops4; oops4.length = nkb; foreach(ref item; oops4) { item = -1; }
        int[] nkcs = getJSONintarray(json_data, "nkcs", oops4);
        bool sf_flag = getJSONbool(json_data, "shock_fitting", false);
        this(nib, njb, nkb, ids, niv, njv, nkv, nics, njcs, nkcs, sf_flag);
    }

    override string toString()
    {
        string result = format("FBArray(nib=%d, njb=%d, nkb=%d, "~
                               "blockIds=%s, blockArray=%s, "~
                               "niv=%d, njv=%d, nkv=%d, "~
                               "nics=%s, njcs=%s, nkcs=%s, shock_fitting=%s)",
                               nib, njb, nkb, blockIds, blockArray,
                               niv, njv, nkv, nics, njcs, nkcs, shock_fitting);
        return result;
    }

    void read_velocity_weights(string filename)
    {
        // File written by function write_shock_fitting_helper_files() in output.lua.
        velocity_weights.length = niv;
        foreach (i; 0 .. niv) {
            velocity_weights[i].length = njv;
            foreach (j; 0 .. njv) {
                velocity_weights[i][j].length = nkv;
            }
        }
        auto f = File(filename, "r");
        auto line = f.readln(); // Discard the comment line.
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    line = f.readln().strip();
                    formattedRead(line, "%e", &(velocity_weights[i][j][k]));
                }
            }
        }
    } // end read_velocity_weights()

    void read_rails_file(string filename)
    {
        // File written by function write_shock_fitting_helper_files() in output.lua.
        p_west.length = njv;
        p_east.length = njv;
        foreach (j; 0 .. njv) {
            p_west[j].length = nkv;
            p_east[j].length = nkv;
        }
        auto f = File(filename, "r");
        auto line = f.readln().strip(); // Discard the comment line.
        double pwx, pwy, pwz, pex, pey, pez;
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                line = f.readln().strip();
                formattedRead(line, "%e %e %e %e %e %e",
                              &pwx, &pwy, &pwz, &pex, &pey, &pez);
                p_west[j][k].set(pwx, pwy, pwz);
                p_east[j][k].set(pex, pey, pez);
            }
        }
    } // end read_rails_file()
} // end class FBArray
