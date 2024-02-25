/**
 * fvvertex.d
 * Finite-volume cell-vertex class for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvvertex;

import std.string;
import std.conv;
import geom;
import gas;
import globalconfig;
import flowstate;
import flowgradients;
import lmr.fluidfvcell : FluidFVCell;
import fvinterface;
import lsqinterp;
import ntypes.complex;
import nm.number;

class FVVertex {
public:
    size_t id;
    // Geometry
    Vector3[] pos;  // x,y,z-Coordinates for time-levels, m
    Vector3[] vel;  // vertex velocity for time-levels, m/s
    // Derivatives of primary-cell variables.
    FlowGradients* grad;
    number radial_pos_norm = 0;
    Vector3*[] cloud_pos; // Positions of flow points for derivative calculation.
    FlowState*[] cloud_fs; // References to flow states at those points.
    FluidFVCell[] cell_cloud; // for the MLP limiter we need access to the gradients within each cell
    WLSQGradWorkspace* ws_grad;
    LSQInterpGradients* gradients; // needed for the MLP limiter

    @disable this();
    
    this(LocalConfig myConfig,
         bool allocate_spatial_deriv_lsq_workspace,
         int id_init=-1)
    {
        id = id_init;
        pos.length = myConfig.n_grid_time_levels;
        vel.length = myConfig.n_grid_time_levels;
        foreach (ref v; vel) v.set(0.0, 0.0, 0.0);
        grad = new FlowGradients(myConfig);
        if (allocate_spatial_deriv_lsq_workspace) {
            ws_grad = new WLSQGradWorkspace();
        }
    }

    this(FVVertex other) // not const; see note below
    {
        id = other.id;
        pos = other.pos.dup;
        vel = other.vel.dup;
        grad = new FlowGradients(*(other.grad));
        if (other.ws_grad) {
            ws_grad = new WLSQGradWorkspace(*(other.ws_grad));
        }
        // Because we copy the following pointers,
        // we cannot have const (or "in") qualifier on other.
        cloud_pos = other.cloud_pos.dup();
        cloud_fs = other.cloud_fs.dup();
    }

    @nogc
    void copy_values_from(ref const(FVVertex) other)
    {
        if (!(this is other)) {
            id = other.id;
            foreach (i; 0 .. pos.length) { pos[i].set(other.pos[i]); }
            foreach (i; 0 .. vel.length) { vel[i].set(other.vel[i]); }
            grad.copy_values_from(*(other.grad));
            // omit ws_grad
        }
    } // end copy_values_from()

    @nogc
    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
        pos[to_level] = pos[from_level];
        vel[to_level] = vel[from_level];
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "FVVertex(";
        repr ~= "id=" ~ to!string(id);
        repr ~= ", pos=" ~ to!string(pos);
        repr ~= ", vel=" ~ to!string(vel);
        repr ~= ", grad=" ~ to!string(grad);
        repr ~= ", cloud_pos=" ~ to!string(cloud_pos);
        repr ~= ", cloud_fs=" ~ to!string(cloud_fs);
        // omit ws_grad
        repr ~= ")";
        return to!string(repr);
    }

} // end class FVVertex
