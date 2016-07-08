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
import fvcore;
import geom;
import gas;
import globalconfig;
import flowstate;
import flowgradients;

class FVVertex {
public:
    size_t id;  // allows us to work out where, in the block, the vertex is
    // Geometry
    Vector3[] pos;  // x,y,z-Coordinates for time-levels, m
    Vector3[] vel;  // vertex velocity for time-levels, m/s
    // Derivatives of primary-cell variables.
    FlowGradients grad;
    Vector3*[] cloud_pos; // Positions of flow points for derivative calculation.
    FlowState[] cloud_fs; // References to flow states at those points.
    double[] cloud_weights; // Weights used in the least-squares gradient calculation.
    
    this(LocalConfig myConfig, size_t id_init=0)
    {
	id = id_init;
	pos.length = myConfig.n_grid_time_levels;
	vel.length = myConfig.n_grid_time_levels;
	grad = new FlowGradients(myConfig.gmodel.n_species);
    }

    this(FVVertex other)
    {
	id = other.id;
	pos = other.pos.dup;
	vel = other.vel.dup;
	grad = new FlowGradients(other.grad);
	// Because we copy the following pointers and references,
	// we cannot have const (or "in") qualifier on other.
	cloud_pos = other.cloud_pos.dup();
	cloud_fs = other.cloud_fs.dup();
    }

    @nogc 
    void copy_values_from(ref const(FVVertex) other)
    {
	if (!(this is other)) {
	    id = other.id;
	    foreach (i; 0 .. pos.length) {
	    	pos[i].refx = other.pos[i].x;
	    	pos[i].refy = other.pos[i].y;
	    	pos[i].refz = other.pos[i].z;
	    }
	    foreach (i; 0 .. vel.length) {
		vel[i].refx = other.vel[i].x;
		vel[i].refy = other.vel[i].y;
		vel[i].refz = other.vel[i].z;
	    }
	    grad.copy_values_from(other.grad);
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
	repr ~= ")";
	return to!string(repr);
    }

} // end class FVVertex
