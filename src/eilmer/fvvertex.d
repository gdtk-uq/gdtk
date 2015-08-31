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
import flowstate;

class FVVertex {
public:
    size_t id;  // allows us to work out where, in the block, the vertex is
    // Geometry
    Vector3[] pos;  // x,y,z-Coordinates for time-levels, m
    Vector3[] vel;  // vertex velocity for time-levels, m/s
    // Derivatives of primary-cell variables.
    double[][] grad_vel; // velocity derivatives stored as a second-order tensor
                         // [[du/dx du/dy du/dz]
                         //  [dv/dx dv/dy dv/dz]
                         //  [dw/dx dw/dy dw/dz]]
    Vector3 grad_T;      // Temperature derivatives (static temperature only)
    Vector3 grad_tke;    // turbulence kinetic energy
    Vector3 grad_omega;  // pseudo vorticity for k-omega turbulence
    Vector3[] grad_f;    // mass fraction derivatives
    Vector3 grad_pe;     // electron pressure derivatives
    Vector3[] cloud_pos; // Positions of flow points for derivative calculation.
    FlowState[] cloud_fs; // References to flow states at those points.

    this(in GasModel gm, size_t id_init=0)
    {
	id = id_init;
	pos.length = n_time_levels;
	vel.length = n_time_levels;
	grad_vel.length = 3;
	foreach(ref e; grad_vel) e.length = 3;
	grad_f.length = gm.n_species;
    }

    this(in FVVertex other)
    {
	id = other.id;
	pos = other.pos.dup;
	vel = other.vel.dup;
	grad_vel.length = 3;
	foreach(i; 0 .. 3) grad_vel[i] = other.grad_vel[i].dup; 
	grad_T = other.grad_T;
	grad_tke = other.grad_tke;
	grad_omega = other.grad_omega;
	foreach(i; 0 .. grad_f.length) grad_f[i] = other.grad_f[i];
	grad_pe = other.grad_pe;
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
	    foreach(i; 0 .. grad_vel.length) {
		foreach(j; 0 .. grad_vel[i].length)
		    grad_vel[i][j] = other.grad_vel[i][j];
	    }
	    grad_T.refx = other.grad_T.x;
	    grad_T.refy = other.grad_T.y;
	    grad_T.refz = other.grad_T.z;
	    grad_tke.refx = other.grad_tke.x;
	    grad_tke.refy = other.grad_tke.y;
	    grad_tke.refz = other.grad_tke.z;
	    grad_omega.refx = other.grad_omega.x;
	    grad_omega.refy = other.grad_omega.y;
	    grad_omega.refz = other.grad_omega.z;
	    foreach(i; 0 .. grad_f.length) {
		grad_f[i].refx = other.grad_f[i].x;
		grad_f[i].refy = other.grad_f[i].y;
		grad_f[i].refz = other.grad_f[i].z;
	    }
	    grad_pe.refx = other.grad_pe.x;
	    grad_pe.refy = other.grad_pe.y;
	    grad_pe.refz = other.grad_pe.z;
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
	repr ~= ", grad_vel=" ~ to!string(grad_vel);
	repr ~= ", grad_T=" ~ to!string(grad_T);
	repr ~= ", grad_tke=" ~ to!string(grad_tke);
	repr ~= ", grad_omega=" ~ to!string(grad_omega);
	repr ~= ", grad_f=" ~ to!string(grad_f);
	repr ~= ", grad_pe=" ~ to!string(grad_pe);
	repr ~= ")";
	return to!string(repr);
    }

/+ [TODO]
    int copy_grid_level_to_level(size_t from_level, size_t to_level);
+/
} // end class FVVertex
