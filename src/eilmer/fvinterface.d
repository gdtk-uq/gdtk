/**
 * fvinterface.d
 * Finite-volume cell-interface class, for use in the CFD codes.
 * Fluxes of conserved quantities are transported (between cells) across cell interfaces.

 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 *          2015-02-13: Keep an eye on the future of the moving_grid option.
 *          2015-05-04: keep references to adjoining cells and defining vertices.
 */

module fvinterface;

import std.conv;
import geom;
import gas;
import fvcore;
import fvvertex;
import fvcell;
import flowstate;
import conservedquantities;

class FVInterface {
public:
    size_t id;  // allows us to work out where, in the block, the interface is
    // Geometry
    Vector3 pos;           // position of the (approx) midpoint
    Vector3 gvel;          // grid velocity at interface, m/s
    double Ybar;           // Y-coordinate of the mid-point
    double length;         // Interface length in the x,y-plane
    double[] area;         // Area m**2 for each time-level.
                           // Area per radian in axisymmetric geometry
    Vector3 n;             // Direction cosines for unit normal
    Vector3 t1;            // tangent vector 1 (aka p)
    Vector3 t2;            // tangent vector 2 (aka q)
    FVVertex[] vtx;        // references to vertices for line (2D) and quadrilateral (3D) faces
    FVCell left_cell;      // interface normal points out of this adjoining cell
    FVCell right_cell;     // interface normal points into this adjoining cell
    // Flow
    FlowState fs;          // Flow properties
    ConservedQuantities F; // Flux conserved quantity per unit area

    this(GasModel gm, size_t id_init=0)
    {
	id = id_init;
	area.length = n_time_levels;
	gvel = Vector3(0.0,0.0,0.0); // default to fixed grid
	fs = new FlowState(gm, 100.0e3, [300.0,], Vector3(0.0,0.0,0.0));
	F = new ConservedQuantities(gm.n_species, gm.n_modes);
    }

    this(in FVInterface other, GasModel gm)
    {
	id = other.id;
	pos = other.pos;
	gvel = other.gvel;
	Ybar = other.Ybar;
	length = other.length;
	area = other.area.dup;
	n = other.n;
	t1 = other.t1;
	t2 = other.t2;
	fs = new FlowState(other.fs, gm);
	F = new ConservedQuantities(other.F);
    }

    @nogc
    void copy_values_from(in FVInterface other, uint type_of_copy)
    {
	switch (type_of_copy) {
	case CopyDataOption.minimal_flow:
	case CopyDataOption.all_flow:
	    fs.copy_values_from(other.fs);
	    F.copy_values_from(other.F);
	    break;
	case CopyDataOption.grid:
	    pos.refx = other.pos.x; pos.refy = other.pos.y; pos.refz = other.pos.z;
	    gvel.refx = other.gvel.x; gvel.refy = other.gvel.y; gvel.refz = other.gvel.z;
	    Ybar = other.Ybar;
	    length = other.length;
	    area[] = other.area[];
	    n.refx = other.n.x; n.refy = other.n.y; n.refz = other.n.z;
	    t1.refx = other.t1.x; t1.refy = other.t1.y; t1.refz = other.t1.z;
	    t2.refx = other.t2.x; t2.refy = other.t2.y; t2.refz = other.t2.z;
	    break;
	case CopyDataOption.all: 
	default:
	    id = other.id;
	    pos.refx = other.pos.x; pos.refy = other.pos.y; pos.refz = other.pos.z;
	    gvel.refx = other.gvel.x; gvel.refy = other.gvel.y; gvel.refz = other.gvel.z;
	    Ybar = other.Ybar;
	    length = other.length;
	    area[] = other.area[];
	    n.refx = other.n.x; n.refy = other.n.y; n.refz = other.n.z;
	    t1.refx = other.t1.x; t1.refy = other.t1.y; t1.refz = other.t1.z;
	    t2.refx = other.t2.x; t2.refy = other.t2.y; t2.refz = other.t2.z;
	    fs.copy_values_from(other.fs);
	    F.copy_values_from(other.F);
	} // end switch
    }

    @nogc
    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
	area[to_level] = area[from_level];
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FVInterface(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", pos=" ~ to!string(pos);
	repr ~= ", gvel=" ~ to!string(gvel);
	repr ~= ", Ybar=" ~ to!string(Ybar);
	repr ~= ", length=" ~ to!string(length);
	repr ~= ", area=" ~ to!string(area);
	repr ~= ", n=" ~ to!string(n);
	repr ~= ", t1=" ~ to!string(t1);
	repr ~= ", t2=" ~ to!string(2);
	repr ~= ", fs=" ~ to!string(fs);
	repr ~= ", F=" ~ to!string(F);
	repr ~= ")";
	return to!string(repr);
    }
} // end of class FV_Interface
