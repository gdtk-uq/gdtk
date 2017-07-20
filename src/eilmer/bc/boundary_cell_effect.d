/**
 * Boundary cell effects are confined to the layer of
 * cells against a boundary.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2017-07-20
 */

import std.json;
import std.string;

import grid;
import block;
import globaldata;
import json_helper;

void make_BCE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    

}

class BoundaryCellEffect {
public:
    Block blk;
    int which_boundary;
    string type;

    this(int id, int boundary, string _type)
    {
	blk = gasBlocks[id];
	which_boundary = boundary;
	type = _type;
    }
    // Most boundary cell effects will not need to do anything
    // special after construction.
    // However, the user-defined ghost cells bc need some
    // extra work done to set-up the Lua_state after all
    // of the blocks and bcs have been constructed.
    void post_bc_construction() {}
    override string toString() const
    {
	return "BoundaryCellEffect()";
    }
    void apply(double t, int gtl, int ftl)
    {
	final switch (blk.grid_type) {
	case Grid_t.unstructured_grid: 
	    apply_unstructured_grid(t, gtl, ftl);
	    break;
	case Grid_t.structured_grid:
	    apply_structured_grid(t, gtl, ftl);
	}
    }
    abstract void apply_unstructured_grid(double t, int gtl, int ftl);
    abstract void apply_structured_grid(double t, int gtl, int ftl);
}
