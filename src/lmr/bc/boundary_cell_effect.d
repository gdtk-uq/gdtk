/**
 * Boundary cell effects are confined to the layer of
 * cells against a boundary.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2017-07-20
 */

module bc.boundary_cell_effect;

import std.json;
import std.string;

import globalconfig: FlowSolverException;
import lmr.fluidfvcell;
import fvinterface;
import geom;
import fluidblock;
import sfluidblock: SFluidBlock;
import globaldata;
import util.json_helper;
import bc;

BoundaryCellEffect make_BCE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bceType = jsonData["type"].str;
    BoundaryCellEffect newBCE;
    switch (bceType) {
    case "wall_function_cell_effect":
        newBCE = new BCE_WallFunction(blk_id, boundary);
        break;
    default:
        string errMsg = format("ERROR: The BoundaryCellEffect type: '%s' is unknown.", bceType);
        throw new FlowSolverException(errMsg);
    }
    return newBCE;
}

class BoundaryCellEffect {
public:
    FluidBlock blk;
    int which_boundary;
    string desc;

    this(int id, int boundary, string description)
    {
        blk = cast(FluidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        which_boundary = boundary;
        desc = description;
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

/**
 * The BCE_WallFunction object must be called AFTER the BIE_WallFunction object
 * because it relies on updated values for tke and omega supplied at the interface.
 *
 * In simcore.d, the application of boundary effects respects this ordering:
 * boundary interface actions are called first, then boundary cell actions.
 */

class BCE_WallFunction : BoundaryCellEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "WallFunction_CellEffect");
        _cells_need_to_be_flagged = true;
    } // end constructor

    override string toString() const
    {
        return "WallFunction_CellEffect()";
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new FlowSolverException("WallFunction_CellEffect bc not implemented for unstructured grids.");
    } // end apply_unstructured_grid()

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        // Apply the cell flags, just once.
        if (_cells_need_to_be_flagged) {
            foreach (i, f; bc.faces) {
                FluidFVCell c = (bc.outsigns[i] == 1) ? f.left_cells[0] : f.right_cells[0];
                c.allow_k_omega_update = false;
            }
            _cells_need_to_be_flagged = false;
        }
        // Do some real work.
        foreach (i, f; bc.faces) {
            FluidFVCell c = (bc.outsigns[i] == 1) ? f.left_cells[0] : f.right_cells[0];
            version(turbulence) {
                c.fs.turb[0] = f.fs.turb[0];
                c.fs.turb[1] = f.fs.turb[1];
            }
        }
    } // end apply_structured_grid()

private:
    bool _cells_need_to_be_flagged = true;
} // end class BCE_WallFunction
