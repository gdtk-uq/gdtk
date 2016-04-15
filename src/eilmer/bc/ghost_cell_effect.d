// ghost_cell_effect.d
//
// RG & PJ 2015-12-03 : first hack

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;

import geom;
import sgrid;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
import fvcore;
import fvinterface;
import fvcell;
import block;
import sblock;
import gas;
import user_defined_effects;
import bc;

@nogc
void reflect_normal_velocity(ref FlowState fs, in FVInterface IFace)
// Reflects normal velocity with respect to the supplied interface.
//
// The process is to rotate the velocity vector into the local frame of
// the interface, negate the normal (local x-component) velocity and
// rotate back to the global frame.
{
    fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.vel.refx = -(fs.vel.x);
    fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

@nogc
void reflect_normal_magnetic_field(ref FlowState fs, in FVInterface IFace)
{
    fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.B.refx = -(fs.B.x);
    fs.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

GhostCellEffect make_GCE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string gceType = jsonData["type"].str;
    // At the point at which we call this function, we may be inside the block-constructor.
    // Don't attempt the use the block-owned gas model.
    auto gmodel = GlobalConfig.gmodel_master; 
    GhostCellEffect newGCE;
    switch (gceType) {
    case "internal_copy_then_reflect":
	newGCE = new GhostCellInternalCopyThenReflect(blk_id, boundary);
	break;
    case "flowstate_copy":
	auto flowstate = new FlowState(jsonData["flowstate"], gmodel);
	newGCE = new GhostCellFlowStateCopy(blk_id, boundary, flowstate);
	break;
    case "extrapolate_copy":
	int xOrder = getJSONint(jsonData, "x_order", 0);
	newGCE = new GhostCellExtrapolateCopy(blk_id, boundary, xOrder);
	break;
    case "fixed_pressure":
	double p_outside = getJSONdouble(jsonData, "p_outside", 1.0e5);
	newGCE = new GhostCellFixedP(blk_id, boundary, p_outside);
	break;
    case "fixed_pressure_temperature":
	double p_outside = getJSONdouble(jsonData, "p_outside", 1.0e5);
	double T_outside = getJSONdouble(jsonData, "T_outside", 300.0);
	newGCE = new GhostCellFixedPT(blk_id, boundary, p_outside, T_outside);
	break;
    case "from_stagnation_condition":
	auto stagnation_condition = new FlowState(jsonData["stagnation_condition"], gmodel);
	string direction_type = getJSONstring(jsonData, "direction_type", "normal");
	double direction_x = getJSONdouble(jsonData, "direction_x", 1.0);
	double direction_y = getJSONdouble(jsonData, "direction_y", 0.0);
	double direction_z = getJSONdouble(jsonData, "direction_z", 0.0);
	double alpha = getJSONdouble(jsonData, "alpha", 0.0);
	double beta = getJSONdouble(jsonData, "beta", 0.0);
	double mass_flux = getJSONdouble(jsonData, "mass_flux", 0.0);
	double relax_factor = getJSONdouble(jsonData, "relax_factor", 0.1);
	newGCE = new GhostCellFromStagnation(blk_id, boundary, stagnation_condition,
					     direction_type, 
					     Vector3(direction_x, direction_y, direction_z),
					     alpha, beta,
					     mass_flux, relax_factor);
	break;
    case "full_face_exchange_copy":
	int otherBlock = getJSONint(jsonData, "other_block", -1);
	string otherFaceName = getJSONstring(jsonData, "other_face", "none");
	int neighbourOrientation = getJSONint(jsonData, "orientation", 0);
	bool rvq = getJSONbool(jsonData, "reorient_vector_quantities", false);
	double[] Rmatrix = getJSONdoublearray(jsonData, "Rmatrix",
					      [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]);
	newGCE = new GhostCellFullFaceExchangeCopy(blk_id, boundary,
						   otherBlock, face_index(otherFaceName),
						   neighbourOrientation,
						   rvq, Rmatrix);
	break;
    case "mapped_cell_exchange_copy":
	bool transform_pos = getJSONbool(jsonData, "transform_position", false);
	Vector3 c0 = getJSONVector3(jsonData, "c0", Vector3(0.0,0.0,0.0));
	Vector3 n = getJSONVector3(jsonData, "n", Vector3(0.0,0.0,1.0));
	double alpha = getJSONdouble(jsonData, "alpha", 0.0);
	Vector3 delta = getJSONVector3(jsonData, "delta", Vector3(0.0,0.0,0.0));
	bool lmc = getJSONbool(jsonData, "list_mapped_cells", false);
	bool rvq = getJSONbool(jsonData, "reorient_vector_quantities", false);
	double[] Rmatrix = getJSONdoublearray(jsonData, "Rmatrix",
					      [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]);
	newGCE = new GhostCellMappedCellExchangeCopy(blk_id, boundary,
						     transform_pos, c0, n, alpha, delta, lmc,
						     rvq, Rmatrix);
	break;
    case "user_defined":
	string fname = getJSONstring(jsonData, "filename", "none");
	newGCE = new UserDefinedGhostCell(blk_id, boundary, fname);
	break;
    default:
	string errMsg = format("ERROR: The GhostCellEffect type: '%s' is unknown.", gceType);
	throw new Exception(errMsg);
    }
    return newGCE;
}

class GhostCellEffect {
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
    void post_bc_construction() {}
    override string toString() const
    {
	return "GhostCellEffect()";
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
} // end class GhostCellEffect

class GhostCellInternalCopyThenReflect : GhostCellEffect {
public:
    
    this(int id, int boundary)
    {
	super(id, boundary, "InternalCopyThenReflect");
    }

    override string toString() const 
    {
	return "InternalCopyThenReflect()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	// In contrast with the structured-grid implementation, below,
	// we use a single source cell because we're not sure that
	// the next cell in the interior list (if there is one) will
	// be positioned as a mirror image of the second ghost cell.
	// If the grid is clustered toward the boundary, the error
	// introduced by this zero-order reconstruction will be mitigated.
	// PJ 2016-04-12
	BasicCell src_cell, ghost0, ghost1;
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		src_cell = f.left_cells[0];
		ghost0 = f.right_cells[0];
		ghost1 = f.right_cells[1];
	    } else {
		src_cell = f.right_cells[0];
		ghost0 = f.left_cells[0];
		ghost1 = f.left_cells[1];
	    }
	    ghost0.fs.copy_values_from(src_cell.fs);
	    reflect_normal_velocity(ghost0.fs, f);
	    ghost1.fs.copy_values_from(src_cell.fs);
	    reflect_normal_velocity(ghost1.fs, f);
	    if (blk.myConfig.MHD) {
		reflect_normal_magnetic_field(ghost0.fs, f);
		reflect_normal_magnetic_field(ghost1.fs, f);
	    }
	} // end foreach face
    } // end apply_unstructured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVInterface IFace;
	auto copy_opt = CopyDataOption.minimal_flow;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.north];
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j-1,k);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.east];
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i-1,j,k);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.south];
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j+1,k);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.west];
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i+1,j,k);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.top];
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j,k-1);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.bottom];
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j,k+1);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell.fs, IFace);
		    if (blk.myConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell.fs, IFace);
		    }
		} // end j loop
	    } // for i
	    break;
	} // end switch which_boundary
    } // end apply_structured_grid()
} // end class GhostCellInternalCopyThenReflect

class GhostCellFlowStateCopy : GhostCellEffect {
public:
    FlowState fstate;

    this(int id, int boundary, in FlowState _fstate)
    {
	super(id, boundary, "flowStateCopy");
	fstate = new FlowState(_fstate);
    }

    override string toString() const
    {
	return "flowStateCopy(fstate=" ~ to!string(fstate) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BasicCell ghost0, ghost1;
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		ghost0 = f.right_cells[0];
		ghost1 = f.right_cells[1];
	    } else {
		ghost0 = f.left_cells[0];
		ghost1 = f.left_cells[1];
	    }
	    ghost0.fs.copy_values_from(fstate);
	    ghost1.fs.copy_values_from(fstate);
	} // end foreach face
    } // end apply_unstructured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
    	// Fill ghost cells with data from just inside the boundary
	// using zero-order extrapolation (i.e. just copy the data).
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVInterface dest_face;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_structured_grid()

} // end class GhostCellFlowStateCopy

class GhostCellExtrapolateCopy : GhostCellEffect {
public:
    int xOrder;
    
    this(int id, int boundary, int x_order)
    {
	super(id, boundary, "ExtrapolateCopy");
	xOrder = x_order;
    }

    override string toString() const 
    {
	return "ExtrapolateCopy(x_order=" ~ to!string(xOrder) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BasicCell src_cell, ghost0, ghost1;
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		src_cell = f.left_cells[0];
		ghost0 = f.right_cells[0];
		ghost1 = f.right_cells[1];
	    } else {
		src_cell = f.right_cells[0];
		ghost0 = f.left_cells[0];
		ghost1 = f.left_cells[1];
	    }
	    if (xOrder == 1) {
		throw new Error("First order extrapolation not implemented.");
	    } else {
		// Zero-order extrapolation.
		ghost0.fs.copy_values_from(src_cell.fs);
		ghost1.fs.copy_values_from(src_cell.fs);
	    }
	} // end foreach face
    } // end apply_unstructured_grid()
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	// Fill ghost cells with data from just inside the boundary
	// using zero-order extrapolation (i.e. just copy the data).
	// We assume that this boundary is an outflow boundary.
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVCell cell_1, cell_2;
	auto gmodel = blk.myConfig.gmodel;
	size_t nsp = gmodel.n_species;
	size_t nmodes = gmodel.n_modes;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (j-1)        (j)           (j+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j-1,k);
			dest_cell = blk.get_cell(i,j+1,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (j)        (j+1)       (j+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j+2,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    }
		    else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j+1,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j+2,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    } 
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (i-1)        (i)           (i+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i-1,j,k);
			dest_cell = blk.get_cell(i+1,j,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (i)        (i+1)       (i+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i+2,j,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    }
		    else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i+1,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i+2,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (j+1)        (j)           (j-1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j+1,k);
			dest_cell = blk.get_cell(i,j-1,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (j)        (j-1)       (j-2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j-2,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j-1,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j-2,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  ---[ghost cell 2]---|--- [dest] ---|||--- [1] ---|---[2]----
			//      (i-2)                 (i-1)           (i)       (i+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i+1,j,k);
			dest_cell = blk.get_cell(i-1,j,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[dest]---|---[1]---|||---[2]---|------|
			//       (i-2)       (i-1)       (i)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i-2,j,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i-1,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i-2,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (k-1)        (k)           (k+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j,k-1);
			dest_cell = blk.get_cell(i,j,k+1);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (k)        (k+1)       (k+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j,k+2);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j,k+1);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j,k+2);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (k+1)        (k)           (k-1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j,k+2);
			dest_cell = blk.get_cell(i,j,k-1);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (k)        (k-1)       (k-2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j,k-2);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j,k-1);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j,k-2);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_structured_grid()
} // end class GhostCellExtrapolateCopy

class GhostCellFixedP : GhostCellEffect {
public:
    double p_outside;
    
    this(int id, int boundary, double p_outside)
    {
	super(id, boundary, "FixedP");
	this.p_outside = p_outside;
    }

    override string toString() const 
    {
	return "FixedP(p_outside=" ~ to!string(p_outside) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BasicCell src_cell, ghost0, ghost1;
	BoundaryCondition bc = blk.bc[which_boundary];
	auto gmodel = blk.myConfig.gmodel;
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		src_cell = f.left_cells[0];
		ghost0 = f.right_cells[0];
		ghost1 = f.right_cells[1];
	    } else {
		src_cell = f.right_cells[0];
		ghost0 = f.left_cells[0];
		ghost1 = f.left_cells[1];
	    }
	    ghost0.fs.copy_values_from(src_cell.fs);
	    ghost0.fs.gas.p = p_outside;
	    gmodel.update_thermo_from_pT(ghost0.fs.gas);
	    ghost1.fs.copy_values_from(src_cell.fs);
	    ghost1.fs.gas.p = p_outside;
	    gmodel.update_thermo_from_pT(ghost1.fs.gas);
	} // end foreach face
    } // end apply_unstructured_grid()
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_structured_grid()
} // end class GhostCellFixedP

class GhostCellFixedPT : GhostCellEffect {
public:
    double p_outside;
    double T_outside;
    
    this(int id, int boundary, double p_outside, double T_outside)
    {
	super(id, boundary, "FixedPT");
	this.p_outside = p_outside;
	this.T_outside = T_outside;
    }

    override string toString() const 
    {
	return "FixedPT(p_outside=" ~ to!string(p_outside) ~ ", T_outside=" ~ to!string(T_outside) ~")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BasicCell src_cell, ghost0, ghost1;
	BoundaryCondition bc = blk.bc[which_boundary];
	auto gmodel = blk.myConfig.gmodel;
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		src_cell = f.left_cells[0];
		ghost0 = f.right_cells[0];
		ghost1 = f.right_cells[1];
	    } else {
		src_cell = f.right_cells[0];
		ghost0 = f.left_cells[0];
		ghost1 = f.left_cells[1];
	    }
	    ghost0.fs.copy_values_from(src_cell.fs);
	    ghost0.fs.gas.p = p_outside;
	    foreach(ref elem; ghost0.fs.gas.T) elem = T_outside; 
	    gmodel.update_thermo_from_pT(ghost0.fs.gas);
	    ghost1.fs.copy_values_from(src_cell.fs);
	    ghost1.fs.gas.p = p_outside;
	    foreach(ref elem; ghost1.fs.gas.T) elem = T_outside; 
	    gmodel.update_thermo_from_pT(ghost1.fs.gas);
	} // end foreach face
    } // end apply_unstructured_grid()
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = p_outside;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = T_outside; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_structured_grid()
} // end class GhostCellFixedPT

class GhostCellFromStagnation : GhostCellEffect {
public:
    FlowState stagnation_condition;
    string direction_type;
    Vector3 direction_vector;
    double alpha;
    double beta;
    double mass_flux;
    double relax_factor;
private:
    FlowState inflow_condition;
    double stagnation_entropy;
    double stagnation_enthalpy;
    double p0_min;
    double p0_max;

public:    
    this(int id, int boundary, in FlowState stagnation_condition,
	 string direction_type, in Vector3 vec, double alpha, double beta,
	 double mass_flux, double relax_factor)
    {
	super(id, boundary, "FromStagnation");
	this.stagnation_condition = new FlowState(stagnation_condition);
	this.direction_type = direction_type;
	this.direction_vector = vec;
	this.direction_vector.normalize();
	this.alpha = alpha;
	this.beta = beta;
	this.mass_flux = mass_flux;
	this.relax_factor = relax_factor;
	auto gmodel = blk.myConfig.gmodel;
	stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
	stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
	inflow_condition = new FlowState(stagnation_condition);
	// The following limits are used when adjusting the stagnation pressure
	// to achieve a fixed mass flow per unit area.
	// The initially-specified stagnation condition needs to be a reasonable guess
	// because T0 will be held fixed and p0 adjusted within the following bounds.
	p0_min = 0.1 * stagnation_condition.gas.p;
	p0_max = 10.0 * stagnation_condition.gas.p;
    }

    override string toString() const 
    {
	return "FromStagnation(stagnation_condition=" ~ to!string(stagnation_condition) ~ 
	    ", direction_type=" ~ direction_type ~ 
	    ", direction_vector=" ~ to!string(direction_vector) ~
	    ", alpha=" ~ to!string(alpha) ~ ", beta=" ~ to!string(beta) ~ 
	    ", mass_flux=" ~ to!string(mass_flux) ~
	    ", relax_factor=" ~ to!string(relax_factor) ~ ")";
    }

    void set_velocity_components(ref Vector3 vel, double speed, ref FVInterface face)
    {
	switch ( direction_type ) {
	case "uniform":
	    // Given a flow direction.
	    vel.refx = speed * direction_vector.x;
   	    vel.refy = speed * direction_vector.y;
	    vel.refz = speed * direction_vector.z;
	    break;
	case "axial":
	    // Axial-flow through a presumably circular surface.
	    // [TODO] 27-Feb-2014 through 16-Oct-2015 PJ:
	    // check that this fall-through is OK.
	case "radial":
	    // For turbo inflow, beta sets the axial flow.
	    double vz = speed * sin(beta);
	    // ...and alpha sets the angle of the flow in the plane of rotation.
	    double vt = speed * cos(beta) * sin(alpha);
	    double vr = -speed * cos(beta) * cos(alpha);
	    double x = face.pos.x;
	    double y = face.pos.y;
	    double rxy = sqrt(x*x + y*y);
	    vel.refx = vr * x/rxy - vt * y/rxy;
   	    vel.refy = vt * x/rxy + vr * y/rxy;
	    vel.refz = vz;
	    break;
	case "normal":
	default:
	    vel.refx = speed * face.n.x;
   	    vel.refy = speed * face.n.y;
	    vel.refz = speed * face.n.z;
	}
    } // end set_velocity_components()

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("GhostCellFromStagnation.apply_unstructured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVInterface face;
	auto gmodel = blk.myConfig.gmodel;

	double p_stag = 0.0;
	double T_stag = 0.0; // temporary

	final switch (which_boundary) {
	case Face.north:
	    throw new Error("GhostCellFromStagnation.apply_structured_grid() " ~
			    "not yet implemented for north boundary");
	    // j = blk.jmax;
	    // for (k = blk.kmin; k <= blk.kmax; ++k) {
	    // 	for (i = blk.imin; i <= blk.imax; ++i) {
	    // 	    src_cell = blk.get_cell(i,j,k);
	    // 	    dest_cell = blk.get_cell(i,j+1,k);
	    // 	    dest_cell.fs.gas.p = p_stag; // [TODO] FIX-ME
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	    dest_cell = blk.get_cell(i,j+2,k);
	    // 	    dest_cell.fs.gas.p = p_stag;
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	} // end i loop
	    // } // for k
	    // break;
	case Face.east:
	    throw new Error("GhostCellFromStagnation.apply_structured_grid() " ~
			    "not yet implemented for east boundary");
	    // i = blk.imax;
	    // for (k = blk.kmin; k <= blk.kmax; ++k) {
	    // 	for (j = blk.jmin; j <= blk.jmax; ++j) {
	    // 	    src_cell = blk.get_cell(i,j,k);
	    // 	    dest_cell = blk.get_cell(i+1,j,k);
	    // 	    dest_cell.fs.gas.p = p_stag; // [TODO] FIX-ME
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	    dest_cell = blk.get_cell(i+2,j,k);
	    // 	    dest_cell.fs.gas.p = p_stag;
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	} // end j loop
	    // } // for k
	    // break;
	case Face.south:
	    throw new Error("GhostCellFromStagnation.apply_structured_grid() " ~
			    "not yet implemented for south boundary");
	    // j = blk.jmin;
	    // for (k = blk.kmin; k <= blk.kmax; ++k) {
	    // 	for (i = blk.imin; i <= blk.imax; ++i) {
	    // 	    src_cell = blk.get_cell(i,j,k);
	    // 	    dest_cell = blk.get_cell(i,j-1,k);
	    // 	    dest_cell.fs.gas.p = p_stag; // [TODO] FIX-ME
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	    dest_cell = blk.get_cell(i,j-2,k);
	    // 	    dest_cell.fs.gas.p = p_stag;
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	} // end i loop
	    // } // for k
	    // break;
	case Face.west:
	    i = blk.imin;
	    // First, estimate the current bulk inflow condition.
	    double area = 0.0;
	    double rhoUA = 0.0; // current mass_flux through boundary
	    double rhovxA = 0.0; // mass-weighted x-velocity
	    double rhovyA = 0.0;
	    double rhovzA = 0.0;
	    double rhoA = 0.0;
	    double pA = 0.0;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    auto cell = blk.get_cell(i,j,k);
		    face = cell.iface[Face.west];
		    area += face.area[0];
		    double local_rhoA = cell.fs.gas.rho * face.area[0];
		    rhoA += local_rhoA;
		    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
		    rhovxA += local_rhoA * cell.fs.vel.x;
		    rhovyA += local_rhoA * cell.fs.vel.y;
		    rhovzA += local_rhoA * cell.fs.vel.z;
		    pA += cell.fs.gas.p * face.area[0];
		} // end j loop
	    } // end k loop
	    if ( mass_flux > 0.0 && ftl == 0 ) {
		// Adjust the stagnation pressure to better achieve the specified mass flux.
		// Note that we only do this adjustment once, at the start of a
		// multi-level gas-dynamic update.
		double p = pA / area;
		double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
		    (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
		double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
		new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
		stagnation_condition.gas.p = new_p0;
		gmodel.update_thermo_from_pT(stagnation_condition.gas);
		stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
		stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
	    }
	    double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
	    // Assume an isentropic process from a known total enthalpy.
	    double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
	    gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
	    // Now, apply the ghost-cell conditions
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    face = src_cell.iface[Face.west];
		    // Velocity components may vary with position on the block face.
		    set_velocity_components(inflow_condition.vel, bulk_speed, face);
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.fs.copy_values_from(inflow_condition);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.fs.copy_values_from(inflow_condition);
		} // end j loop
	    } // end k loop
	    break;
	case Face.top:
	    throw new Error("GhostCellFromStagnation.apply_structured_grid() " ~
			    "not yet implemented for top boundary");
	    // k = blk.kmax;
	    // for (i = blk.imin; i <= blk.imax; ++i) {
	    // 	for (j = blk.jmin; j <= blk.jmax; ++j) {
	    // 	    src_cell = blk.get_cell(i,j,k);
	    // 	    dest_cell = blk.get_cell(i,j,k+1);
	    // 	    dest_cell.fs.gas.p = p_stag; // [TODO] FIX-ME
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	    dest_cell = blk.get_cell(i,j,k+2);
	    // 	    dest_cell.fs.gas.p = p_stag;
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	} // end j loop
	    // } // for i
	    // break;
	case Face.bottom:
	    throw new Error("GhostCellFromStagnation.apply_structured_grid() " ~
			    "not yet implemented for bottom boundary");
	    // k = blk.kmin;
	    // for (i = blk.imin; i <= blk.imax; ++i) {
	    // 	for (j = blk.jmin; j <= blk.jmax; ++j) {
	    // 	    src_cell = blk.get_cell(i,j,k);
	    // 	    dest_cell = blk.get_cell(i,j,k-1);
	    // 	    dest_cell.fs.gas.p = p_stag; // [TODO] FIX-ME
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	    dest_cell = blk.get_cell(i,j,k-2);
	    // 	    dest_cell.fs.gas.p = p_stag;
	    // 	    foreach(ref elem; dest_cell.fs.gas.T) elem = T_stag; 
	    // 	    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
	    // 	} // end j loop
	    // } // for i
	    // break;
	} // end switch
    } // end apply()
} // end class GhostCellFixedStagnationPT

class GhostCellFullFaceExchangeCopy : GhostCellEffect {
public:
    Block neighbourBlock;
    int neighbourFace;
    int neighbourOrientation;
    bool reorient_vector_quantities;
    double[] Rmatrix;

    this(int id, int boundary,
	 int otherBlock, int otherFace, int orient,
	 bool reorient_vector_quantities,
	 ref const(double[]) Rmatrix)
    {
	super(id, boundary, "FullFaceExchangeCopy");
	neighbourBlock = gasBlocks[otherBlock];
	neighbourFace = otherFace;
	neighbourOrientation = orient;
	this.reorient_vector_quantities = reorient_vector_quantities;
	this.Rmatrix = Rmatrix.dup();
    }

    override string toString() const
    { 
	string str = "FullFaceExchangeCopy(otherBlock=" ~ to!string(neighbourBlock.id) ~ 
	    ", otherFace=" ~ to!string(neighbourFace) ~ 
	    ", orient=" ~ to!string(neighbourOrientation) ~
	    ", reorient_vector_quantities=" ~ to!string(reorient_vector_quantities) ~
	    ", Rmatrix=[";
	foreach(i, v; Rmatrix) {
	    str ~= to!string(v);
	    str ~= (i < Rmatrix.length-1) ? ", " : "]";
	}
	str ~= ")";
	return str;
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("GhostCellFullFaceExchangeCopy.apply_unstructured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	blk.copy_into_ghost_cells(which_boundary, 
				  neighbourBlock,
				  neighbourFace, neighbourOrientation,
				  CopyDataOption.all, true,
				  reorient_vector_quantities, Rmatrix);
    }
} // end class GhostCellFullFaceExchangeCopy

class GhostCellMappedCellExchangeCopy : GhostCellEffect {
public:
    // For each ghost cell associated with the boundary,
    // we will have a corresponding "mapped cell" from which we will copy
    // the flow conditions.
    BasicCell[] ghost_cells;
    BasicCell[] mapped_cells;
    // Parameters for the calculation of the mapped-cell location.
    bool transform_position;
    Vector3 c0 = Vector3(0.0, 0.0, 0.0); // default origin
    Vector3 n = Vector3(0.0, 0.0, 1.0); // z-axis
    double alpha = 0.0; // rotation angle (radians) about specified axis vector
    Vector3 delta = Vector3(0.0, 0.0, 0.0); // default zero translation
    bool list_mapped_cells;
    // Parameters for the optional rotation of copied vector data.
    bool reorient_vector_quantities;
    double[] Rmatrix;

    this(int id, int boundary,
	 bool transform_pos,
	 ref const(Vector3) c0, ref const(Vector3) n, double alpha,
	 ref const(Vector3) delta,
	 bool list_mapped_cells,
	 bool reorient_vector_quantities,
	 ref const(double[]) Rmatrix)
    {
	super(id, boundary, "MappedCellExchangeCopy");
	this.transform_position = transform_pos;
	this.c0 = c0;
	this.n = n; this.n.normalize();
	this.alpha = alpha;
	this.delta = delta;
	this.list_mapped_cells = list_mapped_cells;
	this.reorient_vector_quantities = reorient_vector_quantities;
	this.Rmatrix = Rmatrix.dup();
    }

    override string toString() const
    { 
	string str = "MappedCellExchangeCopy(" ~
	    "transform_position=" ~ to!string(transform_position) ~
	    ", c0=" ~ to!string(c0) ~ 
	    ", n=" ~ to!string(n) ~ 
	    ", alpha=" ~ to!string(alpha) ~
	    ", delta=" ~ to!string(delta) ~
	    ", list_mapped_cells=" ~ to!string(list_mapped_cells) ~
	    ", reorient_vector_quantities=" ~ to!string(reorient_vector_quantities) ~
	    ", Rmatrix=[";
	foreach(i, v; Rmatrix) {
	    str ~= to!string(v);
	    str ~= (i < Rmatrix.length-1) ? ", " : "]";
	}
	str ~= ")";
	return str;
    }

    void set_up_cell_mapping()
    {
	// Stage-2 construction for this boundary condition.
	// Needs to be called after the cell geometries have been computed.
	write("GhostCellMappedCellExchangeCopy: ");
	final switch (blk.grid_type) {
	case Grid_t.unstructured_grid: 
	    writeln("Set up mapping to unstructured-grid ghost cells");
	    BoundaryCondition bc = blk.bc[which_boundary];
	    foreach (i, f; bc.faces) {
		if (bc.outsigns[i] == 1) {
		    ghost_cells ~= f.right_cells[0];
		    ghost_cells ~= f.right_cells[1];
		} else {
		    ghost_cells ~= f.left_cells[0];
		    ghost_cells ~= f.left_cells[1];
		}
	    } // end foreach face
	    break;
	case Grid_t.structured_grid:
	    writeln("Set up mapping to structured-grid ghost cells");
	    size_t i, j, k;
	    final switch (which_boundary) {
	    case Face.north:
		j = blk.jmax;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (i = blk.imin; i <= blk.imax; ++i) {
			ghost_cells ~= blk.get_cell(i,j+1,k);
			ghost_cells ~= blk.get_cell(i,j+2,k);
		    } // end i loop
		} // for k
		break;
	    case Face.east:
		i = blk.imax;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			ghost_cells ~= blk.get_cell(i+1,j,k);
			ghost_cells ~= blk.get_cell(i+2,j,k);
		    } // end j loop
		} // for k
		break;
	    case Face.south:
		j = blk.jmin;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (i = blk.imin; i <= blk.imax; ++i) {
			ghost_cells ~= blk.get_cell(i,j-1,k);
			ghost_cells ~= blk.get_cell(i,j-2,k);
		    } // end i loop
		} // for k
		break;
	    case Face.west:
		i = blk.imin;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			ghost_cells ~= blk.get_cell(i-1,j,k);
			ghost_cells ~= blk.get_cell(i-2,j,k);
		    } // end j loop
		} // for k
		break;
	    case Face.top:
		k = blk.kmax;
		for (i = blk.imin; i <= blk.imax; ++i) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			ghost_cells ~= blk.get_cell(i,j,k+1);
			ghost_cells ~= blk.get_cell(i,j,k+2);
		    } // end j loop
		} // for i
		break;
	    case Face.bottom:
		k = blk.kmin;
		for (i = blk.imin; i <= blk.imax; ++i) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			ghost_cells ~= blk.get_cell(i,j,k-1);
			ghost_cells ~= blk.get_cell(i,j,k-2);
		    } // end j loop
		} // for i
		break;
	    } // end switch
	} // end switch blk.grid_type
	// Now that we have a collection of the local ghost cells,
	// locate the corresponding active cell so that we can later
	// copy that cell's flow state.
	if (list_mapped_cells) {
	    writefln("Mapped cells for block[%d] boundary[%d]:", blk.id, which_boundary);
	}
	foreach (mygc; ghost_cells) {
	    Vector3 ghostpos = mygc.pos[0];
	    Vector3 mypos = ghostpos;
	    if (transform_position) {
		Vector3 c1 = c0 + dot(n, (ghostpos - c0)) * n;
		Vector3 t1 = (ghostpos - c1);
		t1.normalize();
		Vector3 t2 = cross(n, t1);
		mypos = c1 + cos(alpha) * t1 + sin(alpha) * t2;
		mypos += delta;
	    }
	    // Because we need to access all of the gas blocks in the following search,
	    // we have to run this set_up_cell_mapping function from a serial loop.
	    // In parallel code, threads other than the main thread get uninitialized
	    // gersions of the gasBlocks array.
	    BasicCell closest_cell = gasBlocks[0].cells[0];
	    Vector3 cellpos = closest_cell.pos[0];
	    double min_distance = abs(cellpos - mypos);
	    foreach (b; gasBlocks) {
		foreach (cell; b.cells) {
		    double distance = abs(cell.pos[0] - mypos);
		    if (distance < min_distance) {
			closest_cell = cell;
			min_distance = distance;
		    }
		}
	    }
	    mapped_cells ~= closest_cell;
	    if (list_mapped_cells) {
		writeln("    ghost-cell-pos=", to!string(mygc.pos[0]), 
			" mapped-cell-pos=", to!string(closest_cell.pos[0]));
	    }
	} // end foreach mygc
    } // end set_up_cell_mapping()

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	foreach (i; 0 .. ghost_cells.length) {
	    ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
	    if (reorient_vector_quantities) {
		ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
	    }
	}
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	foreach (i; 0 .. ghost_cells.length) {
	    ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
	    if (reorient_vector_quantities) {
		ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
	    }
	}
    }
} // end class GhostCellMappedCellExchangeCopy

