// ghost_cell_effect.d
//
// RG & PJ 2015-12-03 : first hack

import std.json;
import std.string;
import std.conv;
import std.stdio;

import geom;
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

@nogc
void reflect_normal_velocity(ref FVCell cell, in FVInterface IFace)
// Reflects normal velocity with respect to the cell interface.
//
// The process is to rotate the velocity vector into the local frame of
// the interface, negate the normal (local x-component) velocity and
// rotate back to the global frame.
{
    cell.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    cell.fs.vel.refx = -(cell.fs.vel.x);
    cell.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

@nogc
void reflect_normal_magnetic_field(ref FVCell cell, in FVInterface IFace)
{
    cell.fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    cell.fs.B.refx = -(cell.fs.B.x);
    cell.fs.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
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
    case "fixed_pressure_temperature":
	double POut = getJSONdouble(jsonData, "p_out", 1.0e5);
	double TOut = getJSONdouble(jsonData, "T_out", 300.0);
	newGCE = new GhostCellFixedPT(blk_id, boundary, POut, TOut);
	break;
    case "full_face_exchange_copy":
	int otherBlock = getJSONint(jsonData, "other_block", -1);
	string otherFaceName = getJSONstring(jsonData, "other_face", "none");
	
	int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
	newGCE = new GhostCellFullFaceExchangeCopy(blk_id, boundary,
						   otherBlock, face_index(otherFaceName),
						   neighbourOrientation);
	break;
    case "mapped_cell_exchange_copy":
	int[][] mapped_cells;
	newGCE = new GhostCellMappedCellExchangeCopy(blk_id, boundary, mapped_cells);
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
    SBlock blk;
    int which_boundary;
    string type;

    this(int id, int boundary, string _type)
    {
	blk = gasBlocks[id];
	which_boundary = boundary;
	type = _type;
    }
    override string toString() const
    {
	return "GhostCellEffect()";
    }
    abstract void apply(double t, int gtl, int ftl);
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

    override void apply(double t, int gtl, int ftl)
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
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j-1,k);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
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
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i-1,j,k);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
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
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j+1,k);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
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
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i+1,j,k);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
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
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j,k-1);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
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
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j,k+1);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end j loop
	    } // for i
	    break;
	} // end switch which_boundary
    } // end apply()
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

    override void apply(double t, int gtl, int ftl)
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
    } // end apply()

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
    
    override void apply(double t, int gtl, int ftl)
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
    } // end apply()
} // end class GhostCellExtrapolateCopy

class GhostCellFixedPT : GhostCellEffect {
public:
    double Pout;
    double Tout;
    
    this(int id, int boundary, double Pout, double Tout)
    {
	super(id, boundary, "FixedPT");
	this.Pout = Pout;
	this.Tout = Tout;
    }

    override string toString() const 
    {
	return "FixedPT(p_out=" ~ to!string(Pout) ~ ", T_out=" ~ to!string(Tout) ~")";
    }
    
    override void apply(double t, int gtl, int ftl)
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
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
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
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
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
		    dest_cell.fs.gas.p = Pout;
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
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
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
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
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
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
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.fs.gas.p = Pout;
		    foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply()
} // end class GhostCellFixedPT

class GhostCellFullFaceExchangeCopy : GhostCellEffect {
public:
    SBlock neighbourBlock;
    int neighbourFace;
    int neighbourOrientation;

    this(int id, int boundary,
	 int otherBlock, int otherFace, int orient)
    {
	super(id, boundary, "FullFaceExchangeCopy");
	neighbourBlock = gasBlocks[otherBlock];
	neighbourFace = otherFace;
	neighbourOrientation = orient;
    }

    override string toString() const
    { 
	return "FullFaceExchangeCopy(otherBlock=" ~ to!string(neighbourBlock.id) ~ 
	    ", otherFace=" ~ to!string(neighbourFace) ~ 
	    ", orient=" ~ to!string(neighbourOrientation) ~ ")";
    }

    override void apply(double t, int gtl, int ftl)
    {
	// TODO Check me!  This is a work-around.
	// We should be able to directly reference the BCs block as blk.
	blk.copy_into_ghost_cells(which_boundary, 
				  neighbourBlock,
				  neighbourFace, neighbourOrientation,
				  CopyDataOption.all, true);
    }
} // end class GhostCellFullFaceExchangeCopy

class GhostCellMappedCellExchangeCopy : GhostCellEffect {
public:
    // For each ghost cell associated with the boundary.
    // the following data structure stores 4 integers,
    // specifying the mapped-cell block and its ijk-indices,
    int[][] mapped_cells;
    // These data are (i,j,k)-triples indexed by [other_block][other_face]
    // and are used by the distributed-memory mapped-cell copying functions
    // to marshall and send the requested data.
    int[][][] incoming_mapped_cells; 
    int[][][] outgoing_mapped_cells;

    this(int id, int boundary, int[][] mapped_cells)
    {
	super(id, boundary, "MappedCellExchangeCopy");
	assert(false, "Not implemented yet");
    }

    override string toString() const
    { 
	return "MappedCellExchangeCopy(" ~ ")";
    }

    override void apply(double t, int gtl, int ftl)
    {
    }
} // end class GhostCellMappedCellExchangeCopy

