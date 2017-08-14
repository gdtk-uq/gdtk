// boundary_interface_effect.d
//
// Effects needed to compute viscous fluxes and the like.
//
// PJ and RG, 2015-04-28, initial code mainly from 
//    the break-up of the Fixed_T boundary condition.
//

module boundary_interface_effect;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;

import geom;
import grid;
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
import solidfvcell;
import solidfvinterface;
import gas_solid_interface;
import bc;

BoundaryInterfaceEffect make_BIE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    // At the point at which we call this function, we may be inside the block-constructor.
    // Don't attempt the use the block-owned gas model.
    auto gmodel = GlobalConfig.gmodel_master; 
    // If we need access to a gas model in here, 
    // be sure to use GlobalConfig.gmodel_master.
    BoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "copy_cell_data":
	newBIE = new BIE_CopyCellData(blk_id, boundary);
	break;
    case "flow_state_copy_to_interface":
	auto flowstate = new FlowState(jsonData["flowstate"], gmodel);
	newBIE = new BIE_FlowStateCopy(blk_id, boundary, flowstate);
	break;
    case "flow_state_copy_from_profile_to_interface":
	string fname = getJSONstring(jsonData, "filename", "");
	string match = getJSONstring(jsonData, "match", "xyz");
	newBIE = new BIE_FlowStateCopyFromProfile(blk_id, boundary, fname, match);
	break;
    case "zero_velocity":
	newBIE = new BIE_ZeroVelocity(blk_id, boundary);
	break;
    case "translating_surface":
	Vector3 v_trans = getJSONVector3(jsonData, "v_trans", Vector3(0.0,0.0,0.0));
	newBIE = new BIE_TranslatingSurface(blk_id, boundary, v_trans);
	break;
    case "rotating_surface":
	Vector3 r_omega = getJSONVector3(jsonData, "r_omega", Vector3(0.0,0.0,0.0));
	Vector3 centre = getJSONVector3(jsonData, "centre", Vector3(0.0,0.0,0.0));
	newBIE = new BIE_RotatingSurface(blk_id, boundary, r_omega, centre);
	break;
    case "fixed_temperature":
	double Twall = getJSONdouble(jsonData, "Twall", 300.0);
	newBIE = new BIE_FixedT(blk_id, boundary, Twall);
	break;
    case "fixed_composition":
	double[] massfAtWall = getJSONdoublearray(jsonData, "wall_massf_composition", [1.0,]);
	newBIE = new BIE_FixedComposition(blk_id, boundary, massfAtWall);
	break;
    case "update_thermo_trans_coeffs":
	newBIE = new BIE_UpdateThermoTransCoeffs(blk_id, boundary);
	break;
    case "wall_k_omega":
	newBIE = new BIE_WallKOmega(blk_id, boundary);
	break;
    case "wall_function_interface_effect":
	string thermalCond = getJSONstring(jsonData, "thermal_condition", "FIXED_T");
	thermalCond = toUpper(thermalCond);
	newBIE = new BIE_WallFunction(blk_id, boundary, thermalCond);
	break;
    case "temperature_from_gas_solid_interface":
	int otherBlock = getJSONint(jsonData, "other_block", -1);
	string otherFaceName = getJSONstring(jsonData, "other_face", "none");
	int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
	newBIE = new BIE_TemperatureFromGasSolidInterface(blk_id, boundary,
							  otherBlock, face_index(otherFaceName),
							  neighbourOrientation);
	break;
    case "user_defined":
     	string fname = getJSONstring(jsonData, "filename", "none");
	newBIE = new BIE_UserDefined(blk_id, boundary, fname);
	break;
    default:
	string errMsg = format("ERROR: The BoundaryInterfaceEffect type: '%s' is unknown.", bieType);
	throw new FlowSolverException(errMsg);
    }
    return newBIE;
}


class BoundaryInterfaceEffect {
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
	return "BoundaryInterfaceEffect()";
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
} // end class BoundaryInterfaceEffect


class BIE_CopyCellData : BoundaryInterfaceEffect {
    this(int id, int boundary, double Twall=300.0)
    {
	super(id, boundary, "CopyCellData");
    }

    override string toString() const 
    {
	return "CopyCellData()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		f.fs.copy_values_from(f.left_cell.fs);
	    } else {
		f.fs.copy_values_from(f.right_cell.fs);
	    }
	} // end foreach face
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_CopyCellData

class BIE_FlowStateCopy : BoundaryInterfaceEffect {

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
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    f.fs.copy_values_from(fstate);
	} // end foreach face
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(fstate);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(fstate);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(fstate);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(fstate);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(fstate);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(fstate);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()

private:
    FlowState fstate;

} // end class BIE_FlowStateCopy


class BIE_FlowStateCopyFromProfile : BoundaryInterfaceEffect {
public:
    this(int id, int boundary, string fileName, string match)
    {
	super(id, boundary, "flowStateCopyFromProfile");
	fprofile = new FlowProfile(fileName, match);
    }

    override string toString() const 
    {
	return format("flowStateCopyFromProfile(filename=\"%s\", match=\"%s\")",
		      fprofile.fileName, fprofile.posMatch);
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
	    fprofile.adjust_velocity(f.fs, f.pos);
	}
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface f;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    f = cell.iface[Face.north];
		    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
		    fprofile.adjust_velocity(f.fs, f.pos);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    f = cell.iface[Face.east];
		    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
		    fprofile.adjust_velocity(f.fs, f.pos);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    f = cell.iface[Face.south];
		    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
		    fprofile.adjust_velocity(f.fs, f.pos);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    f = cell.iface[Face.west];
		    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
		    fprofile.adjust_velocity(f.fs, f.pos);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    f = cell.iface[Face.top];
		    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
		    fprofile.adjust_velocity(f.fs, f.pos);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    f = cell.iface[Face.bottom];
		    f.fs.copy_values_from(fprofile.get_flowstate(f.id, f.pos));
		    fprofile.adjust_velocity(f.fs, f.pos);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()

private:
    FlowProfile fprofile;

} // end class BIE_FlowStateCopyFromProfile


class BIE_ZeroVelocity : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
	super(id, boundary, "ZeroVelocity");
    }

    override string toString() const 
    {
	return "ZeroVelocity()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
      	auto gmodel = blk.myConfig.gmodel;	
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    FlowState fs = f.fs;
	    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
	} // end foreach face
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_ZeroVelocity


class BIE_TranslatingSurface : BoundaryInterfaceEffect {
    // The boundary surface is translating with fixed velocity v_trans.
    Vector3 v_trans;

    this(int id, int boundary, Vector3 v_trans)
    {
	super(id, boundary, "TranslatingSurface");
	this.v_trans = v_trans;
    }

    override string toString() const 
    {
	return "TranslatingSurface(v_trans=" ~ to!string(v_trans) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("BIE_TranslatingSurface.apply_unstructured_grid() not implemented yet");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = v_trans.x; fs.vel.refy = v_trans.y; fs.vel.refz = v_trans.z;
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = v_trans.x; fs.vel.refy = v_trans.y; fs.vel.refz = v_trans.z;
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = v_trans.x; fs.vel.refy = v_trans.y; fs.vel.refz = v_trans.z;
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = v_trans.x; fs.vel.refy = v_trans.y; fs.vel.refz = v_trans.z;
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = v_trans.x; fs.vel.refy = v_trans.y; fs.vel.refz = v_trans.z;
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = v_trans.x; fs.vel.refy = v_trans.y; fs.vel.refz = v_trans.z;
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_TranslatingSurface


class BIE_RotatingSurface : BoundaryInterfaceEffect {
    // The boundary surface is rotating with fixed angular velocity r_omega
    // about centre.
    Vector3 r_omega;
    Vector3 centre;

    this(int id, int boundary, Vector3 r_omega, Vector3 centre)
    {
	super(id, boundary, "RotatingSurface");
	this.r_omega = r_omega;
	this.centre = centre;
    }

    override string toString() const 
    {
	return "RotatingSurface(r_omega=" ~ to!string(r_omega) ~ 
	    ", centre=" ~ to!string(centre) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("BIE_RotatingSurface.apply_unstructured_grid() not implemented yet");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.vel = cross(r_omega, IFace.pos-centre);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.vel = cross(r_omega, IFace.pos-centre);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.vel = cross(r_omega, IFace.pos-centre);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.vel = cross(r_omega, IFace.pos-centre);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.vel = cross(r_omega, IFace.pos-centre);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.vel = cross(r_omega, IFace.pos-centre);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_RotatingSurface


class BIE_FixedT : BoundaryInterfaceEffect {
public:
    double Twall;

    this(int id, int boundary, double Twall)
    {
	super(id, boundary, "FixedT");
	this.Twall = Twall;
    }

    override string toString() const 
    {
	return "FixedT(Twall=" ~ to!string(Twall) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	auto gmodel = blk.myConfig.gmodel;	
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    FlowState fs = f.fs;
	    fs.gas.Ttr = Twall;
	    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
	} // end foreach face
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.gas.Ttr = Twall;
		    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.gas.Ttr = Twall;
		    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.gas.Ttr = Twall;
		    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.gas.Ttr = Twall;
		    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.gas.Ttr = Twall;
		    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.gas.Ttr = Twall;
		    foreach(ref elem; fs.gas.T_modes) { elem = Twall; }
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_FixedT

class BIE_FixedComposition : BoundaryInterfaceEffect {
public:    
    double[] massfAtWall;
    
    this(int id, int boundary, double[] massfAtWall)
    {
	super(id, boundary, "FixedComposition");
	this.massfAtWall = massfAtWall;
    }

    override string toString() const 
    {
	return "FixedComposition(massfAtWall=" ~ to!string(massfAtWall) ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	auto gmodel = blk.myConfig.gmodel;	
	int nsp = gmodel.n_species;
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    FlowState fs = f.fs;
	    for(int isp; isp<nsp; isp++) {
                fs.gas.massf[isp] = massfAtWall[isp];	
	    }	
	} // end foreach face
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;
	int nsp = gmodel.n_species;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    for(int isp; isp<nsp; isp++) {
                	fs.gas.massf[isp] = massfAtWall[isp];	
	    	    }
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    for(int isp; isp<nsp; isp++) {
			fs.gas.massf[isp] = massfAtWall[isp];	
		    }
                } // en for j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    for(int isp; isp<nsp; isp++) {
                	fs.gas.massf[isp] = massfAtWall[isp];	
	    	    }
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    for(int isp; isp<nsp; isp++) {
                	fs.gas.massf[isp] = massfAtWall[isp];	
	    	    }
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    for(int isp; isp<nsp; isp++) {
                	fs.gas.massf[isp] = massfAtWall[isp];	
	    	    }
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    for(int isp; isp<nsp; isp++) {
                	fs.gas.massf[isp] = massfAtWall[isp];	
	    	    }
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_FixedComposition


class BIE_UpdateThermoTransCoeffs : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
	super(id, boundary, "UpdateThermoTransCoeffs");
    }

    override string toString() const 
    {
	return "UpdateThermoTransCoeffs()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BoundaryCondition bc = blk.bc[which_boundary];
	auto gmodel = blk.myConfig.gmodel;
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		FlowState fs = f.fs;
		gmodel.update_thermo_from_pT(fs.gas);
		gmodel.update_trans_coeffs(fs.gas);
	    } else {
		FlowState fs = f.fs;
		gmodel.update_thermo_from_pT(fs.gas);
		gmodel.update_trans_coeffs(fs.gas);
	    }
	} // end foreach face
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_UpdateThermoTransCoeffs

class BIE_WallKOmega : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
	super(id, boundary, "WallKOmega");
    }

    override string toString() const 
    {
	return "WallKOmega()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		f.fs.omega = ideal_omega_at_wall(f.left_cell);
	    } else {
		f.fs.omega = ideal_omega_at_wall(f.right_cell);
	    }
	    f.fs.tke = 0.0;
	} // end foreach face
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()

    @nogc
    double ideal_omega_at_wall(in FVCell cell)
    // As recommended by Wilson Chan, we use Menter's correction
    // for omega values at the wall. This appears as Eqn A12 in 
    // Menter's paper.
    // Reference:
    // Menter (1994)
    // Two-Equation Eddy-Viscosity Turbulence Models for
    // Engineering Applications.
    // AIAA Journal, 32:8, pp. 1598--1605
    {
	auto wall_gas = cell.cell_at_nearest_wall.fs.gas;
	double d0 = cell.half_cell_width_at_wall;
	double nu = wall_gas.mu / wall_gas.rho;
	double beta1 = 0.075;
	return 10 * (6 * nu) / (beta1 * d0 * d0);
    }

    @nogc
    double ideal_omega(in FVCell cell)
    {
	double d0 = cell.half_cell_width_at_wall;
	double d = cell.distance_to_nearest_wall;
	return ideal_omega_at_wall(cell) * (d0 * d0) / ((d0 + d) * (d0 + d));
    }
} // end class BIE_WallKOmega

class BIE_WallFunction : BoundaryInterfaceEffect {
    this(int id, int boundary, string thermalCond)
    {
	super(id, boundary, "WallFunction_InterfaceEffect");
	_isFixedTWall = (thermalCond == "FIXED_T") ? true : false;
	_faces_need_to_be_flagged = true;
    }

    override string toString() const 
    {
	return "WallFunction_InterfaceEffect()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new FlowSolverException("WallFunction_InterfaceEffect bc not implemented for unstructured grids.");
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	if (_faces_need_to_be_flagged) {
	    // Flag faces, just once.
	    final switch (which_boundary) {
	    case Face.north:
		j = blk.jmax;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (i = blk.imin; i <= blk.imax; ++i) {
			blk.get_cell(i,j,k).iface[Face.north].use_wall_function_shear_and_heat_flux = true;
		    }
		}
		break;
	    case Face.east:
		i = blk.imax;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			blk.get_cell(i,j,k).iface[Face.east].use_wall_function_shear_and_heat_flux = true;
		    }
		}
		break;
	    case Face.south:
		j = blk.jmin;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (i = blk.imin; i <= blk.imax; ++i) {
			blk.get_cell(i,j,k).iface[Face.south].use_wall_function_shear_and_heat_flux = true;
		    }
		}
		break;
	    case Face.west:
		i = blk.imin;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			blk.get_cell(i,j,k).iface[Face.west].use_wall_function_shear_and_heat_flux = true;
		    }
		}
		break;
	    case Face.top:
		k = blk.kmax;
		for (i = blk.imin; i <= blk.imax; ++i) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			blk.get_cell(i,j,k).iface[Face.top].use_wall_function_shear_and_heat_flux = true;
		    }
		}
		break;
	    case Face.bottom:
		k = blk.kmin;
		for (i = blk.imin; i <= blk.imax; ++i) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			blk.get_cell(i,j,k).iface[Face.bottom].use_wall_function_shear_and_heat_flux = true;
		    }
		}
		break;
	    } // end switch which_boundary
	    _faces_need_to_be_flagged = false;
	} // end if _faces_need_to_be_flagged
	//
	// Do some real work.
	//
	FVCell cell, second_cell_from_wall;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    second_cell_from_wall = blk.get_cell(i,j-1,k);
		    IFace = cell.iface[Face.north];
	            wall_function(cell, second_cell_from_wall, IFace); 
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
                    second_cell_from_wall = blk.get_cell(i-1,j,k);
		    IFace = cell.iface[Face.east];
                    wall_function(cell, second_cell_from_wall, IFace);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
                    second_cell_from_wall = blk.get_cell(i,j+1,k);
		    IFace = cell.iface[Face.south];
                    wall_function(cell, second_cell_from_wall, IFace);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
                    second_cell_from_wall = blk.get_cell(i+1,j,k);
		    IFace = cell.iface[Face.west];
                    wall_function(cell, second_cell_from_wall, IFace);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
                    second_cell_from_wall = blk.get_cell(i,j,k-1);
		    IFace = cell.iface[Face.top];
                    wall_function(cell, second_cell_from_wall, IFace);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
                    second_cell_from_wall = blk.get_cell(i,j,k-1);
		    IFace = cell.iface[Face.bottom];
                    wall_function(cell, second_cell_from_wall, IFace);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()

    void wall_function(const FVCell cell, const FVCell second_cell_from_wall, FVInterface IFace)
    // Implement Nichols' and Nelson's wall function boundary condition
    // Reference:
    //  Nichols RH & Nelson CC (2004)
    //  Wall Function Boundary Conditions Inclding Heat Transfer
    //  and Compressibility.
    //  AIAA Journal, 42:6, pp. 1107--1114
    // NOTE: IFace.fs will receive updated values of tke and omega for later copying
    //       to boundary cells.
    {
	auto gmodel = blk.myConfig.gmodel; 
	//
	// Compute recovery factor
	double cp = gmodel.Cp(cell.fs.gas); 
	double Pr = cell.fs.gas.mu * cp / cell.fs.gas.k;
	double gas_constant = gmodel.R(cell.fs.gas);
	double recovery = pow(Pr, (1.0/3.0));
	// Compute tangent velocity at nearest interior point and wall interface
	Vector3 cellVel = cell.fs.vel;
	Vector3 faceVel = IFace.fs.vel;
        cellVel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2); 
	faceVel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	double cell_tangent = sqrt( pow(cellVel.y, 2.0) + pow(cellVel.z, 2.0) );
	double face_tangent = sqrt( pow(faceVel.y, 2.0) + pow(faceVel.z, 2.0) );
	//TODO 3D formulation (not implemented)
	//double cell_tangent0 = cellVel.y;
	//double cell_tangent1 = cellVel.z;
	//double face_tangent0 = faceVel.y;
	//double face_tangent1 = faceVel.z;
	//double vt1_2_angle = atan2(fabs(cellVel.z - faceVel.z), fabs(cellVel.y - faceVel.y));
	double du = fabs(cell_tangent - face_tangent);
	//TODO 3D formulation (not implemented)
	//double du = sqrt( pow((cell_tangent0-face_tangent0),2.0) + pow((cell_tangent1-face_tangent1),2.0) );
	// Compute wall gas properties from either ...
	double T_wall, rho_wall;
	if ( _isFixedTWall ) {
	    // ... user-specified wall temperature, or ... 
	    T_wall = IFace.fs.gas.Ttr; 
	    rho_wall = IFace.fs.gas.rho;
	} else {
	    // ... using Crocco-Busemann relation (Eq 11) 
	    T_wall = cell.fs.gas.Ttr + recovery * du * du / (2.0 * cp);
	    // Update gas properties at the wall, assuming static pressure
	    // at the wall is the same as that in the first wall cell
	    IFace.fs.gas.Ttr = T_wall;
	    IFace.fs.gas.p = cell.fs.gas.p;
	    gmodel.update_thermo_from_pT(IFace.fs.gas);
	    gmodel.update_trans_coeffs(IFace.fs.gas);
	    rho_wall = IFace.fs.gas.rho;
	}
	// Compute wall shear stess (and heat flux for fixed temperature wall) 
	// using the surface stress tensor. This provides initial values to solve
	// for tau_wall and q_wall iteratively
	double wall_dist = cell.half_cell_width_at_wall;
	double dudy = du / wall_dist;
	double mu_lam_wall = IFace.fs.gas.mu;
	double mu_lam = cell.fs.gas.mu;
	double tau_wall_old = mu_lam_wall * dudy;
	double dT, dTdy, k_lam_wall, q_wall_old;
	if ( _isFixedTWall ) {
	    dT = fabs( cell.fs.gas.Ttr - IFace.fs.gas.Ttr );
	    dTdy = dT / wall_dist;
	    k_lam_wall = IFace.fs.gas.k;
	    q_wall_old = k_lam_wall * dTdy;
	}
	// Constants from Spalding's Law of the Wall theory and 
	// Nichols' wall function implementation
	double kappa = 0.4;
	double B = 5.5;
	double C_mu = 0.09;
	size_t counter = 0; 
	double tolerance = 1.0e-10;
	double diff_tau = 100.0;
	double diff_q;
	if ( _isFixedTWall ) {
	    diff_q = 100.0;
	} else {
	    // Make diff_q smaller than tolerance so that while loop
	    // will only check for diff_tau for adiabatic wall cases
	    diff_q = tolerance / 10.0; 
	}
	double tau_wall = 0.0; double q_wall = 0.0;
	double u_tau = 0.0; double u_plus = 0.0;
	double Gam = 0.0; double Beta = 0.0; double Q = 0.0; double Phi = 0.0;
	double y_plus_white = 0.0; double y_plus = 0.0;
	//
	// Iteratively solve for the wall-function corrected shear stress
	// (and heat flux for fixed temperature walls)
	while ( diff_tau > tolerance || diff_q > tolerance) {
	    // Friction velocity and u+ (Eq 2)
	    u_tau = sqrt( tau_wall_old / rho_wall );
	    u_plus = du / u_tau;
	    // Gamma, Beta, Qm and Phi (Eq 7)
	    Gam = recovery * u_tau * u_tau / (2.0 * cp * T_wall); 
	    if (_isFixedTWall) {
		Beta = q_wall_old * mu_lam_wall / (rho_wall*T_wall*k_lam_wall*u_tau);
	    } else {
		Beta = 0.0;
	    }
	    Q = sqrt(Beta*Beta + 4.0*Gam);
	    Phi = asin(-1.0 * Beta / Q);
	    // y+ defined by White and Christoph (Eq 9)
	    y_plus_white = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*u_plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);
	    // Spalding's unified form for defining y+ (Eq 8)
	    y_plus = u_plus + y_plus_white - exp(-1.0*kappa*B) * ( 1.0 + kappa*u_plus
								       + pow((kappa*u_plus), 2.0) / 2.0
								       + pow((kappa*u_plus), 3.0) / 6.0 );
	    // Calculate an updated value for the wall shear stress and heat flux 
	    tau_wall = 1.0/rho_wall * pow(y_plus*mu_lam_wall/wall_dist, 2.0);
	    if (_isFixedTWall) {
		Beta = (cell.fs.gas.Ttr/T_wall - 1.0 + Gam*u_plus*u_plus) / u_plus;
		q_wall = Beta * (rho_wall*T_wall*k_lam_wall*u_tau) / mu_lam_wall;
	    }
	    // Difference between old and new tau_wall and q_wall. Update old value
	    diff_tau = fabs(tau_wall - tau_wall_old);
	    tau_wall_old += 0.25 * (tau_wall - tau_wall_old);
	    if ( _isFixedTWall ) {
		diff_q = fabs( q_wall - q_wall_old );
		q_wall_old += 0.25 * (q_wall - q_wall_old);
	    }
	    // Set counter to limit number of times we iterate this loop
	    counter++;
	    if (counter > 500) break;
	}
	//
	// Store wall shear stress and heat flux to be used later to replace viscous stress 
	// in flux calculations. Also, for wall shear stress, transform value back to the
	// global frame of reference.
	double reverse_flag = 1.0;
	if ( face_tangent > cell_tangent ) reverse_flag = -1.0;
	//TODO 3D formulation (not implemented)
	//double reverse_flag0 = 1.0;
	//double reverse_flag1 = 1.0;
	//if ( face_tangent0  > cell_tangent0 ) reverse_flag0 = -1.0;\
	//if ( face_tangent1  > cell_tangent1 ) reverse_flag1 = -1.0;
	double tau_wall_x, tau_wall_y, tau_wall_z;
	if ( IFace.bc_id == Face.north || IFace.bc_id == Face.east || IFace.bc_id == Face.top ) {
	    // TODO 2D formulation (to be integrated with 3D)
	    IFace.tau_wall_x = -1.0 * reverse_flag * tau_wall * IFace.n.y;
	    IFace.tau_wall_y = -1.0 * reverse_flag * tau_wall * IFace.n.x;
	    IFace.tau_wall_z = 0.0;
	    //TODO 3D formulation (not implemented)
	    //local_tau_wall_x = 0.0
	    //local_tau_wall_y = -1.0 * reverse_flag0 * tau_wall * cos(vt1_2_angle);
	    //local_tau_wall_z = -1.0 * reverse_flag1 * tau_wall * sin(vt1_2_angle);
	    if (_isFixedTWall) {
		IFace.q = -1.0 * q_wall;
	    } else {
		IFace.q = 0.0;
	    }
	} else { // South, West and Bottom
	    // 2D formulation (to be integrated with 3D)
	    IFace.tau_wall_x = reverse_flag * tau_wall * IFace.n.y;
	    IFace.tau_wall_y = reverse_flag * tau_wall * IFace.n.x;
            IFace.tau_wall_z = 0.0;
	    //TODO 3D formulation (not implemented)
	    //local_tau_wall_x = 0.0
	    //local_tau_wall_y = reverse_flag0 * tau_wall * cos(vt1_2_angle);
	    //local_tau_wall_z = reverse_flag1 * tau_wall * sin(vt1_2_angle);
	    if ( _isFixedTWall ) {
		IFace.q = q_wall;
            } else {
		IFace.q = 0.0;
            }
	}
	//TODO 3D formulation (not implemented)
	//IFace.tau.transform_to_global(IFace.n, IFace.t1, IFace.t2);
	//
	// Turbulence model boundary conditions (Eq 15 & 14)
	double y_white_y_plus = 2.0 * y_plus_white * kappa*sqrt(Gam)/Q
		* pow((1.0 - pow(2.0*Gam*u_plus-Beta,2.0)/(Q*Q)), 0.5);
	double mu_coeff = 1.0 + y_white_y_plus
		- kappa*exp(-1.0*kappa*B) * (1.0 + kappa*u_plus + kappa*u_plus*kappa*u_plus/2.0)
		- mu_lam/mu_lam_wall;
	double mu_t = mu_lam_wall * mu_coeff;
	// omega (Eq 19 - 21)
	double omega_i = 6.0*mu_lam_wall / (0.075*rho_wall*wall_dist*wall_dist);
	double omega_o = u_tau / (sqrt(C_mu)*kappa*wall_dist);
	double omega = sqrt(omega_i*omega_i + omega_o*omega_o);
	// tke (Eq 22)
	double tke = omega * mu_t / cell.fs.gas.rho;
	// tke may be negative at initial stages. When this happens, we just use a value of
	// tke that is scaled based off the squared of the ratio of distance between the
	// first and second cells from the wall <- this was implemented by Jason Qin in E3.
	if ( tke < 0.0 || isNaN(tke) != 0 ) { 
	    double wall_dist2 = second_cell_from_wall.distance_to_nearest_wall;
	    tke = wall_dist*wall_dist / (wall_dist2*wall_dist2) * second_cell_from_wall.fs.tke;
	}
	IFace.fs.tke = tke;
	IFace.fs.omega = omega;
	return;
    } // end wall_function()

private:
    bool _isFixedTWall;
    bool _faces_need_to_be_flagged = true;
} // end class BIE_WallFunction


// NOTE: This GAS DOMAIN boundary effect has a large
//       and important side-effect:
//       IT ALSO SETS THE FLUX IN THE ADJACENT SOLID DOMAIN
//       AT THE TIME IT IS CALLED.
// TODO: We need to work out a way to coordinate this 
//       interface effect (ie. the setting of temperature)
//       with the flux effect. Ideally, we only want to compute
//       the temperature/flux once per update. This will require
//       some storage at the block level, or in the in the
//       gas/solid interface module since we can't share information
//       (easily) between different types of boundary condition
//       objects. We need to store the energy flux somewhere where it
//       so that we can use it again in the boundary flux effect.
//       It's no good storing the flux in the interface object since
//       that will be changed during the diffusive flux calculation.
 
class BIE_TemperatureFromGasSolidInterface : BoundaryInterfaceEffect {
public:
    int neighbourSolidBlk;
    int neighbourSolidFace;
    int neighbourOrientation;

    this(int id, int boundary, 
	 int otherBlock, int otherFace, int orient)
    {
	super(id, boundary, "TemperatureFromGasSolidInterface");
	neighbourSolidBlk = otherBlock;
	neighbourSolidFace = otherFace;
	neighbourOrientation = orient;
    }

    override string toString() const 
    {
	return "TemperatureFromGasSolidInterface()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	throw new Error("BIE_TemperatureFromGasSolidInterface.apply_unstructured_grid() not implemented yet");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	double kS = solidBlocks[neighbourSolidBlk].sp.k;
	computeFluxesAndTemperatures(ftl, kS,
				     _gasCells, _gasIFaces,
				     _solidCells, _solidIFaces);
    }

private:
    // Some private working arrays.
    // We'll pack data into these can pass out
    // to a routine that can compute the flux and
    // temperatures that balance at the interface.
    FVCell[] _gasCells;
    FVInterface[] _gasIFaces;
    SolidFVCell[] _solidCells;
    SolidFVInterface[] _solidIFaces;

public:
    void initSolidCellsAndIFaces()
    {
	size_t i, j, k;
	auto blk = solidBlocks[neighbourSolidBlk];
	switch ( neighbourSolidFace ) {
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    _solidCells ~= blk.getCell(i, j, k);
		    _solidIFaces ~= _solidCells[$-1].iface[Face.south];
		}
	    }
	    break;
	default:
	    throw new Error("initSolidCellsAndIFaces() only implemented for SOUTH face.");
	}
    }

    void initGasCellsAndIFaces()
    {
	size_t i, j, k;
	switch ( which_boundary ) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    _gasCells ~= blk.get_cell(i, j, k);
		    _gasIFaces ~= _gasCells[$-1].iface[Face.north];
		}
	    }
	    break;
	default:
	    throw new Error("initGasCellsAndIFaces() only implemented for NORTH gas face.");
	}
    }


} // end class BIE_TemperatureFromGasSolidInterface
