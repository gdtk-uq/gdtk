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
    case "update_thermo_trans_coeffs":
	newBIE = new BIE_UpdateThermoTransCoeffs(blk_id, boundary);
	break;
    case "wall_k_omega":
	newBIE = new BIE_WallKOmega(blk_id, boundary);
	break;
    case "wall_function":
	newBIE = new BIE_WallFunction(blk_id, boundary);
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
	throw new Exception(errMsg);
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
    this(int id, int boundary)
    {
	super(id, boundary, "WallFunction");
    }

    override string toString() const 
    {
	return "WallFunction()";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    // TODO: Wilson please fill me in.
	    //       In fact, you could leave the unstructured version unimplemented for the moment.
	    throw new FlowSolverException("WallFunction bc not implemented for unstructured grids.");
	} // end foreach face
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto gmodel = blk.myConfig.gmodel;

	// TODO: Wilson please fill in details for each boundary.
	throw new FlowSolverException("WallFunction bc not implemented for structured grids.");

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;

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

		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()

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
