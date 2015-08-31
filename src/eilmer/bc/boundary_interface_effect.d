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

BoundaryInterfaceEffect make_BIE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    // If we need access to a gas model in here, 
    // be sure to use GlobalConfig.gmodel_master.
    BoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "copy_cell_data":
	newBIE = new BIE_CopyCellData(blk_id, boundary);
	break;
    case "zero_velocity":
	newBIE = new BIE_ZeroVelocity(blk_id, boundary);
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
	return "BoundaryInterfaceEffect()";
    }
    abstract void apply(double t, int gtl, int ftl);
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

    override void apply(double t, int gtl, int ftl)
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


class BIE_ZeroVelocity : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
	super(id, boundary, "ZeroVelocity");
    }

    override string toString() const 
    {
	return "ZeroVelocity()";
    }

    override void apply(double t, int gtl, int ftl)
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

    override void apply(double t, int gtl, int ftl)
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
		    foreach(ref elem; fs.gas.T) elem = Twall;
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
		    foreach(ref elem; fs.gas.T) elem = Twall;
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
		    foreach(ref elem; fs.gas.T) elem = Twall;
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
		    foreach(ref elem; fs.gas.T) elem = Twall;
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
		    foreach(ref elem; fs.gas.T) elem = Twall;
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
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_FixedT


class BIE_UpdateThermoTransCoeffs : BoundaryInterfaceEffect {
    this(int id, int boundary, double Twall=300.0)
    {
	super(id, boundary, "UpdateThermoTransCoeffs");
    }

    override string toString() const 
    {
	return "UpdateThermoTransCoeffs()";
    }

    override void apply(double t, int gtl, int ftl)
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
    // Menter's slightly-rough-surface boundary condition is described
    // in Wilcox's 2006 text, eqn 7.36.
    // For low-resolution grids, the k-omega model is reported to over-estimate
    // the magnitude of omega, well out into the boundary layer so,
    // to get reasonable values for omega close to the wall, we propagate
    // the 1/y**2 form of the omega data out a few cells from the wall.
    this(int id, int boundary)
    {
	super(id, boundary, "WallKOmega");
    }

    override string toString() const 
    {
	return "WallKOmega()";
    }

    override void apply(double t, int gtl, int ftl)
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
    {
	auto wall_gas = cell.cell_at_nearest_wall.fs.gas;
	double d0 = cell.half_cell_width_at_wall;
	return 400.0 * wall_gas.mu / wall_gas.rho / (d0 * d0);
    }

    @nogc
    double ideal_omega(in FVCell cell)
    {
	double d0 = cell.half_cell_width_at_wall;
	double d = cell.distance_to_nearest_wall;
	return ideal_omega_at_wall(cell) * (d0 * d0) / ((d0 + d) * (d0 + d));
    }
} // end class BIE_WallKOmega

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

    override void apply(double t, int gtl, int ftl)
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
