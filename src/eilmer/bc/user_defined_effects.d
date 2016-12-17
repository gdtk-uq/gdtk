// user_defined_effects.d
//
// Authors: RG & PJ
// Date: 2015-03-14

import std.string;
import std.stdio;
import util.lua;
import util.lua_service;
import gas.gas_model;
import gas.luagas_model;

import geom;
import simcore;
import flowstate;
import fvcore;
import fvcell;
import fvinterface;
import globalconfig;
import globaldata;
import ghost_cell_effect;
import boundary_interface_effect;
import luaflowstate;
import lua_helper;
import bc;

class UserDefinedGhostCell : GhostCellEffect {
public:
    string luafname;
    this(int id, int boundary, string fname)
    {
	super(id, boundary, "UserDefined");
	luafname = fname;
    }
    override void post_bc_construction()
    {
	if (blk.bc[which_boundary].myL == null) {
	    blk.bc[which_boundary].myL = luaL_newstate();
	    auto L = blk.bc[which_boundary].myL;
	    luaL_openlibs(L);
	    lua_pushinteger(L, blk.id); lua_setglobal(L, "blkId");
	    registerGasModel(L, LUA_GLOBALSINDEX);
	    pushObj!(GasModel, GasModelMT)(L, blk.myConfig.gmodel);
	    lua_setglobal(L, "gmodel");
	    lua_pushinteger(L, blk.myConfig.gmodel.n_species);
	    lua_setglobal(L, "n_species");
	    lua_pushinteger(L, blk.myConfig.gmodel.n_modes);
	    lua_setglobal(L, "n_modes");
	    lua_pushinteger(L, blk.nicell); lua_setglobal(L, "nicell");
	    lua_pushinteger(L, blk.njcell); lua_setglobal(L, "njcell");
	    lua_pushinteger(L, blk.nkcell); lua_setglobal(L, "nkcell");
	    lua_pushinteger(L, Face.north); lua_setglobal(L, "north");
	    lua_pushinteger(L, Face.east); lua_setglobal(L, "east");
	    lua_pushinteger(L, Face.south); lua_setglobal(L, "south");
	    lua_pushinteger(L, Face.west); lua_setglobal(L, "west");
	    lua_pushinteger(L, Face.top); lua_setglobal(L, "top");
	    lua_pushinteger(L, Face.bottom); lua_setglobal(L, "bottom");
	    lua_pushcfunction(L, &luafn_sampleFlow); lua_setglobal(L, "sampleFlow");
	}
	luaL_dofile(blk.bc[which_boundary].myL, luafname.toStringz);
    }
    override string toString() const
    {
	return "UserDefinedGhostCellEffect(fname=" ~ luafname ~ ")";
    }

    override ref FVCell get_mapped_cell(size_t i)
    {
	assert(0, "not implemented for this ghost_cell_effect");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	size_t j = 0, k = 0;
	FVCell ghost0, ghost1;
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    if (bc.outsigns[i] == 1) {
		ghost0 = f.right_cell;
	    } else {
		ghost0 = f.left_cell;
	    }
	    callGhostCellUDF(t, gtl, ftl, i, j, k, f, ghost0, ghost1);
	} // end foreach face
    }  // end apply_unstructured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell ghostCell0, ghostCell1;
	FVInterface IFace;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k)  {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    // ghostCell0 is closest to domain
		    // ghostCell1 is one layer out.
		    ghostCell0 = blk.get_cell(i,j+1,k);
		    ghostCell1 = blk.get_cell(i,j+2,k);
		    IFace = ghostCell0.iface[Face.south];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end i loop
	    } // end k loop
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i+1,j,k);
		    ghostCell1 = blk.get_cell(i+2,j,k);
		    IFace = ghostCell0.iface[Face.west];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end k loop
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i=blk.imin; i <= blk.imax; ++i) {
		    ghostCell0 = blk.get_cell(i,j-1,k);
		    ghostCell1 = blk.get_cell(i,j-2,k);
		    IFace = ghostCell0.iface[Face.north];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end i loop
	    } // end j loop
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i-1,j,k);
		    ghostCell1 = blk.get_cell(i-2,j,k);
		    IFace = ghostCell0.iface[Face.east];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end k loop
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i= blk.imin; i <= blk.imax; ++i) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i,j,k+1);
		    ghostCell1 = blk.get_cell(i,j,k+2);
		    IFace = ghostCell0.iface[Face.bottom];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end i loop
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i,j,k-1);
		    ghostCell1 = blk.get_cell(i,j,k-2);
		    IFace = ghostCell0.iface[Face.top];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end i loop
	    break;
	} // end switch which boundary
    }
			
private:
    void putFlowStateIntoGhostCell(lua_State* L, int tblIdx, FVCell ghostCell)
    {
	auto gmodel = blk.myConfig.gmodel;
	try {
	    ghostCell.fs.gas.p = getDouble(L, tblIdx, "p");
	    ghostCell.fs.gas.Ttr = getDouble(L, tblIdx, "T");
	    if (gmodel.n_modes > 0) {
		getArrayOfDoubles(L, tblIdx, "T_modes", ghostCell.fs.gas.T_modes);
	    }
	    lua_getfield(L, tblIdx, "massf");
	    if ( lua_istable(L, -1) ) {
		int massfIdx = lua_gettop(L);
		getSpeciesValsFromTable(L, gmodel, massfIdx, ghostCell.fs.gas.massf, "massf");
	    }
	    else {
		if ( gmodel.n_species() == 1 ) {
		    ghostCell.fs.gas.massf[0] = 1.0;
		}
		else {
		    // There's no clear choice for multi-species.
		    // Maybe best to set everything to zero to 
		    // trigger some bad behaviour rather than
		    // one value to 1.0 and have the calculation
		    // proceed but not follow the users' intent.
		    ghostCell.fs.gas.massf[] = 0.0;
		}
	    }
	    lua_pop(L, 1);
	}
	catch (Exception e) {
	    string errMsg = "There was an error trying to read p, T or massf in user-supplied table.\n";
	    errMsg ~= "The error message from the lua state follows.\n";
	    errMsg ~= e.toString();
	    throw new Exception(errMsg);
	}
	gmodel.update_thermo_from_pT(ghostCell.fs.gas);
	gmodel.update_sound_speed(ghostCell.fs.gas);
	ghostCell.fs.vel.refx = getNumberFromTable(L, tblIdx, "velx", false, 0.0);
	ghostCell.fs.vel.refy = getNumberFromTable(L, tblIdx, "vely", false, 0.0);
	ghostCell.fs.vel.refz = getNumberFromTable(L, tblIdx, "velz", false, 0.0);
	ghostCell.fs.tke = getNumberFromTable(L, tblIdx, "tke", false, 0.0);
	ghostCell.fs.omega = getNumberFromTable(L, tblIdx, "omega", false, 0.0);
    }

    void callGhostCellUDF(double t, int gtl, int ftl, size_t i, size_t j, size_t k,
			  in FVInterface IFace, FVCell ghostCell0, FVCell ghostCell1)
    {
	// 1. Set up for calling function
	auto L = blk.bc[which_boundary].myL;
	// 1a. Place function to call at TOS
	lua_getglobal(L, "ghostCells");
	// 1b. Then put arguments (as single table) at TOS
	lua_newtable(L);
	lua_pushnumber(L, t); lua_setfield(L, -2, "t");
	lua_pushnumber(L, dt_global); lua_setfield(L, -2, "dt");
	lua_pushinteger(L, step); lua_setfield(L, -2, "timeStep");
	lua_pushinteger(L, gtl); lua_setfield(L, -2, "gridTimeLevel");
	lua_pushinteger(L, ftl); lua_setfield(L, -2, "flowTimeLevel");
	lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "boundaryId");
	// Geometric information for the face.
	lua_pushnumber(L, IFace.pos.x); lua_setfield(L, -2, "x");
	lua_pushnumber(L, IFace.pos.y); lua_setfield(L, -2, "y");
	lua_pushnumber(L, IFace.pos.z); lua_setfield(L, -2, "z");
	lua_pushnumber(L, IFace.n.x); lua_setfield(L, -2, "csX"); // unit-normal, cosine x-coordinate
	lua_pushnumber(L, IFace.n.y); lua_setfield(L, -2, "csY");
	lua_pushnumber(L, IFace.n.z); lua_setfield(L, -2, "csZ");
	lua_pushnumber(L, IFace.t1.x); lua_setfield(L, -2, "csX1");
	lua_pushnumber(L, IFace.t1.y); lua_setfield(L, -2, "csY1");
	lua_pushnumber(L, IFace.t1.z); lua_setfield(L, -2, "csZ1");
	lua_pushnumber(L, IFace.t2.x); lua_setfield(L, -2, "csX2");
	lua_pushnumber(L, IFace.t2.y); lua_setfield(L, -2, "csY2");
	lua_pushnumber(L, IFace.t2.z); lua_setfield(L, -2, "csZ2");
	// Structured-grid indices for the cell just inside the boundary.
	lua_pushinteger(L, i); lua_setfield(L, -2, "i");
	lua_pushinteger(L, j); lua_setfield(L, -2, "j");
	lua_pushinteger(L, k); lua_setfield(L, -2, "k");
	// Geometric information for the ghost cells (just outside the boundary).
	lua_pushnumber(L, ghostCell0.pos[0].x); lua_setfield(L, -2, "gc0x"); // ghostcell0 x-coordinate
	lua_pushnumber(L, ghostCell0.pos[0].y); lua_setfield(L, -2, "gc0y");
	lua_pushnumber(L, ghostCell0.pos[0].z); lua_setfield(L, -2, "gc0z");
	if (ghostCell1) {
	    lua_pushnumber(L, ghostCell1.pos[0].x); lua_setfield(L, -2, "gc1x");
	    lua_pushnumber(L, ghostCell1.pos[0].y); lua_setfield(L, -2, "gc1y");
	    lua_pushnumber(L, ghostCell1.pos[0].z); lua_setfield(L, -2, "gc1z");
	}

	// 2. Call LuaFunction and expect two tables of ghost cell flow state
	int number_args = 1;
	int number_results = 1; // default to expecting 1 ghostCell table
	if (ghostCell1) {
	    number_results = 2; // expecting two table of ghostCells
	}
	if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	    luaL_error(L, "error running user-defined b.c. ghostCell function on boundaryId %d: %s\n",
		       which_boundary, lua_tostring(L, -1));
	}

	// 3. Grab Flowstate data from table and populate ghost cell
	if (ghostCell1) {
	    // Stack positions for two ghost cells:
	    //    -2 :: ghostCell0
	    //    -1 :: ghostCell1
	    putFlowStateIntoGhostCell(L, -2, ghostCell0);
	    putFlowStateIntoGhostCell(L, -1, ghostCell1);
	} else {
	    // Just the first ghost cell.
	    putFlowStateIntoGhostCell(L, -1, ghostCell0);
	}

	// 4. Clear stack
	lua_settop(L, 0);
    }
} // end class UserDefinedGhostCell


class BIE_UserDefined : BoundaryInterfaceEffect {
public:
    string luafname;
    this(int id, int boundary, string fname)
    {
	super(id, boundary, "UserDefined");
	luafname = fname;
    }
    override void post_bc_construction()
    {
	if (blk.bc[which_boundary].myL == null) {
	    blk.bc[which_boundary].myL = luaL_newstate();
	    auto L = blk.bc[which_boundary].myL;
	    luaL_openlibs(L);
	    lua_pushinteger(L, blk.id); lua_setglobal(L, "blkId");
	    registerGasModel(L, LUA_GLOBALSINDEX);
	    pushObj!(GasModel, GasModelMT)(L, blk.myConfig.gmodel);
	    lua_setglobal(L, "gmodel");
	    lua_pushinteger(L, blk.myConfig.gmodel.n_species);
	    lua_setglobal(L, "n_species");
	    lua_pushinteger(L, blk.myConfig.gmodel.n_modes);
	    lua_setglobal(L, "n_modes");
	    lua_pushinteger(L, blk.nicell); lua_setglobal(L, "nicell");
	    lua_pushinteger(L, blk.njcell); lua_setglobal(L, "njcell");
	    lua_pushinteger(L, blk.nkcell); lua_setglobal(L, "nkcell");
	    lua_pushinteger(L, Face.north); lua_setglobal(L, "north");
	    lua_pushinteger(L, Face.east); lua_setglobal(L, "east");
	    lua_pushinteger(L, Face.south); lua_setglobal(L, "south");
	    lua_pushinteger(L, Face.west); lua_setglobal(L, "west");
	    lua_pushinteger(L, Face.top); lua_setglobal(L, "top");
	    lua_pushinteger(L, Face.bottom); lua_setglobal(L, "bottom");
	    lua_pushcfunction(L, &luafn_sampleFlow); lua_setglobal(L, "sampleFlow");
	}
	luaL_dofile(blk.bc[which_boundary].myL, luafname.toStringz);
    }

    override string toString() const
    {
	return "UserDefined(fname=" ~ luafname ~ ")";
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
	size_t j = 0, k = 0;
	BoundaryCondition bc = blk.bc[which_boundary];
	foreach (i, f; bc.faces) {
	    callInterfaceUDF(t, gtl, ftl, i, j, k, f);
	} // end foreach face
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k)  {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end i loop
	    } // end k loop
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end k loop
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i=blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end i loop
	    } // end j loop
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end k loop
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i= blk.imin; i <= blk.imax; ++i) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end i loop
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end i loop
	    break;
	} // end switch which boundary
    }
private:
    void putFlowStateIntoInterface(lua_State* L, int tblIdx, FVInterface iface)
    {
	// Now the user might only set some of the flowstate
	// since they might be relying on another boundary
	// effect to do some work.
	// So we need to test every possibility and only set
	// the non-nil values.
	auto gmodel = blk.myConfig.gmodel;
	FlowState fs = iface.fs;
	
	lua_getfield(L, tblIdx, "p");
	if ( !lua_isnil(L, -1) ) {
	    fs.gas.p = getDouble(L, tblIdx, "p");
	}
	lua_pop(L, 1);
	
	lua_getfield(L, tblIdx, "T");
	if ( !lua_isnil(L, -1) ) {
	    fs.gas.Ttr = getDouble(L, tblIdx, "T");
	}
	lua_pop(L, 1);

	if (gmodel.n_modes > 0) {
	    lua_getfield(L, tblIdx, "T_modes");
	    if ( !lua_isnil(L, -1) ) {
		// Interior-modes temperatures should be provided as an array.
		getArrayOfDoubles(L, tblIdx, "T_modes", fs.gas.T_modes);
	    }
	    lua_pop(L, 1);
	}

	lua_getfield(L, tblIdx, "massf");
	if ( !lua_isnil(L, -1) ) {
	    int massfIdx = lua_gettop(L);
	    getSpeciesValsFromTable(L, gmodel, massfIdx, fs.gas.massf, "massf");
	}
	lua_pop(L, 1);

	lua_getfield(L, tblIdx, "velx");
	if ( !lua_isnil(L, -1) ) {
	    fs.vel.refx = getDouble(L, tblIdx, "velx");
	}
	lua_pop(L, 1);

	lua_getfield(L, tblIdx, "vely");
	if ( !lua_isnil(L, -1) ) {
	    fs.vel.refy = getDouble(L, tblIdx, "vely");
	}
	lua_pop(L, 1);

	lua_getfield(L, tblIdx, "velz");
	if ( !lua_isnil(L, -1) ) {
	    fs.vel.refz = getDouble(L, tblIdx, "velz");
	}
	lua_pop(L, 1);

	lua_getfield(L, tblIdx, "tke");
	if ( !lua_isnil(L, -1) ) {
	    fs.tke = getDouble(L, tblIdx, "tke");
	}
	lua_pop(L, 1);

	lua_getfield(L, tblIdx, "omega");
	if ( !lua_isnil(L, -1) ) {
	    fs.omega = getDouble(L, tblIdx, "omega");
	}
	lua_pop(L, 1);
    }
	    
    void callInterfaceUDF(double t, int gtl, int ftl, size_t i, size_t j, size_t k,
			  FVInterface IFace)
    {
	// 1. Set up for calling function
	auto L = blk.bc[which_boundary].myL;
	// 1a. Place function to call at TOS
	lua_getglobal(L, "interface");
	// 1b. Then put arguments (as single table) at TOS
	lua_newtable(L);
	lua_pushnumber(L, t); lua_setfield(L, -2, "t");
	lua_pushnumber(L, dt_global); lua_setfield(L, -2, "dt");
	lua_pushinteger(L, step); lua_setfield(L, -2, "timeStep");
	lua_pushinteger(L, gtl); lua_setfield(L, -2, "gridTimeLevel");
	lua_pushinteger(L, ftl); lua_setfield(L, -2, "flowTimeLevel");
	lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "boundaryId");
	lua_pushnumber(L, IFace.pos.x); lua_setfield(L, -2, "x");
	lua_pushnumber(L, IFace.pos.y); lua_setfield(L, -2, "y");
	lua_pushnumber(L, IFace.pos.z); lua_setfield(L, -2, "z");
	lua_pushnumber(L, IFace.n.x); lua_setfield(L, -2, "csX");
	lua_pushnumber(L, IFace.n.y); lua_setfield(L, -2, "csY");
	lua_pushnumber(L, IFace.n.z); lua_setfield(L, -2, "csZ");
	lua_pushnumber(L, IFace.t1.x); lua_setfield(L, -2, "csX1");
	lua_pushnumber(L, IFace.t1.y); lua_setfield(L, -2, "csY1");
	lua_pushnumber(L, IFace.t1.z); lua_setfield(L, -2, "csZ1");
	lua_pushnumber(L, IFace.t2.x); lua_setfield(L, -2, "csX2");
	lua_pushnumber(L, IFace.t2.y); lua_setfield(L, -2, "csY2");
	lua_pushnumber(L, IFace.t2.z); lua_setfield(L, -2, "csZ2");
	lua_pushinteger(L, i); lua_setfield(L, -2, "i");
	lua_pushinteger(L, j); lua_setfield(L, -2, "j");
	lua_pushinteger(L, k); lua_setfield(L, -2, "k");

	// 2. Call LuaFunction and expect a table of values for flow state
	int number_args = 1;
	int number_results = 1;
	if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	    luaL_error(L, "error running user-defined b.c. interface function on boundaryId %d: %s\n",
		       which_boundary, lua_tostring(L, -1));
	}

	// 3. Grab Flowstate data from table and populate interface
	int tblIdx = lua_gettop(L);
	putFlowStateIntoInterface(L, tblIdx, IFace);

	// 4. Clear stack
	lua_settop(L, 0);
    }
} // end class BIEUserDefined

