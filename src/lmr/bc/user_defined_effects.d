// user_defined_effects.d
//
// Authors: RG & PJ
// Date: 2015-03-14
//       2017-05-09 clean up Lua interpreter initialization.
//
// Notes:
// There is one Lua interpreter for the BoundaryCondition object
// and each effect may initialize the interpreter if it is sees
// that the interpreter has not yet been initialized.
// The consequence of this initialize when needed approach is that
// only one of the effects to do the initialization,
// using the file name that it possesses.  Thus, you need to have
// all of your user-defined effects within the one file.

module bc.user_defined_effects;

import std.conv;
import std.string;
import std.stdio;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas.gas_model;
import gas.luagas_model;
import gas.gas_state;

import geom;
import simcore;
import flowstate;
import lmr.fluidfvcell;
import fvinterface;
import sfluidblock: SFluidBlock;
import sfluidblock : cell_index_to_logical_coordinates;
import globalconfig;
import globaldata;
import luaflowstate;
import bc;

class UserDefinedGhostCell : GhostCellEffect {
public:
    string luaFileName;

    this(int id, int boundary, string fname)
    {
        super(id, boundary, "UserDefined");
        luaFileName = fname;
    }
    override void post_bc_construction()
    {
        if ((luaFileName.length > 0) && (blk.bc[which_boundary].myL == null)) {
            blk.bc[which_boundary].init_lua_State(luaFileName);
        }
    }
    override string toString() const
    {
        return "UserDefinedGhostCellEffect(luaFileName=" ~ luaFileName ~ ")";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        size_t j = 0, k = 0;
        FluidFVCell[1] ghostCells;
        BoundaryCondition bc = blk.bc[which_boundary];
	if (bc.outsigns[f.i_bndry] == 1) {
	    ghostCells[0] = f.right_cell;
	} else {
	    ghostCells[0] = f.left_cell;
	}
	callGhostCellUDF(t, gtl, ftl, f.i_bndry, j, k, f, ghostCells);
	lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }  // end apply_for_interface_unstructured_grid()

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        size_t j = 0, k = 0;
        FluidFVCell[1] ghostCells;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            if (bc.outsigns[i] == 1) {
                ghostCells[0] = f.right_cell;
            } else {
                ghostCells[0] = f.left_cell;
            }
            callGhostCellUDF(t, gtl, ftl, i, j, k, f, ghostCells);
        } // end foreach face
        lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }  // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        size_t[3] ijk;
        if (bc.outsigns[f.i_bndry] == 1) {
            // Indices passed in are for the cell just inside the boundary.
            ijk = cell_index_to_logical_coordinates(f.left_cells[0].id, blk.nic, blk.njc);
            callGhostCellUDF(t, gtl, ftl, ijk[0], ijk[1], ijk[2], f, f.right_cells);
	} else {
            // Indices passed in are for the cell just inside the boundary.
            ijk = cell_index_to_logical_coordinates(f.right_cells[0].id, blk.nic, blk.njc);
            callGhostCellUDF(t, gtl, ftl, ijk[0], ijk[1], ijk[2], f, f.left_cells);
	}
    }  // end apply_for_interface_structured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        //
        final switch (which_boundary) {
        case Face.north:
            size_t j = blk.njc;
            foreach (k; 0 .. blk.nkc)  {
                foreach (i; 0 .. blk.nic) {
                    auto f = blk.get_ifj(i,j,k);
                    // Indices passed in are for the cell just inside the boundary.
                    callGhostCellUDF(t, gtl, ftl, i, j-1, k, f, f.right_cells);
                }
            }
            break;
        case Face.east:
            size_t i = blk.nic;
            foreach (k; 0 .. blk.nkc) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifi(i,j,k);
                    callGhostCellUDF(t, gtl, ftl, i-1, j, k, f, f.right_cells);
                }
            }
            break;
        case Face.south:
            size_t j = 0;
            foreach (k; 0 .. blk.nkc)  {
                foreach (i; 0 .. blk.nic) {
                    auto f = blk.get_ifj(i,j,k);
                    callGhostCellUDF(t, gtl, ftl, i, j, k, f, f.left_cells);
                }
            }
            break;
        case Face.west:
            size_t i = 0;
            foreach (k; 0 .. blk.nkc) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifi(i,j,k);
                    callGhostCellUDF(t, gtl, ftl, i, j, k, f, f.left_cells);
                }
            }
            break;
        case Face.top:
            size_t k = blk.nkc;
            foreach (i; 0 .. blk.nic) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifk(i,j,k);
                    callGhostCellUDF(t, gtl, ftl, i, j, k-1, f, f.right_cells);
                }
            }
            break;
        case Face.bottom:
            size_t k = 0;
            foreach (i; 0 .. blk.nic) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifk(i,j,k);
                    callGhostCellUDF(t, gtl, ftl, i, j, k, f, f.left_cells);
                }
            }
            break;
        } // end switch which boundary
        lua_gc(blk.bc[which_boundary].myL, LUA_GCCOLLECT, 0);
    }

private:
    // not @nogc
    void putFlowStateIntoGhostCell(lua_State* L, int tblIdx, FluidFVCell ghostCell)
    {
        int old_top = lua_gettop(L);
        auto gmodel = blk.myConfig.gmodel;
        GasState* gs = &(ghostCell.fs.gas);
        try {
            gs.p = getDouble(L, tblIdx, "p");
            gs.T = getDouble(L, tblIdx, "T");
            version(multi_T_gas) {
                if (gmodel.n_modes > 0) {
                    getArrayOfDoubles(L, tblIdx, "T_modes", gs.T_modes);
                }
            }
            version(multi_species_gas) {
                lua_getfield(L, tblIdx, "massf");
                if ( lua_istable(L, -1) ) {
                    int massfIdx = lua_gettop(L);
                    getSpeciesValsFromTable(L, gmodel, massfIdx, gs.massf, "massf");
                } else {
                    if ( gmodel.n_species() == 1 ) {
                        gs.massf[0] = 1.0;
                    } else {
                        // There's no clear choice for multi-species.
                        // Maybe best to set everything to zero to
                        // trigger some bad behaviour rather than
                        // one value to 1.0 and have the calculation
                        // proceed but not follow the users' intent.
                        foreach (ref mf; gs.massf) { mf = 0.0; }
                    }
                }
                lua_pop(L, 1);
            }
        }
        catch (Exception e) {
            string errMsg = "There was an error trying to read p, T or massf in user-supplied table.\n";
            errMsg ~= "The error message from the lua state follows.\n";
            errMsg ~= e.toString();
            throw new Exception(errMsg);
        }
        gmodel.update_thermo_from_pT(*gs);
        gmodel.update_sound_speed(*gs);
        // Once we have mass fractions and density, we need to construct values for the species densities
        // because these individual densities are now used by the core gas dynamics (since 20 June 2023).
        foreach (i; 0 .. gmodel.n_species()) { gs.rho_s[i] = gs.massf[i] * gs.rho; }
        // For UserDefinedGhostCellBC, the following call to update_trans_coeffs() is done
        // a little later via an action in the preSpatialDerivActionAtBndryFaces list.
        // gmodel.update_trans_coeffs(*gs);

        ghostCell.fs.vel.x = getNumberFromTable(L, tblIdx, "velx", false, 0.0);
        ghostCell.fs.vel.y = getNumberFromTable(L, tblIdx, "vely", false, 0.0);
        ghostCell.fs.vel.z = getNumberFromTable(L, tblIdx, "velz", false, 0.0);

        version(turbulence) {
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                string tname = blk.myConfig.turb_model.primitive_variable_name(it);
                ghostCell.fs.turb[it] = getNumberFromTable(L, -1, tname, false, 0.0);
            }
        }
        warnOnStackChange(L, old_top);
    } // end putFlowStateIntoGhostCell()

    // not @nogc because of Lua functions
    void callGhostCellUDF(double t, int gtl, int ftl, size_t i, size_t j, size_t k,
                          in FVInterface IFace, FluidFVCell[] ghostCells)
    {
        // 1. Set up for calling function
        auto L = blk.bc[which_boundary].myL;
        lua_settop(L, 0);
        bool nghost3 = (blk.n_ghost_cell_layers == 3);
        // 1a. Place function to call at TOS
        lua_getglobal(L, "ghostCells");
        // 1b. Then put arguments (as single table) at TOS
        lua_newtable(L);
        int old_top = lua_gettop(L);
        lua_pushnumber(L, t); lua_setfield(L, -2, "t");
        lua_pushnumber(L, SimState.dt_global); lua_setfield(L, -2, "dt");
        lua_pushinteger(L, SimState.step); lua_setfield(L, -2, "timeStep");
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
        foreach (ig; 0 .. blk.n_ghost_cell_layers) {
            lua_pushnumber(L, ghostCells[ig].pos[0].x); lua_setfield(L, -2, format("gc%dx", ig).toStringz);
            lua_pushnumber(L, ghostCells[ig].pos[0].y); lua_setfield(L, -2, format("gc%dy", ig).toStringz);
            lua_pushnumber(L, ghostCells[ig].pos[0].z); lua_setfield(L, -2, format("gc%dz", ig).toStringz);
        }
        warnOnStackChange(L, old_top);

        // 2. Call LuaFunction and expect two tables of ghost cell flow state
        int number_args = 1;
        int number_results = to!int(blk.n_ghost_cell_layers);
        if (lua_pcall(L, number_args, number_results, 0) != 0) {
            luaL_error(L, "error running user-defined b.c. ghostCell function on boundaryId %d: %s\n",
                       which_boundary, lua_tostring(L, -1));
        }

        // 3. Grab Flowstate data from tables from stack and populate ghost-cell FlowStates.
        // We assume that the user-defined function returns the correct number of tables
        // representing the flow states in the ghost cells.
        // Stack positions for ghost cells:
        //    -3 :: ghostCell0  -2 :: ghostCell0  -1 :: ghostCell0
        //    -2 :: ghostCell1  -1 :: ghostCell1
        //    -1 :: ghostCell2
        foreach (ig; 0 .. blk.n_ghost_cell_layers) {
            int stack_location = -(to!int(blk.n_ghost_cell_layers-ig));
            if (lua_isnil(L, stack_location)) {
                writefln("Oops, expected a table for the ghost-cell flow-state ig=%d, but got nil.", ig);
                break;
            }
            if (lua_istable(L, stack_location) && !tableEmpty(L, stack_location)) {
                putFlowStateIntoGhostCell(L, stack_location, ghostCells[ig]);
            }
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
            blk.bc[which_boundary].init_lua_State(luafname);
        }
    }

    override string toString() const
    {
        return "UserDefined(fname=" ~ luafname ~ ")";
    }

    // not @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        size_t j = 0, k = 0;
        BoundaryCondition bc = blk.bc[which_boundary];
	callInterfaceUDF(t, gtl, ftl, f.i_bndry, j, k, f);
	lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }

    // not @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        size_t j = 0, k = 0;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            callInterfaceUDF(t, gtl, ftl, i, j, k, f);
        }
        lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        size_t j = 0, k = 0;
        BoundaryCondition bc = blk.bc[which_boundary];
	callInterfaceUDF(t, gtl, ftl, f.i_bndry, j, k, f);
	lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        //
        final switch (which_boundary) {
        case Face.north:
            size_t j = blk.njc;
            foreach (k; 0 .. blk.nkc)  {
                foreach (i; 0 .. blk.nic) {
                    auto f = blk.get_ifj(i,j,k);
                    // Indices passed in are for the cell just inside the boundary.
                    callInterfaceUDF(t, gtl, ftl, i, j-1, k, f);
                }
            }
            break;
        case Face.east:
            size_t i = blk.nic;
            foreach (k; 0 .. blk.nkc) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifi(i,j,k);
                    callInterfaceUDF(t, gtl, ftl, i-1, j, k, f);
                }
            }
            break;
        case Face.south:
            size_t j = 0;
            foreach (k; 0 .. blk.nkc)  {
                foreach (i; 0 .. blk.nic) {
                    auto f = blk.get_ifj(i,j,k);
                    callInterfaceUDF(t, gtl, ftl, i, j, k, f);
                }
            }
            break;
        case Face.west:
            size_t i = 0;
            foreach (k; 0 .. blk.nkc) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifi(i,j,k);
                    callInterfaceUDF(t, gtl, ftl, i, j, k, f);
                }
            }
            break;
        case Face.top:
            size_t k = blk.nkc;
            foreach (i; 0 .. blk.nic) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifk(i,j,k);
                    callInterfaceUDF(t, gtl, ftl, i, j, k-1, f);
                }
            }
            break;
        case Face.bottom:
            size_t k = blk.nkc;
            foreach (i; 0 .. blk.nic) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifk(i,j,k);
                    callInterfaceUDF(t, gtl, ftl, i, j, k, f);
                }
            }
            break;
        } // end switch which boundary
        lua_gc(blk.bc[which_boundary].myL, LUA_GCCOLLECT, 0);
    }
private:
    // not @nogc because of Lua function calls
    void putFlowStateIntoInterface(lua_State* L, int tblIdx, FVInterface iface)
    {
        // Now the user might only set some of the flowstate
        // since they might be relying on another boundary
        // effect to do some work.
        // So we need to test every possibility and only set
        // the non-nil values.
        auto gmodel = blk.myConfig.gmodel;
        FlowState* fs = iface.fs;

        lua_getfield(L, tblIdx, "p");
        if ( !lua_isnil(L, -1) ) {
            fs.gas.p = getDouble(L, tblIdx, "p");
        }
        lua_pop(L, 1);

        lua_getfield(L, tblIdx, "T");
        if ( !lua_isnil(L, -1) ) {
            fs.gas.T = getDouble(L, tblIdx, "T");
        }
        lua_pop(L, 1);

        version(multi_T_gas) {
            if (gmodel.n_modes > 0) {
                lua_getfield(L, tblIdx, "T_modes");
                if ( !lua_isnil(L, -1) ) {
                    // Interior-modes temperatures should be provided as an array.
                    getArrayOfDoubles(L, tblIdx, "T_modes", fs.gas.T_modes);
                }
                lua_pop(L, 1);
            }
        }

        version(multi_species_gas) {
            // If we only have one species, set it to one for safety
            if ( gmodel.n_species() == 1 ) {
                fs.gas.massf[0] = 1.0;
            } else {
                // If the user provided a table, pull it in, otherwise do nothing
                lua_getfield(L, tblIdx, "massf");
                if ( lua_istable(L, -1) ) {
                    int massfIdx = lua_gettop(L);
                    getSpeciesValsFromTable(L, gmodel, massfIdx, fs.gas.massf, "massf");
                }
                lua_pop(L, 1);
            }
        }

        lua_getfield(L, tblIdx, "velx");
        if ( !lua_isnil(L, -1) ) {
            fs.vel.x = getDouble(L, tblIdx, "velx");
        }
        lua_pop(L, 1);

        lua_getfield(L, tblIdx, "vely");
        if ( !lua_isnil(L, -1) ) {
            fs.vel.y = getDouble(L, tblIdx, "vely");
        }
        lua_pop(L, 1);

        lua_getfield(L, tblIdx, "velz");
        if ( !lua_isnil(L, -1) ) {
            fs.vel.z = getDouble(L, tblIdx, "velz");
        }
        lua_pop(L, 1);

        version(turbulence) {
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                string tname = blk.myConfig.turb_model.primitive_variable_name(it);
                fs.turb[it] = getNumberFromTable(L, -1, tname, false, 0.0);
            }
        }

        lua_getfield(L, tblIdx, "mu_t");
        if ( !lua_isnil(L, -1) ) {
            fs.mu_t = getDouble(L, tblIdx, "mu_t");
        }
        lua_pop(L, 1);

        lua_getfield(L, tblIdx, "k_t");
        if ( !lua_isnil(L, -1) ) {
            fs.k_t = getDouble(L, tblIdx, "k_t");
        }

        gmodel.update_thermo_from_pT(fs.gas);
        // Once we have mass fractions and density, we need to construct values for the species densities
        // because these individual densities are now used by the core gas dynamics (since 20 June 2023).
        foreach (i; 0 .. gmodel.n_species()) { fs.gas.rho_s[i] = fs.gas.massf[i] * fs.gas.rho; }
        gmodel.update_trans_coeffs(fs.gas);
        gmodel.update_sound_speed(fs.gas);

        lua_pop(L, 1);
    } // end putFlowStateIntoInterface()

    // not @nogc
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
        lua_pushnumber(L, SimState.dt_global); lua_setfield(L, -2, "dt");
        lua_pushinteger(L, SimState.step); lua_setfield(L, -2, "timeStep");
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
        // Note that the following indices are for the cell just inside the boundary.
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
        if (!tableEmpty(L, tblIdx)) { putFlowStateIntoInterface(L, tblIdx, IFace); }

        // 4. Clear stack
        lua_settop(L, 0);
    }
} // end class BIEUserDefined

class BFE_UserDefined : BoundaryFluxEffect {
public:
    string luafname;
    string luaFnName;
    this(int id, int boundary, string fname, string funcName)
    {
        super(id, boundary, "UserDefinedFluxEffect");
        luafname = fname;
        luaFnName = funcName;
    }
    override void post_bc_construction()
    {
        if (blk.bc[which_boundary].myL == null) {
            blk.bc[which_boundary].init_lua_State(luafname);
        }
        if (blk.myConfig.n_species>1){
            massflux.length = blk.myConfig.n_species;
        }
        if (blk.myConfig.n_modes>0){
            modesflux.length = blk.myConfig.n_modes;
        }
   }
    override string toString() const
    {
        return "UserDefinedFluxEffect(fname=" ~ luafname ~ ", luaFnName=" ~ luaFnName ~ ")";
    }

    // not @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        size_t j = 0, k = 0;
        FluidFVCell ghost0, ghost1;
        BoundaryCondition bc = blk.bc[which_boundary];
        int outsign = bc.outsigns[f.i_bndry];
        callFluxUDF(t, gtl, ftl, f.i_bndry, j, k, f, outsign);
        lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }  // end apply_unstructured_grid()

    // not @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        size_t j = 0, k = 0;
        FluidFVCell ghost0, ghost1;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            int outsign = bc.outsigns[f.i_bndry];
            callFluxUDF(t, gtl, ftl, i, j, k, f, outsign);
        }
        lua_gc(bc.myL, LUA_GCCOLLECT, 0);
    }  // end apply_unstructured_grid()

    //@nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        size_t[3] ijk;
        if (bc.outsigns[f.i_bndry] == 1) {
            // Indices passed in are for the cell just inside the boundary.
            ijk = cell_index_to_logical_coordinates(f.left_cells[0].id, blk.nic, blk.njc);
	} else {
            // Indices passed in are for the cell just inside the boundary.
            ijk = cell_index_to_logical_coordinates(f.right_cells[0].id, blk.nic, blk.njc);
	}
        int outsign = bc.outsigns[f.i_bndry];
        callFluxUDF(t, gtl, ftl, ijk[0], ijk[1], ijk[2], f, outsign);
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        //
        final switch (which_boundary) {
        case Face.north:
            size_t j = blk.njc;
            foreach (k; 0 .. blk.nkc)  {
                foreach (i; 0 .. blk.nic) {
                    auto f = blk.get_ifj(i,j,k);
                    int outsign = bc.outsigns[f.i_bndry];
                    // Indices passed in are for the cell just inside the boundary.
                    callFluxUDF(t, gtl, ftl, i, j-1, k, f, outsign);
                }
            }
            break;
        case Face.east:
            size_t i = blk.nic;
            foreach (k; 0 .. blk.nkc) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifi(i,j,k);
                    int outsign = bc.outsigns[f.i_bndry];
                    callFluxUDF(t, gtl, ftl, i-1, j, k, f, outsign);
                }
            }
            break;
        case Face.south:
            size_t j = 0;
            foreach (k; 0 .. blk.nkc)  {
                foreach (i; 0 .. blk.nic) {
                    auto f = blk.get_ifj(i,j,k);
                    int outsign = bc.outsigns[f.i_bndry];
                    callFluxUDF(t, gtl, ftl, i, j, k, f, outsign);
                }
            }
            break;
        case Face.west:
            size_t i = 0;
            foreach (k; 0 .. blk.nkc) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifi(i,j,k);
                    int outsign = bc.outsigns[f.i_bndry];
                    callFluxUDF(t, gtl, ftl, i, j, k, f, outsign);
                }
            }
            break;
        case Face.top:
            size_t k = blk.nkc;
            foreach (i; 0 .. blk.nic) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifk(i,j,k);
                    int outsign = bc.outsigns[f.i_bndry];
                    callFluxUDF(t, gtl, ftl, i, j, k-1, f, outsign);
                }
            }
            break;
        case Face.bottom:
            size_t k = blk.nkc;
            foreach (i; 0 .. blk.nic) {
                foreach (j; 0 .. blk.njc) {
                    auto f = blk.get_ifk(i,j,k);
                    int outsign = bc.outsigns[f.i_bndry];
                    callFluxUDF(t, gtl, ftl, i, j, k, f, outsign);
                }
            }
            break;
        } // end switch which boundary
        lua_gc(blk.bc[which_boundary].myL, LUA_GCCOLLECT, 0);
    } // end apply_structured_grid()

private:
    double[] massflux, modesflux;
    // not @nogc
    void putFluxIntoInterface(lua_State* L, int tblIdx, FVInterface iface)
    {
        // Now the user might only set some of the fluxes
        // since they might be relying on another boundary
        // effect to do some work.
        // So we need to test every possibility and only set
        // the non-nil values.
        auto gmodel = blk.myConfig.gmodel;
        auto cqi = blk.myConfig.cqi;

        lua_getfield(L, tblIdx, "mass");
        if ( !lua_isnil(L, -1) ) {
            if (cqi.mass==0) iface.F[cqi.mass] += getDouble(L, tblIdx, "mass");
        }
        lua_pop(L, 1);

        lua_getfield(L, tblIdx, "momentum_x");
        if ( !lua_isnil(L, -1) ) {
            iface.F[cqi.xMom] += getDouble(L, tblIdx, "momentum_x");
        }
        lua_pop(L, 1);

        lua_getfield(L, tblIdx, "momentum_y");
        if ( !lua_isnil(L, -1) ) {
            iface.F[cqi.yMom] += getDouble(L, tblIdx, "momentum_y");
        }
        lua_pop(L, 1);

        if (cqi.threeD) {
            lua_getfield(L, tblIdx, "momentum_z");
            if ( !lua_isnil(L, -1) ) {
                iface.F[cqi.zMom] += getDouble(L, tblIdx, "momentum_z");
            }
            lua_pop(L, 1);
        }

        lua_getfield(L, tblIdx, "total_energy");
        if ( !lua_isnil(L, -1) ) {
            iface.F[cqi.totEnergy] += getDouble(L, tblIdx, "total_energy");
        }
        lua_pop(L, 1);

        version(multi_T_gas) {
            if (gmodel.n_modes > 0) {
                lua_getfield(L, tblIdx, "T_modes");
                if ( !lua_isnil(L, -1) ) {
                    // Interior-modes fluxes should be provided as an array.
                    getArrayOfDoubles(L, tblIdx, "T_modes", modesflux);
                    foreach (i; 0 .. cqi.n_modes) {
                        iface.F[cqi.modes+i] += modesflux[i];
                    }
                }
                lua_pop(L, 1);
            }
        }

        version(multi_species_gas) {
            if (cqi.n_species>1){
                lua_getfield(L, tblIdx, "species");
                if (!lua_isnil(L, -1)) {
                    int massfIdx = lua_gettop(L);
                    getSpeciesValsFromTable(L, gmodel, massfIdx, massflux, "species");

                    version(steady_state) {
                        foreach (i; 0 .. cqi.n_species) {
                            iface.F[cqi.species+i] += massflux[i];
                        }
                    } else { // The transient code uses mass fractions instead
                        double totalmassflux = 0.0;
                        foreach (i; 0 .. cqi.n_species) {
                            totalmassflux += massflux[i];
                        }
                        if (cqi.mass==0) iface.F[cqi.mass] += totalmassflux;

                        foreach (i; 0 .. cqi.n_species) {
                            iface.F[cqi.species+i] += massflux[i]/totalmassflux;
                        }
                    }
                }
                lua_pop(L, 1);

            }
        }

        version(turbulence) {
            // TODO: Think on tke and omega fluxes.
            // foreach(it; 0 .. blk.myConfig.turb_model.nturb){
            //     string tname = blk.myConfig.turb_model.primitive_variable_name(it);
            //     fs.turb[it] = getNumberFromTable(L, -1, tname, false, 0.0);
            // }
        }
    } // end putFluxIntoInterface()

    // not @nogc
    void callFluxUDF(double t, int gtl, int ftl, size_t i, size_t j, size_t k,
                     FVInterface IFace, int outsign)
    {
        // 1. Set up for calling function
        auto L = blk.bc[which_boundary].myL;
        // 1a. Place function to call at TOS
        lua_getglobal(L, luaFnName.toStringz);
        // 1b. Then put arguments (as single table) at TOS
        lua_newtable(L);
        lua_pushnumber(L, t); lua_setfield(L, -2, "t");
        lua_pushnumber(L, SimState.dt_global); lua_setfield(L, -2, "dt");
        lua_pushinteger(L, SimState.step); lua_setfield(L, -2, "timeStep");
        lua_pushinteger(L, gtl); lua_setfield(L, -2, "gridTimeLevel");
        lua_pushinteger(L, ftl); lua_setfield(L, -2, "flowTimeLevel");
        lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "boundaryId");
        lua_pushinteger(L, outsign); lua_setfield(L, -2, "outsign");
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
        // Note that the following indices are for the cell just inside the boundary.
        lua_pushinteger(L, i); lua_setfield(L, -2, "i");
        lua_pushinteger(L, j); lua_setfield(L, -2, "j");
        lua_pushinteger(L, k); lua_setfield(L, -2, "k");

        // 2. Call LuaFunction and expect a table of flux values
        int number_args = 1;
        int number_results = 1;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            luaL_error(L, "error running user-defined b.c. interface function on boundaryId %d: %s\n",
                       which_boundary, lua_tostring(L, -1));
        }

        // 3. Grab flux data from table and populate interface fluxes
        int tblIdx = lua_gettop(L);
        if (!tableEmpty(L, tblIdx)) { putFluxIntoInterface(L, tblIdx, IFace); }

        // 4. Clear stack
        lua_settop(L, 0);
    }

} // end class BFE_UserDefined
