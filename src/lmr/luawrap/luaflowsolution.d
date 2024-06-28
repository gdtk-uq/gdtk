/**
 * A lua wrapper for the D FlowSolution class.
 *
 * The uses cases for this in the Lua layer are:
 *   1. For reading in a flow solution during preparation stage.
 *      Typically, one might want to read a previously-generated
 *      flow field solution and use this as part of a initialisation
 *      condition. We sometimes call this warm starting.
 *
 *   2. For reading in a flow solution and doing actions such as
 *      manipulation or query during post-processing.
 *
 * Authors: PAJ, RJG, NNG and KAD
 * Date: 2023-12-22
 *
 * History: This code is a modified version of the Eilmer4 code
 *          that began life on 2015-10-20. The main difference is
 *          to remove dependencies on a tindx and times file.
 */

module luaflowsolution;

import std.algorithm;
import std.array;
import std.format;
import std.stdio;
import std.string;
import std.conv;
import std.traits;
import core.memory : GC;
import util.lua;
import util.lua_service;
import globalconfig;
import lmrconfig;
import cmdhelper : determineAvailableSnapshots;
import geom;
import geom.luawrap;
import flowsolution;

/// name for FlowSolution object in Lua scripts.
immutable string FlowSolutionMT = "FlowSolution";

static const(FlowSolution)[] flowSolutionStore;

FlowSolution checkFlowSolution(lua_State* L, int index)
{
    return checkObj!(FlowSolution, FlowSolutionMT)(L, index);
}

/**
 * This function implements our constructor for the Lua interface.
 *
 * Construction of a FlowSolution object from in Lua will accept:
 * fs = FlowSolution:new{dir=<string>, snapshot=<int|string>, nBlocks=<int>, make_kdtree=<bool>}
 */
extern(C) int newFlowSolution(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument "this".
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " A table is expected as first (and only) argument.";
        luaL_error(L, errMsg.toStringz);
    }
    string[] allowedNames = ["dir", "snapshot", "nBlocks", "sim_time", "tag", "make_kdtree"];
    if ( !checkAllowedNames(L, 1, allowedNames) ) {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " The table contains unexpected names.";
        errMsg ~= format(" Allowed names: %s ", allowedNames);
        luaL_error(L, errMsg.toStringz);
    }


    string dir;
    lua_getfield(L, 1, "dir");
    if ( lua_isnil(L, -1) ) {
        dir = ".";
    } else if ( lua_isstring(L, -1) ) {
        dir = to!string(luaL_checkstring(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " A field for 'dir' was found, but the content was not valid.";
        errMsg ~= " The 'dir' field, if given, should be a string.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    string snapshotDir = (dir == ".") ? lmrCfg.snapshotDir : dir ~ "/" ~ lmrCfg.snapshotDir;
    string[] availSnapshots = determineAvailableSnapshots(snapshotDir);

    int snapshot;
    lua_getfield(L, 1, "snapshot");
    if ( lua_isnil(L, -1) ) {
        snapshot = 0;
    } else if ( lua_isnumber(L, -1) ) {
        snapshot = to!int(luaL_checknumber(L, -1));
	if (snapshot >= availSnapshots.length) {
	    string errMsg = "Error in call to FlowSolution:new.\n";
	    errMsg ~= format("An integer value was passed to the field 'snapshot': snapshot=%d\n", snapshot);
	    errMsg ~= format("That snapshot is not available. Total number of snapshots found is: %d\n", availSnapshots.length);
	    errMsg ~= "(And snapshots are counted from 0.)\n";
	    luaL_error(L, errMsg.toStringz);
	}
    } else if ( lua_isstring(L, -1) ) {
        string snapStr = to!string(luaL_checkstring(L, -1));
        if ( snapStr == "last" || snapStr == "final" ) {
            snapshot = to!int(availSnapshots[$-1]);
        } else {
            string errMsg = "Error in call to FlowSolution:new.\n";
            errMsg ~= " A string value was passed to the field 'snapshot', but the content was not valid.\n";
            errMsg ~= " The only valid strings are 'last' and 'final'.\n";
            luaL_error(L, errMsg.toStringz);
        }
    } else {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " A value in field 'snapshot' was found, but the content was not valid.";
        errMsg ~= " The snapshot field, if given, should be an integer or the string 'last'.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int nBlocks;
    lua_getfield(L, 1, "nBlocks");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " A field for 'nBlocks' was not found.";
        luaL_error(L, errMsg.toStringz);
    } else if ( lua_isnumber(L, -1) ) {
        nBlocks = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " A field for 'nBlocks' was found, but the content was not valid.";
        errMsg ~= " The nBlocks field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    bool make_kdtree;
    lua_getfield(L, 1, "make_kdtree");
    if ( lua_isnil(L, -1) ) {
        make_kdtree = true;
    } else if ( lua_isboolean(L, -1) ) {
        make_kdtree = to!bool(lua_toboolean(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:new.";
        errMsg ~= " A field for 'make_kdtree' was found, but the content was not valid.";
        errMsg ~= " The value should be given as a boolean: true or false.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    double simTime;
    lua_getfield(L, 1, "sim_time");
    if (lua_isnil(L, -1)) {
        simTime = -1.0;
    } else if (lua_isnumber(L, -1)) {
        simTime = luaL_checknumber(L, -1);
    } else {
        string errMsg = "Error in call to FlowSolution:new.\n";
        errMsg ~= "  A field for 'sim_time' was found, but the content was not valid.\n";
        errMsg ~= "  The value should be given as a number.\n";
        luaL_error(L, errMsg.toStringz);
    }

    auto fsol = new FlowSolution(snapshot, nBlocks, simTime, dir, make_kdtree);
    flowSolutionStore ~= pushObj!(FlowSolution, FlowSolutionMT)(L, fsol);
    return 1;
} // end newFlowSolution()

// Do NOT try to use a FlowSolution object after calling close() on it.
extern(C) int closeFlowSolution(lua_State *L)
{
    auto fsol = checkFlowSolution(L, 1);
    fsol.releaseMemory();
    GC.collect();
    return 0;
}

extern(C) int find_enclosing_cell_from_lua(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
        string errMsg = "Error in call to FlowSolution:find_enclosing_cell.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }

    double x;
    lua_getfield(L, 2, "x");
    if ( lua_isnil(L, -1) ) {
        x = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
        x = to!double(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_enclosing_cell.";
        errMsg ~= " A field for x was found, but the content was not valid.";
        errMsg ~= " The x field, if given, should be a double.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1); // dispose of x item
    //
    double y;
    lua_getfield(L, 2, "y");
    if ( lua_isnil(L, -1) ) {
        y = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
        y = to!double(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_enclosing_cell.";
        errMsg ~= " A field for y was found, but the content was not valid.";
        errMsg ~= " The y field, if given, should be a double.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1); // dispose of y item
    //
    double z;
    lua_getfield(L, 2, "z");
    if ( lua_isnil(L, -1) ) {
        z = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
        z = to!double(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_enclosing_cell.";
        errMsg ~= " A field for z was found, but the content was not valid.";
        errMsg ~= " The z field, if given, should be a double.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1); // dispose of z item

    auto cell_identity = fsol.find_enclosing_cell(x, y, z);

    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    if (cell_identity[2] != 0) {
        lua_pushnumber(L, cell_identity[0]); lua_setfield(L, tblIndx, "ib");
        lua_pushnumber(L, cell_identity[1]); lua_setfield(L, tblIndx, "i");
    } else {
        lua_pushnil(L); lua_setfield(L, tblIndx, "ib");
        lua_pushnil(L); lua_setfield(L, tblIndx, "i");
    }
    return 1; // Just the table of labelled indices is left on the stack.
} // end find_enclosing_cell_from_lua()


extern(C) int find_enclosing_cells_along_line_from_lua(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
        string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }

    Vector3 p0;
    lua_getfield(L, 2, "p0");
    if (lua_istable(L, -1)) {
        // Seems that we have been given a table for the p0 value.
        double x;
        lua_getfield(L, -1, "x");
        if ( lua_isnil(L, -1) ) {
            x = 0.0;
        } else if ( lua_isnumber(L, -1) ) {
            x = to!double(luaL_checknumber(L, -1));
        } else {
            string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
            errMsg ~= " A field for p0.x was found, but the content was not valid.";
            errMsg ~= " The x field, if given, should be a double.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of x item
        //
        double y;
        lua_getfield(L, -1, "y");
        if ( lua_isnil(L, -1) ) {
            y = 0.0;
        } else if ( lua_isnumber(L, -1) ) {
            y = to!double(luaL_checknumber(L, -1));
        } else {
            string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
            errMsg ~= " A field for p0.y was found, but the content was not valid.";
            errMsg ~= " The y field, if given, should be a double.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of y item
        //
        double z;
        lua_getfield(L, -1, "z");
        if ( lua_isnil(L, -1) ) {
            z = 0.0;
        } else if ( lua_isnumber(L, -1) ) {
            z = to!double(luaL_checknumber(L, -1));
        } else {
            string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
            errMsg ~= " A field for p0.z was found, but the content was not valid.";
            errMsg ~= " The z field, if given, should be a double.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of z item
        p0 = Vector3(x, y, z);
    } else {
        // Was not a table, so try extracting a Vector3 pointer.
        p0 = *checkVector3(L, -1);
    }
    lua_pop(L, 1); // dispose of p0 item

    Vector3 p1;
    lua_getfield(L, 2, "p1");
    if (lua_istable(L, -1)) {
        // Seems that we have been given a table for the p1 value.
        double x;
        lua_getfield(L, -1, "x");
        if ( lua_isnil(L, -1) ) {
            x = 0.0;
        } else if ( lua_isnumber(L, -1) ) {
            x = to!double(luaL_checknumber(L, -1));
        } else {
            string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
            errMsg ~= " A field for p1.x was found, but the content was not valid.";
            errMsg ~= " The x field, if given, should be a double.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of x item
        //
        double y;
        lua_getfield(L, -1, "y");
        if ( lua_isnil(L, -1) ) {
            y = 0.0;
        } else if ( lua_isnumber(L, -1) ) {
            y = to!double(luaL_checknumber(L, -1));
        } else {
            string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
            errMsg ~= " A field for p1.y was found, but the content was not valid.";
            errMsg ~= " The y field, if given, should be a double.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of y item
        //
        double z;
        lua_getfield(L, -1, "z");
        if ( lua_isnil(L, -1) ) {
            z = 0.0;
        } else if ( lua_isnumber(L, -1) ) {
            z = to!double(luaL_checknumber(L, -1));
        } else {
            string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
            errMsg ~= " A field for p1.z was found, but the content was not valid.";
            errMsg ~= " The z field, if given, should be a double.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of z item
        p1 = Vector3(x, y, z);
    } else {
        // Was not a table, so try extracting a Vector3 pointer.
        p1 = *checkVector3(L, -1);
    }
    lua_pop(L, 1); // dispose of p1 item

    size_t n = 0;
    lua_getfield(L, -1, "n");
    if ( lua_isnil(L, -1) ) {
        n = 0;
    } else if ( lua_isnumber(L, -1) ) {
        n = to!size_t(lua_tointeger(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_enclosing_cells_along_line.";
        errMsg ~= " A field for n was found, but the content was not valid.";
        errMsg ~= " The n field, if given, should be a nonnegative integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1); // dispose of z item

    size_t[2][] cells_found;
    auto count = fsol.find_enclosing_cells_along_line(p0, p1, n, cells_found);

    lua_settop(L, 0); // clear stack
    lua_newtable(L); // start anonymous table { } for the array
    foreach(j, cell_identity; cells_found) {
        lua_newtable(L); // start anonymous table {} for the pair indices representing cellid
        auto tblIndx = lua_gettop(L);
        lua_pushnumber(L, cell_identity[0]); lua_setfield(L, tblIndx, "ib");
        lua_pushnumber(L, cell_identity[1]); lua_setfield(L, tblIndx, "i");
        lua_rawseti(L, 1, to!int(j+1)); // push the cellid pair into the array
    }
    return 1; // Just the array of tables of labelled indices is left on the stack.
} // end find_enclosing_cells_along_line_from_lua()


extern(C) int find_nearest_cell_centre_from_lua(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
        string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }

    double x;
    lua_getfield(L, 2, "x");
    if ( lua_isnil(L, -1) ) {
        x = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
        x = to!double(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
        errMsg ~= " A field for x was found, but the content was not valid.";
        errMsg ~= " The x field, if given, should be a double.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    double y;
    lua_getfield(L, 2, "y");
    if ( lua_isnil(L, -1) ) {
        y = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
        y = to!double(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
        errMsg ~= " A field for y was found, but the content was not valid.";
        errMsg ~= " The y field, if given, should be a double.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    double z;
    lua_getfield(L, 2, "z");
    if ( lua_isnil(L, -1) ) {
        z = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
        z = to!double(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
        errMsg ~= " A field for z was found, but the content was not valid.";
        errMsg ~= " The z field, if given, should be a double.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    auto indices = fsol.find_nearest_cell_centre(x, y, z);

    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    lua_pushnumber(L, indices[0]);
    lua_setfield(L, tblIndx, "ib");
    lua_pushnumber(L, indices[1]);
    lua_setfield(L, tblIndx, "i");
    return 1; // Just the table of indices is left on the stack.
} // end find_nearest_cell_centre_from_lua()

extern(C) int get_sim_time(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.simTime);
    return 1;
}

extern(C) int get_nic(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    int ib = to!int(luaL_checknumber(L, 2));
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.flowBlocks[ib].nic);
    return 1;
}

extern(C) int get_njc(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    int ib = to!int(luaL_checknumber(L, 2));
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.flowBlocks[ib].njc);
    return 1;
}

extern(C) int get_nkc(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    int ib = to!int(luaL_checknumber(L, 2));
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.flowBlocks[ib].nkc);
    return 1;
}


extern(C) int get_vtx_pos(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
        string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }

    int ib;
    lua_getfield(L, 2, "ib");
    if ( lua_isnil(L, -1) ) {
        ib = 0;
    } else if ( lua_isnumber(L, -1) ) {
        ib = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
        errMsg ~= " A field for ib was found, but the content was not valid.";
        errMsg ~= " The ib field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int i;
    lua_getfield(L, 2, "i");
    if ( lua_isnil(L, -1) ) {
        i = 0;
    } else if ( lua_isnumber(L, -1) ) {
        i = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
        errMsg ~= " A field for i was found, but the content was not valid.";
        errMsg ~= " The i field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int j;
    lua_getfield(L, 2, "j");
    if ( lua_isnil(L, -1) ) {
        j = 0;
    } else if ( lua_isnumber(L, -1) ) {
        j = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
        errMsg ~= " A field for j was found, but the content was not valid.";
        errMsg ~= " The j field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int k;
    lua_getfield(L, 2, "k");
    if ( lua_isnil(L, -1) ) {
        k = 0;
    } else if ( lua_isnumber(L, -1) ) {
        k = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
        errMsg ~= " A field for k was found, but the content was not valid.";
        errMsg ~= " The k field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    Vector3 vtx = *(fsol.gridBlocks[ib][i, j, k]);
    lua_settop(L, 0);
    return pushVector3(L, vtx);
} // end get_vtx_pos()


extern(C) int get_var_names(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    foreach (ivar; 0 .. fsol.flowBlocks[0].variableNames.length) {
        lua_pushstring(L, fsol.flowBlocks[0].variableNames[ivar].toStringz);
        lua_rawseti(L, tblIndx, to!int(ivar+1)); // Lua table indexing starts from 1
    }
    return 1; // Just the table of indices is left on the stack.
} // end get_var_names()


void plottingTableToFlowStateTable(lua_State *L)
{
    // Assume that table is at top-of-stack.
    // We are going to massage the contents of the table
    // such that it is sutiable to pass to
    // luaFlowState:fromTable()
    int tblIdx = lua_gettop(L);
    auto managedGasModel = GlobalConfig.gmodel_master;
    auto n_species = managedGasModel.n_species;
    auto n_modes = managedGasModel.n_modes;
    auto n_turb  = to!int(GlobalConfig.turb_model.nturb);

    // 0. Set a type string so that we may later identify this table as
    // having all the relevant data for making a FlowState object.
    lua_pushstring(L, toStringz("CellData"));
    lua_setfield(L, tblIdx, toStringz("myType"));

    // 1. Convert velocities
    lua_getfield(L, tblIdx, "vel.x");
    lua_setfield(L, tblIdx, "velx");
    lua_getfield(L, tblIdx, "vel.y");
    lua_setfield(L, tblIdx, "vely");
    lua_getfield(L, tblIdx, "vel.z");
    lua_setfield(L, tblIdx, "velz");

    // 2. Convert magnetic field components
    lua_getfield(L, tblIdx, "B.x");
    lua_setfield(L, tblIdx, "Bx");
    lua_getfield(L, tblIdx, "B.y");
    lua_setfield(L, tblIdx, "By");
    lua_getfield(L, tblIdx, "B.z");
    lua_setfield(L, tblIdx, "Bz");

    // 3. Convert temperatures
    lua_newtable(L);
    foreach ( imode; 0 .. n_modes ) {
        string modeName = managedGasModel.energy_mode_name(imode);
        string key = format("T-%s", modeName);
        lua_getfield(L, tblIdx, toStringz(key));
        lua_rawseti(L, -2, imode+1);
    }
    lua_setfield(L, tblIdx, "T_modes");

    // 4. Convert mass fractions
    lua_newtable(L);
    foreach ( isp; 0 .. n_species ) {
        string spName = managedGasModel.species_name(isp);
        string key = format("massf-%s", spName);
        lua_getfield(L, tblIdx, toStringz(key));
        lua_setfield(L, -2, toStringz(spName));
    }
    lua_setfield(L, tblIdx, "massf");

    // 5. Convert turbulence variables
    lua_newtable(L);
    foreach ( iturb; 0 .. n_turb ) {
        string turbName = GlobalConfig.turb_model.primitive_variable_name(iturb);
        string key = format("tq-%s", turbName);
        lua_getfield(L, tblIdx, toStringz(key));
        lua_rawseti(L, -2, iturb+1);
    }
    lua_setfield(L, tblIdx, "turb");

}

extern(C) int get_cell_data(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
        string errMsg = "Error in call to FlowSolution:get_cell_data.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }

    string fmt;
    lua_getfield(L, 2, "fmt");
    if ( lua_isnil(L, -1) ) {
        fmt = "Plotting";
    }
    else if ( lua_isstring(L, -1) ) {
        fmt = to!string(luaL_checkstring(L, -1));
        if ( (fmt != "Plotting") && (fmt != "FlowState") ) {
            string errMsg = "Error in call to FlowSolution:get_cell_data.";
            errMsg ~= " The fmt field should be one of 'Plotting' or 'FlowState'.";
            errMsg ~= " The fmt string received was: " ~ fmt;
            luaL_error(L, errMsg.toStringz);
        }
    }
    else {
        string errMsg = "Error in call to FlowSolution:get_cell_data.";
        errMsg ~= " A field for fmt was found, but the content was not valid.";
        errMsg ~= " A string was expected. Either 'Plotting' or 'FlowState' are acceptable strings.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int ib;
    lua_getfield(L, 2, "ib");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to FlowSolution:get_cell_data.\n";
        errMsg ~= " No field for ib was found.\n";
        errMsg ~= " A block number in field ib should be provided.\n";
        luaL_error(L, errMsg.toStringz);
    } else if ( lua_isnumber(L, -1) ) {
        ib = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_cell_data.";
        errMsg ~= " A field for ib was found, but the content was not valid.";
        errMsg ~= " The ib field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int i;
    lua_getfield(L, 2, "i");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to FlowSolution:get_cell_data.\n";
        errMsg ~= " No field for i was found.\n";
        errMsg ~= " A cell index number in field i should be provided.\n";
        luaL_error(L, errMsg.toStringz);
    } else if ( lua_isnumber(L, -1) ) {
        i = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_cell_data.";
        errMsg ~= " A field for i was found, but the content was not valid.";
        errMsg ~= " The i field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    // We make a decision at this point. If the "j" field is nil,
    // then we assume that the caller is using single-index (unstructured)
    // access to the cell data. In this case, we delegate the work,
    // return the result, and clean-up and leave.
    //
    // Otherwise we go onto gather the "j" and possibly "k" fields.

    int j;
    lua_getfield(L, 2, "j");
    if ( lua_isnil(L, -1) ) {
        lua_pop(L, 1);
        lua_settop(L, 0); // clear stack
        lua_newtable(L); // anonymous table { }
        auto tblIndx = lua_gettop(L);
        foreach (varName; fsol.flowBlocks[ib].variableNames) {
            lua_pushnumber(L, fsol.flowBlocks[ib][varName,i]);
            lua_setfield(L, tblIndx, varName.toStringz);
        }
        if ( fmt == "FlowState" ) {
            plottingTableToFlowStateTable(L);
        }
        return 1; // Just the table of indices is left on the stack.
    }

    // otherwise, assume structured-grid type access.
    if ( lua_isnumber(L, -1) ) {
        j = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_cell_data.";
        errMsg ~= " A field for j was found, but the content was not valid.";
        errMsg ~= " The j field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    int k;
    lua_getfield(L, 2, "k");
    if ( lua_isnil(L, -1) ) {
        k = 0;
    } else if ( lua_isnumber(L, -1) ) {
        k = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_cell_data.";
        errMsg ~= " A field for k was found, but the content was not valid.";
        errMsg ~= " The k field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    foreach (varName; fsol.flowBlocks[ib].variableNames) {
        lua_pushnumber(L, fsol.flowBlocks[ib][varName,i,j,k]);
        lua_setfield(L, tblIndx, varName.toStringz);
    }
    if ( fmt == "FlowState" ) {
        plottingTableToFlowStateTable(L);
    }

    return 1; // Just the table of indices is left on the stack.
} // end get_cell_data()


extern(C) int get_sgrid(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
        string errMsg = "Error in call to FlowSolution:get_sgrid.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }

    int ib;
    lua_getfield(L, 2, "ib");
    if ( lua_isnil(L, -1) ) {
        ib = 0;
    } else if ( lua_isnumber(L, -1) ) {
        ib = to!int(luaL_checknumber(L, -1));
    } else {
        string errMsg = "Error in call to FlowSolution:get_grid.";
        errMsg ~= " A field for ib was found, but the content was not valid.";
        errMsg ~= " The ib field, if given, should be an integer.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);

    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, cast(StructuredGrid) fsol.gridBlocks[ib]);
    return 1;

}

extern(C) int add_aux_variables(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if (!lua_istable(L, 2)) {
        string errMsg = "Error in call to FlowSolution:add_aux_variables.";
        errMsg ~= " A table is expected as first (and only) argument to the method.";
        luaL_error(L, errMsg.toStringz);
    }
    string[] names;
    int nEntries = to!int(lua_objlen(L, 2));
    foreach (i; 0 .. nEntries) {
        lua_rawgeti(L, 2, i+1);
        names ~= to!string(luaL_checkstring(L, -1));
        lua_pop(L, 1);
    }
    fsol.add_aux_variables(names);
    lua_settop(L, 0);
    return 0;
}

extern(C) int write_vtk_files(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    string plotDir = to!string(luaL_checkstring(L, 2));
    string plotName = to!string(luaL_checkstring(L, 3));
    fsol.write_vtk_files(plotDir, plotName, false);
    lua_settop(L, 0);
    return 0;
}

extern(C) int subtract_flow_solution(lua_State* L)
{
    auto fsol1 = checkFlowSolution(L, 1);
    auto fsol2 = checkFlowSolution(L, 2);
    fsol1.subtract(fsol2);
    lua_settop(L, 0);
    return 0;
}


void registerFlowSolution(lua_State* L)
{
    luaL_newmetatable(L, FlowSolutionMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newFlowSolution);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(FlowSolution, FlowSolutionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &closeFlowSolution);
    lua_setfield(L, -2, "close");
    lua_pushcfunction(L, &find_enclosing_cell_from_lua);
    lua_setfield(L, -2, "find_enclosing_cell");
    lua_pushcfunction(L, &find_enclosing_cells_along_line_from_lua);
    lua_setfield(L, -2, "find_enclosing_cells_along_line");
    lua_pushcfunction(L, &find_nearest_cell_centre_from_lua);
    lua_setfield(L, -2, "find_nearest_cell_centre");
    lua_pushcfunction(L, &get_sim_time);
    lua_setfield(L, -2, "sim_time");
    lua_pushcfunction(L, &get_nic);
    lua_setfield(L, -2, "get_nic");
    // Alias get_nic == nic
    lua_pushcfunction(L, &get_nic);
    lua_setfield(L, -2, "nic");
    lua_pushcfunction(L, &get_njc);
    lua_setfield(L, -2, "get_njc");
    // Alias get_njc == njc
    lua_pushcfunction(L, &get_njc);
    lua_setfield(L, -2, "njc");
    lua_pushcfunction(L, &get_nkc);
    lua_setfield(L, -2, "get_nkc");
    // Alias get_nkc == nkc
    lua_pushcfunction(L, &get_nkc);
    lua_setfield(L, -2, "nkc");
    lua_pushcfunction(L, &get_vtx_pos);
    lua_setfield(L, -2, "get_vtx");
    // Alias get_vtx == vtx
    lua_pushcfunction(L, &get_vtx_pos);
    lua_setfield(L, -2, "vtx");
    lua_pushcfunction(L, &get_var_names);
    lua_setfield(L, -2, "get_var_names");
    // Alias get_var_names == var_names
    lua_pushcfunction(L, &get_var_names);
    lua_setfield(L, -2, "var_names");
    lua_pushcfunction(L, &get_cell_data);
    lua_setfield(L, -2, "get_cell_data");
    // Alias get_cell_data == cell_data
    lua_pushcfunction(L, &get_cell_data);
    lua_setfield(L, -2, "cell_data");
    lua_pushcfunction(L, &get_sgrid);
    lua_setfield(L, -2, "get_sgrid");
    // Alias get_sgrid == sgrid
    lua_pushcfunction(L, &get_sgrid);
    lua_setfield(L, -2, "sgrid");
    lua_pushcfunction(L, &add_aux_variables);
    lua_setfield(L, -2, "add_aux_variables");
    lua_pushcfunction(L, &write_vtk_files);
    lua_setfield(L, -2, "write_vtk_files");
    lua_pushcfunction(L, &subtract_flow_solution);
    lua_setfield(L, -2, "subtract");
    // Make class visible
    lua_setglobal(L, FlowSolutionMT.toStringz);
} // end registerFlowSolution()
