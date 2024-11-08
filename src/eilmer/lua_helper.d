// lua_helper.d
//
// A place to put some frequently used functions
// for interacting with the Lua stack. Some of the functions
// are for interacting in the D code. Other functions
// (marked extern(C)) are functions made available in the
// Lua script.
//
// RG & PJ 2015-03-17 -- First hack (with Guiness in hand)
// 2022-04-13 -- Here, hold my beer while I add stuff for the SolidCell.

import std.stdio;
import std.conv;
import std.string;
import std.algorithm;

import util.lua;
import ntypes.complex;
import nm.number;
import geom.luawrap;
import gas;
import gas.luagas_model;
import nm.luabbla;
import geom: gridTypeName, Grid_t;

import fvcell;
import fvinterface;
import luaflowstate;
import luaflowsolution;
import gasdyn.luaidealgasflow;
import gasdyn.luagasflow;
import globalconfig;
import globaldata;
import solidfvcell;
import ssolidblock;
import sfluidblock;
import ufluidblock;
import fluidblock;

struct LuaEnvOptions {
    bool withGlobalConfig = true;
    bool withGeom = true;
    bool withGas = true;
    bool withFlow = true;
    bool withNumerics = true;
};

lua_State* initLuaStateForPrep()
{
    auto L = luaL_newstate();
    luaL_openlibs(L);
    LuaEnvOptions luaOpts;
    luaOpts.withGlobalConfig = true;
    luaOpts.withGeom = true;
    luaOpts.withGas = true;
    luaOpts.withFlow = true;
    luaOpts.withNumerics = true;
    registerLuaEnvironment(L, luaOpts);
    return L;
}

// ------------------------------------------------------------------
void registerLuaEnvironment(lua_State *L, LuaEnvOptions opt)
{
    if (opt.withGlobalConfig) registerGlobalConfig(L);
    if (opt.withGeom) {
        registerVector3(L);
        registerGeomNomenclature(L);
        registerPaths(L);
        registerGpathUtils(L);
        registerSurfaces(L);
        registerVolumes(L);
        registerUnivariateFunctions(L);
        registerStructuredGrid(L);
        registerUnstructuredGrid(L);
        registerSketch(L);
    }
    if (opt.withGas) {
        registerGasModel(L);
    }
    if (opt.withFlow) {
        registerFlowSolution(L);
        registerFlowState(L);
        registeridealgasflowFunctions(L);
        registergasflowFunctions(L);
    }
    if (opt.withNumerics) {
        registerBBLA(L);
    }
}

// -----------------------------------------------------
// Functions to synchronise an array in the Dlang domain with a table in Lua

void push_array_to_Lua(T)(lua_State *L, const(T[]) array_in_dlang, string name_in_Lua)
{
    lua_getglobal(L, name_in_Lua.toStringz);
    if (!lua_istable(L, -1)) {
        // TOS is not a table, so dispose of it and make a fresh table.
        lua_pop(L, 1);
        lua_newtable(L);
        lua_setglobal(L, name_in_Lua.toStringz);
        lua_getglobal(L, name_in_Lua.toStringz);
    }
    assert(lua_istable(L, -1), format("Did not find Lua table %s", name_in_Lua));
    foreach (i, elem; array_in_dlang) {
        lua_pushnumber(L, elem);
        lua_rawseti(L, -2, to!int(i+1));
    }
    lua_pop(L, 1); // dismiss the table
}

void fill_array_from_Lua_table_at_locn(T)(lua_State *L, T[] array_in_dlang, int locn)
{
    assert(lua_istable(L, locn), "Did not find Lua table at stack location %d");
    foreach (i; 0 .. array_in_dlang.length) {
        lua_rawgeti(L, locn, to!int(i+1)); // get an item to top of stack
        array_in_dlang[i] = (lua_isnumber(L, -1)) ? to!double(lua_tonumber(L, -1)) : 0.0;
        lua_pop(L, 1); // discard item
    }
    // We leave the table sitting at its original location on the stack.
}

void fill_array_from_Lua(T)(lua_State *L, T[] array_in_dlang, string name_in_Lua)
{
    lua_getglobal(L, name_in_Lua.toStringz);
    assert(lua_istable(L, -1), format("Did not find Lua table %s", name_in_Lua));
    fill_array_from_Lua_table_at_locn(L, array_in_dlang, -1);
    lua_pop(L, 1); // dismiss the table
}

// -----------------------------------------------------
// Convenience functions for user's Lua script

void setSampleHelperFunctions(lua_State *L)
{
    lua_pushcfunction(L, &luafn_infoFluidBlock);
    lua_setglobal(L, "infoFluidBlock");
    lua_pushcfunction(L, &luafn_sampleFluidFace);
    lua_setglobal(L, "sampleFluidFace");
    lua_pushcfunction(L, &luafn_sampleFluidFace);
    lua_setglobal(L, "sampleFace"); // alias for sampleFluidFace; [TODO] remove eventually
    lua_pushcfunction(L, &luafn_sampleFluidCell);
    lua_setglobal(L, "sampleFluidCell");
    lua_pushcfunction(L, &luafn_sampleFluidCell);
    lua_setglobal(L, "sampleFlow"); // alias for sampleFluidCell; [TODO] remove eventually
    lua_pushcfunction(L, &luafn_setBxyzInFluidCell);
    lua_setglobal(L, "setBxyzInFluidCell");
    lua_pushcfunction(L, &luafn_runTimeLoads);
    lua_setglobal(L, "getRunTimeLoads");
    //
    lua_pushcfunction(L, &luafn_infoSolidBlock);
    lua_setglobal(L, "infoSolidBlock");
    lua_pushcfunction(L, &luafn_sampleSolidCell);
    lua_setglobal(L, "sampleSolidCell");
}

extern(C) int luafn_infoFluidBlock(lua_State *L)
{
    // Expect FluidBlock index on the lua_stack.
    auto blkId = lua_tointeger(L, 1);
    auto blk = cast(FluidBlock) globalBlocks[blkId];
    assert(blk !is null, "Oops, this should be a FluidBlock object.");
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    lua_pushinteger(L, GlobalConfig.dimensions); lua_setfield(L, tblIdx, "dimensions");
    lua_pushstring(L, blk.label.toStringz); lua_setfield(L, tblIdx, "label");
    lua_pushstring(L, gridTypeName(blk.grid_type).toStringz); lua_setfield(L, tblIdx, "grid_type");
    if (blk.grid_type == Grid_t.structured_grid) {
        auto sblk = cast(SFluidBlock) blk;
        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
        // For a structured_grid
        lua_pushinteger(L, sblk.nic); lua_setfield(L, tblIdx, "nicell");
        lua_pushinteger(L, sblk.njc); lua_setfield(L, tblIdx, "njcell");
        lua_pushinteger(L, sblk.nkc); lua_setfield(L, tblIdx, "nkcell");
        lua_pushinteger(L, 0); lua_setfield(L, tblIdx, "imin");
        lua_pushinteger(L, 0); lua_setfield(L, tblIdx, "jmin");
        lua_pushinteger(L, 0); lua_setfield(L, tblIdx, "kmin");
        lua_pushinteger(L, sblk.nic-1); lua_setfield(L, tblIdx, "imax");
        lua_pushinteger(L, sblk.njc-1); lua_setfield(L, tblIdx, "jmax");
        lua_pushinteger(L, sblk.nkc-1); lua_setfield(L, tblIdx, "kmax");
        string[] corner_names;
        if (GlobalConfig.dimensions == 3) {
            corner_names = ["p000","p100","p110","p010","p001","p101","p111","p011"];
        } else {
            corner_names = ["p00","p10","p11","p01"];
        }
        foreach (i; 0 .. corner_names.length) {
            lua_newtable(L);
            lua_pushnumber(L, sblk.corner_coords[i*3+0]); lua_setfield(L, -2, "x");
            lua_pushnumber(L, sblk.corner_coords[i*3+1]); lua_setfield(L, -2, "y");
            lua_pushnumber(L, sblk.corner_coords[i*3+2]); lua_setfield(L, -2, "z");
            lua_setfield(L, tblIdx, corner_names[i].toStringz);
        }
    }
    // For an unstructured_grid or structured_grid
    lua_pushinteger(L, blk.cells.length); lua_setfield(L, tblIdx, "ncells");
    lua_pushinteger(L, blk.faces.length); lua_setfield(L, tblIdx, "nfaces");
    lua_pushinteger(L, blk.vertices.length); lua_setfield(L, tblIdx, "nvertices");
    return 1;
} // end luafn_infoFluidBlock()

extern(C) int luafn_sampleFluidCell(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
        string msg = format("Block id %d is not local to process.", blkId);
        luaL_error(L, msg.toStringz);
    }
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);
    //
    // Grab the appropriate cell
    auto sblk = cast(SFluidBlock) globalBlocks[blkId];
    FVCell cell;
    if (sblk) {
        try {
            cell = sblk.get_cell(i, j, k);
        } catch (Exception e) {
            string msg = format("Failed to locate cell[%d,%d,%d] in block %d.", i, j, k, blkId);
            luaL_error(L, msg.toStringz);
        }
    } else {
        string msg = "Not implemented.";
        msg ~= " You have asked for an ijk-index cell in an unstructured-grid block.";
        luaL_error(L, msg.toStringz);
    }
    //
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    auto blk = cast(FluidBlock) globalBlocks[blkId];
    assert(blk !is null, "Oops, this should be a FluidBlock object.");
    pushFluidCellToTable(L, tblIdx, cell, 0, blk.myConfig);
    return 1;
} // end luafn_sampleFluidCell()

extern(C) int luafn_setBxyzInFluidCell(lua_State *L)
// In Lua: setBxyzInFluidCell(blkid, i, j, k, Bx, By, Bz)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
        string msg = format("Block id %d is not local to process.", blkId);
        luaL_error(L, msg.toStringz);
    }
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);
    double Bx = lua_tonumber(L, 5);
    double By = lua_tonumber(L, 6);
    double Bz = lua_tonumber(L, 7);
    //
    version(MHD) {
        // Actually set the magnetic field components.
        // Grab the appropriate cell
        auto sblk = cast(SFluidBlock) globalBlocks[blkId];
        FVCell cell;
        if (sblk) {
            try {
                cell = sblk.get_cell(i, j, k);
            } catch (Exception e) {
                string msg = format("Failed to locate cell[%d,%d,%d] in block %d.", i, j, k, blkId);
                luaL_error(L, msg.toStringz);
            }
        } else {
            string msg = "Not implemented.";
            msg ~= " You have asked for an ijk-index cell in an unstructured-grid block.";
            luaL_error(L, msg.toStringz);
        }
        cell.fs.B.set(Bx, By, Bz);
    }
    return 0;
} // end luafn_setBxyzInFluidCell()

extern(C) int luafn_sampleFluidFace(lua_State *L)
{
    // Get arguments from lua_stack
    string which_face = to!string(lua_tostring(L, 1));
    auto blkId = lua_tointeger(L, 2);
    if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
        string msg = format("Block id %d is not local to process.", blkId);
        luaL_error(L, msg.toStringz);
    }
    auto i = lua_tointeger(L, 3);
    auto j = lua_tointeger(L, 4);
    auto k = lua_tointeger(L, 5);

    FVInterface face;
    // Grab the appropriate face
    auto sblk = cast(SFluidBlock) globalBlocks[blkId];
    auto ublk = cast(UFluidBlock) globalBlocks[blkId];
    try {
        switch (which_face) {
        case "i": face = sblk.get_ifi(i, j, k); break;
        case "j": face = sblk.get_ifj(i, j, k); break;
        case "k": face = sblk.get_ifk(i, j, k); break;
        case "u": face = ublk.faces[i]; break; // unstructured grid
        default:
            string msg = "You have asked for an unknown type of face.";
            luaL_error(L, msg.toStringz);
        }
    } catch (Exception e) {
        string msg = format("Failed to locate face[%d,%d,%d] in block %d.", i, j, k, blkId);
        luaL_error(L, msg.toStringz);
    }
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    auto blk = cast(FluidBlock) globalBlocks[blkId];
    assert(blk !is null, "Oops, this should be a FluidBlock object.");
    pushFluidFaceToTable(L, tblIdx, face, 0, blk.myConfig.gmodel);
    return 1;
} // end luafn_sampleFluidFace()

extern(C) int luafn_runTimeLoads(lua_State *L)
{
    string loadsGroup = to!string(lua_tostring(L, 1));
    size_t grpIdx;
    size_t* grpIdxPtr = (loadsGroup in runTimeLoadsByName);
    if (grpIdxPtr !is null) {
        grpIdx = *grpIdxPtr;
    }
    else {
        string msg = "You have asked for an unknown loads group: "~loadsGroup;
        luaL_error(L, msg.toStringz);
    }
    // Set force as table {x=.., y=..., z=...}
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    lua_pushnumber(L, runTimeLoads[grpIdx].resultantForce.x);
    lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, runTimeLoads[grpIdx].resultantForce.y);
    lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, runTimeLoads[grpIdx].resultantForce.z);
    lua_setfield(L, tblIdx, "z");
     // Set moment as table {x=.., y=..., z=...}
    lua_newtable(L);
    tblIdx = lua_gettop(L);
    lua_pushnumber(L, runTimeLoads[grpIdx].resultantMoment.x);
    lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, runTimeLoads[grpIdx].resultantMoment.y);
    lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, runTimeLoads[grpIdx].resultantMoment.z);
    lua_setfield(L, tblIdx, "z");

    return 2;
} // end luafn_runTimeLoads()


// -----------------------------------------------------
// D code functions

/**
 * Push the interesting data from a FVCell and FVInterface to a Lua table
 *
 */
void pushFluidCellToTable(lua_State* L, int tblIdx, ref const(FVCell) cell,
                          size_t gtl, LocalConfig myConfig)
{
    lua_pushnumber(L, cell.pos[gtl].x); lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, cell.pos[gtl].y); lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, cell.pos[gtl].z); lua_setfield(L, tblIdx, "z");
    lua_pushnumber(L, cell.volume[gtl]); lua_setfield(L, tblIdx, "vol");
    lua_pushnumber(L, cell.iLength); lua_setfield(L, tblIdx, "iLength");
    lua_pushnumber(L, cell.jLength); lua_setfield(L, tblIdx, "jLength");
    lua_pushnumber(L, cell.kLength); lua_setfield(L, tblIdx, "kLength");
    pushFlowStateToTable(L, tblIdx, cell.fs, myConfig.gmodel);
    // For Carrie Xie 2022-05-24, we want access to the user-defined energy source term
    // when we sample the Fluid cell during the UDF evaluation for the corresponding solid cell.
    auto cqi = myConfig.cqi;
    lua_pushnumber(L, cell.Qudf[cqi.totEnergy]); lua_setfield(L, tblIdx, "Qudf_totEnergy");
} // end pushFluidCellToTable()

void pushFluidFaceToTable(lua_State* L, int tblIdx, ref const(FVInterface) face,
                          size_t gtl, GasModel gmodel)
{
    lua_pushnumber(L, face.pos.x); lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, face.pos.y); lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, face.pos.z); lua_setfield(L, tblIdx, "z");
    lua_pushnumber(L, face.area[gtl]); lua_setfield(L, tblIdx, "area");
    lua_pushnumber(L, face.n.x); lua_setfield(L, tblIdx, "nx");
    lua_pushnumber(L, face.n.y); lua_setfield(L, tblIdx, "ny");
    lua_pushnumber(L, face.n.z); lua_setfield(L, tblIdx, "nz");
    lua_pushnumber(L, face.t1.x); lua_setfield(L, tblIdx, "t1x");
    lua_pushnumber(L, face.t1.y); lua_setfield(L, tblIdx, "t1y");
    lua_pushnumber(L, face.t1.z); lua_setfield(L, tblIdx, "t1z");
    lua_pushnumber(L, face.t2.x); lua_setfield(L, tblIdx, "t2x");
    lua_pushnumber(L, face.t2.y); lua_setfield(L, tblIdx, "t2y");
    lua_pushnumber(L, face.t2.z); lua_setfield(L, tblIdx, "t2z");
    lua_pushnumber(L, face.Ybar); lua_setfield(L, tblIdx, "Ybar");
    lua_pushnumber(L, face.gvel.x); lua_setfield(L, tblIdx, "gvelx");
    lua_pushnumber(L, face.gvel.y); lua_setfield(L, tblIdx, "gvely");
    lua_pushnumber(L, face.gvel.z); lua_setfield(L, tblIdx, "gvelz");
    pushFlowStateToTable(L, tblIdx, face.fs, gmodel);
} // end pushFluidFaceToTable()

// ----------------------------------------------------------------------
// Functions related to solid domains

extern(C) int luafn_infoSolidBlock(lua_State *L)
{
    // Expect SSolidBlock index on the lua_stack.
    auto blkId = lua_tointeger(L, 1);
    auto blk = cast(SSolidBlock) globalBlocks[blkId];
    assert(blk !is null, "Oops, this should be a SSolidBlock object.");
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    lua_pushinteger(L, GlobalConfig.dimensions); lua_setfield(L, tblIdx, "dimensions");
    lua_pushstring(L, blk.label.toStringz); lua_setfield(L, tblIdx, "label");
    // Always a structured_grid
    lua_pushinteger(L, blk.nicell); lua_setfield(L, tblIdx, "nicell");
    lua_pushinteger(L, blk.njcell); lua_setfield(L, tblIdx, "njcell");
    lua_pushinteger(L, blk.nkcell); lua_setfield(L, tblIdx, "nkcell");
    lua_pushinteger(L, blk.imin); lua_setfield(L, tblIdx, "imin");
    lua_pushinteger(L, blk.jmin); lua_setfield(L, tblIdx, "jmin");
    lua_pushinteger(L, blk.kmin); lua_setfield(L, tblIdx, "kmin");
    lua_pushinteger(L, blk.imax); lua_setfield(L, tblIdx, "imax");
    lua_pushinteger(L, blk.jmax); lua_setfield(L, tblIdx, "jmax");
    lua_pushinteger(L, blk.kmax); lua_setfield(L, tblIdx, "kmax");
    lua_pushinteger(L, blk.cells.length); lua_setfield(L, tblIdx, "ncells");
    return 1;
} // end luafn_infoSolidBlock()

extern(C) int luafn_sampleSolidCell(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    if (!canFind(GlobalConfig.localSolidBlockIds, blkId)) {
        string msg = format("Block id %d is not local to process.", blkId);
        luaL_error(L, msg.toStringz);
    }
    auto blk = cast(SSolidBlock) globalBlocks[blkId];
    assert(blk !is null, "Oops, this should be a SSolidBlock object.");
    auto i = lua_tointeger(L, 2);
    assert(i < blk.nicell, "i-index out of range.");
    auto j = lua_tointeger(L, 3);
    assert(j < blk.njcell, "j-index out of range.");
    size_t k = 0;
    if (GlobalConfig.dimensions == 3) {
        k = lua_tointeger(L, 4);
        assert(k < blk.nkcell, "k-index out of range.");
    }
    //
    // Grab the appropriate cell.
    // Note that we want the indexing in the Lua domain to look like that for SFluidBlocks.
    SolidFVCell cell;
    try {
        if (GlobalConfig.dimensions == 3) {
            cell = blk.getCell(i+blk.imin, j+blk.jmin, k+blk.kmin);
        } else {
            cell = blk.getCell(i+blk.imin, j+blk.jmin);
        }
    } catch (Exception e) {
        string msg = format("Failed to locate cell[%d,%d,%d] in block %d.", i, j, k, blkId);
        luaL_error(L, msg.toStringz);
    }
    //
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushSolidCellToTable(L, tblIdx, cell);
    return 1;
} // end luafn_sampleSolidCell()

void pushSolidCellToTable(lua_State* L, int tblIdx, ref const(SolidFVCell) cell)
{
    lua_pushnumber(L, cell.pos.x); lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, cell.pos.y); lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, cell.pos.z); lua_setfield(L, tblIdx, "z");
    lua_pushnumber(L, cell.volume); lua_setfield(L, tblIdx, "vol");
    lua_pushnumber(L, cell.T); lua_setfield(L, tblIdx, "T");
    lua_pushnumber(L, cell.sp.rho); lua_setfield(L, tblIdx, "rho");
    lua_pushnumber(L, cell.sp.Cp); lua_setfield(L, tblIdx, "Cp");
    lua_pushnumber(L, cell.sp.k); lua_setfield(L, tblIdx, "k");

    lua_pushnumber(L, cell.sp.k11); lua_setfield(L, tblIdx, "k11");
    lua_pushnumber(L, cell.sp.k12); lua_setfield(L, tblIdx, "k12");
    lua_pushnumber(L, cell.sp.k13); lua_setfield(L, tblIdx, "k13");
    lua_pushnumber(L, cell.sp.k21); lua_setfield(L, tblIdx, "k21");
    lua_pushnumber(L, cell.sp.k22); lua_setfield(L, tblIdx, "k22");
    lua_pushnumber(L, cell.sp.k23); lua_setfield(L, tblIdx, "k23");
    lua_pushnumber(L, cell.sp.k31); lua_setfield(L, tblIdx, "k31");
    lua_pushnumber(L, cell.sp.k32); lua_setfield(L, tblIdx, "k32");
    lua_pushnumber(L, cell.sp.k33); lua_setfield(L, tblIdx, "k33");

} // end pushSolidCellToTable()
