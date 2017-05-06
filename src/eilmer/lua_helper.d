// lua_helper.d
//
// A place to put some frequently used functions
// for interacting with the Lua stack. Some of the functions
// are for interacting in the D code. Other functions
// (marked extern(C)) are functions made available in the
// Lua script.
//
// RG & PJ 2015-03-17 -- First hack (with Guiness in hand)
import std.stdio;
import std.conv;
import util.lua;

import gas;
import fvcell;
import fvinterface;
import luaflowstate;
import globaldata;

// -----------------------------------------------------
// Convenience functions for user's Lua script
// [TODO] [FIXME] the following functions won't work in parallel loops
// because the gasBlocks array probably won't be initialized correctly
// for any thread other than the main thread.


void setSampleHelperFunctions(lua_State *L)
{
    lua_pushcfunction(L, &luafn_sampleFace);
    lua_setglobal(L, "sampleFace");
    lua_pushcfunction(L, &luafn_sampleFlow);
    lua_setglobal(L, "sampleCell");
    lua_pushcfunction(L, &luafn_sampleFlow);
    lua_setglobal(L, "sampleFlow"); // alias for sampleCell
}

extern(C) int luafn_sampleFlow(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);

    // Grab the appropriate cell
    auto cell = gasBlocks[blkId].get_cell(i, j, k);
    
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushCellToTable(L, tblIdx, cell, 0, gasBlocks[blkId].myConfig.gmodel);
    return 1;
} // end luafn_sampleFlow()

extern(C) int luafn_sampleFace(lua_State *L)
{
    // Get arguments from lua_stack
    string which_face = to!string(lua_tostring(L, 1));
    auto blkId = lua_tointeger(L, 2);
    auto i = lua_tointeger(L, 3);
    auto j = lua_tointeger(L, 4);
    auto k = lua_tointeger(L, 5);

    FVInterface face;
    // Grab the appropriate face
    switch (which_face) {
    case "i": face = gasBlocks[blkId].get_ifi(i, j, k); break;
    case "j": face = gasBlocks[blkId].get_ifj(i, j, k); break;
    case "k": face = gasBlocks[blkId].get_ifk(i, j, k); break;
    case "u": face = gasBlocks[blkId].get_ifi(i, j, k); break; // unstructured grid
    default:  face = gasBlocks[blkId].get_ifi(i, j, k);
    }
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushFaceToTable(L, tblIdx, face, 0, gasBlocks[blkId].myConfig.gmodel);
    return 1;
} // end luafn_sampleIFace()

// -----------------------------------------------------
// D code functions

/**
 * Push the interesting data from a cell to a Lua table
 *
 */
void pushCellToTable(lua_State* L, int tblIdx, ref const(FVCell) cell, 
		     size_t gtl, GasModel gmodel)
{
    lua_pushnumber(L, cell.pos[gtl].x); lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, cell.pos[gtl].y); lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, cell.pos[gtl].z); lua_setfield(L, tblIdx, "z");
    lua_pushnumber(L, cell.volume[gtl]); lua_setfield(L, tblIdx, "vol");
    pushFlowStateToTable(L, tblIdx, cell.fs, gmodel);
}

void pushFaceToTable(lua_State* L, int tblIdx, ref const(FVInterface) face, 
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
}






