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
import util.lua;
import fvcell;
import luaflowstate;
import globaldata;

// -----------------------------------------------------
// Convenience functions for user's Lua script

// [TODO] [FIXME] the following function won't work in parallel loops
// because the gasBlocks array probably won't be initialized correctly
// for any thread other than the main thread.
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
    pushCellToTable(L, tblIdx, cell, 0);
    return 1;
}

// -----------------------------------------------------
// D code functions

/**
 * Push the interesting data from a cell to a Lua table
 *
 */
void pushCellToTable(lua_State* L, int tblIdx, in FVCell cell, size_t gtl)
{
    lua_pushnumber(L, cell.pos[gtl].x); lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, cell.pos[gtl].y); lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, cell.pos[gtl].z); lua_setfield(L, tblIdx, "z");
    lua_pushnumber(L, cell.volume[gtl]); lua_setfield(L, tblIdx, "vol");
    pushFlowStateToTable(L, tblIdx, cell.fs);
}






