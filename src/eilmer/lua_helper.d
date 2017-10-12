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
import std.string;
import util.lua;

import grid;
import gas;
import fvcell;
import fvinterface;
import luaflowstate;
import globalconfig;
import globaldata;
import solidfvcell;
import sblock: SBlock;

// -----------------------------------------------------
// Convenience functions for user's Lua script
// [TODO] [FIXME] the following functions won't work in parallel loops
// because the gasBlocks array probably won't be initialized correctly
// for any thread other than the main thread.


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
}

extern(C) int luafn_infoFluidBlock(lua_State *L)
{
    // Expect FluidBlock index on the lua_stack.
    auto blkId = lua_tointeger(L, 1);
    auto blk = gasBlocks[blkId];
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    lua_pushinteger(L, GlobalConfig.dimensions); lua_setfield(L, tblIdx, "dimensions");
    lua_pushstring(L, blk.label.toStringz); lua_setfield(L, tblIdx, "label");
    lua_pushstring(L, gridTypeName(blk.grid_type).toStringz); lua_setfield(L, tblIdx, "grid_type");
    if (blk.grid_type == Grid_t.structured_grid) {
	SBlock sblk = cast(SBlock) blk;
	assert(sblk !is null, "Oops, this should be an SBlock object.");
	// For a structured_grid
	lua_pushinteger(L, sblk.nicell); lua_setfield(L, tblIdx, "nicell");
	lua_pushinteger(L, sblk.njcell); lua_setfield(L, tblIdx, "njcell");
	lua_pushinteger(L, sblk.nkcell); lua_setfield(L, tblIdx, "nkcell");
	lua_pushinteger(L, sblk.imin); lua_setfield(L, tblIdx, "imin");
	lua_pushinteger(L, sblk.jmin); lua_setfield(L, tblIdx, "jmin");
	lua_pushinteger(L, sblk.kmin); lua_setfield(L, tblIdx, "kmin");
	lua_pushinteger(L, sblk.imax); lua_setfield(L, tblIdx, "imax");
	lua_pushinteger(L, sblk.jmax); lua_setfield(L, tblIdx, "jmax");
	lua_pushinteger(L, sblk.kmax); lua_setfield(L, tblIdx, "kmax");
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
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);

    // Grab the appropriate cell
    auto cell = gasBlocks[blkId].get_cell(i, j, k);
    
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushFluidCellToTable(L, tblIdx, cell, 0, gasBlocks[blkId].myConfig.gmodel);
    return 1;
} // end luafn_sampleFluidCell()

extern(C) int luafn_sampleFluidFace(lua_State *L)
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
    pushFluidFaceToTable(L, tblIdx, face, 0, gasBlocks[blkId].myConfig.gmodel);
    return 1;
} // end luafn_sampleFluidFace()

// -----------------------------------------------------
// D code functions

/**
 * Push the interesting data from a FVCell and FVInterface to a Lua table
 *
 */
void pushFluidCellToTable(lua_State* L, int tblIdx, ref const(FVCell) cell, 
			  size_t gtl, GasModel gmodel)
{
    lua_pushnumber(L, cell.pos[gtl].x); lua_setfield(L, tblIdx, "x");
    lua_pushnumber(L, cell.pos[gtl].y); lua_setfield(L, tblIdx, "y");
    lua_pushnumber(L, cell.pos[gtl].z); lua_setfield(L, tblIdx, "z");
    lua_pushnumber(L, cell.volume[gtl]); lua_setfield(L, tblIdx, "vol");
    pushFlowStateToTable(L, tblIdx, cell.fs, gmodel);
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
extern(C) int luafn_sampleSolidCell(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);

    // Grab the appropriate cell
    auto cell = solidBlocks[blkId].getCell(i, j, k);
    
    // Return the interesting bits as a table.
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushSolidCellToTable(L, tblIdx, cell);
    return 1;
} // end luafn_sampleFluidCell()

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



