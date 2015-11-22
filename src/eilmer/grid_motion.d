// Authors: RG, PJ & KD
// Date: 2015-11-20

import std.string;

import util.lua;
import util.lua_service;
import lua_helper;
import globalconfig;
import globaldata;
import luageom;
import sgrid;
import block;

void init_master_lua_State(string fname)
{
    GlobalConfig.master_lua_State = init_lua_State(fname);
    // Give me a conveniently-named pointer for use in this function.
    auto L = GlobalConfig.master_lua_State;
    // Set some globally available constants for the
    // Lua state.
    lua_pushnumber(L, GlobalConfig.nBlocks);
    lua_setglobal(L, "nBlocks");
    lua_pushnumber(L, nghost);
    lua_setglobal(L, "nGhost");
    // Give the user a table that holds information about
    // all of the blocks
    lua_newtable(L);
    foreach ( int i, blk; gasBlocks ) {
	lua_newtable(L);
	lua_pushnumber(L, blk.cells.length);
	lua_setfield(L, -2, "nCells");
	lua_pushnumber(L, blk.vertices.length);
	lua_setfield(L, -2, "nVertices");
	if ( blk.grid_type == Grid_t.structured_grid ) {
	    lua_pushnumber(L, blk.nicell);
	    lua_setfield(L, -2, "niCells");
	    lua_pushnumber(L, blk.njcell);
	    lua_setfield(L, -2, "njCells");
	    lua_pushnumber(L, blk.nkcell);
	    lua_setfield(L, -2, "nkCells");
	    lua_pushnumber(L, blk.imin);
	    lua_setfield(L, -2, "vtxImin");
	    lua_pushnumber(L, blk.imax+1);
	    lua_setfield(L, -2, "vtxImax");
	    lua_pushnumber(L, blk.jmin);
	    lua_setfield(L, -2, "vtxJmin");
	    lua_pushnumber(L, blk.jmax+1);
	    lua_setfield(L, -2, "vtxJmax");
	    lua_pushnumber(L, blk.kmin);
	    lua_setfield(L, -2, "vtxKmin");
	    lua_pushnumber(L, blk.kmax+1);
	    lua_setfield(L, -2, "vtxKmax");
	}
	lua_rawseti(L, -2, i);
    }
    lua_setglobal(L, "blockData");

    // Now set some helper functions
    lua_pushcfunction(L, &luafn_sampleFlow);
    lua_setglobal(L, "sampleFlow");
    lua_pushcfunction(L, &luafn_setVtxVelocityForDomain);
    lua_setglobal(L, "setVtxVelocityForDomain");
    lua_pushcfunction(L, &luafn_setVtxVelocityForBlock);
    lua_setglobal(L, "setVtxVelocityForBlock");
    lua_pushcfunction(L, &luafn_setVtxVelocity);
    lua_setglobal(L, "setVtxVelocity");

}

extern(C) int luafn_setVtxVelocityForDomain(lua_State* L)
{
    // Expect a single argument: a Vector3 object
    auto vel = checkVector3(L, 1);

    foreach ( blk; gasBlocks ) {
	foreach ( vtx; blk.vertices ) {
	    /* We assume that we'll only update grid positions
	       at the start of the increment. This should work
	       well except in the most critical cases of time
	       accuracy.
	    */
	    vtx.vel[0] = *vel;
	}
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

extern(C) int luafn_setVtxVelocityForBlock(lua_State* L)
{
    // Expect two arguments: 1. a Vector3 object
    //                       2. a block id
    auto vel = checkVector3(L, 1);
    auto blkId = lua_tointeger(L, 2);

    foreach ( vtx; gasBlocks[blkId].vertices ) {
	/* We assume that we'll only update grid positions
	   at the start of the increment. This should work
	   well except in the most critical cases of time
	   accuracy.
	*/
	vtx.vel[0] = *vel;
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

/**
 * Sets the velocity of an individual vertex.
 *
 * This function can be called for structured 
 * or unstructured grids. We'll determine what
 * type grid is meant by the number of arguments
 * supplied. The following calls are allowed:
 *
 * setVtxVelocity(vel, blkId, vtxId)
 *   Sets the velocity vector for vertex vtxId in
 *   block blkId. This works for both structured
 *   and unstructured grids.
 *
 * setVtxVelocity(vel, blkId, i, j)
 *   Sets the velocity vector for vertex "i,j" in
 *   block blkId in a two-dimensional structured grid.
 *
 * setVtxVelocity(vel, blkId, i, j, k)
 *   Set the velocity vector for vertex "i,j,k" in
 *   block blkId in a three-dimensional structured grid.
 */
extern(C) int luafn_setVtxVelocity(lua_State* L)
{
    int narg = lua_gettop(L);
    auto vel = checkVector3(L, 1);
    auto blkId = lua_tointeger(L, 2);

    if ( narg == 3 ) {
	auto vtxId = lua_tointeger(L, 3);
	gasBlocks[blkId].vertices[vtxId].vel[0] = *vel;
    }
    else if ( narg == 4 ) {
	auto i = lua_tointeger(L, 3);
	auto j = lua_tointeger(L, 4);
	gasBlocks[blkId].get_vtx(i,j).vel[0] = *vel;
    }
    else if ( narg >= 5 ) {
	auto i = lua_tointeger(L, 3);
	auto j = lua_tointeger(L, 4);
	auto k = lua_tointeger(L, 5);
	gasBlocks[blkId].get_vtx(i,j,k).vel[0] = *vel;
    }
    else {
	string errMsg = "ERROR: Too few arguments passed to luafn: setVtxVelocity()\n";
	luaL_error(L, errMsg.toStringz);
    }
    lua_settop(L, 0);
    return 0;
}

void assign_vertex_velocities_via_udf(double sim_time)
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "assignVtxVelocities");
    lua_pushnumber(L, sim_time);
    int number_args = 1;
    int number_results = 0;

    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	luaL_error(L, "error running user-defined function assignVtxVelocities()");
    }
}
