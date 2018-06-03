// Authors: RG, PJ, KD & IJ
// Date: 2015-11-20

module grid_motion;

import std.string;
import std.conv;

import util.lua;
import util.lua_service;
import nm.luabbla;
import lua_helper;
import fvcore;
import globalconfig;
import globaldata;
import geom;
import geom.luawrap;
import fluidblock;
import sfluidblock;
import std.stdio;

void setGridMotionHelperFunctions(lua_State *L)
{
    lua_pushcfunction(L, &luafn_getVtxPosition);
    lua_setglobal(L, "getVtxPosition");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForDomain);
    lua_setglobal(L, "setVtxVelocitiesForDomain");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForBlock);
    lua_setglobal(L, "setVtxVelocitiesForBlock");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForRotatingBlock);
    lua_setglobal(L, "setVtxVelocitiesForRotatingBlock");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesByCorners);
    lua_setglobal(L, "setVtxVelocitiesByCorners");
    lua_pushcfunction(L, &luafn_setVtxVelocity);
    lua_setglobal(L, "setVtxVelocity");
}

extern(C) int luafn_getVtxPosition(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);

    // Grab the appropriate vtx
    auto vtx = globalFluidBlocks[blkId].get_vtx(i, j, k);
    
    // Return the interesting bits as a table with entries x, y, z.
    lua_newtable(L);
    lua_pushnumber(L, vtx.pos[0].x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, vtx.pos[0].y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, vtx.pos[0].z); lua_setfield(L, -2, "z");
    return 1;
}

extern(C) int luafn_setVtxVelocitiesForDomain(lua_State* L)
{
    // Expect a single argument: a Vector3 object
    auto vel = checkVector3(L, 1);

    foreach ( blk; localFluidBlocks ) {
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

extern(C) int luafn_setVtxVelocitiesForBlock(lua_State* L)
{
    // Expect two arguments: 1. a Vector3 object
    //                       2. a block id
    auto vel = checkVector3(L, 1);
    auto blkId = lua_tointeger(L, 2);

    foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
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
 * Sets the velocity of an entire block, for the case
 * that the block is rotating about an axis with direction 
 * (0 0 1) located at a point defined by vector (x y 0).
 *
 * setVtxVelocitiesRotatingBlock(blkId, omega, vector3)
 *      Sets rotational speed omega (rad/s) for rotation about 
 *      (0 0 1) axis defined by Vector3.
 *
 * setVtxVelocitiesRotatingBlock(blkId, omega)
 *      Sets rotational speed omega (rad/s) for rotation about 
 *      Z-axis.
 */
extern(C) int luafn_setVtxVelocitiesForRotatingBlock(lua_State* L)
{
    // Expect two/three arguments: 1. a block id
    //                             2. a float object
    //                             3. a vector (optional) 
    int narg = lua_gettop(L);
    auto blkId = lua_tointeger(L, 1);
    double omega = lua_tonumber(L, 2);
    double velx, vely;

    if ( narg == 2 ) {
        // assume rotation about Z-axis 
        foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
            velx = - omega * vtx.pos[0].y;
            vely =   omega * vtx.pos[0].x;
            vtx.vel[0] = Vector3(velx, vely, 0.);
        }
    }
    else if ( narg == 3 ) {
        auto axis = checkVector3(L, 3);
        foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
            velx = - omega * (vtx.pos[0].y - axis.y);
            vely =   omega * (vtx.pos[0].x - axis.x);
            vtx.vel[0] = Vector3(velx, vely, 0.);
        }
    }
    else {
        string errMsg = "ERROR: Too few arguments passed to luafn: setVtxVelocitiesRotatingBlock()\n";
        luaL_error(L, errMsg.toStringz);
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

/**
 * Sets the velocity of vertices in a block based on
 * specified corner velocities. The velocity of any cell 
 * is estimated using interpolation based on cell indices. 
 * Note, for clustered blocks or blocks with curved edges 
 * this can lead to distortion.
 *
 * This function can be called for structured 
 * grids only. The following calls are allowed:
 *
 * setVtxVelocitiesByCorners(blkId, p00vel, p10vel, p01vel, p11vel)
 *   Sets the velocity vectors for block with corner velocities 
 *   specified by four corner velocities. This works for both 
 *   2-D and 3- meshes. In 3-D a uniform velocity is applied 
 *   in k direction.
 * 
 * setVtxVelocitiesByCorners(blkId, p000vel, p100vel, p010vel, p110vel, 
 *                       p001vel, p101vel, p011vel, p111vel)
 * (eqivalent to)
 * setVtxVelocitiesByCorners(blkId, p0vel, p1vel, p3vel, p2vel, 
 *                       p4vel, p5vel, p7vel, p6vel)
 *   As above but suitable for 3-D meshes with eight specified 
 *   velocities.
 */
extern(C) int luafn_setVtxVelocitiesByCorners(lua_State* L)
{
    // Expect five/nine arguments: 1. a block id
    //                             2-5. corner velocities for 2-D motion
    //                             6-9. corner velocities for full 3-D motion 
    //                                  (optional)
    int narg = lua_gettop(L);
    double u, v;
    auto blkId = lua_tointeger(L, 1);
    size_t i, j, k;
    Vector3 velw, vele, veln, vels, vel;
    auto blk = cast(SFluidBlock) globalFluidBlocks[blkId];
    // get corner velocities
    auto p00vel = checkVector3(L, 2);
    auto p10vel = checkVector3(L, 3);
    auto p01vel = checkVector3(L, 4);
    auto p11vel = checkVector3(L, 5);  

    if ( narg == 5 ) {
        if (blk.myConfig.dimensions == 2) {
            // deal with 2-D meshes
            k = blk.kmin;
            for (j = blk.jmin; j <= blk.jmax+1; ++j) {
                // find velocity along west and east edge
                v = to!double(j-blk.jmin) / to!double(blk.jmax+1-blk.jmin); 
                velw = v * *p01vel + (1 - v) * *p00vel;
                vele = v * *p11vel + (1 - v) * *p10vel;

                for (i = blk.imin; i <= blk.imax+1; ++i) {
                    //// interpolate in i direction
                    u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin); 
                    vel = u * vele + (1-u) * velw;
                    //// set velocity
                    blk.get_vtx(i,j,k).vel[0] = vel;

                    // transfinite interpolation is yielding same result, but omitted as more expensive.
                    //u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin);
                    //vels = u * *p10vel + (1 - u) * *p00vel;
                    //veln = u * *p11vel + (1 - u) * *p01vel;
                    //// do transfinite interpolation
                    //vel = (1-v)*vels + v*veln + (1-u)*velw + u*vele
                    //    - ( (1-u)*(1-v)* *p00vel + u*v* *p11vel + u*(1-v)* *p10vel + (1-u)*v* *p01vel );
                    //// set velocity
                    //blk.get_vtx(i,j,k).vel[0] = vel;
                }
            }
        }
        else { // deal with 3-D meshesv (assume constant properties wrt k index)
            for (j = blk.jmin; j <= blk.jmax+1; ++j) {
                // find velocity along west and east edge
                v = to!double(j-blk.jmin) / to!double(blk.jmax+1-blk.jmin); 
                velw = v * *p01vel + (1 - v) * *p00vel;
                vele = v * *p11vel + (1 - v) * *p10vel;


                for (i = blk.imin; i <= blk.imax+1; ++i) {
                    // interpolate in i direction
                    u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin); 
                    vel = u * vele + (1-u) * velw;

                    // set velocity for all k
                    for (k = blk.kmin; k <= blk.kmax+1; ++i) {
                        blk.get_vtx(i,j,k).vel[0] = vel;
                    }
                }
            }
        }
    }
    else if ( narg == 9 ) {
        // assume all blocks are 3-D
        writeln("setVtxVelocitiesByCorners not verified for 3-D. 
                Proceed at own peril. See /src/eilmer/grid_motion.d");
        double w;
        Vector3 velwt, velet, velt;
        // get corner velocities
        auto p001vel = checkVector3(L, 6);
        auto p101vel = checkVector3(L, 7);
        auto p011vel = checkVector3(L, 8);
        auto p111vel = checkVector3(L, 9);  
        for (j = blk.jmin; j <= blk.jmax+1; ++j) {
            // find velocity along west and east edge (four times)
            v = to!double(j-blk.jmin) / to!double(blk.jmax+1-blk.jmin); 
            velw = v * *p01vel + (1 - v) * *p00vel;
            vele = v * *p11vel + (1 - v) * *p10vel;
            velwt = v * *p011vel + (1 - v) * *p001vel;
            velet = v * *p111vel + (1 - v) * *p101vel;

            for (i = blk.imin; i <= blk.imax+1; ++i) {
                // interpolate in i direction (twice)
                u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin); 
                vel = u * vele + (1-u) * velw;
                velt = u * velet + (1-u) * velwt;

                // set velocity by interpolating in k.
                for (k = blk.kmin; k <= blk.kmax+1; ++i) {
                    w = to!double(k-blk.kmin) / to!double(blk.kmax+1-blk.kmin); 
                    blk.get_vtx(i,j,k).vel[0] = w * velt + (1 - w) * vel;
                }
            }
        }
    }
    else {
        string errMsg = "ERROR: Wrong number of arguments passed to luafn: setVtxVelocitiesByCorners()\n";
        luaL_error(L, errMsg.toStringz);
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
        globalFluidBlocks[blkId].vertices[vtxId].vel[0] = *vel;
    }
    else if ( narg == 4 ) {
        auto i = lua_tointeger(L, 3);
        auto j = lua_tointeger(L, 4);
        globalFluidBlocks[blkId].get_vtx(i,j).vel[0] = *vel;
    }
    else if ( narg >= 5 ) {
        auto i = lua_tointeger(L, 3);
        auto j = lua_tointeger(L, 4);
        auto k = lua_tointeger(L, 5);
        globalFluidBlocks[blkId].get_vtx(i,j,k).vel[0] = *vel;
    }
    else {
        string errMsg = "ERROR: Too few arguments passed to luafn: setVtxVelocity()\n";
        luaL_error(L, errMsg.toStringz);
    }
    lua_settop(L, 0);
    return 0;
}

void assign_vertex_velocities_via_udf(double sim_time, double dt)
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "assignVtxVelocities");
    lua_pushnumber(L, sim_time);
    lua_pushnumber(L, dt);
    int number_args = 2;
    int number_results = 0;

    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
        string errMsg = "ERROR: while running user-defined function assignVtxVelocities()\n";
        errMsg ~= to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
}
