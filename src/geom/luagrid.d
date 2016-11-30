/**
 * A Lua interface for the D grid module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2016-11-30
 *
 * We need a base class to collect the common wrapper
 * implementations between structured and unstructured grids.
 */

module luagrid;

import std.conv;

import util.lua;
import util.lua_service;
import grid;
import sgrid;
import usgrid;


extern(C) int get_dimensions(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.dimensions);
    return 1;
}

extern(C) int get_nvertices(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.nvertices);
    return 1;
}

extern(C) int get_ncells(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.ncells);
    return 1;
}

extern(C) int cellVolume(T, string MTname)(lua_State *L)
{
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t indx = to!size_t(luaL_checkint(L, 2));
    double vol = grid.cell_volume(indx);
    lua_pushnumber(L, vol);
    return 1;
}

extern(C) int write_to_vtk_file(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(T, MTname)(L, 1);
    auto fileName = to!string(luaL_checkstring(L, 2));
    grid.write_to_vtk_file(fileName);
    return 0;
}

extern(C) int write_to_gzip_file(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(T, MTname)(L, 1);
    auto fileName = to!string(luaL_checkstring(L, 2));
    grid.write_to_gzip_file(fileName);
    return 0;
}
