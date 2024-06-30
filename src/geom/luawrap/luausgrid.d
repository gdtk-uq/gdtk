/**
 * A Lua interface for the D usgrid (UnstructuredGrid) module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2015-11-06 adapted from luasgrid.d
 */

module geom.luawrap.luausgrid;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import geom.misc.univariatefunctions;
import geom;
import geom.luawrap.luagrid;
import geom.luawrap.luasgrid;
import geom.luawrap.luageom;
// Might put the following back later when we generate our own
// UnstructuredGrid objects from parametric surfaces and volumes.
// import luasurface;
// import luavolume;
// import luaunifunction;

// Name of metatables
immutable string UnstructuredGridMT = "UnstructuredGrid";

static const(UnstructuredGrid)[] unstructuredGridStore;

UnstructuredGrid checkUnstructuredGrid(lua_State* L, int index) {
    if ( isObjType(L, index, UnstructuredGridMT) ) {
        return checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, index);
    }
    // all else fails
    return null;
}

extern(C) int copyUnstructuredGrid(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a grid.
    auto grid = checkObj!(T, MTname)(L, 1);
    unstructuredGridStore ~= pushObj!(T, MTname)(L, grid);
    return 1;
}

extern(C) int get_vtx(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < niv
    Vector3 vtx = grid.vertices[i];
    return pushVector3(L, vtx);
}

extern(C) int get_vtx_id_list_for_cell(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t iCell = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < ncells
    auto vtxList = grid.get_vtx_id_list_for_cell(iCell, 0);
    lua_newtable(L);
    foreach (i, id; vtxList) {
        lua_pushinteger(L, id);
        lua_rawseti(L, -2, cast(int)(i+1));
    }
    return 1;
}

extern(C) int get_connected_cell_ids_for_cell(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t iCell = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < ncells
    auto cellList = grid.get_connected_cell_ids_for_cell(iCell);
    lua_newtable(L);
    foreach (i, id; cellList) {
        lua_pushinteger(L, id);
        lua_rawseti(L, -2, cast(int)(i+1));
    }
    return 1;
}

extern(C) int get_nfaces(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.nfaces);
    return 1;
}

extern(C) int get_nboundaries(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.nboundaries);
    return 1;
}

extern(C) int get_boundaryset_tag(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < nboundaries
    if (i < grid.boundaries.length) {
        lua_pushstring(L, grid.boundaries[i].tag.toStringz());
    } else {
        lua_pushnil(L);
    }
    return 1;
}

extern(C) int set_boundaryset_tag(lua_State* L)
{
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < nboundaries
    string tag = to!string(luaL_checkstring(L, 3));
    if (i < grid.boundaries.length) {
        grid.boundaries[i].tag = tag;
    }
    return 0;
}

extern(C) int is_boundaryset_empty(lua_State* L)
{
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < nboundaries
    if (i < grid.boundaries.length) {
        if (grid.boundaries[i].face_id_list.length == 0)
            lua_pushboolean(L, 1); // 1 is true, yes boundary is empty
        else
            lua_pushboolean(L, 0); // zero is false, boundar is NOT empty
    }
    else {
        // If there's no boundary, then it must be empty
        lua_pushboolean(L, 1);
    }
    return 1;
}

extern(C) int add_boundaryset_faces_to_table(lua_State* L)
{
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < nboundaries
    int idx = to!int(lua_objlen(L, 3)); // Assuming table at location 3
    foreach (faceId; grid.boundaries[i].face_id_list) {
        lua_pushnumber(L, faceId);
        lua_rawseti(L, 3, idx+1);
        idx++;
    }
    return 0;
}

extern(C) int find_nearest_cell_centre_usg(lua_State *L)
{
    double x, y, z;
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    lua_getfield(L, 2, "x");
    if ( !lua_isnil(L, -1) ) {
        x = luaL_checknumber(L, -1);
    } else {
        x = 0.0;
    }
    lua_pop(L, 1);
    lua_getfield(L, 2, "y");
    if ( !lua_isnil(L, -1) ) {
        y = luaL_checknumber(L, -1);
    } else {
        y = 0.0;
    }
    lua_pop(L, 1);
    lua_getfield(L, 2, "z");
    if ( !lua_isnil(L, -1) ) {
        z = luaL_checknumber(L, -1);
    } else {
        z = 0.0;
    }
    lua_pop(L, 1);
    size_t indx;
    double dist;
    Vector3 p = Vector3(x, y, z);
    grid.find_nearest_cell_centre(p, indx, dist);
    lua_pushinteger(L, indx);
    lua_pushnumber(L, dist);
    return 2;
} // end find_nearest_cell_centre_usg()

extern(C) int splitTriangleCell(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t id = to!size_t(luaL_checkint(L, 2));
    grid.splitTriangleCell(id);
    return 0;
}

/**
 * The Lua constructor for a UnstructuredGrid.
 *
 * Example construction in Lua:
 * grid = UnstructuredGrid:new{sgrid=sgrid_object, new_label="a-new-label"}
 * grid = UnstructuredGrid:new{filename="a-file-name", new_label="a-new-label",
 *                             fmt="gziptext", scale=1.0, expect_gmsh_order_for_wedges=true}
 */
extern(C) int newUnstructuredGrid(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected UnstructuredGrid:new{}; ";
        errMsg ~= "maybe you tried UnstructuredGrid.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to UnstructuredGrid:new{}.; ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["sgrid","new_label","filename","fileName","scale","fmt",
                                  "expect_gmsh_order_for_wedges", "label"])) {
        string errMsg = "Error in call to UnstructuredGrid:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    UnstructuredGrid usgrid;
    string label = "";
    lua_getfield(L, 1, "label".toStringz);
    if (lua_isstring(L, -1)) {
        label = to!string(luaL_checkstring(L, -1));
    }
    lua_pop(L, 1); // dispose of label entry, even if it is a nil
    //
    // First, look for a StructuredGrid field.
    lua_getfield(L, 1, "sgrid".toStringz);
    if ( !lua_isnil(L, -1) ) {
        StructuredGrid sgrid = checkStructuredGrid(L, -1);
        if (!sgrid) {
            string errMsg = "Error in UnstructuredGrid:new{}. sgrid not a StructuredGrid.";
            luaL_error(L, errMsg.toStringz);
        }
        usgrid = new UnstructuredGrid(sgrid, label);
    }
    lua_pop(L, 1);
    //
    if (!usgrid) {
        // We didn't find a sgrid entry, so try for a file name.
        lua_getfield(L, 1, "filename".toStringz);
        if (lua_isnil(L, -1)) {
            lua_pop(L, 1);
            lua_getfield(L, 1, "fileName".toStringz);
        }
        if ( !lua_isnil(L, -1) ) {
            string filename = to!string(luaL_checkstring(L, -1));
            double scale = 1.0;
            lua_getfield(L, 1, "scale".toStringz);
            if ( !lua_isnil(L, -1) ) { scale = to!double(luaL_checknumber(L, -1)); }
            lua_pop(L, 1); // dispose of scale item
            string fmt = "gziptext";
            lua_getfield(L, 1, "fmt".toStringz);
            if ( !lua_isnil(L, -1) ) { fmt = to!string(luaL_checkstring(L, -1)); }
            lua_pop(L, 1); // dispose of fmt item
            bool true_centroids = false;
            lua_getfield(L, 1, "true_centroids".toStringz);
            if ( !lua_isnil(L, -1) ) {
                true_centroids = to!bool(lua_toboolean(L, -1));
            }
            lua_pop(L, 1); // dispose of true_centroids item
            bool expect_gmsh_order_for_wedges = false;
            lua_getfield(L, 1, "expect_gmsh_order_for_wedges".toStringz);
            if ( !lua_isnil(L, -1) ) {
                expect_gmsh_order_for_wedges = to!bool(lua_toboolean(L, -1));
            }
            lua_pop(L, 1); // dispose of expect_gmsh_order_for_wedges item
            if (filename.length > 0) {
                usgrid = new UnstructuredGrid(filename, fmt, true_centroids, scale,
                                              expect_gmsh_order_for_wedges,
                                              label);
            }
        } else {
            string errMsg = "Error in UnstructuredGrid:new{}. expected a string for filename.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1); // dispose of filename item
    }
    //
    if (usgrid) {
        unstructuredGridStore ~= pushObj!(UnstructuredGrid, UnstructuredGridMT)(L, usgrid);
        return 1;
    } else {
        string errMsg = "Error in UnstructuredGrid:new{}. Failed to construct a grid.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
} // end newUnstructuredGrid()


extern(C) int usg_joinGrid(lua_State* L)
{
    int narg = lua_gettop(L);
    UnstructuredGrid master = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    if (!master) {
        string errMsg = "Error in UnstructuredGrid:joinGrid(). master grid not an UnstructuredGrid.";
        luaL_error(L, errMsg.toStringz);
    }
    UnstructuredGrid other = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 2);
    if (!other) {
        string errMsg = "Error in UnstructuredGrid:joinGrid(). other grid not an UnstructuredGrid.";
        luaL_error(L, errMsg.toStringz);
    }
    double relTol = 1.0e-6;
    double absTol = 1.0e-9;
    bool checkBoundariesOnly = true;
    bool openFoam = false;
    int dimensions = 3;
    if (narg > 2) { relTol = to!double(luaL_checknumber(L, 3)); }
    if (narg > 3) { absTol = to!double(luaL_checknumber(L, 4)); }
    if (narg > 4) { checkBoundariesOnly = to!bool(lua_toboolean(L, 5)); }
    if (narg > 5) { openFoam = to!bool(lua_toboolean(L, 6)); }
    if (narg > 6) { dimensions = luaL_checkint(L, 7); }
    //
    master.joinGrid(other, relTol, absTol, checkBoundariesOnly, openFoam, dimensions);
    lua_settop(L, 1); // We leave the master grid at stack position 1
    return 1;
} // end usg_joinGrid()

extern(C) int usg_writeStats(lua_State* L)
{
    int narg = lua_gettop(L);
    UnstructuredGrid usgrid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    if (!usgrid) {
        string errMsg = "Error in UnstructuredGrid:writeStats(). grid not an UnstructuredGrid.";
        luaL_error(L, errMsg.toStringz);
    }
    usgrid.writeStats();
    lua_settop(L, 0);
    return 0;
} // end usg_writeStats()

extern(C) int usg_writeOpenFoamPolyMesh(lua_State* L)
{
    int narg = lua_gettop(L);
    UnstructuredGrid usgrid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    if (!usgrid) {
        string errMsg = "Error in UnstructuredGrid:writeOpenFoamPolyMesh(). grid not an UnstructuredGrid.";
        luaL_error(L, errMsg.toStringz);
    }
    string topLevelDir = ".";
    if (lua_isstring(L, 2)) { topLevelDir = to!string(luaL_checkstring(L, 2)); }
    usgrid.write_openFoam_polyMesh(topLevelDir);
    lua_settop(L, 0);
    return 0;
} // end usg_writeOpenFoamPolyMesh()


void registerUnstructuredGrid(lua_State* L)
{
    // Register the UnstructuredGrid object
    luaL_newmetatable(L, UnstructuredGridMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newUnstructuredGrid);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnstructuredGrid!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &get_dimensions!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_dimensions");
    lua_pushcfunction(L, &get_type!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_type");
    lua_pushcfunction(L, &get_nvertices!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_nvertices");
    lua_pushcfunction(L, &get_ncells!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_ncells");
    lua_pushcfunction(L, &get_nfaces);
    lua_setfield(L, -2, "get_nfaces");
    lua_pushcfunction(L, &get_nboundaries);
    lua_setfield(L, -2, "get_nboundaries");
    lua_pushcfunction(L, &get_boundaryset_tag);
    lua_setfield(L, -2, "get_boundaryset_tag");
    lua_pushcfunction(L, &set_boundaryset_tag);
    lua_setfield(L, -2, "set_boundaryset_tag");
    lua_pushcfunction(L, &is_boundaryset_empty);
    lua_setfield(L, -2, "is_boundaryset_empty");
    lua_pushcfunction(L, &add_boundaryset_faces_to_table);
    lua_setfield(L, -2, "add_boundaryset_faces_to_table");
    lua_pushcfunction(L, &get_vtx!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &get_vtx_id_list_for_cell);
    lua_setfield(L, -2, "get_vtx_id_list_for_cell");
    lua_pushcfunction(L, &get_connected_cell_ids_for_cell);
    lua_setfield(L, -2, "get_connected_cell_ids_for_cell");
    lua_pushcfunction(L, &cellVolume!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "cellVolume");
    lua_pushcfunction(L, &cellCentroid!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "cellCentroid");
    lua_pushcfunction(L, &write_to_gzip_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_gzip_file");
    lua_pushcfunction(L, &write_to_raw_binary_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_raw_binary_file");
    lua_pushcfunction(L, &write_to_vtk_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_vtk_file");
    lua_pushcfunction(L, &write_to_stl_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_stl_file");
    lua_pushcfunction(L, &write_to_su2_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_su2_file");
    lua_pushcfunction(L, &rotate!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "rotate");
    lua_pushcfunction(L, &find_nearest_cell_centre_usg);
    lua_setfield(L, -2, "find_nearest_cell_centre");
    lua_pushcfunction(L, &splitTriangleCell);
    lua_setfield(L, -2, "splitTriangleCell");
    lua_pushcfunction(L, &usg_joinGrid);
    lua_setfield(L, -2, "joinGrid");
    lua_pushcfunction(L, &usg_writeStats);
    lua_setfield(L, -2, "writeStats");
    lua_pushcfunction(L, &usg_writeOpenFoamPolyMesh);
    lua_setfield(L, -2, "writeOpenFoamPolyMesh");

    lua_setglobal(L, UnstructuredGridMT.toStringz);

} // end registerUnstructuredGrid()

