/**
 * A Lua interface for the D sgrid (StructuredGrid) module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-26
 */

module geom.luawrap.luasgrid;

import std.stdio;
import std.string;
import std.conv;
import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import geom.misc.univariatefunctions;
import geom;
import geom.luawrap;

// Name of metatables
immutable string StructuredGridMT = "StructuredGrid";

static const(StructuredGrid)[] structuredGridStore;

StructuredGrid checkStructuredGrid(lua_State* L, int index) {
    if ( isObjType(L, index, StructuredGridMT) ) {
        return checkObj!(StructuredGrid, StructuredGridMT)(L, index);
    }
    return null; // on fail
}

extern(C) int copyStructuredGrid(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a grid.
    auto grid = checkObj!(T, MTname)(L, 1);
    structuredGridStore ~= pushObj!(T, MTname)(L, grid.dup());
    return 1;
}

extern(C) int get_niv(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushinteger(L, grid.niv);
    return 1;
}

extern(C) int get_njv(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushinteger(L, grid.njv);
    return 1;
}

extern(C) int get_nkv(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushinteger(L, grid.nkv);
    return 1;
}

extern(C) int get_vtx(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2));
    // Note that we expect 0 <= i < niv
    size_t j = 0; if (narg > 2) { j = to!size_t(luaL_checkint(L, 3)); }
    size_t k = 0; if (narg > 3) { k = to!size_t(luaL_checkint(L, 4)); }
    Vector3* vtx = grid[i,j,k];
    return pushVector3(L, *vtx);
}

extern(C) int set_vtx(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2));
    // Note that we expect 0 <= i < niv
    size_t j = 0; if (narg > 3) { j = to!size_t(luaL_checkint(L, 3)); }
    size_t k = 0; if (narg > 4) { k = to!size_t(luaL_checkint(L, 4)); }
    Vector3 vtx = toVector3(L, narg);
    Vector3* grid_vtx = grid[i,j,k];
    grid_vtx.set(vtx);
    return 0;
}

extern(C) int get_corner_vtx(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    string corner_tag = to!string(lua_tostring(L, 2));
    Vector3* vtx;
    if (grid.dimensions == 2) {
        switch (corner_tag) {
        case "00": vtx = grid[0, 0]; break;
        case "10": vtx = grid[grid.niv-1, 0]; break;
        case "11": vtx = grid[grid.niv-1, grid.njv-1]; break;
        case "01": vtx = grid[0, grid.njv-1]; break;
        default:
            string errMsg = format("Error in StructuredGrid:get_corner_vtx(). \"%s\" is not a valid corner name.",
                                   corner_tag);
            luaL_error(L, errMsg.toStringz);
        }
    } else if (grid.dimensions == 3) {
        switch (corner_tag) {
        case "000": vtx = grid[0, 0, 0]; break;
        case "100": vtx = grid[grid.niv-1, 0, 0]; break;
        case "110": vtx = grid[grid.niv-1, grid.njv-1, 0]; break;
        case "010": vtx = grid[0, grid.njv-1, 0]; break;
        case "001": vtx = grid[0, 0, grid.nkv-1]; break;
        case "101": vtx = grid[grid.niv-1, 0, grid.nkv-1]; break;
        case "111": vtx = grid[grid.niv-1, grid.njv-1, grid.nkv-1]; break;
        case "011": vtx = grid[0, grid.njv-1, grid.nkv-1]; break;
        default:
            string errMsg = format("Error in StructuredGrid:get_corner_vtx(). \"%s\" is not a valid corner name.",
                                   corner_tag);
            luaL_error(L, errMsg.toStringz);
        }
    } else {
        // Assume grid dimensions == 1.
        if (grid.dimensions != 1) {
            string errMsg = format("Error in StructuredGrid:get_corner_vtx(). \"%d\" invalid number of dimensions.",
                                   grid.dimensions);
            luaL_error(L, errMsg.toStringz);
        }
        switch (corner_tag) {
        case "0": vtx = grid[0]; break;
        case "1": vtx = grid[grid.niv-1]; break;
        default:
            string errMsg = format("Error in StructuredGrid:get_corner_vtx(). \"%s\" is not a valid corner name.",
                                   corner_tag);
            luaL_error(L, errMsg.toStringz);
        }
    }
    return pushVector3(L, *vtx);
} // end get_corner_vtx()

/*
 * Set one or more string tags for the boundaries.
 * my_sgrid:set_tags{north="my-north-tag", south="special-axis"}
 */
extern(C) int set_boundary_tags(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    if (lua_istable(L, 2)) {
        lua_getfield(L, 2, "north");
        if (!lua_isnil(L, -1)) {
            string tag = to!string(lua_tostring(L, -1));
            grid.set_tag(Face.north, tag);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "east");
        if (!lua_isnil(L, -1)) {
            string tag = to!string(lua_tostring(L, -1));
            grid.set_tag(Face.east, tag);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "south");
        if (!lua_isnil(L, -1)) {
            string tag = to!string(lua_tostring(L, -1));
            grid.set_tag(Face.south, tag);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "west");
        if (!lua_isnil(L, -1)) {
            string tag = to!string(lua_tostring(L, -1));
            grid.set_tag(Face.west, tag);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "top");
        if (!lua_isnil(L, -1)) {
            string tag = to!string(lua_tostring(L, -1));
            grid.set_tag(Face.top, tag);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "bottom");
        if (!lua_isnil(L, -1)) {
            string tag = to!string(lua_tostring(L, -1));
            grid.set_tag(Face.bottom, tag);
        }
        lua_pop(L, 1);
    } else {
        string errMsg = "myGrid:set_tags{} was expecting a table of named tag strings.";
        luaL_error(L, errMsg.toStringz);
    }
    return 0;
}

extern(C) int get_boundary_tag(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    string tag = ""; // default value of empty string
    if (lua_isstring(L, 2)) {
        string bName = to!string(lua_tostring(L, 2));
        int bIndx = face_index(bName);
        if (bIndx < grid.bcTags.length) { tag = grid.bcTags[bIndx]; }
    }
    if (lua_isnumber(L, 2)) {
        int bIndx = to!int(lua_tointeger(L, 2));
        if (bIndx < grid.bcTags.length) { tag = grid.bcTags[bIndx]; }
    }
    lua_pushstring(L, tag.toStringz);
    return 1;
}

extern(C) int subgrid(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    // Note that we expect 0 <= i0 < niv
    size_t i0 = to!size_t(luaL_checkint(L, 2));
    size_t ni = to!size_t(luaL_checkint(L, 3));
    size_t j0 = 0; if (narg > 3) { j0 = to!size_t(luaL_checkint(L, 4)); }
    size_t nj = 1; if (narg > 4) { nj = to!size_t(luaL_checkint(L, 5)); }
    size_t k0 = 0; if (narg > 5) { k0 = to!size_t(luaL_checkint(L, 6)); }
    size_t nk = 1; if (narg > 6) { nk = to!size_t(luaL_checkint(L, 7)); }
    auto new_subgrid = grid.subgrid(i0,ni,j0,nj,k0,nk);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, new_subgrid);
    return 1;
} // end subgrid()

extern(C) int get_boundary_grid(lua_State *L)
{
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    size_t boundary_indx = to!size_t(luaL_checkint(L, 2));
    StructuredGrid extracted_grid = cast(StructuredGrid) grid.get_boundary_grid(boundary_indx);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, extracted_grid);
    return 1;
}

extern(C) int joinGrid(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 3;
    auto selfGrid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    auto otherGrid = checkObj!(StructuredGrid, StructuredGridMT)(L, 2);
    auto joinLocation = to!string(luaL_checkstring(L, 3));
    selfGrid.joinGrid(otherGrid, joinLocation);
    return 0;
}

extern(C) int rotateGrid(lua_State* L)
{
    // Flynn Hack 2023-03-23
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_remove(L, 1); // remove grid from stack
    if (narg == 2 && lua_istable(L, 1)) {
        lua_getfield(L, 1, "q0");
        double q0 = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "q1");
        double q1 = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "q2");
        double q2 = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "q3");
        double q3 = lua_tonumber(L, -1);
        lua_pop(L, 1);
        grid.rotateGrid(q0, q1, q2, q3); // rotate our grid
        return 0;
    } else {
        string errMsg = "Expected grid:rotateGrid{q0=number, q1=number,
                                                    q2=number, q3=number} ";
        luaL_error(L, errMsg.toStringz);
        return 1;
    }
} // end rotateGrid

extern(C) int find_nearest_cell_centre_sg(lua_State *L)
{
    double x, y, z;
    int narg = lua_gettop(L); // assume narg == 2;
    if (narg != 2) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected myGrid:find_nearest_cell_centre_sg{}; ";
        errMsg ~= "maybe you tried myGrid:find_nearest_cell_centre_sg{}.";
        luaL_error(L, errMsg.toStringz);
    }
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
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
} // end find_nearest_cell_centre_sg()

/**
 * The Lua constructor for a StructuredGrid.
 *
 * Example construction in Lua:
 * grid1D = StructuredGrid:new{path=somePath, niv=10, cf=my_cf01}
 * grid2D = StructuredGrid:new{psurface=someParametricSurface, niv=10, njv=10,
 *                             cfList={north=cfn, east=cfe, south=cfs, west=cfw},
 *                             interpolation="linear",
 *                             label="A-2D-Grid"}
 * grid3D = StructuredGrid:new{pvolume=someParametricVolume,
 *                             niv=11, njv=21, nkv=11,
 *                             cfList={edge01=cf01, edge32=cf32, edge45=cf45, edge76=cf76},
 *                             interpolation="linear",
 *                             label="A-3D-Grid"}
 *
 * importedGrid = StructuredGrid:new{filename="myFileName.dat", fmt="vtk"}
 * The (default) format "vtk" is the legacy text format for the VTK system.
 * Other format values could be "gziptext", which is the Eilmer4 native file
 * or "text", which is essentially the Eilmer3 native file uncompressed.
 * There is also "rawbinary".
 */
extern(C) int newStructuredGrid(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected StructuredGrid:new{}; ";
        errMsg ~= "maybe you tried StructuredGrid.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to StructuredGrid:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["path","psurface","pvolume","niv","njv","nkv",
                                  "cf","cfList","label","filename","fileName","fmt",
                                  "r_grid","s_grid","interpolation"])) {
        string errMsg = "Error in call to StructuredGrid:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Check if label has been supplied
    string label = "";
    lua_getfield(L, 1, "label".toStringz);
    if (!lua_isnil(L, -1)) {
	label = to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);

    // Let's first see if we have been given a file that contains the grid.
    lua_getfield(L, 1, "filename".toStringz);
    if (lua_isnil(L, -1)) {
        lua_pop(L, 1);
        lua_getfield(L, 1, "fileName".toStringz);
    }
    if ( lua_isstring(L, -1) ) {
        string fileName = to!string(lua_tostring(L, -1));
        lua_pop(L, 1);
        string fmt = "vtk";
        lua_getfield(L, 1, "fmt".toStringz);
        if ( lua_isstring(L, -1) ) { fmt = to!string(lua_tostring(L, -1)); }
        lua_pop(L, 1);
        StructuredGrid importedGrid = new StructuredGrid(fileName, fmt, label);
        structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, importedGrid);
        return 1; // our work is done
    }
    lua_pop(L, 1);
    //
    // If we haven't already obtained a grid from a file,
    // it's time to try to construct one on top of
    // a geometric description of the region.
    Path mypath;
    ParametricSurface psurface;
    ParametricVolume pvolume;
    int dimensions = 0;
    // First, look for a Path field, for a 1D grid.
    lua_getfield(L, 1, "path".toStringz);
    if ( !lua_isnil(L, -1) ) {
        dimensions = 1;
        mypath = checkPath(L, -1);
        if (!mypath) {
            string errMsg = "Error in StructuredGrid:new{}. path not a Path.";
            luaL_error(L, errMsg.toStringz);
        }
    }
    lua_pop(L, 1);
    if (dimensions == 0) {
        // Next, look for a ParametricSurface field, for a 2D grid.
        lua_getfield(L, 1, "psurface".toStringz);
        if ( !lua_isnil(L, -1) ) {
            dimensions = 2;
            psurface = checkSurface(L, -1);
            if (!psurface) {
                string errMsg = "Error in StructuredGrid:new{}. psurface not a ParametricSurface.";
                luaL_error(L, errMsg.toStringz);
            }
        }
        lua_pop(L, 1);
    }
    if (dimensions == 0) {
        // No path or surface, so try for a pvolume, for a 3D grid.
        lua_getfield(L, 1, "pvolume".toStringz);
        if ( !lua_isnil(L, -1) ) {
            dimensions = 3;
            pvolume = checkVolume(L, -1);
            if (!pvolume) {
                string errMsg = "Error in StructuredGrid:new{}. pvolume not a ParametricVolume.";
                luaL_error(L, errMsg.toStringz);
            }
        }
        lua_pop(L, 1);
    }
    if (dimensions == 0) {
        string errMsg = "Error in StructuredGrid:new{}. no path, psurface or pvolume found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    // Get clustering functions, filling in nil or invalid entries with LinearFunction.
    UnivariateFunction cf;
    UnivariateFunction[] cfList;
    if (dimensions == 1) {
        lua_getfield(L, 1, "cf".toStringz);
        cf = checkUnivariateFunction(L, -1);
        if (!cf) { cf = new LinearFunction(); }
        lua_pop(L, 1);
    } else {
        // There will be a list of cluster functions for 2D or 3D.
        int number_of_edges = (dimensions == 2) ? 4 : 12;
        lua_getfield(L, 1, "cfList".toStringz);
        if ( lua_istable(L, -1) ) {
            // Before going on, check the caller gave us named parameters.
            if ( lua_objlen(L, -1) != 0 ) {
                // It appears that the caller has tried to set arguments as an array
                string errMsg = "Error in StructuredGrid:new{}. " ~
                    "A table of named parameters is expected for the cfList. " ~
                    "An array style of parameters was found.";
                luaL_error(L, errMsg.toStringz);
            }

            string[] edges;
            // Extract the cluster functions from the table at top of stack.
            if (dimensions == 2) {
                // The new face order, in 2023.  See ../geom/elements/nomenclature.d.
                edges = ["west", "east", "south", "north"];
            } else {
                // This 3D edge order needs to be kept consistent with
                // the order expected by make_grid_from_volume() over in ../grid/sgrid.d.
                edges = ["edge01", "edge12", "edge32", "edge03",
                         "edge45", "edge56", "edge76", "edge47",
                         "edge04", "edge15", "edge26", "edge37"];
            }
            assert(number_of_edges == edges.length, "Mismatch in edge list.");
            foreach (edge_name; edges) {
                lua_getfield(L, -1, edge_name.toStringz);
                if ( lua_isnil(L, -1) ) {
                    cfList ~= new LinearFunction();
                } else {
                    auto mycf = checkUnivariateFunction(L, -1);
                    if ( mycf ) {
                        cfList ~= mycf;
                    } else {
                        cfList ~= new LinearFunction();
                    }
                }
                lua_pop(L, 1);
            }
        } else {
            // Didn't find a table of cluster functions.
            foreach (i; 0 .. number_of_edges) {
                cfList ~= new LinearFunction();
            }
        }
        lua_pop(L, 1);
    }
    //
    string interpolation = "simple";
    lua_getfield(L, 1, "interpolation".toStringz);
    if ( lua_isstring(L, -1) ) {
        interpolation = to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    //
    string errMsgTmplt = "Error in StructuredGrid:new{}. " ~
        "A valid value for '%s' was not found in list of arguments. " ~
        "The value, if present, should be a number.";
    int niv = getIntegerFromTable(L, 1, "niv", true, 0, true, format(errMsgTmplt, "niv"));
    int njv = 1;
    if (dimensions > 1) {
        njv = getIntegerFromTable(L, 1, "njv", true, 0, true, format(errMsgTmplt, "njv"));
    }
    int nkv = 1;
    if (dimensions == 3) {
        nkv = getIntegerFromTable(L, 1, "nkv", true, 0, true, format(errMsgTmplt, "nkv"));
    }
    //
    // Look for an r_grid and s_grid. If present, call special construction.
    bool useSpecial2DConstructor = false;
    double[4][4] r_grid, s_grid;
    lua_getfield(L, 1, "r_grid");
    if (!lua_isnil(L, -1)) {
	useSpecial2DConstructor = true;
        getMatrixOfDoubles!(4)(L, -1, r_grid);
    }
    lua_pop(L, 1);
    //
    lua_getfield(L, 1, "s_grid");
    if (!lua_isnil(L, -1)) { getMatrixOfDoubles!(4)(L, -1, s_grid); }
    lua_pop(L, 1);
    //
    StructuredGrid grid;
    switch (dimensions) {
    case 1: grid = new StructuredGrid(mypath, niv, cf, label); break;
    case 2:
	grid = useSpecial2DConstructor ?
	    new StructuredGrid(psurface, niv, njv, cfList, r_grid, s_grid, label)
            : new StructuredGrid(psurface, niv, njv, cfList, interpolation, label);
	break;
    case 3: grid = new StructuredGrid(pvolume, niv, njv, nkv, cfList, interpolation, label); break;
    default: throw new GeometryException("invalid number of dimensions");
    }
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid);
    return 1;
} // end newStructuredGrid()


extern(C) int get_measure_of_badness(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.measure_of_badness());
    return 1;
}

extern(C) int determine_rs_grids(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2
    if (narg != 2) {
        string errMsg = "Error in StructuredGrid:determine_rs_grids(). " ~
            "Expected call of form grid:determine_rs_grids(psurface).";
        luaL_error(L, errMsg.toStringz);
    }
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    if (!grid) {
        string errMsg = "Error in StructuredGrid:determine_rs_grids(). " ~
            "Argument 1 should be a 2D grid.";
        luaL_error(L, errMsg.toStringz);
    }
    auto psurface = checkSurface(L, 2);
    if (!psurface) {
        string errMsg = "Error in StructuredGrid:determine_rs_grids(). " ~
            "Argument 2 should be a ParametricSurface.";
        luaL_error(L, errMsg.toStringz);
    }
    auto cf = [new LinearFunction(), new LinearFunction(),
               new LinearFunction(), new LinearFunction()];
    double[4][4] r_grid = [[0.0, 1.0/3, 2.0/3, 1.0],
                           [0.0, 1.0/3, 2.0/3, 1.0],
                           [0.0, 1.0/3, 2.0/3, 1.0],
                           [0.0, 1.0/3, 2.0/3, 1.0]];
    double[4][4] s_grid = [[0.0, 0.0, 0.0, 0.0],
                           [1.0/3, 1.0/3, 1.0/3, 1.0/3],
                           [2.0/3, 2.0/3, 2.0/3, 2.0/3],
                           [1.0, 1.0, 1.0, 1.0]];
    grid.determine_rs_grids(psurface, cf, r_grid, s_grid);
    pushMatrixOfDoubles!4(L, r_grid);
    pushMatrixOfDoubles!4(L, s_grid);
    return 2;
}


extern(C) int makeSlabGrid(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg != 2) {
        string errMsg = "Error in StructuredGrid:makeSlabGrid{}. ";
        errMsg ~= "Expected a 2D StructuredGrid object. ";
        errMsg ~= "Be sure to use myGrid:makeSlabGrid{} rather than myGrid.makeSlabGrid{}.";
        luaL_error(L, errMsg.toStringz);
    }
    auto grid2D = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    if (grid2D is null || grid2D.nkv != 1) {
        string errMsg = "Error in StructuredGrid:makeSlabGrid{}. ";
        errMsg ~= "Expected a 2D StructuredGrid object.";
        luaL_error(L, errMsg.toStringz);
    }
    if ( narg < 2 || !lua_istable(L, 2) ) {
        string errMsg = "Error in call to StructuredGrid:makeSlabGrid{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 2, ["dz","symmetric","label"])) {
        string errMsg = "Error in call to StructuredGrid:makeSlabGrid{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    double dz;
    lua_getfield(L, 2, "dz");
    if (!lua_isnil(L, -1)) {
        dz = luaL_checknumber(L, -1);
    } else {
        dz = 0.0;
        string errMsg = "Error in call to StructuredGrid:makeSlabGrid{}. expected number for dz.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    bool symmetric = true;
    lua_getfield(L, 2, "symmetric");
    if (!lua_isnil(L, -1)) {
        symmetric = to!bool(lua_toboolean(L, -1));
    }
    lua_pop(L, 1);
    string label = grid2D.label; // carry over original label, by default
    lua_getfield(L, 2, "label");
    if (!lua_isnil(L, -1)) {
        label = to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    //
    lua_settop(L, 0); // clear stack
    auto grid3D = grid2D.makeSlabGrid(dz, symmetric, label);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid3D);
    return 1;
} // end makeSlabGrid()

extern(C) int makeWedgeGrid(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg != 2) {
        string errMsg = "Error in StructuredGrid:makeWedgeGrid{}. ";
        errMsg ~= "Expected a 2D StructuredGrid object and a table. ";
        errMsg ~= "Be sure to use myGrid:makeWedgeGrid{} rather than myGrid.makeWedgeGrid{}.";
        luaL_error(L, errMsg.toStringz);
    }
    auto grid2D = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    if (grid2D is null || grid2D.nkv != 1) {
        string errMsg = "Error in StructuredGrid:makeWedgeGrid{}. ";
        errMsg ~= "Expected a 2D StructuredGrid object.";
        luaL_error(L, errMsg.toStringz);
    }
    if ( narg < 2 || !lua_istable(L, 2) ) {
        string errMsg = "Error in call to StructuredGrid:makeWedgeGrid{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 2, ["dtheta","symmetric","label","nkv_new"])) {
        string errMsg = "Error in call to StructuredGrid:makeWedgeGrid{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    double dtheta;
    lua_getfield(L, 2, "dtheta");
    if (!lua_isnil(L, -1)) {
        dtheta = luaL_checknumber(L, -1);
    } else {
        dtheta = 0.0;
        string errMsg = "Error in call to StructuredGrid:makeWedgeGrid{}. expected number for dtheta.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    bool symmetric = true;
    lua_getfield(L, 2, "symmetric");
    if (!lua_isnil(L, -1)) {
        symmetric = to!bool(lua_toboolean(L, -1));
    }
    lua_pop(L, 1);
    string label = grid2D.label; // carry over original label, by default
    lua_getfield(L, 2, "label");
    if (!lua_isnil(L, -1)) {
        label = to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    int nkv_new = 2;
    lua_getfield(L, 2, "nkv_new");
    if (!lua_isnil(L, -1)) {
        nkv_new = luaL_checkint(L, -1);
    }
    lua_pop(L, 1);
    //
    lua_settop(L, 0); // clear stack
    auto grid3D = grid2D.makeWedgeGrid(dtheta, symmetric, label, nkv_new);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid3D);
    return 1;
} // end makeWedgeGrid()

extern(C) int importGridproGrid(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg == 0 ) {
        string errMsg = "Error in call to importGridproGrid(). " ~
            "At least one argument is required: the name of the Gridpro file.";
        luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, 1));
    double scale = 1.0;
    if ( narg >= 2 ) {
        scale = luaL_checknumber(L, 2);
    }
    auto sgrids = import_gridpro_grid(fname, scale);
    lua_newtable(L);
    foreach (i, grid; sgrids) {
        structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid);
        lua_rawseti(L, -2, to!int(i)+1);
    }
    return 1;
} // end importGridproGrid()

extern(C) int importPlot3DGrid(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg < 2 ) {
        string errMsg = "Error in call to importPlot3DGrid(). " ~
            "At least two arguments are required: the name of the Plot3D file and the dimensionality (2 or 3).";
        luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, 1));

    int dim = luaL_checkint(L, 2);

    double scale = 1.0;
    if ( narg >= 3 ) {
        scale = luaL_checknumber(L, 3);
    }
    auto sgrids = geom.grid.sgrid.importPlot3DGrid(fname, dim, scale);
    lua_newtable(L);
    foreach (i, grid; sgrids) {
        structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid);
        lua_rawseti(L, -2, to!int(i)+1);
    }
    return 1;
} // end importPlot3DGrid()


extern(C) int writeGridsAsPlot3D(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 3 ) {
        string errMsg = "Error in call to writeGridsAsPlot3D(). " ~
            "3 arguments are required: name of file, array of grids, and dimension/";
        luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, 1));
    if (!lua_istable(L, 2)) {
        string errMsg = "Error in call to writeGridsAsPlot3D(). " ~
            "The second argument should be an array of grids.";
        luaL_error(L, errMsg.toStringz);
    }
    int n = to!int(lua_objlen(L, 2));
    StructuredGrid[] grids;
    foreach (i; 1 .. n+1) {
        lua_rawgeti(L, 2, i);
        grids ~= checkObj!(StructuredGrid, StructuredGridMT)(L, -1);
        lua_pop(L, 1);
    }
    int dim = to!int(luaL_checkinteger(L, 3));
    geom.grid.sgrid.writeGridsAsPlot3D(fname, grids, dim);
    return 0;
}

extern(C) int rotateGridproBlocks(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg != 3) {
        string errMsg = "Error in call to rotateGridproBlocks(). " ~
            "3 arguments are required: grid, rotateSouthToThis, rotateWestToThis";
        luaL_error(L, errMsg.toStringz);
    }
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    auto rotateSouthToThis = to!string(luaL_checkstring(L, 2));
    auto rotateWestToThis = to!string(luaL_checkstring(L, 3));
    auto newGrid = rotate_gridpro_blocks(grid, rotateSouthToThis, rotateWestToThis);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, newGrid);
    return 1;
}

extern(C) int gridFaceSwap(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg != 3) {
        string errMsg = "Error in call to gridFaceSwap(). " ~
            "3 arguments are required: grid, swapNorthToSouth, swapEastToWest";
        luaL_error(L, errMsg.toStringz);
    }
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    auto swapNorthToSouth = to!bool(lua_toboolean(L, 2));
    auto swapEastToWest = to!bool(lua_toboolean(L, 3));
    auto newGrid = grid_faceswap(grid, swapNorthToSouth, swapEastToWest);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, newGrid);
    return 1;
}

void registerStructuredGrid(lua_State* L)
{
    // Register the StructuredGrid object
    luaL_newmetatable(L, StructuredGridMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newStructuredGrid);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &get_measure_of_badness);
    lua_setfield(L, -2, "measure_of_badness");
    lua_pushcfunction(L, &determine_rs_grids);
    lua_setfield(L, -2, "determine_rs_grids");
    lua_pushcfunction(L, &toStringObj!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyStructuredGrid!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &get_dimensions!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_dimensions");
    lua_pushcfunction(L, &get_type!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_type");
    lua_pushcfunction(L, &get_ncells!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_ncells");
    lua_pushcfunction(L, &get_niv);
    lua_setfield(L, -2, "get_niv");
    lua_pushcfunction(L, &get_njv);
    lua_setfield(L, -2, "get_njv");
    lua_pushcfunction(L, &get_nkv);
    lua_setfield(L, -2, "get_nkv");
    lua_pushcfunction(L, &get_vtx!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &set_vtx);
    lua_setfield(L, -2, "set_vtx");
    lua_pushcfunction(L, &get_corner_vtx);
    lua_setfield(L, -2, "get_corner_vtx");
    lua_pushcfunction(L, &set_boundary_tags);
    lua_setfield(L, -2, "set_tags");
    lua_pushcfunction(L, &get_boundary_tag);
    lua_setfield(L, -2, "get_tag");
    lua_pushcfunction(L, &subgrid);
    lua_setfield(L, -2, "subgrid");
    lua_pushcfunction(L, &get_boundary_grid);
    lua_setfield(L, -2, "get_boundary_grid");
    lua_pushcfunction(L, &cellVolume!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "cellVolume");
    lua_pushcfunction(L, &cellCentroid!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "cellCentroid");
    lua_pushcfunction(L, &write_to_vtk_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_vtk_file");
    lua_pushcfunction(L, &write_to_stl_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_stl_file");
    lua_pushcfunction(L, &write_to_gzip_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_gzip_file");
    lua_pushcfunction(L, &write_to_raw_binary_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_raw_binary_file");
    lua_pushcfunction(L, &rotate!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "rotate");
    lua_pushcfunction(L, &joinGrid);
    lua_setfield(L, -2, "joinGrid");
    lua_pushcfunction(L, &rotateGrid);
    lua_setfield(L, -2, "rotateGrid");
    lua_pushcfunction(L, &find_nearest_cell_centre_sg);
    lua_setfield(L, -2, "find_nearest_cell_centre");
    lua_pushcfunction(L, &makeSlabGrid);
    lua_setfield(L, -2, "makeSlabGrid");
    lua_pushcfunction(L, &makeWedgeGrid);
    lua_setfield(L, -2, "makeWedgeGrid");

    lua_setglobal(L, StructuredGridMT.toStringz);

    // Global functions available for use
    lua_pushcfunction(L, &importGridproGrid);
    lua_setglobal(L, "importGridproGrid");
    lua_pushcfunction(L, &importPlot3DGrid);
    lua_setglobal(L, "importPlot3DGrid");
    lua_pushcfunction(L, &writeGridsAsPlot3D);
    lua_setglobal(L, "writeGridsAsPlot3D");
    lua_pushcfunction(L, &rotateGridproBlocks);
    lua_setglobal(L, "rotateGridproBlocks");
    lua_pushcfunction(L, &gridFaceSwap);
    lua_setglobal(L, "gridFaceSwap");

} // end registerStructuredGrid()
