/**
 * A Lua interface to D linear algebra package.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2017-04-03
 */

module nm.luabbla;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import nm.bbla;

immutable string MatrixMT = "Matrix"; // Name of Matrix metatable

// A place to hang on to references to objects that are pushed into the Lua domain.
// We don't want the D garbage collector to get rid of them too early.
static const(Matrix!double)[] matrixStore;

Matrix!double checkMatrix(lua_State *L, int index)
{
    return checkObj!(Matrix!double, MatrixMT)(L, index);
}

/**
 * This function implements the constructor for a Matrix!double
 * in the Lua interface.
 *
 * Construction of a Matrix in Lua will look like any of the following:
 * ----------------
 * a = Matrix:new{n=4} -- for a square matrix, 4x4
 * b = Matrix:new{nrows=4, ncols=9} -- for a 4x9 matrix
 * c = Matrix:new{other=a} -- for construction from an existing matrix
 * v = {7, 8, 15, 9} 
 * d = Matrix:new{vec=v, orient='column'} -- for taking a vector and
 *                                        -- making a column or row matrix.
 * ----------------
 * 
 * The preference for construction is as per the order above.
 * We will look for the keyword arguments in the supplied table
 * and act on a complete set when found.
 */

extern(C) int newMatrix(lua_State *L)
{
    lua_remove(L, 1); // remove first argument "this".
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in Matrix:new{} constructor.\n";
        errMsg ~= "A table of keyword arguments is expected.\n";
        throw new Error(errMsg);
    }
    
    // Look for square-matrix constructor: Matrix:new{n=...}
    lua_getfield(L, 1, "n");
    if ( !lua_isnil(L, -1) ) {
        int n = luaL_checkint(L, -1);
        lua_pop(L, 1);
        auto mat = new Matrix!double(n);
        matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, mat);
        return 1;
    }
    lua_pop(L, 1);
    
    // Look for matrix constructor: Matrix:new{nrows= ..., ncols= ....}
    lua_getfield(L, 1, "nrows");
    if ( !lua_isnil(L, -1) ) {
        int nrows = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_getfield(L, 1, "ncols");
        if ( !lua_isnumber(L, -1) ) {
            string errMsg = "Error in Matrix:new{} constructor.\n";
            errMsg ~= "You have supplied an 'nrows' value but not a corresponding 'ncols' value.\n";
            errMsg ~= "The form of this constructor is Matrix:new{nrows=..., ncols=...}\n";
            throw new Error(errMsg);
        }
        int ncols = luaL_checkint(L, -1);
        lua_pop(L, 1);

        auto mat = new Matrix!double(nrows, ncols);
        matrixStore ~=  pushObj!(Matrix!double, MatrixMT)(L, mat);
        return 1;
    }
    lua_pop(L, 1);

    // Look for matrix constructor: Matrix:new{other=...}
    lua_getfield(L, 1, "other");
    if ( !lua_isnil(L, -1) ) {
        Matrix!double a = checkMatrix(L, -1);
        if (a is null) {
            string errMsg = "Error in Matrix:new{} constructor.\n";
            errMsg ~= "You have used the 'other' keyword argument but have not supplied a valid Matrix object.\n";
            throw new Error(errMsg);
        }
        auto mat = new Matrix!double(a);
        matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, mat);
        return 1;
    }
    lua_pop(L, 1);

    // Look for matrix constructor: Matrix:new{vec=..., orient=".."}
    lua_getfield(L, 1, "vec");
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in Matrix:new{} constructor.\n";
        errMsg ~= "You have not supplied a valid array to 'vec' keyword.\n";
        throw new Error(errMsg);
    }
    auto n = to!int(lua_objlen(L, -1));
    double[] vec;
    foreach (i; 1..n+1) {
        lua_rawgeti(L, -1, i);
        if ( lua_isnumber(L, -1) ) {
            vec ~= lua_tonumber(L, -1);
        }
        else {
            string errMsg = "Error in Matrix:new{} constructor.\n";
            errMsg ~= "One of the values in your 'vec' array is not a valid number.\n";
            throw new Error(errMsg);
        }
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "orient");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in Matrix:new{} constructor.\n";
        errMsg ~= "You supplied a 'vec' keyword argument, but no valid 'orient' argument.\n";
        throw new Error(errMsg);
    }
    string orient = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);

    auto mat = new Matrix!double(vec, orient);
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, mat);
    return 1;
}

extern(C) int nrows(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    lua_pushinteger(L, M.nrows);
    return 1;
}

extern(C) int ncols(lua_State *L)
{
    auto mat = checkMatrix(L, 1);
    lua_pushinteger(L, mat.nrows);
    return 1;
}

extern(C) int get(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    int i = luaL_checkint(L, 2);
    int j = luaL_checkint(L, 3);
    auto val = M[i,j];
    lua_pushnumber(L, val);
    return 1;
}

extern(C) int set(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    int i = luaL_checkint(L, 2);
    int j = luaL_checkint(L, 3);
    auto val = lua_tonumber(L, 4);
    M[i,j] = val;
    return 0;
}

extern(C) int add(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    auto rhs = checkMatrix(L, 2);
    auto result = M + rhs;
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, result);
    return 1;
}

extern(C) int subtract(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    auto rhs = checkMatrix(L, 2);
    auto result = M - rhs;
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, result);
    return 1;
}

extern(C) int tostring(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    lua_pushstring(L, toStringz(M.toString()));
    return 1;
}


extern(C) int zerosMatrix(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    M.zeros();
    return 0;
}

extern(C) int eyeMatrix(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    M.eye();
    return 0;
}

extern(C) int zeros(lua_State *L)
{
    size_t nrows = luaL_checkint(L, 1);
    size_t ncols = luaL_checkint(L, 2);
    auto M = nm.bbla.zeros!double(nrows, ncols);
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, M);
    return 1;

}

extern(C) int eye(lua_State *L)
{
    size_t n = luaL_checkint(L, 1);
    auto I = nm.bbla.eye!double(n);
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, I);
    return 1;

}

extern(C) int transpose(lua_State *L)
{
    auto other = checkMatrix(L, 1);
    auto M = nm.bbla.transpose(other);
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, M);
    return 1;
}

extern(C) int solve(lua_State *L)
{
    // Unlike the D implementation, we do this in one pass.
    // The D implementation will be much more efficient for
    // many uses of the decomposition. In this case, we supply
    // a simple result to the equation A x = b.
    auto A = checkMatrix(L, 1);
    auto b = checkMatrix(L, 2);
    auto c = new Matrix!double(A);
    auto perm = decomp!double(c);
    auto x = new Matrix!double(b);
    nm.bbla.solve!double(c, x, perm);
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, x);
    return 1;
}

void registerBBLA(lua_State *L)
{
    luaL_newmetatable(L, MatrixMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newMatrix);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &nrows);
    lua_setfield(L, -2, "nrows");
    lua_pushcfunction(L, &ncols);
    lua_setfield(L, -2, "ncols");
    lua_pushcfunction(L, &get);
    lua_setfield(L, -2, "get");
    lua_pushcfunction(L, &set);
    lua_setfield(L, -2, "set");
    lua_pushcfunction(L, &add);
    lua_setfield(L, -2, "__add");
    lua_pushcfunction(L, &subtract);
    lua_setfield(L, -2, "__sub");
    lua_pushcfunction(L, &tostring);
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &zerosMatrix);
    lua_setfield(L, -2, "zeros");
    lua_pushcfunction(L, &eyeMatrix);
    lua_setfield(L, -2, "eye");

    lua_setglobal(L, MatrixMT.toStringz);

    lua_pushcfunction(L, &zeros);
    lua_setglobal(L, "zeros");
    lua_pushcfunction(L, &eye);
    lua_setglobal(L, "eye");
    lua_pushcfunction(L, &transpose);
    lua_setglobal(L, "transpose");
    lua_pushcfunction(L, &solve);
    lua_setglobal(L, "solve");

}

version(luabbla_test) {
    import util.msg_service;
    int main() {
        auto L = luaL_newstate();
        luaL_openlibs(L);
        registerBBLA(L);
        string testCode = `
a = Matrix:new{n=4}
assert(a:nrows() == 4)
assert(a:ncols() == 4)
a:zeros()
a:set(0,0, 5.0)
assert(a:get(0,0) == 5.0)
b = Matrix:new{nrows=4, ncols=4}
b:zeros()
c = b + a
assert(c:get(0,0) == 5.0)
d = b - a
assert(d:get(0,0) == -5.0)
a:eye()
assert(a:get(1,1) == 1.0)
e = zeros(4, 2)
assert(e:get(3,1) == 0.0)
f = eye(5)
assert(f:get(2,2) == 1.0)
A = Matrix:new{n=4}
A:set(0,0, 0.0); A:set(0,1, 2.0); A:set(0,2, 0.0); A:set(0,3, 1.0); 
A:set(1,0, 2.0); A:set(1,1, 2.0); A:set(1,2, 3.0); A:set(1,3, 2.0); 
A:set(2,0, 4.0); A:set(2,1, -3.0); A:set(2,2, 0.0); A:set(2,3, 1.0); 
A:set(3,0, 6.0); A:set(3,1, 1.0); A:set(3,2, -6.0); A:set(3,3, -5.0); 
b = Matrix:new{vec={0.0, -2.0, -7.0, 6.0}, orient='column'}
x = solve(A, b)
assert(math.abs(x:get(0,0) - -0.5) < 1.0e-6)
assert(math.abs(x:get(1,0) - 1.0) < 1.0e-6)
assert(math.abs(x:get(2,0) - 0.333333) < 1.0e-6)
assert(math.abs(x:get(3,0) - -2.0) < 1.0e-6)


        `;
        assert(luaL_dostring(L, toStringz(testCode)) == 0, failedUnitTest());
        return 0;
    }
}

