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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Matrix:new{}; ";
        errMsg ~= "maybe you tried Matrix.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this".
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in Matrix:new{} constructor.\n";
        errMsg ~= "A table of keyword arguments is expected.\n";
        luaL_error(L, errMsg.toStringz);
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
            luaL_error(L, errMsg.toStringz);
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
            luaL_error(L, errMsg.toStringz);
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
            luaL_error(L, errMsg.toStringz);
        }
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "orient");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in Matrix:new{} constructor.\n";
        errMsg ~= "You supplied a 'vec' keyword argument, but no valid 'orient' argument.\n";
        luaL_error(L, errMsg.toStringz);
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
    lua_pushinteger(L, mat.ncols);
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

extern(C) int getSlice(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    int row0 = luaL_checkint(L, 2);
    int row1 = luaL_checkint(L, 3);
    int col0 = luaL_checkint(L, 4);
    int col1 = luaL_checkint(L, 5);
    if (row0 >= 0 && row0 < row1 && row1 <= M.nrows &&
        col0 >= 0 && col0 < col1 && col1 <= M.ncols) {
        auto mat = new Matrix!double(row1-row0, col1-col0);
        foreach (i; row0 .. row1) {
            foreach (j; col0 .. col1) {
                mat[i-row0, j-col0] = M[i, j];
            }
        }
        matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, mat);
    } else {
        lua_pushnil(L);
        string errMsg = "Error in Matrix:getSlice() method.\n";
        errMsg ~= "Something is wrong with the specified ranges.\n";
        luaL_error(L, errMsg.toStringz);
    }
    return 1;
} // end getSlice()

extern(C) int setSlice(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    int row0 = luaL_checkint(L, 2);
    int row1 = luaL_checkint(L, 3);
    int col0 = luaL_checkint(L, 4);
    int col1 = luaL_checkint(L, 5);
    auto mat = checkMatrix(L, 6);
    if ( (row1-row0) == mat.nrows && (col1-col0) == mat.ncols ) {
        if (row0 >= 0 && row0 < row1 && row1 <= M.nrows &&
            col0 >= 0 && col0 < col1 && col1 <= M.ncols) {
            foreach (i; row0 .. row1) {
                foreach (j; col0 .. col1) {
                    M[i, j] = mat[i-row0, j-col0];
                }
            }
            matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, mat);
        } else {
            lua_pushnil(L);
            string errMsg = "Error in Matrix:setSlice() method.\n";
            errMsg ~= "Something is wrong with the specified ranges.\n";
            errMsg ~= "Specified ranges don't match Matrix.\n";
            luaL_error(L, errMsg.toStringz);
        }
    } else {
        lua_pushnil(L);
        string errMsg = "Error in Matrix:setSlice() method.\n";
        errMsg ~= "Something is wrong with the specified ranges.\n";
        errMsg ~= "Specified ranges don't match Slice.\n";
        luaL_error(L, errMsg.toStringz);
    }

    return 1;
} // end setSlice()

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

extern(C) int multiply(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    auto rhs = checkMatrix(L, 2);
    auto result = M * rhs;
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, result);
    return 1;
}

extern(C) int divide(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    auto rhs = checkMatrix(L, 2);
    auto result = M / rhs;
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, result);
    return 1;
}

extern(C) int scale(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    double rhs = luaL_checknumber(L, 2);
    auto result = M * rhs;
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, result);
    return 1;
}

extern(C) int scale_in_place(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    double rhs = luaL_checknumber(L, 2);
    M.scale(rhs);
    return 0;
}

extern(C) int add_in_place(lua_State *L)
{
    // Use as M:dot(other), where other is a matrix object of matching size.
    // Euqivalent to M = M + other
    auto M = checkMatrix(L, 1);
    auto other = checkMatrix(L, 2);
    M.add(other);
    return 0;
}

extern(C) int dotProduct(lua_State *L)
{
    // Use as result = M:dot(other)
    auto M = checkMatrix(L, 1); // this
    auto other = checkMatrix(L, 2);
    if (M.ncols == other.nrows) {
        // We have consistent dimensions and can compute the dot product.
        auto mat = new Matrix!double(M.nrows, other.ncols);
        foreach (i; 0 .. M.nrows) {
            foreach (j; 0 .. other.ncols) {
                double v = 0.0;
                foreach (k; 0 .. M.ncols) { v += M[i,k] * other[k,j]; }
                mat[i, j] = v;
            }
        }
        matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, mat);
    } else {
        lua_pushnil(L);
        string errMsg = "Error in Matrix:dot() method. ";
        errMsg ~= "Inconsistent dimensions.\n";
        luaL_error(L, errMsg.toStringz);
    }
    return 1;
} // end dotProduct()

extern(C) int dotProduct2(lua_State *L)
{
    // Use as M:dot2(other, result) or as M.dot2(M, other, result)
    auto M = checkMatrix(L, 1); // this
    auto other = checkMatrix(L, 2);
    auto result = checkMatrix(L, 3); // result matrix must exist already
    if ((M.ncols == other.nrows) &&
        (result.nrows == M.nrows) &&
        (result.ncols == other.ncols)) {
        // We have consistent dimensions and can compute the dot product.
        dot(M, other, result);
    } else {
        string errMsg = "Error in Matrix:dot2() method. ";
        errMsg ~= "Inconsistent dimensions.\n";
        luaL_error(L, errMsg.toStringz);
    }
    return 0;
} // end dotProduct2()

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

extern(C) int onesMatrix(lua_State *L)
{
    auto M = checkMatrix(L, 1);
    M.ones();
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

extern(C) int ones(lua_State *L)
{
    size_t nrows = luaL_checkint(L, 1);
    size_t ncols = luaL_checkint(L, 2);
    auto M = nm.bbla.ones!double(nrows, ncols);
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

extern(C) int hstack(lua_State *L)
{
    // We expect a single argument that is a table (array) of Matrix objects.
    if (!lua_istable(L, 1)) {
        string errMsg = "Error in call to hstack().; " ~
            "An array of Matrix objects is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Matrix objects at array positions within that table.
    Matrix!double[] mList;
    int position = 1;
    while (true) {
        lua_rawgeti(L, 1, position);
        if (lua_isnil(L, -1)) { lua_pop(L, 1); break; }
        auto m = checkMatrix(L, -1);
        lua_pop(L, 1);
        if (m is null) break;
        mList ~= m;
        ++position;
    }
    lua_pop(L, 1); // dispose of original table
    if (mList.length == 0) {
        string errMsg = "Error in call to hstack(). No valid Matrix objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto result = nm.bbla.hstack!double(mList);
    matrixStore ~= pushObj!(Matrix!double, MatrixMT)(L, result);
    return 1;
} // end hstack()

extern(C) int gaussJordan(lua_State *L)
{
    auto Ab = checkMatrix(L, 1);
    nm.bbla.gaussJordanElimination!double(Ab);
    return 0;
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
    lua_pushcfunction(L, &getSlice);
    lua_setfield(L, -2, "getSlice");
    lua_pushcfunction(L, &setSlice);
    lua_setfield(L, -2, "setSlice");
    lua_pushcfunction(L, &add);
    lua_setfield(L, -2, "__add");
    lua_pushcfunction(L, &subtract);
    lua_setfield(L, -2, "__sub");
    lua_pushcfunction(L, &multiply);
    lua_setfield(L, -2, "__mul");
    lua_pushcfunction(L, &divide);
    lua_setfield(L, -2, "__div");
    lua_pushcfunction(L, &scale);
    lua_setfield(L, -2, "scale");
    lua_pushcfunction(L, &scale_in_place);
    lua_setfield(L, -2, "scale_in_place");
    lua_pushcfunction(L, &add_in_place);
    lua_setfield(L, -2, "add_in_place");
    lua_pushcfunction(L, &dotProduct);
    lua_setfield(L, -2, "dot");
    lua_pushcfunction(L, &dotProduct2);
    lua_setfield(L, -2, "dot2");
    lua_pushcfunction(L, &tostring);
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &zerosMatrix);
    lua_setfield(L, -2, "zeros");
    lua_pushcfunction(L, &onesMatrix);
    lua_setfield(L, -2, "ones");
    lua_pushcfunction(L, &eyeMatrix);
    lua_setfield(L, -2, "eye");

    lua_setglobal(L, MatrixMT.toStringz);

    lua_pushcfunction(L, &zeros);
    lua_setglobal(L, "zeros");
    lua_pushcfunction(L, &ones);
    lua_setglobal(L, "ones");
    lua_pushcfunction(L, &eye);
    lua_setglobal(L, "eye");
    lua_pushcfunction(L, &transpose);
    lua_setglobal(L, "transpose");
    lua_pushcfunction(L, &solve);
    lua_setglobal(L, "solve");
    lua_pushcfunction(L, &hstack);
    lua_setglobal(L, "hstack");
    lua_pushcfunction(L, &gaussJordan);
    lua_setglobal(L, "gaussJordan");

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
b2 = zeros(4, 1)
A:dot2(x, b2)
assert(math.abs(b2:get(0,0) - b:get(0,0)) < 1.0e-6)
assert(math.abs(b2:get(1,0) - b:get(1,0)) < 1.0e-6)
assert(math.abs(b2:get(2,0) - b:get(2,0)) < 1.0e-6)
assert(math.abs(b2:get(3,0) - b:get(3,0)) < 1.0e-6)
        `;
        assert(luaL_dostring(L, toStringz(testCode)) == 0, failedUnitTest());
        return 0;
    }
}
