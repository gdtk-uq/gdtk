module util.lua;


// A.M.M - 01/12/2023
// Create a trait that allows us to check if a type can be cast to another.
// We require this, as the traditional "isFloatingPoint" cannot be extended
// to custom types, such as Complex and Dual. 
public template canCastTo(T, U) {
    import std.traits : ReturnType;
    // Check if T can be implicitly or explicitly cast to U
    // A more complete version will hijack the `to!` checks
    static if ( is(T : U) || is(ReturnType!(T.opCast!U) == U) ) {
        enum canCastTo = true;
    } else {
        enum canCastTo = false;
    }
}

@("canCastTo Trait")
unittest {
    struct MyType(T) {
        T inner;

        R opCast(R : T)() const {
            return inner;
        }
    }

    auto value = MyType!(double)(1.0);
    auto casted = cast(double)(value);
    assert( is(typeof(value) == MyType!double) );
    assert( is(typeof(casted) == double) );

    static assert(canCastTo!(MyType!double, double));
    static assert(canCastTo!(MyType!float, float));
    static assert(canCastTo!(MyType!double, float));
    static assert(canCastTo!(MyType!float, double));

    static assert(canCastTo!(double, double));
    static assert(canCastTo!(float, double));
    static assert(canCastTo!(double, float));

    int i1 = 1;
    double d1 = cast(double) i1;
    assert( is(typeof(d1) == double) );
    assert( d1 == 1.0 );
    
    static assert(canCastTo!(int, double));
    static assert(canCastTo!(char, double));
    static assert(canCastTo!(bool, double));
    static assert(!canCastTo!(string, double));
}

extern (C):

const LUAI_MAXSTACK = 1_000_000;


const LUA_REGISTRYINDEX = -LUAI_MAXSTACK - 1000;
const LUA_RIDX_GLOBALS = 2;
const LUA_MULTRET = -1;

const LUA_TNIL = 0;
const LUA_TBOOLEAN = 1;
const LUA_TTABLE = 5;
const LUA_TFUNCTION = 6;

const LUA_GCSTOP = 0;
const LUA_GCRESTART = 1;
const LUA_GCCOLLECT = 2;
const LUA_GCCOUNT = 3;
const LUA_GCCOUNTB = 4;
const LUA_GCSTEP = 5;
const LUA_GCSETPAUSE = 6;
const LUA_GCSETSTEPMUL = 7;

struct lua_State {}
alias lua_CFunction = int function(lua_State* L);
alias lua_KFunction = int function(lua_State *L, int status, lua_KContext ctx);

alias lua_Number = double;
alias lua_KContext = int;
alias lua_Integer = ptrdiff_t;

void lua_close(lua_State* L);
int lua_gettop(lua_State* L) nothrow;
void lua_settop(lua_State *L, int idx) nothrow;
void lua_pushvalue(lua_State* L, int idx) nothrow;

void lua_pop(lua_State* L, int n) nothrow { lua_settop(L, -(n)-1); }
void lua_rotate (lua_State *L, int idx, int n);
void lua_remove(lua_State* L, int idx) {lua_rotate(L, idx, -1); lua_pop(L, 1);}
void lua_replace(lua_State* L, int idx);

void lua_pushnil(lua_State* L);
void lua_pushnumber(lua_State* L, lua_Number n);
// Templated version, to allow for Complex and Dual numbers without coupling packages
void lua_pushnumber(T)(lua_State* L, T n) if (canCastTo!(T, lua_Number)) { lua_pushnumber(L, cast(lua_Number) n); }
void lua_pushinteger(lua_State* L, lua_Integer n);
void lua_pushstring(lua_State* L, const(char)* s);
void lua_pushcclosure(lua_State* L, lua_CFunction fn, int n);
void lua_pushboolean(lua_State* L, int b);

void lua_getfield(lua_State* L, int idx, const(char)* k);
void lua_pushglobaltable(lua_State* L) {lua_rawgeti(L, LUA_REGISTRYINDEX, LUA_RIDX_GLOBALS);}

void lua_setfield(lua_State* L, int idx, const(char)* k);
void lua_rawseti(lua_State* L, int idx, int n);
void lua_rawgeti(lua_State* L, int idx, int n) nothrow;
void  lua_createtable(lua_State *L, int narr, int nrec);
void *lua_newuserdatauv (lua_State *L, size_t size, int nuvalue);
void* lua_newuserdata(lua_State* L, size_t sz) {return lua_newuserdatauv(L,sz,1);}

int lua_getmetatable(lua_State* L, int objindex);
const(char)* lua_tolstring(lua_State *L, int idx, size_t *len);
void lua_setglobal(lua_State* L, const(char)* s);
void lua_getglobal(lua_State* L, const(char)* s);
int lua_pcallk(lua_State *L, int nargs, int nresults, int errfunc, lua_KContext ctx, lua_KFunction k);
int lua_pcall(lua_State *L, int nargs, int nresults, int errfunc) {return lua_pcallk(L, nargs, nresults, errfunc, 0, null);}
int lua_gc(lua_State *L, int what, int data);

int lua_setmetatable(lua_State* L, int idx);
int lua_next(lua_State *L, int idx);

int lua_isnumber(lua_State* L, int idx);
int lua_isinteger(lua_State* L, int idx);
int lua_isstring(lua_State* L, int idx);
int lua_isuserdata(lua_State* L, int idx);
int lua_type(lua_State* L, int idx) nothrow;
bool lua_rawequal(lua_State* L, int idx1, int idx2);

bool lua_istable(lua_State* L, int n) { return lua_type(L, n) == LUA_TTABLE; }
bool lua_isnil(lua_State* L, int n) { return lua_type(L, n) == LUA_TNIL; }
bool lua_isboolean(lua_State* L, int n) { return lua_type(L, n) == LUA_TBOOLEAN; }

void lua_newtable(lua_State* L) { lua_createtable(L, 0, 0); }
void lua_pushcfunction(lua_State* L, lua_CFunction f) { lua_pushcclosure(L, f, 0); }
bool lua_isfunction(lua_State* L, int n) { return lua_type(L, n) == LUA_TFUNCTION; }

lua_Number lua_tonumberx(lua_State *L, int idx, int *isnum);
lua_Number lua_tonumber(lua_State *L, int idx) {return lua_tonumberx(L,idx,null);}
lua_Integer lua_tointegerx(lua_State* L, int idx, int *pisnum);
lua_Integer lua_tointeger(lua_State* L, int idx){return lua_tointegerx(L,idx,null);}
bool lua_toboolean(lua_State* L, int idx);
const(char)* lua_tostring(lua_State* L, int i) { return lua_tolstring(L, i, null); }
size_t lua_rawlen (lua_State *L, int idx);
size_t lua_objlen(lua_State* L, int idx) {return lua_rawlen(L,idx);}
void* lua_touserdata(lua_State* L, int idx);

const(char)* luaL_checklstring(lua_State *L, int numArg, size_t *l);
lua_Number luaL_checknumber(lua_State* L, int numArg);
lua_Integer luaL_checkinteger(lua_State* L, int numArg);
int luaL_newmetatable(lua_State* L, const(char)* tname);
void* luaL_checkudata(lua_State* L, int ud, const(char)* tname);
int luaL_error(lua_State* L, const(char)* fmt, ...);

lua_State* luaL_newstate();
const(char)* luaL_checkstring(lua_State* L, int n) { return luaL_checklstring(L, n, null); }
int luaL_checkint(lua_State* L, int n) { return cast(int) luaL_checkinteger(L, n); }

void luaL_openlibs(lua_State* L);
int luaL_loadfilex(lua_State *L, const char *filename, const char *mode);
int luaL_loadfile(lua_State* L, const(char)* filename){return luaL_loadfilex(L,filename,null);}
int luaL_loadstring(lua_State* L, const(char)* s);
int luaL_dofile(lua_State* L, const(char)* fn) { return luaL_loadfile(L, fn) || lua_pcall(L, 0, LUA_MULTRET, 0); }
int luaL_dostring(lua_State* L, const(char)* s) { return luaL_loadstring(L, s) || lua_pcall(L, 0, LUA_MULTRET, 0); }
void luaL_getmetatable(lua_State* L, const(char)* s) { lua_getfield(L, LUA_REGISTRYINDEX, s); }
