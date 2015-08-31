module util.lua;

extern (C):
const LUA_REGISTRYINDEX = -10000;
const LUA_GLOBALSINDEX = -10002;
const LUA_MULTRET = -1;

const LUA_TNIL = 0;
const LUA_TBOOLEAN = 1;
const LUA_TTABLE = 5;
const LUA_TFUNCTION = 6;

struct lua_State {}
alias int function(lua_State* L) lua_CFunction;
alias double lua_Number;
alias ptrdiff_t lua_Integer;
void lua_close(lua_State* L);
int lua_gettop(lua_State* L) nothrow;
void lua_settop(lua_State *L, int idx) nothrow;
void lua_pushvalue(lua_State* L, int idx) nothrow;
void lua_remove(lua_State* L, int idx);
void lua_replace(lua_State* L, int idx);

void lua_pushnil(lua_State* L);
void lua_pushnumber(lua_State* L, lua_Number n);
void lua_pushinteger(lua_State* L, lua_Integer n);
void lua_pushstring(lua_State* L, const(char)* s);
void lua_pushcclosure(lua_State* L, lua_CFunction fn, int n);
void lua_pushboolean(lua_State* L, int b);

void lua_getfield(lua_State* L, int idx, const(char)* k);

void lua_setfield(lua_State* L, int idx, const(char)* k);
void lua_rawseti(lua_State* L, int idx, int n);
void lua_rawgeti(lua_State* L, int idx, int n) nothrow;
void  lua_createtable(lua_State *L, int narr, int nrec);
void* lua_newuserdata(lua_State* L, size_t sz);
int lua_getmetatable(lua_State* L, int objindex);
const(char)* lua_tolstring(lua_State *L, int idx, size_t *len);
void lua_setglobal(lua_State* L, const(char)* s) { lua_setfield(L, LUA_GLOBALSINDEX, s); }
void lua_getglobal(lua_State* L, const(char)* s) { lua_getfield(L, LUA_GLOBALSINDEX, s); }
int lua_pcall(lua_State *L, int nargs, int nresults, int errfunc);


int lua_setmetatable(lua_State* L, int idx);

int lua_isnumber(lua_State* L, int idx);
int lua_isstring(lua_State* L, int idx);
int lua_isuserdata(lua_State* L, int idx);
int lua_type(lua_State* L, int idx) nothrow;
bool lua_rawequal(lua_State* L, int idx1, int idx2);

bool lua_istable(lua_State* L, int n) { return lua_type(L, n) == LUA_TTABLE; }
bool lua_isnil(lua_State* L, int n) { return lua_type(L, n) == LUA_TNIL; }
bool lua_isboolean(lua_State* L, int n) { return lua_type(L, n) == LUA_TBOOLEAN; }

void lua_pop(lua_State* L, int n) nothrow { lua_settop(L, -(n)-1); }
void lua_newtable(lua_State* L) { lua_createtable(L, 0, 0); }
void lua_pushcfunction(lua_State* L, lua_CFunction f) { lua_pushcclosure(L, f, 0); }
bool lua_isfunction(lua_State* L, int n) { return lua_type(L, n) == LUA_TFUNCTION; }

lua_Number lua_tonumber(lua_State *L, int idx);
lua_Integer lua_tointeger(lua_State* L, int idx);
bool lua_toboolean(lua_State* L, int idx);
const(char)* lua_tostring(lua_State* L, int i) { return lua_tolstring(L, i, null); }
size_t lua_objlen(lua_State* L, int idx);
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
int luaL_loadfile(lua_State* L, const(char)* filename);
int luaL_loadstring(lua_State* L, const(char)* s);
int luaL_dofile(lua_State* L, const(char)* fn) { return luaL_loadfile(L, fn) || lua_pcall(L, 0, LUA_MULTRET, 0); }
int luaL_dostring(lua_State* L, const(char)* s) { return luaL_loadstring(L, s) || lua_pcall(L, 0, LUA_MULTRET, 0); }
void luaL_getmetatable(lua_State* L, const(char)* s) { lua_getfield(L, LUA_REGISTRYINDEX, s); }





