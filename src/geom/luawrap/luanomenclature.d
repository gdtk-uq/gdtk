/**
 * Expose the Face name mapping to Lua interpreter.
 *
 * This could be shoved in one of the other wrapped geometry functions.
 * however it seemed cleaner to keep the parallel pattern of:
 * modulename.d --> luamodulename.d
 *
 * Authors: Rowan J. Gollan and Peter A. Jacobs
 * Date: 2023-07-06
 */

module geom.luawrap.luanomenclature;

import util.lua;
import geom.elements.nomenclature;

void registerGeomNomenclature(lua_State *L)
{
    lua_newtable(L);
    lua_pushinteger(L, face_index("west")); lua_setfield(L, -2, "west");
    lua_pushinteger(L, face_index("east")); lua_setfield(L, -2, "east");
    lua_pushinteger(L, face_index("south")); lua_setfield(L, -2, "south");
    lua_pushinteger(L, face_index("north")); lua_setfield(L, -2, "north");
    lua_pushinteger(L, face_index("bottom")); lua_setfield(L, -2, "bottom");
    lua_pushinteger(L, face_index("top")); lua_setfield(L, -2, "top");
    lua_pushinteger(L, face_index("iminus")); lua_setfield(L, -2, "iminus");
    lua_pushinteger(L, face_index("iplus")); lua_setfield(L, -2, "iplus");
    lua_pushinteger(L, face_index("jminus")); lua_setfield(L, -2, "jminus");
    lua_pushinteger(L, face_index("jplus")); lua_setfield(L, -2, "jplus");
    lua_pushinteger(L, face_index("kminus")); lua_setfield(L, -2, "kminus");
    lua_pushinteger(L, face_index("kplus")); lua_setfield(L, -2, "kplus");
    lua_setglobal(L, "Face");
}

