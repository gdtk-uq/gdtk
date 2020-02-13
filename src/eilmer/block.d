// block.d
// Base class for blocks of cells, for use within the Eilmer flow solver.
// Kyle A. Damm 2020-02-11 first cut.

module block;

import globalconfig;
import util.lua;

class Block {
public:
    int id;      // block identifier: assumed to be the same as the block number.
    string label;
    bool active; // if true, block participates in the time integration
                 // The active flag is used principally for the block-marching calculation,
                 // where we want to integrate a few blocks at a time.
    lua_State* myL;
    LocalConfig myConfig;
    
    this(int id, string label)
    {
        this.id = id;
        this.label = label;

        // Lua interpreter for the block.
        myL = luaL_newstate();
        luaL_openlibs(myL);
        lua_pushinteger(myL, id);
        lua_setglobal(myL, "blkId");
    }

} // end class Block
