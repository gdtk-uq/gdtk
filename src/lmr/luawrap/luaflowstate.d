/**
 * luaflowstate.d
 * Lua interface to the FlowState class for users of the prep program.
 *
 * Authors: Rowan G. and Peter J.
 * Version: Initial cut.
 */

module luaflowstate;

import std.algorithm;
import std.array;
import std.format;
import std.stdio;
import std.string;
import std.conv;
import std.traits;
import gzip;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas;
import gas.luagas_model;
import flowstate;
import geom;
import geom.luawrap;
import globalconfig;

// Name for FlowState objects in Lua scripts are a pure Lua play,
// however, we wish to provide a few FlowState-related functions
// that have been coded in Dlang.
// These functions are provided in a Lua table called _FlowState.

immutable string[] validFlowStateFields = ["p", "T", "T_modes", "p_e",
                                           "quality", "massf",
                                           "mu", "k",
                                           "velx", "vely", "velz",
                                           "Bx", "By", "Bz", "psi", "divB",
                                           "turb", "mu_t", "k_t", "S"];
static const(FlowState*)[] flowStateStore;

FlowState* checkFlowState(lua_State* L, int index)
{
    size_t fs_id = to!size_t(luaL_checknumber(L, index));
    if (fs_id < flowStateStore.length) {
        return cast(FlowState*)flowStateStore[fs_id];
    } else {
        throw new Error(format("id number to FlowState is not valid: %d", fs_id));
    }
}

extern(C) int toStringFlowState(lua_State* L)
{
    FlowState* a = checkFlowState(L, 1);
    lua_pushstring(L, a.toString.toStringz);
    return 1;
}

/**
 * This function provides a way of making Dlanf FlowState structs from the Lua interface.
 *
 * Construction of a _FlowState struct from in Lua will accept:
 * -----------------
 * fs_id = _FlowState.new{p=1.0e5, T=300.0, velx=1000.0, vely=200.0, massf={spName=1.0}}
 * fs_id = _FlowState.new{p=1.0e7, T=300.0}
 * fs_id = _FlowState.new()
 * fs_id = _FlowState.new{}
 * -----------------
 * Missing velocity components are set to 0.0.
 * Missing mass fraction list is set to {1.0}.
 * For one-temperature gas models, single value for T is OK.
 * Temperature will always be accepted as an array.
 * For all other missing components, the values
 * are the defaults as given by the first constructor
 * in flowstate.d
 * The empty constructors forward through to PJ's
 * constructor that accepts a GasModel argument only.
 *
 * 2017-12-08
 * Note that the user is not expected to use this constructor
 * directly in their input script, however, we want it available
 * so that we can call the toJSONString function when writing
 * the config file.
 */
extern(C) int newFlowState(lua_State* L)
{
    auto managedGasModel = GlobalConfig.gmodel_master;
    if (managedGasModel is null) {
        string errMsg = `Error in call to FlowState:new.
It appears that you have not yet set the GasModel.
Be sure to call setGasModel(fname) before using a FlowState object.`;
        luaL_error(L, errMsg.toStringz);
    }
    FlowState* fs;
    int narg = lua_gettop(L);
    if (narg == 0) {
        // Make an empty FlowState in the context of a GasModel.
        fs = new FlowState(managedGasModel, GlobalConfig.turb_model.nturb);
    } else {
        // Presume that the supplied argument is a table representing the FlowState.
        fs = makeFlowStateFromTable(L, 1);
    }
    flowStateStore ~= fs;
    lua_pushinteger(L, flowStateStore.length-1);
    return 1;
}

FlowState* makeFlowStateFromTable(lua_State* L, int tblindx)
{
    string errMsg;
    auto managedGasModel = GlobalConfig.gmodel_master;
    if (managedGasModel is null) {
        errMsg = `Error in call to makeFlowStateFromTable.
It appears that you have not yet set the GasModel.
Be sure to call setGasModel(fname) before using a FlowState object.`;
        luaL_error(L, errMsg.toStringz);
    }
    if (!lua_istable(L, tblindx)) {
        errMsg = "Error in call to makeFlowStateFromTable. A table is expected as first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    // At this point we have a table at idx=1.
    //
    // If we have received a table that happens to be a Lua FlowState or CellData,
    // we should be able to trust the content of the table and so will not check fields.
    // Otherwise, we don't trust the table contents and we will check
    // that all fields in the table are valid.
    bool allFieldsAreValid = true;
    lua_getfield(L, tblindx, "myType");
    string myType = "";
    if (lua_isstring(L, -1)) {
        myType = to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    if (myType != "FlowState" && myType != "CellData") {
        // Proceed to check all entries.
        lua_pushnil(L);
        while (lua_next(L, tblindx) != 0) {
            string key = to!string(lua_tostring(L, -2));
            if (find(validFlowStateFields, key).empty) {
                allFieldsAreValid = false;
                errMsg ~= format("ERROR: '%s' is not a valid input in makeFlowStateFromTable\n", key);
            }
            lua_pop(L, 1);
        }
    }
    if (!allFieldsAreValid) {
        luaL_error(L, errMsg.toStringz);
    }
    // Now we are committed to using the first constructor
    // in class _FlowState. So we have to find at least
    // a pressure and temperature(s).
    errMsg = `Error in call to makeFlowStateFromTable.
A valid pressure value 'p' is not found in arguments.
The value should be a number.`;
    double p = getNumberFromTable(L, tblindx, "p", true, double.init, true, errMsg);

    errMsg = `Error in call to makeFlowStateFromTable.
A valid pressure value 'T' is not found in arguments.
The value should be a number.`;
    double T = getNumberFromTable(L, tblindx, "T", true, double.init, true, errMsg);

    // If we have a multi-temperature gas, then we must also find an entry for
    // T_modes.
    double[] T_modes;
    if (managedGasModel.n_modes >= 1) {
        lua_getfield(L, tblindx, "T_modes");
        if (lua_isnil(L, -1)) {
            errMsg = "Error in call to makeFlowStateFromTable.\n";
            errMsg ~= "T_modes is not set in table, but we have a multi-temperature gas.\n";
            throw new LuaInputException(errMsg);
        }

        // Next test for T_modes and see if it is a scalar or an array.
        if (lua_isnumber(L, -1)) {
            double Tval = lua_tonumber(L, -1);
            foreach (i; 0 .. managedGasModel.n_modes) { T_modes ~= Tval; }
        } else if (lua_istable(L, -1)) {
            getArrayOfDoubles(L, tblindx, "T_modes", T_modes);
            if ( T_modes.length != managedGasModel.n_modes ) {
                errMsg = "Error in call to makeFlowStateFromTable.";
                errMsg ~= "Length of T_modes vector does not match number of modes in gas model.";
                errMsg ~= format("T_modes.length= %d; n_modes= %d\n", T_modes.length, managedGasModel.n_modes);
                throw new LuaInputException(errMsg);
            }
        } else {
            foreach (i; 0 .. managedGasModel.n_modes) { T_modes ~= T; }
        }
        lua_pop(L, 1);
    }

    // Now everything else is optional. If it has been set, then we will
    // ensure that it can be retrieved correctly, or signal the user.

    // Values related to velocity.
    double velx = 0.0;
    double vely = 0.0;
    double velz = 0.0;
    string errMsgTmplt = "Error in call to makeFlowStateFromTable.\n";
    errMsgTmplt ~= "A valid value for '%s' is not found in arguments.\n";
    errMsgTmplt ~= "The value, if present, should be a number.";
    velx = getNumberFromTable(L, tblindx, "velx", false, 0.0, true, format(errMsgTmplt, "velx"));
    vely = getNumberFromTable(L, tblindx, "vely", false, 0.0, true, format(errMsgTmplt, "vely"));
    velz = getNumberFromTable(L, tblindx, "velz", false, 0.0, true, format(errMsgTmplt, "velz"));
    auto vel = Vector3(velx, vely, velz);

    // Values related to mass fractions.
    double[] massf;
    auto nsp = managedGasModel.n_species();
    massf.length = nsp;
    lua_getfield(L, tblindx, "massf");
    if (lua_isnil(L, -1)) {
        if (nsp == 1) {
            massf[0] = 1.0;
        } else {
            errMsg = "ERROR: in call to makeFlowStateFromTable.\n";
            errMsg ~= format("You are using a multi-component gas with n_species= %d\n", nsp);
            errMsg ~= "However, you have not set any mass fraction values.\n";
            throw new LuaInputException(errMsg);
        }
    } else if (lua_istable(L, -1)) {
        int massfIdx = lua_gettop(L);
        getSpeciesValsFromTable(L, managedGasModel, massfIdx, massf, "massf");
    } else {
        errMsg = "Error in call to makeFlowStateFromTable.\n";
        errMsg ~= "A field for mass fractions was found, but the contents are not valid.";
        errMsg ~= "The mass fraction should be given as a table of key-value pairs { speciesName=val }.";
        throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    // Value for quality
    double quality = getNumberFromTable(L, tblindx, "quality", false, 1.0, true, format(errMsgTmplt, "quality"));

    // Values for B (magnetic field)
    double Bx = 0.0;
    double By = 0.0;
    double Bz = 0.0;
    Bx = getNumberFromTable(L, tblindx, "Bx", false, 0.0, true, format(errMsgTmplt, "Bx"));
    By = getNumberFromTable(L, tblindx, "By", false, 0.0, true, format(errMsgTmplt, "By"));
    Bz = getNumberFromTable(L, tblindx, "Bz", false, 0.0, true, format(errMsgTmplt, "Bz"));
    auto B = Vector3(Bx, By, Bz);

    //Divergence of the magnetic field
    double divB = getNumberFromTable(L, tblindx, "divB", false, 0.0, true, format(errMsgTmplt, "divB"));
    //Divergence cleaning parameter psi for MHD
    double psi = getNumberFromTable(L, tblindx, "psi", false, 0.0, true, format(errMsgTmplt, "psi"));

    // Values related to turbulence modelling.
    double[] turb;
    lua_getfield(L, tblindx, "turb");
    if (lua_isnil(L, -1)) {
        auto tm = GlobalConfig.turb_model;
        turb.length = tm.nturb;
        foreach(it; 0 .. tm.nturb){
            string tvname = tm.primitive_variable_name(it);
            double tv = getNumberFromTable(L, tblindx, tvname, false, 0.0, true, format(errMsgTmplt, tvname));
            turb[it] = tv;
        }
    } else if (lua_istable(L, -1)) {
        // TODO: Consider making LUA always store tvariables by name
        lua_pop(L, 1); // get turb off the stack, getArrayOfDoubles will make its own copy
        getArrayOfDoubles(L, tblindx, "turb", turb);
    } else {
        lua_pop(L, 1);
        errMsg = "Error in call to makeFlowStateFromTable.\n";
        errMsg ~= "turb field not valid";
        throw new LuaInputException(errMsg);
    }
    double[2] turb0 = 0.0;
    foreach(it; 0 .. turb.length) turb0[it] = turb[it];
    double mu_t = getNumberFromTable(L, tblindx, "mu_t", false, 0.0, true, format(errMsgTmplt, "mu_t"));
    double k_t = getNumberFromTable(L, tblindx, "k_t", false, 0.0, true, format(errMsgTmplt, "k_t"));

    // Shock detector value.
    double S = getNumberFromTable(L, tblindx, "S", false, 0.0, true, format(errMsgTmplt, "S"));

    FlowState* fs = new FlowState(managedGasModel, p, T, T_modes, vel, turb0, massf, quality, B,
                                  psi, divB, mu_t, k_t, S);
    return fs;
} // end makeFlowStateFromTable()

/**
 * Provide a peek into the FlowState data as a Lua table.
 *
 * Basically, this gives the user a table to look at the values
 * in a FlowState in a read-only manner. (Well, in truth, the
 * table values can be changed, but they won't be reflected
 * in the FlowState object. This is consistent with the methods
 * of the FlowState object. Presently, there is no automatic
 * update if one fiddles with the gas properties in FlowState
 * object.
 */

string pushGasVar(string var)
{
    return `lua_pushnumber(L, fs.gas.` ~ var ~ `);
lua_setfield(L, tblIdx, "` ~ var ~`");`;
}

string pushGasVar(string var_in_D, string var_in_Lua)
{
    return `lua_pushnumber(L, fs.gas.` ~ var_in_D ~ `);
lua_setfield(L, tblIdx, "` ~ var_in_Lua ~`");`;
}

string pushGasVarArray(string var)
{
    return `lua_newtable(L);
foreach (i, val; fs.gas.` ~ var ~ `) {
    lua_pushnumber(L, val); lua_rawseti(L, -2,to!int(i+1));
}
lua_setfield(L, tblIdx, "` ~ var ~`");`;
}

string pushFSVar(string var)
{
return `lua_pushnumber(L, fs.` ~ var ~ `);
lua_setfield(L, tblIdx, "` ~ var ~`");`;
}

string pushFSVarArray(string var)
{
    return `lua_newtable(L);
foreach (i, val; fs.` ~ var ~ `) {
    lua_pushnumber(L, val); lua_rawseti(L, -2,to!int(i));
}
lua_setfield(L, tblIdx, "` ~ var ~`");`;
}

string pushFSVecVar(string var)
{
return `lua_pushnumber(L, fs.`~var~`.x);
lua_setfield(L, tblIdx, "`~var~`x");
lua_pushnumber(L, fs.`~var~`.y);
lua_setfield(L, tblIdx, "`~var~`y");
lua_pushnumber(L, fs.`~var~`.z);
lua_setfield(L, tblIdx, "`~var~`z");`;
}

/**
 * Push FlowState values to a table at TOS in lua_State.
 */
void pushFlowStateToTable(lua_State* L, int tblIdx, in FlowState fs, GasModel gmodel)
{
    mixin(pushGasVar("p"));
    mixin(pushGasVar("T", "T")); // now same in Lua and Dlang domains, 2017-12-04
    mixin(pushGasVar("u"));
    version(multi_T_gas) {
        mixin(pushGasVarArray("T_modes"));
        mixin(pushGasVarArray("u_modes"));
    }
    mixin(pushGasVar("quality"));
    version(multi_species_gas) {
        // -- massf as key-val table
        lua_newtable(L);
        foreach (isp, mf; fs.gas.massf) {
            lua_pushnumber(L, mf);
            lua_setfield(L, -2, toStringz(gmodel.species_name(isp)));
        }
        lua_setfield(L, tblIdx, "massf");
        // -- done setting massf
    }
    mixin(pushGasVar("a"));
    mixin(pushGasVar("rho"));
    mixin(pushGasVar("mu"));
    mixin(pushGasVar("k", "k"));
    version(multi_T_gas) {
        mixin(pushGasVarArray("k_modes"));
    }
    version(turbulence) {
        mixin(pushFSVarArray("turb"));
    }
    mixin(pushFSVar("mu_t"));
    mixin(pushFSVar("k_t"));
    mixin(pushFSVecVar("vel"));
    version(MHD) {
        mixin(pushFSVecVar("B"));
        mixin(pushFSVar("psi"));
        mixin(pushFSVar("divB"));
    }
}

/**
 * Gives the caller a table populated with FlowState values.
 *
 * Note that the table is flat, and that just a few GasState
 * variables have been unpacked. The fields in the returned table
 * form a superset of those that the user can set.
 */
extern(C) int toTable(lua_State* L)
{
    auto gmodel = GlobalConfig.gmodel_master;
    auto fs = checkFlowState(L, 1);
    lua_newtable(L); // anonymous table { }
    int tblIdx = lua_gettop(L);
    pushFlowStateToTable(L, tblIdx, *fs, gmodel);
    return 1;
}

string checkGasVar(string var)
{
    return `lua_getfield(L, 2, "`~var~`");
if ( !lua_isnil(L, -1) ) {
    fs.gas.`~var~` = luaL_checknumber(L, -1);
}
lua_pop(L, 1);`;
}

string checkGasVar(string var_in_D, string var_in_Lua)
{
    return `lua_getfield(L, 2, "`~var_in_Lua~`");
if ( !lua_isnil(L, -1) ) {
    fs.gas.`~var_in_D~` = luaL_checknumber(L, -1);
}
lua_pop(L, 1);`;
}

string checkGasVarArray(string var)
{
    return `lua_getfield(L, 2, "`~var~`");
if ( lua_istable(L, -1 ) ) {
    fs.gas.`~var~`.length = 0;
    getArrayOfDoubles(L, -2, "`~var~`", fs.gas.`~var~`);
}
lua_pop(L, 1);`;
}

string checkFSVarArray(string var)
{
    return `lua_getfield(L, 2, "`~var~`");
if ( lua_istable(L, -1 ) ) {
    double[] arr;
    getArrayOfDoubles(L, -2, "`~var~`", arr);
    foreach(i; 0 .. arr.length) fs.`~var~`[i] = arr[i];
}
lua_pop(L, 1);`;
}


string checkFSVar(string var)
{
    return `lua_getfield(L, 2, "`~var~`");
if ( !lua_isnil(L, -1) ) {
    fs.`~var~` = luaL_checknumber(L, -1);
}
lua_pop(L, 1);`;
}

extern(C) int fromTable(lua_State* L)
{
    auto managedGasModel = GlobalConfig.gmodel_master;
    auto fs = checkFlowState(L, 1);
    if ( !lua_istable(L, 2) ) {
        return 0;
    }
    // Look for gas variables: "p" and "quality"
    mixin(checkGasVar("p"));
    mixin(checkGasVar("quality"));
    mixin(checkGasVar("T", "T")); // now same name in Lua domain
    version(multi_species_gas) {
        // Look for a table with mass fraction info
        lua_getfield(L, 2, "massf");
        if ( lua_istable(L, -1) ) {
            int massfIdx = lua_gettop(L);
            getSpeciesValsFromTable(L, managedGasModel, massfIdx, fs.gas.massf, "massf");
        }
        lua_pop(L, 1);
    }
    version(multi_T_gas) {
        // Look for an array of internal temperatures.
        mixin(checkGasVarArray("T_modes"));
        if ( fs.gas.T_modes.length != GlobalConfig.gmodel_master.n_modes ) {
            string errMsg = "The temperature array ('T_modes') did not contain"~
                " the correct number of entries.\n";
            errMsg ~= format("T_modes.length= %d; n_modes= %d\n", fs.gas.T_modes.length,
                             GlobalConfig.gmodel_master.n_modes);
            luaL_error(L, errMsg.toStringz);
        }
    }
    // Let's try to find rho and u so that the pT thermo call
    // has a good set of starting values.
    mixin(checkGasVar("rho"));
    mixin(checkGasVar("u"));

    // We should call equation of state to make sure gas state is consistent.
    GlobalConfig.gmodel_master.update_thermo_from_pT(fs.gas);
    GlobalConfig.gmodel_master.update_sound_speed(fs.gas);
    GlobalConfig.gmodel_master.update_trans_coeffs(fs.gas);

    // Look for velocity components: "velx", "vely", "velz"
    lua_getfield(L, 2, "velx");
    if ( !lua_isnil(L, -1 ) ) {
        fs.vel.x = luaL_checknumber(L, -1);
    }
    lua_pop(L, 1);
    lua_getfield(L, 2, "vely");
    if ( !lua_isnil(L, -1 ) ) {
        fs.vel.y = luaL_checknumber(L, -1);
    }
    lua_pop(L, 1);
    lua_getfield(L, 2, "velz");
    if ( !lua_isnil(L, -1 ) ) {
        fs.vel.z = luaL_checknumber(L, -1);
    }
    lua_pop(L, 1);
    version(MHD) {
        // Look for B components: "Bx", "By", "Bz"
        lua_getfield(L, 2, "Bx");
        if ( !lua_isnil(L, -1 ) ) {
            fs.B.x = luaL_checknumber(L, -1);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "By");
        if ( !lua_isnil(L, -1 ) ) {
            fs.B.y = luaL_checknumber(L, -1);
        }
        lua_pop(L, 1);
        lua_getfield(L, 2, "Bz");
        if ( !lua_isnil(L, -1 ) ) {
            fs.B.z = luaL_checknumber(L, -1);
        }
        lua_pop(L, 1);

        // Look for divergence cleaning parameter psi
        mixin(checkFSVar("psi"));
        mixin(checkFSVar("divB"));
    }

    // Now look turbulence quantities
    version(turbulence) {
        mixin(checkFSVarArray("turb"));
    }
    mixin(checkFSVar("mu_t"));
    mixin(checkFSVar("k_t"));
    return 0;
}

extern(C) int toJSONString(lua_State* L)
{
    FlowState* fs = checkFlowState(L, 1);
    lua_pushstring(L, fs.toJSONString().toStringz);
    return 1;
}

void registerFlowState(lua_State* L)
{
    // Register methods for use but put them into a table called _FlowState.
    lua_newtable(L);
    lua_pushcfunction(L, &newFlowState);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringFlowState);
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &toTable);
    lua_setfield(L, -2, "toTable");
    lua_pushcfunction(L, &fromTable);
    lua_setfield(L, -2, "fromTable");
    lua_pushcfunction(L, &toJSONString);
    lua_setfield(L, -2, "toJSONString");
    // Make table visible
    lua_setglobal(L, "_FlowState".toStringz);
}
