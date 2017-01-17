/**
 * luagas_model.d
 * Lua interface to the gas model.
 *
 * Author: Rowan G.
 * Version: Initial cut.
 */

module gas.luagas_model;

import std.math;
import std.algorithm;
import std.stdio;
import std.conv;
import std.string;
import util.lua;
import util.lua_service;
import kinetics.luareaction_mechanism;
import kinetics.luachemistry_update;

import gas.gas_model;
import gas.physical_constants;


// name for GasModel class in Lua scripts
immutable string GasModelMT = "GasModel";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(GasModel)[] GasModelStore;

GasModel checkGasModel(lua_State* L, int index)
{
    return checkObj!(GasModel, GasModelMT)(L, index);
}

/**
 * This function implements the constructor for the GasModel
 * from the Lua interface.
 *
 * Construction of GasModel object in Lua will accept
 * a filename.
 * -------------
 * gmodel = GasModel:new{fname}
 *--------------
 */
extern(C) int newGasModel(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument this

    lua_rawgeti(L, 1, 1);
    string fname = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    GasModel gm = init_gas_model(fname);
    GasModelStore ~= pushObj!(GasModel, GasModelMT)(L, gm);
    return 1;
}

extern(C) int nSpecies(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_pushnumber(L, gm.n_species);
    return 1;
}

extern(C) int nModes(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_pushnumber(L, gm.n_modes);
    return 1;
}

extern(C) int molMasses(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_newtable(L);
    foreach ( int isp, mMass; gm.mol_masses ) {
	lua_pushnumber(L, mMass);
	lua_rawseti(L, -2, isp+1);
    }
    return 1;
}

extern(C) int speciesName(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    int i = to!int(lua_tointeger(L, 2));
    // Use "i-1" since D indexes from 0
    lua_pushstring(L, toStringz(gm.species_name(i-1)));
    return 1;
}

extern(C) int thermoPT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    gm.update_thermo_from_pT(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoRHOE(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    gm.update_thermo_from_rhoe(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoRHOT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    gm.update_thermo_from_rhoT(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoRHOP(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    gm.update_thermo_from_rhop(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoPS(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto s = lua_tonumber(L, 3);
    gm.update_thermo_from_ps(Q, s);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoHS(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto h = lua_tonumber(L, 3);
    auto s = lua_tonumber(L, 4);
    gm.update_thermo_from_hs(Q, h, s);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int soundSpeed(lua_State* L)
{
   auto gm = checkGasModel(L, 1);
   auto Q = new GasState(gm.n_species, gm.n_modes);
   getGasStateFromTable(L, gm, 2, Q); 
   gm.update_sound_speed(Q);
   setGasStateInTable(L, gm, 2, Q);
   return 0;
}

extern(C) int transCoeffs(lua_State* L)
{
   auto gm = checkGasModel(L, 1);
   auto Q = new GasState(gm.n_species, gm.n_modes);
   getGasStateFromTable(L, gm, 2, Q); 
   gm.update_trans_coeffs(Q);
   setGasStateInTable(L, gm, 2, Q);
   return 0;
}


// We don't wrap dudT_const_v as we do Cv in its place.
// We don't wrap dhdT_const_p as we do Cp in its place.

extern(C) int dpdrhoConstT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto dpdrho = gm.dpdrho_const_T(Q);
    lua_pushnumber(L, dpdrho);
    return 1;
}

// We don't wrap gas_constant as we do R in its place.

extern(C) int intEnergy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    int narg = lua_gettop(L);
    auto e = gm.internal_energy(Q);
    lua_pushnumber(L, e);
    return 1;
}

extern(C) int enthalpy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    int narg = lua_gettop(L);
    double h;
    if ( narg >= 3 ) { // Call species-specific version
	int isp = luaL_checkint(L, 3);
	h = gm.enthalpy(Q, isp);
    }
    else { // Call total gas mix version
	h = gm.enthalpy(Q);
    }
    lua_pushnumber(L, h);
    return 1;
}

extern(C) int entropy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    int narg = lua_gettop(L);
    double s;
    if ( narg >= 3 ) { // Call species-specific version
	int isp = luaL_checkint(L, 3);
	s = gm.entropy(Q, isp);
    }
    else {
	s = gm.entropy(Q);
    }
    lua_pushnumber(L, s);
    return 1;
}

extern(C) int Cv(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto Cv = gm.Cv(Q);
    lua_pushnumber(L, Cv);
    return 1;
}

extern(C) int Cp(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto Cp = gm.Cp(Q);
    lua_pushnumber(L, Cp);
    return 1;
}

extern(C) int R(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto R = gm.R(Q);
    lua_pushnumber(L, R);
    return 1;
}

extern(C) int gamma(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto gamma = gm.gamma(Q);
    lua_pushnumber(L, gamma);
    return 1;
}

extern(C) int molMass(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    auto molMass = gm.molecular_mass(Q);
    lua_pushnumber(L, molMass);
    return 1;
}

extern(C) int massf2molef(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    double[] molef; molef.length = gm.n_species;
    gm.massf2molef(Q, molef);
    // Place molef in an array and leave at at
    // top-of-stack as a return to the caller.
    lua_newtable(L);
    foreach ( int i, mf; molef ) {
	lua_pushnumber(L, mf);
	lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

extern(C) int molef2massf(lua_State* L)
{
    int narg = lua_gettop(L);
    auto gm = checkGasModel(L, 1);
    double[] molef;
    molef.length = gm.n_species;
    if ( lua_istable(L, 2) ) {
	getSpeciesValsFromTable(L, gm, 2, molef, "molef");
    }
    else {
	string errMsg = "Error in call to molef2massf():\n";
	errMsg ~= "The value for 'molef' is not a table of key-val pairs.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    auto Q = new GasState(gm.n_species, gm.n_modes);
    if ( narg == 3 ) {
	getGasStateFromTable(L, gm, 3, Q);
    }
    gm.molef2massf(molef, Q);
    // Update table with new mass fractions
    if ( narg == 3 ) {
	setGasStateInTable(L, gm, 3, Q);
    }
    // and return a table to the caller.
    lua_newtable(L);
    foreach ( int i, mf; Q.massf ) {
	lua_pushnumber(L, mf);
	lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

extern(C) int massf2conc(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    double[] conc; conc.length = gm.n_species;
    gm.massf2conc(Q, conc);
    // Place conc in an array and leave at at
    // top-of-stack as a return to the caller.
    lua_newtable(L);
    foreach ( int i, c; conc ) {
	lua_pushnumber(L, c);
	lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

extern(C) int conc2massf(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    double[] conc;
    conc.length = gm.n_species();
    if ( lua_istable(L, 2) ) {
	getSpeciesValsFromTable(L, gm, 2, conc, "conc");
    }
    else {
	string errMsg = "Error in call to conc2massf():\n";
	errMsg ~= "The value for 'conc' is not a table of key-val pairs.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 3, Q);
    gm.conc2massf(conc, Q);
    // Update table with new mass fractions
    setGasStateInTable(L, gm, 3, Q);
    // and return a table to the caller.
    lua_newtable(L);
    foreach ( int i, mf; Q.massf ) {
	lua_pushnumber(L, mf);
	lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

void getSpeciesValsFromTable(lua_State* L, GasModel gm, int idx,
			     double[] vals, string tabName)
{
    // 1. Check all keys are valid species names.
    lua_pushnil(L);
    while ( lua_next(L, idx) != 0 ) {
	string key = to!string(lua_tostring(L, -2));
	auto isp = gm.species_index(key);
	if ( isp == -1 ) {
	    string errMsg = format("Species name used in %s table does not exist: %s\n", tabName, key);
	    lua_pop(L, 1);
	    throw new LuaInputException(errMsg);
	}
	lua_pop(L, 1);
    }
    // 2. Now set all values to 0.0
    //    (then we'll correct that in the next step)
    //    [Or to 1.0 in the case of a n_species = 1]
    if ( gm.n_species() == 1 )
	vals[0] = 1.0;
    else 
	vals[] = 0.0;
    // 3. Now find those values that we have explicitly set
    foreach ( isp; 0 .. gm.n_species() ) {
	lua_getfield(L, -1, toStringz(gm.species_name(isp)));
	if ( lua_isnumber(L, -1) ) {
	    vals[isp] = lua_tonumber(L, -1);
	}
	else if ( lua_isnil(L, -1) ) {
	    vals[isp] = 0.0;
	}
	else {
	    string errMsg = format("The value for species '%s' in the %s table is not a number.\n", gm.species_name(isp), tabName);
	    lua_pop(L, 1);
	    throw new LuaInputException(errMsg);
	}
	lua_pop(L, 1);
    }
}

void checkAndScaleMassFractions(double[] massf, double tol)
{
    auto massfSum = sum(massf);
    if ( fabs(massfSum - 1.0) > tol ) {
	string errMsg = "The given mass fraction values do not sum to 1.0.\n";
	errMsg ~= format("The sum value is: %e\n", massfSum);
	errMsg ~= format("The error in the sum is larger than tolerance= %e \n", tol);
	throw new Error(errMsg);
    }
    // otherwise, do scaling based on normalising by sum
    massf[] /= massfSum;
}

extern(C) int createTableForGasState(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    // For the special case of n_species = 1, it's quite likely
    // that the user might never consider the massf array.
    // In which case, we'll just set that value to 1.0.
    // For all other cases, we take the default action of
    // setting all mass fractions to 0.0.
    // For a multi-component gas, we expect the user to 
    // take care with setting mass fraction values.
    if ( gm.n_species == 1 ) {
	Q.massf[0] = 1.0;
    }
    else {
	Q.massf[] = 0.0;
    }
    lua_newtable(L);
    int idx = lua_gettop(L);
    setGasStateInTable(L, gm, idx, Q);
    // It's convenient in the Lua code to put a reference
    // to the gas model in the GasState table.
    pushObj!(GasModel, GasModelMT)(L, gm);
    lua_setfield(L, idx, "gasmodel");
    return 1;
}

extern(C) int newTableForGasState(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument this
    lua_rawgeti(L, 1, 1);
    auto gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    // For the special case of n_species = 1, it's quite likely
    // that the user might never consider the massf array.
    // In which case, we'll just set that value to 1.0.
    // For all other cases, we take the default action of
    // setting all mass fractions to 0.0.
    // For a multi-component gas, we expect the user to 
    // take care with setting mass fraction values.
    if ( gm.n_species == 1 ) {
	Q.massf[0] = 1.0;
    }
    else {
	Q.massf[] = 0.0;
    }
    lua_newtable(L);
    int idx = lua_gettop(L);
    setGasStateInTable(L, gm, idx, Q);
    // It's convenient in the Lua code to put a reference
    // to the gas model in the GasState table.
    pushObj!(GasModel, GasModelMT)(L, gm);
    lua_setfield(L, idx, "gasmodel");
    return 1;
}

void getGasStateFromTable(lua_State* L, GasModel gm, int idx, GasState Q)
{
    lua_getfield(L, idx, "rho");
    if ( lua_isnumber(L, -1) ) {
	Q.rho = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'rho' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "p");
    if ( lua_isnumber(L, -1) ) {
	Q.p = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'p' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);
    
    lua_getfield(L, idx, "T");
    if ( lua_isnumber(L, -1) ) {
	Q.Ttr = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'T' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "p_e");
    if ( lua_isnumber(L, -1) ) {
	Q.p_e = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'p_e' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);
    
    lua_getfield(L, idx, "a");
    if ( lua_isnumber(L, -1) ) {
	Q.a = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'a' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);
    
    lua_getfield(L, idx, "u");
    if ( lua_isnumber(L, -1) ) {
	Q.u = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'u' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "e_modes");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.e_modes.length ) {
	    string errMsg = format("Wrong number of internal energy values in GasState table.\n Expected: %d, Given: %d\n", Q.e_modes.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i);
	    if ( lua_isnumber(L, -1) ) {
		Q.e_modes[i-1] = lua_tonumber(L, -1);
	    }
	    else {
		string errMsg = format("The value for 'e_modes[%d]' is not a number.\n", i);
		lua_pop(L, 1);
		throw new Error(errMsg);
	    }
	    lua_pop(L, 1);
	}
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'e_modes' is not an array of numbers.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "T_modes");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.T_modes.length ) {
	    string errMsg = format("Wrong number of internal temperature values in GasState table.\n Expected: %d, Given: %d\n", Q.T_modes.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i);
	    if ( lua_isnumber(L, -1) ) {
		Q.T_modes[i-1] = lua_tonumber(L, -1);
	    }
	    else {
		string errMsg = format("The value for 'T_modes[%d]' is not a number.\n", i);
		lua_pop(L, 1);
		throw new Error(errMsg);
	    }
	    lua_pop(L, 1);
	}
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'T' is not an array of numbers.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "mu");
    if ( lua_isnumber(L, -1) ) {
	Q.mu = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'mu' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "k");
    if ( lua_isnumber(L, -1) ) {
	Q.k = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'k' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "k_modes");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.k_modes.length ) {
	    string errMsg = format("Wrong number of internal thermal conductivity values in GasState table.\n Expected: %d, Given: %d\n", Q.k_modes.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i);
	    if ( lua_isnumber(L, -1) ) {
		Q.k_modes[i-1] = lua_tonumber(L, -1);
	    }
	    else {
		string errMsg = format("The value for 'k_modes[%d]' is not a number.\n", i);
		lua_pop(L, 1);
		throw new Error(errMsg);
	    }
	    lua_pop(L, 1);
	}
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'k_modes' is not an array of numbers.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "sigma");
    if ( lua_isnumber(L, -1) ) {
	Q.sigma = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'sigma' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "massf");
    if ( lua_istable(L, -1) ) {
	int massfIdx = lua_gettop(L);
	getSpeciesValsFromTable(L, gm, massfIdx, Q.massf, "massf");
	checkAndScaleMassFractions(Q.massf, MASSF_ERROR_TOL);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'massf' is not a table of key-val pairs.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "quality");
    if ( lua_isnumber(L, -1) ) {
	Q.quality = lua_tonumber(L, -1);
    }
    else if ( lua_isnil(L, -1) ) {
	// leave untouched
    }
    else {
	string errMsg = "The value for 'quality' is not a number.\n";
	lua_pop(L, 1);
	throw new Error(errMsg);
    }
    lua_pop(L, 1);
}

void setGasStateInTable(lua_State* L, GasModel gm, int idx, GasState Q)
{
    lua_pushnumber(L, Q.rho);
    lua_setfield(L, idx, "rho");

    lua_pushnumber(L, Q.p);
    lua_setfield(L, idx, "p");

    lua_pushnumber(L, Q.Ttr);
    lua_setfield(L, idx, "T");

    lua_pushnumber(L, Q.p_e);
    lua_setfield(L, idx, "p_e");

    lua_pushnumber(L, Q.a);
    lua_setfield(L, idx, "a");

    lua_pushnumber(L, Q.u);
    lua_setfield(L, idx, "u");

    lua_newtable(L);
    foreach ( int i, e; Q.e_modes ) {
	lua_pushnumber(L, e); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "e_modes");

    lua_newtable(L);
    foreach ( int i, T; Q.T_modes ) {
	lua_pushnumber(L, T); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "T_modes");
    
    lua_pushnumber(L, Q.mu);
    lua_setfield(L, idx, "mu");
    
    lua_pushnumber(L, Q.k);
    lua_setfield(L, idx, "k");

    lua_newtable(L);
    foreach ( int i, k; Q.k_modes ) {
	lua_pushnumber(L, k); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "k_modes");

    lua_pushnumber(L, Q.sigma);
    lua_setfield(L, idx, "sigma");

    lua_newtable(L);
    foreach ( int i, mf; Q.massf ) {
	lua_pushnumber(L, mf);
	lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    lua_setfield(L, idx, "massf");

    lua_pushnumber(L, Q.quality);
    lua_setfield(L, idx, "quality");
}

extern(C) int printValues(lua_State* L)
{
    lua_getfield(L, 1, "gasmodel");
    auto gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 1, Q);
    writeln(Q.toString);
    return 0;
}

extern(C) int copyValues(lua_State* L)
{
    lua_getfield(L, 1, "gasmodel");
    auto gm = checkGasModel(L, -1);
    lua_pop(L, 1);

    // Source GasState table
    auto Q0 = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 1, Q0);

    // Target GasState table
    // Assume table is available at index 2
    setGasStateInTable(L, gm, 2, Q0);

    return 0;
}

void registerGasModel(lua_State* L, int tblIdx)
{
    luaL_newmetatable(L, GasModelMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newGasModel);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(GasModel, GasModelMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &nSpecies);
    lua_setfield(L, -2, "nSpecies");
    lua_pushcfunction(L, &nModes);
    lua_setfield(L, -2, "nModes");
    lua_pushcfunction(L, &molMasses);
    lua_setfield(L, -2, "molMasses");
    lua_pushcfunction(L, &speciesName);
    lua_setfield(L, -2, "speciesName");
    lua_pushcfunction(L, &createTableForGasState);
    lua_setfield(L, -2, "createGasState");
    lua_pushcfunction(L, &thermoPT);
    lua_setfield(L, -2, "updateThermoFromPT");
    lua_pushcfunction(L, &thermoRHOE);
    lua_setfield(L, -2, "updateThermoFromRHOE");
    lua_pushcfunction(L, &thermoRHOT);
    lua_setfield(L, -2, "updateThermoFromRHOT");
    lua_pushcfunction(L, &thermoPT);
    lua_setfield(L, -2, "updateThermoFromRHOP");
    lua_pushcfunction(L, &thermoPS);
    lua_setfield(L, -2, "updateThermoFromPS");
    lua_pushcfunction(L, &thermoHS);
    lua_setfield(L, -2, "updateThermoFromHS");
    lua_pushcfunction(L, &soundSpeed);
    lua_setfield(L, -2, "updateSoundSpeed");
    lua_pushcfunction(L, &transCoeffs);
    lua_setfield(L, -2, "updateTransCoeffs");
    lua_pushcfunction(L, &dpdrhoConstT);
    lua_setfield(L, -2, "dpdrhoConstT");
    lua_pushcfunction(L, &intEnergy);
    lua_setfield(L, -2, "intEnergy");
    lua_pushcfunction(L, &enthalpy);
    lua_setfield(L, -2, "enthalpy");
    lua_pushcfunction(L, &entropy);
    lua_setfield(L, -2, "entropy");
    lua_pushcfunction(L, &Cv);
    lua_setfield(L, -2, "Cv");
    lua_pushcfunction(L, &Cv);
    lua_setfield(L, -2, "dudTConstV");
    lua_pushcfunction(L, &Cp);
    lua_setfield(L, -2, "Cp");
    lua_pushcfunction(L, &Cp);
    lua_setfield(L, -2, "dhdTConstP");
    lua_pushcfunction(L, &R);
    lua_setfield(L, -2, "R");
    lua_pushcfunction(L, &R);
    lua_setfield(L, -2, "gasConstant");
    lua_pushcfunction(L, &gamma);
    lua_setfield(L, -2, "gamma");
    lua_pushcfunction(L, &molMass);
    lua_setfield(L, -2, "molMass");
    lua_pushcfunction(L, &massf2molef);
    lua_setfield(L, -2, "massf2molef");
    lua_pushcfunction(L, &molef2massf);
    lua_setfield(L, -2, "molef2massf");
    lua_pushcfunction(L, &massf2conc);
    lua_setfield(L, -2, "massf2conc");
    lua_pushcfunction(L, &conc2massf);
    lua_setfield(L, -2, "conc2massf");

    // Make class visible
    lua_setfield(L, tblIdx, GasModelMT.toStringz);

    // Make initialisation of GasState table look like a class constructor
    luaL_newmetatable(L, "GasState");

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newTableForGasState);
    lua_setfield(L, -2, "new");

    // Make GasState constructor visible
    lua_setfield(L, tblIdx, "GasState");
    
    // Set a global function to print values in GasState
    lua_pushcfunction(L, &printValues);
    lua_setglobal(L, "printValues");
    // Set a global function to copy GasState values from one table to another.
    lua_pushcfunction(L, &copyValues);
    lua_setglobal(L, "copyValues");

    // Set some of the physical constants as global
    lua_pushnumber(L, P_atm);
    lua_setglobal(L, "P_atm");
    lua_pushnumber(L, R_universal);
    lua_setglobal(L, "R_universal");

}

version(gas_calc) {
    int main(string[] args) {
	if ( args.length < 2 ) {
	    writeln("ERROR: Wrong number of arguments.");
	    writeln("");
	    writeln("Usage: ");
	    writeln("  gas-calc inputFile (+ script args)");
	    return 1;
	}
	auto L = luaL_newstate();
	luaL_openlibs(L);
	registerGasModel(L, LUA_GLOBALSINDEX);
	registerReactionMechanism(L, LUA_GLOBALSINDEX);
	registerChemistryUpdate(L, LUA_GLOBALSINDEX);
	// Pass on command line args to user's scripts.
	lua_newtable(L);
	int argIdx = 1;
	foreach (i; 2 .. args.length) {
	    lua_pushstring(L, args[i].toStringz);
	    lua_rawseti(L, -2, argIdx);
	    argIdx++;
	}
	lua_setglobal(L, "arg");

	if ( luaL_dofile(L, toStringz(args[1])) != 0 ) {
	    writeln(to!string(lua_tostring(L, -1)));
	    return 1;
	}
	return 0;
    }
}

