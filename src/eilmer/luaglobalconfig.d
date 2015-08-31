// luaglobalconfig.d
// Lua access to the GlobalConfig class data, for use in the preparation script.
//
// Peter J. and Rowan G.
// 2015-03-02: First code adapted from the other lua wrapper modules.

module luaglobalconfig;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;

import gas;
import fvcore;
import globalconfig;

// -------------------------------------------------------------------------------
// Set GlobalConfig fields from a table.

extern(C) int configSetFromTable(lua_State* L)
{
    if (!lua_istable(L, 1)) return 0; // nothing to do
    //
    // Look for fields that may be present.
    lua_getfield(L, 1, "base_file_name");
    if (!lua_isnil(L, -1)) GlobalConfig.base_file_name = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "title");
    if (!lua_isnil(L, -1)) GlobalConfig.title = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "gas_model_file");
    if (!lua_isnil(L, -1)) GlobalConfig.gas_model_file = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "nBlocks");
    if (!lua_isnil(L, -1)) GlobalConfig.nBlocks = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "dimensions");
    if (!lua_isnil(L, -1)) GlobalConfig.dimensions = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "axisymmetric");
    if (!lua_isnil(L, -1)) GlobalConfig.axisymmetric = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    
    lua_getfield(L, 1, "MHD");
    if (!lua_isnil(L, -1)) GlobalConfig.MHD = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "gasdynamic_update_scheme");
    if (!lua_isnil(L, -1)) {
	string name = to!string(luaL_checkstring(L, -1));
	GlobalConfig.gasdynamic_update_scheme = update_scheme_from_name(name);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "separate_update_for_viscous_terms");
    if (!lua_isnil(L, -1)) 
	GlobalConfig.separate_update_for_viscous_terms = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "separate_update_for_k_omega_source");
    if (!lua_isnil(L, -1)) 
	GlobalConfig.separate_update_for_k_omega_source = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "apply_bcs_in_parallel");
    if (!lua_isnil(L, -1)) 
	GlobalConfig.apply_bcs_in_parallel = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "adjust_invalid_cell_data");
    if (!lua_isnil(L, -1))
	GlobalConfig.adjust_invalid_cell_data = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "max_invalid_cells");
    if (!lua_isnil(L, -1)) GlobalConfig.max_invalid_cells = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "dt_reduction_factor");
    if (!lua_isnil(L, -1))
	GlobalConfig.dt_reduction_factor = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "interpolation_order");
    if (!lua_isnil(L, -1)) GlobalConfig.interpolation_order = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "thermo_interpolator");
    if (!lua_isnil(L, -1)) {
	string name = to!string(luaL_checkstring(L, -1));
	GlobalConfig.thermo_interpolator = thermo_interpolator_from_name(name);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "apply_limiter");
    if (!lua_isnil(L, -1)) GlobalConfig.apply_limiter = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "extrema_clipping");
    if (!lua_isnil(L, -1)) GlobalConfig.extrema_clipping = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "interpolate_in_local_frame");
    if (!lua_isnil(L, -1))
	GlobalConfig.interpolate_in_local_frame = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "flux_calculator");
    if (!lua_isnil(L, -1)) {
	string name = to!string(luaL_checkstring(L, -1));
	GlobalConfig.flux_calculator = flux_calculator_from_name(name);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "shear_tolerance");
    if (!lua_isnil(L, -1))
	GlobalConfig.shear_tolerance = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "M_inf");
    if (!lua_isnil(L, -1)) GlobalConfig.M_inf = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "compression_tolerance");
    if (!lua_isnil(L, -1))
	GlobalConfig.compression_tolerance = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "moving_grid");
    if (!lua_isnil(L, -1)) GlobalConfig.moving_grid = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "write_vertex_velocities");
    if (!lua_isnil(L, -1))
	GlobalConfig.write_vertex_velocities = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "viscous");
    if (!lua_isnil(L, -1)) GlobalConfig.viscous = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "viscous_factor");
    if (!lua_isnil(L, -1))
	GlobalConfig.viscous_factor = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "viscous_factor_increment");
    if (!lua_isnil(L, -1))
	GlobalConfig.viscous_factor_increment = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "viscous_delay");
    if (!lua_isnil(L, -1))
	GlobalConfig.viscous_delay = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "viscous_signal_factor");
    if (!lua_isnil(L, -1))
	GlobalConfig.viscous_signal_factor = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "diffusion");
    if (!lua_isnil(L, -1)) GlobalConfig.diffusion = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    // [TODO] other diffusion control parameters


    lua_getfield(L, 1, "turbulence_model");
    if (!lua_isnil(L, -1)) {
	string name = to!string(luaL_checkstring(L, -1));
	GlobalConfig.turbulence_model = turbulence_model_from_name(name);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "turbulence_prandtl_number");
    if (!lua_isnil(L, -1))
	GlobalConfig.turbulence_prandtl_number = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "turbulence_schmidt_number");
    if (!lua_isnil(L, -1))
	GlobalConfig.turbulence_schmidt_number = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "max_mu_t_factor");
    if (!lua_isnil(L, -1))
	GlobalConfig.max_mu_t_factor = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "transient_mu_t_factor");
    if (!lua_isnil(L, -1))
	GlobalConfig.transient_mu_t_factor = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "reacting");
    if (!lua_isnil(L, -1)) GlobalConfig.reacting = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "reactions_file");
    if (!lua_isnil(L, -1)) GlobalConfig.reactions_file = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);


    lua_getfield(L, 1, "max_step");
    if (!lua_isnil(L, -1)) GlobalConfig.max_step = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "halt_now");
    if (!lua_isnil(L, -1)) GlobalConfig.halt_now = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "print_count");
    if (!lua_isnil(L, -1)) GlobalConfig.print_count = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "control_count");
    if (!lua_isnil(L, -1)) GlobalConfig.control_count = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "verbosity_level");
    if (!lua_isnil(L, -1)) GlobalConfig.verbosity_level = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "max_time");
    if (!lua_isnil(L, -1)) GlobalConfig.max_time = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "dt_init");
    if (!lua_isnil(L, -1)) GlobalConfig.dt_init = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "dt_max");
    if (!lua_isnil(L, -1)) GlobalConfig.dt_max = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "cfl_value");
    if (!lua_isnil(L, -1)) GlobalConfig.cfl_value = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "stringent_cfl");
    if (!lua_isnil(L, -1)) GlobalConfig.stringent_cfl = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "cfl_count");
    if (!lua_isnil(L, -1)) GlobalConfig.cfl_count = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "fixed_time_step");
    if (!lua_isnil(L, -1)) GlobalConfig.fixed_time_step = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "write_at_step");
    if (!lua_isnil(L, -1)) GlobalConfig.write_at_step = luaL_checkint(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, 1, "dt_plot");
    if (!lua_isnil(L, -1)) GlobalConfig.dt_plot = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "dt_history");
    if (!lua_isnil(L, -1)) GlobalConfig.dt_history = to!double(luaL_checknumber(L, -1));
    lua_pop(L, 1);
    //
    lua_getfield(L, 1, "udf_source_terms_file");
    if (!lua_isnil(L, -1)) GlobalConfig.udf_source_terms_file = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "udf_source_terms");
    if (!lua_isnil(L, -1)) GlobalConfig.udf_source_terms = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "block_marching");
    if (!lua_isnil(L, -1)) GlobalConfig.block_marching = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "nib");
    if (!lua_isnil(L, -1)) GlobalConfig.nib = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "njb");
    if (!lua_isnil(L, -1)) GlobalConfig.njb = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "nkb");
    if (!lua_isnil(L, -1)) GlobalConfig.nkb = luaL_checkint(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 1, "propagate_inflow_data");
    if (!lua_isnil(L, -1)) GlobalConfig.propagate_inflow_data = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);

    lua_getfield(L, 1, "udf_solid_source_terms_file");
    if (!lua_isnil(L, -1)) GlobalConfig.udfSolidSourceTermsFile = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "udf_solid_source_terms");
    if (!lua_isnil(L, -1)) GlobalConfig.udfSolidSourceTerms = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);


    return 0;
} // end configSetFromTable()

// Get GlobalConfig fields by their string name.
extern(C) int configGet(lua_State* L)
{
    string fieldName = to!string(luaL_checkstring(L, 1));

    switch (fieldName) {
    case "base_file_name":
	lua_pushstring(L, GlobalConfig.base_file_name.toStringz);
	break;
    case "title":
	lua_pushstring(L, GlobalConfig.title.toStringz);
	break;
    case "gas_model_file":
	lua_pushstring(L, GlobalConfig.gas_model_file.toStringz);
	break;

    case "nBlocks":
	lua_pushnumber(L, GlobalConfig.nBlocks);
	break;
    case "dimensions":
	lua_pushnumber(L, GlobalConfig.dimensions);
	break;
    case "axisymmetric":
	lua_pushboolean(L, GlobalConfig.axisymmetric);
	break;
	
	case "MHD":
	lua_pushboolean(L, GlobalConfig.MHD);
	break;

    case "gasdynamic_update_scheme":
	string name = gasdynamic_update_scheme_name(GlobalConfig.gasdynamic_update_scheme);
	lua_pushstring(L, name.toStringz);
	break;
    case "separate_update_for_viscous_terms":
	lua_pushboolean(L, GlobalConfig.separate_update_for_viscous_terms);
	break;
    case "separate_update_for_k_omega_source":
	lua_pushboolean(L, GlobalConfig.separate_update_for_k_omega_source);
	break;
    case "apply_bcs_in_parallel":
	lua_pushboolean(L, GlobalConfig.apply_bcs_in_parallel);
	break;
    case "adjust_invalid_cell_data":
	lua_pushboolean(L, GlobalConfig.adjust_invalid_cell_data);
	break;
    case "max_invalid_cells":
	lua_pushnumber(L, GlobalConfig.max_invalid_cells);
	break;
    case "dt_reduction_factor":
	lua_pushnumber(L, GlobalConfig.dt_reduction_factor);
	break;

    case "interpolation_order":
	lua_pushnumber(L, GlobalConfig.interpolation_order);
	break;
    case "thermo_interpolator":
	string name = thermo_interpolator_name(GlobalConfig.thermo_interpolator);
	lua_pushstring(L, name.toStringz);
	break;
    case "apply_limiter":
	lua_pushboolean(L, GlobalConfig.apply_limiter);
	break;
    case "extrema_clipping":
	lua_pushboolean(L, GlobalConfig.extrema_clipping);
	break;
    case "interpolate_in_local_frame":
	lua_pushboolean(L, GlobalConfig.interpolate_in_local_frame);
	break;

    case "flux_calculator":
	string name = flux_calculator_name(GlobalConfig.flux_calculator);
	lua_pushstring(L, name.toStringz);
	break;
    case "shear_tolerance":
	lua_pushnumber(L, GlobalConfig.shear_tolerance);
	break;
    case "M_inf":
	lua_pushnumber(L, GlobalConfig.M_inf);
	break;
    case "compression_tolerance":
	lua_pushnumber(L, GlobalConfig.compression_tolerance);
	break;

    case "moving_grid":
	lua_pushboolean(L, GlobalConfig.moving_grid);
	break;
    case "write_vertex_velocities":
	lua_pushboolean(L, GlobalConfig.write_vertex_velocities);
	break;

    case "viscous":
	lua_pushboolean(L, GlobalConfig.viscous);
	break;
    case "viscous_factor":
	lua_pushnumber(L, GlobalConfig.viscous_factor);
	break;
    case "viscous_factor_increment":
	lua_pushnumber(L, GlobalConfig.viscous_factor_increment);
	break;
    case "viscous_delay":
	lua_pushnumber(L, GlobalConfig.viscous_delay);
	break;
    case "viscous_signal_factor":
	lua_pushnumber(L, GlobalConfig.viscous_signal_factor);
	break;

    case "diffusion":
	lua_pushboolean(L, GlobalConfig.diffusion);
	break;

    case "turbulence_model":
	string name = turbulence_model_name(GlobalConfig.turbulence_model);
	lua_pushstring(L, name.toStringz);
	break;
    case "turbulence_prandtl_number":
	lua_pushnumber(L, GlobalConfig.turbulence_prandtl_number);
	break;
    case "turbulence_schmidt_number":
	lua_pushnumber(L, GlobalConfig.turbulence_schmidt_number);
	break;
    case "max_mu_t_factor":
	lua_pushnumber(L, GlobalConfig.max_mu_t_factor);
	break;
    case "transient_mu_t_factor":
	lua_pushnumber(L, GlobalConfig.transient_mu_t_factor);
	break;

    case "reacting":
	lua_pushboolean(L, GlobalConfig.reacting);
	break;
    case "reactions_file":
	lua_pushstring(L, GlobalConfig.reactions_file.toStringz);
	break;

    case "max_step":
	lua_pushnumber(L, GlobalConfig.max_step);
	break;
    case "halt_now":
	lua_pushnumber(L, GlobalConfig.halt_now);
	break;
    case "print_count":
	lua_pushnumber(L, GlobalConfig.print_count);
	break;
    case "control_count":
	lua_pushnumber(L, GlobalConfig.control_count);
	break;
    case "verbosity_level":
	lua_pushnumber(L, GlobalConfig.verbosity_level);
	break;
    case "max_time":
	lua_pushnumber(L, GlobalConfig.max_time);
	break;
    case "dt_init":
	lua_pushnumber(L, GlobalConfig.dt_init);
	break;
    case "dt_max":
	lua_pushnumber(L, GlobalConfig.dt_max);
	break;
    case "cfl_value":
	lua_pushnumber(L, GlobalConfig.cfl_value);
	break;
    case "stringent_cfl":
	lua_pushboolean(L, GlobalConfig.stringent_cfl);
	break;
    case "cfl_count":
	lua_pushnumber(L, GlobalConfig.cfl_count);
	break;
    case "fixed_time_step":
	lua_pushboolean(L, GlobalConfig.fixed_time_step);
	break;
    case "write_at_step":
	lua_pushnumber(L, GlobalConfig.write_at_step);
	break;

    case "dt_plot":
	lua_pushnumber(L, GlobalConfig.dt_plot);
	break;
    case "dt_history":
	lua_pushnumber(L, GlobalConfig.dt_history);
	break;

    case "udf_source_terms_file":
	lua_pushstring(L, GlobalConfig.udf_source_terms_file.toStringz);
	break;
    case "udf_source_terms":
	lua_pushboolean(L, GlobalConfig.udf_source_terms);
	break;

    case "block_marching":
	lua_pushboolean(L, GlobalConfig.block_marching);
	break;
    case "nib":
	lua_pushnumber(L, GlobalConfig.nib);
	break;
    case "njb":
	lua_pushnumber(L, GlobalConfig.njb);
	break;
    case "nkb":
	lua_pushnumber(L, GlobalConfig.nkb);
	break;
    case "propagate_inflow_data":
	lua_pushboolean(L, GlobalConfig.propagate_inflow_data);
	break;

    case "udf_solid_source_terms_file":
	lua_pushstring(L, GlobalConfig.udfSolidSourceTermsFile.toStringz);
	break;
    case "udf_solid_source_terms":
	lua_pushboolean(L, GlobalConfig.udfSolidSourceTerms);
	break;
    default:
	lua_pushnil(L);
    }
    return 1;
} // end configGet()

// Interact via __call, __index and __newindex

extern(C) int configSetWithCall(lua_State* L)
{
    // Arguments to __call are: table then call arguments
    // So remove table and delegate to configSetFromTable
    lua_remove(L, 1);
    return configSetFromTable(L);
}

extern(C) int configSetFromValue(lua_State *L)
{
    // Argumnets to __newindex are: table, key, value
    // We aren't interested in the table because we have
    // the GlobalConfig object to use.
    // Let's put the key and value into a table with one entry
    // and delegate to configSetFromTable.
    lua_newtable(L);
    lua_pushvalue(L, 3);
    lua_setfield(L, -2, luaL_checkstring(L, 2));
    // Now set table to position 1 in stack for use in call to configSetFromTable.
    lua_replace(L, 1);
    return configSetFromTable(L);
}

extern(C) int configGetByKey(lua_State* L)
{
    // Arguments to __index are: table, key
    // Just remove table and delegate to configGet.
    lua_remove(L, 1);
    return configGet(L);
} 

//------------------------------------------------------------------------
// Functions related to the managed gas model.

extern(C) int setGasModel(lua_State* L)
{
    string fname = to!string(luaL_checkstring(L, 1));
    GlobalConfig.gas_model_file = fname;
    GlobalConfig.gmodel_master = init_gas_model(fname);
    lua_pushinteger(L, GlobalConfig.gmodel_master.n_species);
    lua_pushinteger(L, GlobalConfig.gmodel_master.n_modes);
    return 2;
    
}

extern(C) int get_nspecies(lua_State* L)
{
    lua_pushinteger(L, GlobalConfig.gmodel_master.n_species);
    return 1;
}

extern(C) int get_nmodes(lua_State* L)
{
    lua_pushinteger(L, GlobalConfig.gmodel_master.n_modes);
    return 1;
}

extern(C) int species_name(lua_State* L)
{
    int i = to!int(luaL_checkinteger(L, 1));
    lua_pushstring(L, GlobalConfig.gmodel_master.species_name(i).toStringz);
    return 1;
}

//-----------------------------------------------------------------------
// Call the following function from the main program to get the
// functions appearing in the Lua interpreter.

void registerGlobalConfig(lua_State* L)
{
    // Register global functions for setting configuration.
    lua_pushcfunction(L, &configSetFromTable);
    lua_setglobal(L, "configSet");
    lua_pushcfunction(L, &configGet);
    lua_setglobal(L, "configGet");

    // Make a 'config' table available
    // First, set its metatable so that the metamethods
    // of __call, __index, and __newindex can do their jobs.
    luaL_newmetatable(L, "config_mt");
    lua_pushcfunction(L, &configSetWithCall);
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &configSetFromValue);
    lua_setfield(L, -2, "__newindex");
    lua_pushcfunction(L, &configGetByKey);
    lua_setfield(L, -2, "__index");
    lua_setglobal(L, "config_mt");
    // Second make a globally available table called 'config'
    lua_newtable(L);
    luaL_getmetatable(L, "config_mt");
    lua_setmetatable(L, -2);
    lua_setglobal(L, "config");

    // Register other global functions related to the managed gas model.
    lua_pushcfunction(L, &setGasModel);
    lua_setglobal(L, "setGasModel");
    lua_pushcfunction(L, &get_nspecies);
    lua_setglobal(L, "get_nspecies");
    lua_pushcfunction(L, &get_nmodes);
    lua_setglobal(L, "get_nmodes");
    lua_pushcfunction(L, &species_name);
    lua_setglobal(L, "species_name");
} // end registerGlobalConfig()
