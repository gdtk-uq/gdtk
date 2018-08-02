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
import gas.luagas_model;
import fvcore;
import globalconfig;
import mass_diffusion;

// -------------------------------------------------------------------------------
// Set GlobalConfig fields from a table.

extern(C) int configSetFromTable(lua_State* L)
{
    if (!lua_istable(L, 1)) return 0; // nothing to do
    //
    // Functions to generate the repetitive code for getting values from
    // the individual fields in the Lua table.
    string get_string_field(string lua_field, string config_field)
    {
        string code = "lua_getfield(L, 1, \""~lua_field~"\");
        if (!lua_isnil(L, -1)) {
            GlobalConfig."~config_field~" = to!string(luaL_checkstring(L, -1));
            lua_pushnil(L); lua_setfield(L, 1, \""~lua_field~"\");
        }
        lua_pop(L, 1);";
        return code;
    }
    string get_enum_field(string lua_field, string config_field, string enum_from_name)
    {
        string code = "lua_getfield(L, 1, \""~lua_field~"\");
        if (!lua_isnil(L, -1)) {
            GlobalConfig."~config_field~" = "~enum_from_name~"(to!string(luaL_checkstring(L, -1)));
            lua_pushnil(L); lua_setfield(L, 1, \""~lua_field~"\");
        }
        lua_pop(L, 1);";
        return code;
    }
    string get_bool_field(string lua_field, string config_field)
    {
        string code = "lua_getfield(L, 1, \""~lua_field~"\");
        if (!lua_isnil(L, -1)) {
            GlobalConfig."~config_field~" = to!bool(lua_toboolean(L, -1));
            lua_pushnil(L); lua_setfield(L, 1, \""~lua_field~"\");
        }
        lua_pop(L, 1);";
        return code;
    }
    string get_int_field(string lua_field, string config_field)
    {
        string code = "lua_getfield(L, 1, \""~lua_field~"\");
        if (!lua_isnil(L, -1)) {
            GlobalConfig."~config_field~" = luaL_checkint(L, -1);
            lua_pushnil(L); lua_setfield(L, 1, \""~lua_field~"\");
        }
        lua_pop(L, 1);";
        return code;
    }
    string get_double_field(string lua_field, string config_field)
    {
        string code = "lua_getfield(L, 1, \""~lua_field~"\");
        if (!lua_isnil(L, -1)) {
            GlobalConfig."~config_field~" = to!double(luaL_checknumber(L, -1));
            lua_pushnil(L); lua_setfield(L, 1, \""~lua_field~"\");
        }
        lua_pop(L, 1);";
        return code;
    }
    //
    // Look for fields that may be present in the Lua table.
    // Once the value is read from the field, remove that field from the table.
    // This will allow us to check for and warn about fields that have been 
    // accidentally set in the Lua table.
    mixin(get_string_field("base_file_name", "base_file_name"));
    mixin(get_string_field("grid_format", "grid_format"));
    mixin(get_string_field("flow_format", "flow_format"));
    mixin(get_string_field("title", "title"));
    mixin(get_string_field("gas_model_file", "gas_model_file"));
    mixin(get_string_field("udf_supervisor_file", "udf_supervisor_file"));
    mixin(get_bool_field("include_quality", "include_quality"));
    //
    mixin(get_int_field("nFluidBlocks","nFluidBlocks"));
    mixin(get_int_field("dimensions","dimensions"));
    mixin(get_bool_field("axisymmetric","axisymmetric"));
    //
    mixin(get_bool_field("MHD","MHD"));
    mixin(get_bool_field("MHD_static_field","MHD_static_field"));
    mixin(get_bool_field("MHD_resistive","MHD_resistive"));
    mixin(get_bool_field("divergence_cleaning","divergence_cleaning"));
    mixin(get_double_field("c_h","c_h"));
    mixin(get_double_field("divB_damping_length","divB_damping_length"));
    //
    mixin(get_enum_field("strang_splitting", "strangSplitting", "strangSplittingModeFromName"));
    mixin(get_enum_field("gasdynamic_update_scheme", "gasdynamic_update_scheme", "update_scheme_from_name"));
    mixin(get_enum_field("coupling_with_solid_domains", "coupling_with_solid_domains", "solidDomainCouplingFromName"));
    mixin(get_bool_field("solid_has_isotropic_properties", "solid_has_isotropic_properties"));
    mixin(get_bool_field("solid_has_homogeneous_properties", "solid_has_homogeneous_properties"));
    mixin(get_bool_field("apply_bcs_in_parallel", "apply_bcs_in_parallel"));
    mixin(get_double_field("flowstate_limits_max_velocity", "flowstate_limits.max_velocity"));
    mixin(get_double_field("flowstate_limits_max_tke", "flowstate_limits.max_tke"));
    mixin(get_double_field("flowstate_limits_min_tke", "flowstate_limits.min_tke"));
    mixin(get_double_field("flowstate_limits_max_temp", "flowstate_limits.max_temp"));
    mixin(get_double_field("flowstate_limits_min_temp", "flowstate_limits.min_temp"));
    mixin(get_bool_field("adjust_invalid_cell_data", "adjust_invalid_cell_data"));
    mixin(get_bool_field("report_invalid_cells", "report_invalid_cells"));
    mixin(get_int_field("max_invalid_cells", "max_invalid_cells"));
    //
    mixin(get_int_field("interpolation_order", "interpolation_order"));
    mixin(get_double_field("interpolation_delay", "interpolation_delay"));
    mixin(get_enum_field("thermo_interpolator", "thermo_interpolator", "thermo_interpolator_from_name"));
    mixin(get_bool_field("allow_reconstruction_for_energy_modes", "allow_reconstruction_for_energy_modes"));
    mixin(get_bool_field("apply_limiter", "apply_limiter"));
    mixin(get_bool_field("extrema_clipping", "extrema_clipping"));
    mixin(get_bool_field("interpolate_in_local_frame", "interpolate_in_local_frame"));
    mixin(get_enum_field("unstructured_limiter", "unstructured_limiter", "unstructured_limiter_from_name"));
    mixin(get_bool_field("use_extended_stencil", "use_extended_stencil"));
    mixin(get_double_field("venkat_K_value", "venkat_K_value"));
    mixin(get_enum_field("flux_calculator", "flux_calculator", "flux_calculator_from_name"));
    mixin(get_double_field("shear_tolerance", "shear_tolerance"));
    mixin(get_double_field("M_inf", "M_inf"));
    mixin(get_double_field("compression_tolerance", "compression_tolerance"));
    mixin(get_bool_field("artificial_compressiblity", "artificial_compressibility"));
    mixin(get_double_field("ac_alpha", "ac_alpha"));
    //
    mixin(get_enum_field("grid_motion", "grid_motion", "grid_motion_from_name"));
    mixin(get_double_field("shock_fitting_delay", "shock_fitting_delay"));
    mixin(get_int_field("shock_fitting_interpolation_order","shock_fitting_interpolation_order"));
    mixin(get_double_field("shock_fitting_scale_factor", "shock_fitting_scale_factor"));
    mixin(get_bool_field("write_vertex_velocities", "write_vertex_velocities"));
    mixin(get_string_field("udf_grid_motion_file", "udf_grid_motion_file"));
    //
    mixin(get_bool_field("viscous", "viscous"));
    mixin(get_bool_field("use_viscosity_from_cells", "use_viscosity_from_cells"));
    mixin(get_bool_field("spatial_deriv_from_many_points", "spatial_deriv_from_many_points"));
    mixin(get_enum_field("spatial_deriv_calc", "spatial_deriv_calc", "spatial_deriv_calc_from_name"));
    mixin(get_enum_field("spatial_deriv_locn", "spatial_deriv_locn", "spatial_deriv_locn_from_name"));
    mixin(get_bool_field("include_ghost_cells_in_spatial_deriv_clouds", "include_ghost_cells_in_spatial_deriv_clouds"));
    mixin(get_bool_field("suppress_reconstruction_at_boundaries", "suppress_reconstruction_at_boundaries"));
    mixin(get_double_field("viscous_factor_increment", "viscous_factor_increment"));
    mixin(get_double_field("viscous_delay", "viscous_delay"));
    mixin(get_double_field("shear_stress_relative_limit", "shear_stress_relative_limit"));
    mixin(get_double_field("viscous_signal_factor", "viscous_signal_factor"));
    mixin(get_double_field("turbulent_signal_factor", "turbulent_signal_factor"));
    mixin(get_enum_field("mass_diffusion_model", "mass_diffusion_model", "massDiffusionModelFromName"));
    mixin(get_bool_field("constant_lewis_number", "constant_lewis_number"));
    mixin(get_double_field("lewis_number", "lewis_number"));
    //
    mixin(get_bool_field("separate_update_for_viscous_terms", "separate_update_for_viscous_terms"));
    mixin(get_bool_field("separate_update_for_k_omega_source", "separate_update_for_k_omega_source"));
    //
    mixin(get_enum_field("turbulence_model", "turbulence_model", "turbulence_model_from_name"));
    mixin(get_double_field("turbulence_prandtl_number", "turbulence_prandtl_number"));
    mixin(get_double_field("turbulence_schmidt_number", "turbulence_schmidt_number"));
    mixin(get_double_field("max_mu_t_factor", "max_mu_t_factor"));
    mixin(get_double_field("transient_mu_t_factor", "transient_mu_t_factor"));
    mixin(get_bool_field("limit_tke_production", "limit_tke_production"));
    mixin(get_double_field("tke_production_limit_in_kelvins", "tke_production_limit_in_kelvins"));
    //
    mixin(get_enum_field("tci_model", "tci_model", "tci_model_from_name"));
    //
    mixin(get_bool_field("reacting", "reacting"));
    mixin(get_string_field("reactions_file", "reactions_file"));
    mixin(get_double_field("reaction_time_delay", "reaction_time_delay"));
    mixin(get_double_field("T_frozen", "T_frozen"));
    mixin(get_double_field("T_frozen_energy", "T_frozen_energy"));
    mixin(get_double_field("ignition_time_start", "ignition_time_start"));
    mixin(get_double_field("ignition_time_stop", "ignition_time_stop"));
    mixin(get_string_field("energy_exchange_file", "energy_exchange_file"));
    //
    mixin(get_int_field("max_step", "max_step"));
    mixin(get_int_field("halt_now", "halt_now"));
    mixin(get_int_field("print_count", "print_count"));
    mixin(get_int_field("control_count", "control_count"));
    mixin(get_int_field("verbosity_level", "verbosity_level"));
    mixin(get_double_field("max_time", "max_time"));
    mixin(get_double_field("dt_init", "dt_init"));
    mixin(get_double_field("dt_max", "dt_max"));
    mixin(get_double_field("cfl_value", "cfl_value"));
    mixin(get_bool_field("stringent_cfl", "stringent_cfl"));
    mixin(get_int_field("cfl_count", "cfl_count"));
    mixin(get_bool_field("fixed_time_step", "fixed_time_step"));
    mixin(get_double_field("dt_plot", "dt_plot"));
    mixin(get_double_field("dt_history", "dt_history"));
    mixin(get_double_field("dt_loads", "dt_loads"));
    mixin(get_string_field("boundary_group_for_loads", "boundary_group_for_loads"));
    mixin(get_bool_field("compute_loads", "compute_loads"));
    //
    mixin(get_bool_field("diffuse_wall_bcs_on_init", "diffuseWallBCsOnInit"));
    mixin(get_int_field("number_init_passes", "nInitPasses")); 
    mixin(get_double_field("wall_temperature_on_init", "initTWall"));
    //
    mixin(get_bool_field("block_marching", "block_marching"));
    mixin(get_int_field("nib", "nib"));
    mixin(get_int_field("njb", "njb"));
    mixin(get_int_field("nkb", "nkb"));
    mixin(get_bool_field("propagate_inflow_data", "propagate_inflow_data"));
    mixin(get_bool_field("save_intermediate_results", "save_intermediate_results"));
    //
    mixin(get_string_field("udf_source_terms_file", "udf_source_terms_file"));
    mixin(get_bool_field("udf_source_terms", "udf_source_terms"));
    //
    mixin(get_string_field("udf_solid_source_terms_file", "udfSolidSourceTermsFile"));
    mixin(get_bool_field("udf_solid_source_terms", "udfSolidSourceTerms"));
    //
    mixin(get_double_field("thermionic_emission_bc_time_delay", "thermionic_emission_bc_time_delay"));

    // Look for unused keys. These are unsupported keys that the user
    // has supplied. Give a warning.
    lua_pushnil(L);
    while ( lua_next(L, 1) != 0 ) {
        string key = to!string(lua_tostring(L, -2));
        writeln("WARNING: -----------------------------------------------------------------------------");
        writeln(format("WARNING: The configuration option '%s' is not supported. It has been ignored.", key));
        writeln("WARNING: -----------------------------------------------------------------------------");
        lua_pop(L, 1);
    }

    return 0;
} // end configSetFromTable()

// Get GlobalConfig fields by their string name.
extern(C) int configGet(lua_State* L)
{
    string fieldName = to!string(luaL_checkstring(L, 1));
    // In the following switch statement, we'll try to keep to one line per configuration item
    // even if that results in fairly long lines.  I do acknowledge that this goes against my
    // long time opinion to keep short lines, however, I think that it's easier to scan.
    // PJ 2016-06-23
    switch (fieldName) {
    case "base_file_name": lua_pushstring(L, GlobalConfig.base_file_name.toStringz); break;
    case "grid_format": lua_pushstring(L, GlobalConfig.grid_format.toStringz); break;
    case "flow_format": lua_pushstring(L, GlobalConfig.flow_format.toStringz); break;
    case "title": lua_pushstring(L, GlobalConfig.title.toStringz); break;
    case "gas_model_file": lua_pushstring(L, GlobalConfig.gas_model_file.toStringz); break;
    case "udf_supervisor_file": lua_pushstring(L, toStringz(GlobalConfig.udf_supervisor_file)); break;
    case "include_quality": lua_pushboolean(L, GlobalConfig.include_quality); break;
        //
    case "nFluidBlocks": lua_pushnumber(L, GlobalConfig.nFluidBlocks); break;
    case "dimensions": lua_pushnumber(L, GlobalConfig.dimensions); break;
    case "axisymmetric": lua_pushboolean(L, GlobalConfig.axisymmetric); break;
        //
    case "MHD": lua_pushboolean(L, GlobalConfig.MHD); break;
    case "MHD_static_field": lua_pushboolean(L, GlobalConfig.MHD_static_field); break;
    case "MHD_resistive": lua_pushboolean(L, GlobalConfig.MHD_resistive); break;

    case "divergence_cleaning": lua_pushboolean(L, GlobalConfig.divergence_cleaning); break;
    case "c_h": lua_pushnumber(L, GlobalConfig.c_h); break;
    case "divB_damping_length": lua_pushnumber(L, GlobalConfig.divB_damping_length); break;
        //
    case "strang_splitting" : lua_pushstring(L, strangSplittingModeName(GlobalConfig.strangSplitting).toStringz); break;    
    case "gasdynamic_update_scheme": lua_pushstring(L, gasdynamic_update_scheme_name(GlobalConfig.gasdynamic_update_scheme).toStringz); break;
    case "coupling_with_solid_domains": lua_pushstring(L, solidDomainCouplingName(GlobalConfig.coupling_with_solid_domains).toStringz); break;
    case "solid_has_isotropic_properties": lua_pushboolean(L, GlobalConfig.solid_has_isotropic_properties); break;
    case "solid_has_homogeneous_properties": lua_pushboolean(L, GlobalConfig.solid_has_homogeneous_properties); break;
    case "apply_bcs_in_parallel": lua_pushboolean(L, GlobalConfig.apply_bcs_in_parallel); break;
    case "flowstate_limits_max_velocity": lua_pushnumber(L, GlobalConfig.flowstate_limits.max_velocity); break;
    case "flowstate_limits_max_tke": lua_pushnumber(L, GlobalConfig.flowstate_limits.max_tke); break;
    case "flowstate_limits_min_tke": lua_pushnumber(L, GlobalConfig.flowstate_limits.min_tke); break;
    case "flowstate_limits_max_temp": lua_pushnumber(L, GlobalConfig.flowstate_limits.max_temp); break;
    case "flowstate_limits_min_temp": lua_pushnumber(L, GlobalConfig.flowstate_limits.min_temp); break;
    case "adjust_invalid_cell_data": lua_pushboolean(L, GlobalConfig.adjust_invalid_cell_data); break;
    case "report_invalid_cells": lua_pushboolean(L, GlobalConfig.report_invalid_cells); break;
    case "max_invalid_cells": lua_pushnumber(L, GlobalConfig.max_invalid_cells); break;
        //
    case "interpolation_order": lua_pushnumber(L, GlobalConfig.interpolation_order); break;
    case "interpolation_delay": lua_pushnumber(L, GlobalConfig.interpolation_delay); break;
    case "thermo_interpolator": lua_pushstring(L, thermo_interpolator_name(GlobalConfig.thermo_interpolator).toStringz); break;
    case "allow_reconstruction_for_energy_modes": lua_pushboolean(L, GlobalConfig.allow_reconstruction_for_energy_modes); break;
    case "apply_limiter": lua_pushboolean(L, GlobalConfig.apply_limiter); break;
    case "extrema_clipping": lua_pushboolean(L, GlobalConfig.extrema_clipping); break;
    case "interpolate_in_local_frame": lua_pushboolean(L, GlobalConfig.interpolate_in_local_frame); break;
    case "unstructured_limiter": lua_pushstring(L, unstructured_limiter_name(GlobalConfig.unstructured_limiter).toStringz); break;
    case "use_extended_stencil": lua_pushboolean(L, GlobalConfig.use_extended_stencil); break;
    case "venkat_K_value": lua_pushnumber(L, GlobalConfig.venkat_K_value); break;
    case "flux_calculator": lua_pushstring(L, flux_calculator_name(GlobalConfig.flux_calculator).toStringz); break;
    case "shear_tolerance": lua_pushnumber(L, GlobalConfig.shear_tolerance); break;
    case "M_inf": lua_pushnumber(L, GlobalConfig.M_inf); break;
    case "compression_tolerance": lua_pushnumber(L, GlobalConfig.compression_tolerance); break;
    case "artificial_compressibility": lua_pushboolean(L, GlobalConfig.artificial_compressibility); break;
    case "ac_alpha": lua_pushnumber(L, GlobalConfig.ac_alpha); break;
        //
    case "grid_motion": lua_pushstring(L, grid_motion_name(GlobalConfig.grid_motion).toStringz); break;
    case "write_vertex_velocities": lua_pushboolean(L, GlobalConfig.write_vertex_velocities); break;
    case "udf_grid_motion_file": lua_pushstring(L, toStringz(GlobalConfig.udf_grid_motion_file)); break;
    case "shock_fitting_delay": lua_pushnumber(L, GlobalConfig.shock_fitting_delay); break;
    case "shock_fitting_interpolation_order": lua_pushinteger(L, GlobalConfig.shock_fitting_interpolation_order); break;
    case "shock_fitting_scale_factor": lua_pushnumber(L, GlobalConfig.shock_fitting_scale_factor); break;
        //
    case "viscous": lua_pushboolean(L, GlobalConfig.viscous); break;
    case "use_viscosity_from_cells": lua_pushboolean(L, GlobalConfig.use_viscosity_from_cells); break;
    case "spatial_deriv_from_many_points": lua_pushboolean(L, GlobalConfig.spatial_deriv_from_many_points); break;
    case "spatial_deriv_calc": lua_pushstring(L, spatial_deriv_calc_name(GlobalConfig.spatial_deriv_calc).toStringz); break;
    case "spatial_deriv_locn": lua_pushstring(L, spatial_deriv_locn_name(GlobalConfig.spatial_deriv_locn).toStringz); break;
    case "include_ghost_cells_in_spatial_deriv_clouds": lua_pushboolean(L, GlobalConfig.include_ghost_cells_in_spatial_deriv_clouds); break;
    case "suppress_reconstruction_at_boundaries": lua_pushboolean(L, GlobalConfig.suppress_reconstruction_at_boundaries); break;
    case "viscous_factor_increment": lua_pushnumber(L, GlobalConfig.viscous_factor_increment); break;
    case "viscous_delay": lua_pushnumber(L, GlobalConfig.viscous_delay); break;
    case "shear_stress_relative_limit": lua_pushnumber(L, GlobalConfig.shear_stress_relative_limit); break;
    case "viscous_signal_factor": lua_pushnumber(L, GlobalConfig.viscous_signal_factor); break;
    case "turbulent_signal_factor": lua_pushnumber(L, GlobalConfig.turbulent_signal_factor); break;
    case "mass_diffusion_model": lua_pushstring(L, massDiffusionModelName(GlobalConfig.mass_diffusion_model).toStringz); break;
    case "constant_lewis_number": lua_pushboolean(L, GlobalConfig.constant_lewis_number); break;
    case "lewis_number": lua_pushnumber(L, GlobalConfig.lewis_number); break;
        //
    case "separate_update_for_viscous_terms": lua_pushboolean(L, GlobalConfig.separate_update_for_viscous_terms); break;
    case "separate_update_for_k_omega_source": lua_pushboolean(L, GlobalConfig.separate_update_for_k_omega_source); break;
        //
    case "turbulence_model": lua_pushstring(L, turbulence_model_name(GlobalConfig.turbulence_model).toStringz); break;
    case "turbulence_prandtl_number": lua_pushnumber(L, GlobalConfig.turbulence_prandtl_number); break;
    case "turbulence_schmidt_number": lua_pushnumber(L, GlobalConfig.turbulence_schmidt_number); break;
    case "max_mu_t_factor": lua_pushnumber(L, GlobalConfig.max_mu_t_factor); break;
    case "transient_mu_t_factor": lua_pushnumber(L, GlobalConfig.transient_mu_t_factor); break;
    case "limit_tke_production": lua_pushboolean(L, GlobalConfig.limit_tke_production); break;
    case "tke_production_limit_in_kelvins": lua_pushnumber(L, GlobalConfig.tke_production_limit_in_kelvins); break;
        //
    case "tci_model": lua_pushstring(L, tci_model_name(GlobalConfig.tci_model).toStringz); break;
        //
    case "reacting": lua_pushboolean(L, GlobalConfig.reacting); break;
    case "reactions_file": lua_pushstring(L, GlobalConfig.reactions_file.toStringz); break;
    case "reaction_time_delay": lua_pushnumber(L, GlobalConfig.reaction_time_delay); break;
    case "T_frozen": lua_pushnumber(L, GlobalConfig.T_frozen); break;
    case "T_frozen_energy": lua_pushnumber(L, GlobalConfig.T_frozen_energy); break;
    case "ignition_time_start": lua_pushnumber(L, GlobalConfig.ignition_time_start); break;
    case "ignition_time_stop": lua_pushnumber(L, GlobalConfig.ignition_time_stop); break;
    case "energy_exchange_file": lua_pushstring(L, GlobalConfig.energy_exchange_file.toStringz); break;
        //
    case "max_step": lua_pushnumber(L, GlobalConfig.max_step); break;
    case "halt_now": lua_pushnumber(L, GlobalConfig.halt_now); break;
    case "print_count": lua_pushnumber(L, GlobalConfig.print_count); break;
    case "control_count": lua_pushnumber(L, GlobalConfig.control_count); break;
    case "verbosity_level": lua_pushnumber(L, GlobalConfig.verbosity_level); break;
    case "max_time": lua_pushnumber(L, GlobalConfig.max_time); break;
    case "dt_init": lua_pushnumber(L, GlobalConfig.dt_init); break;
    case "dt_max": lua_pushnumber(L, GlobalConfig.dt_max); break;
    case "cfl_value": lua_pushnumber(L, GlobalConfig.cfl_value); break;
    case "stringent_cfl": lua_pushboolean(L, GlobalConfig.stringent_cfl); break;
    case "cfl_count": lua_pushnumber(L, GlobalConfig.cfl_count); break;
    case "fixed_time_step": lua_pushboolean(L, GlobalConfig.fixed_time_step); break;
    case "dt_plot": lua_pushnumber(L, GlobalConfig.dt_plot); break;
    case "dt_history": lua_pushnumber(L, GlobalConfig.dt_history); break;
    case "dt_loads": lua_pushnumber(L, GlobalConfig.dt_loads); break;
    case "boundary_group_for_loads": lua_pushstring(L, GlobalConfig.boundary_group_for_loads.toStringz); break;
    case "compute_loads": lua_pushboolean(L, GlobalConfig.compute_loads); break;
        //
    case "diffuse_wall_bcs_on_init": lua_pushboolean(L, GlobalConfig.diffuseWallBCsOnInit); break;
    case "number_init_passes": lua_pushnumber(L, GlobalConfig.nInitPasses); break;
    case "wall_temperature_on_init": lua_pushnumber(L, GlobalConfig.initTWall); break;
        //
    case "block_marching": lua_pushboolean(L, GlobalConfig.block_marching); break;
    case "nib": lua_pushnumber(L, GlobalConfig.nib); break;
    case "njb": lua_pushnumber(L, GlobalConfig.njb); break;
    case "nkb": lua_pushnumber(L, GlobalConfig.nkb); break;
    case "propagate_inflow_data": lua_pushboolean(L, GlobalConfig.propagate_inflow_data); break;
    case "save_intermediate_results": lua_pushboolean(L, GlobalConfig.save_intermediate_results); break;
        //
    case "udf_source_terms_file": lua_pushstring(L, GlobalConfig.udf_source_terms_file.toStringz); break;
    case "udf_source_terms": lua_pushboolean(L, GlobalConfig.udf_source_terms); break;
    case "udf_solid_source_terms_file": lua_pushstring(L, GlobalConfig.udfSolidSourceTermsFile.toStringz); break;
    case "udf_solid_source_terms": lua_pushboolean(L, GlobalConfig.udfSolidSourceTerms); break;
        //
    case "thermionic_emission_bc_time_delay": lua_pushnumber(L, GlobalConfig.thermionic_emission_bc_time_delay); break;
        //       
    default: lua_pushnil(L);
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
    // Let's first see if the value is non-nil. If nil,
    // we'd like to warn the user.
    if ( lua_isnil(L, 3) ) {
        writeln("WARNING: -----------------------------------------------------------------------------");
        writeln(format("WARNING: You tried to set the configuration option '%s' with a nil value."~
                       " It has been ignored.", to!string(luaL_checkstring(L, 2))));
        writeln("WARNING: -----------------------------------------------------------------------------");
    }
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
    GasModelStore ~= pushObj!(GasModel, GasModelMT)(L, GlobalConfig.gmodel_master);
    return 3;
    
}

extern(C) int getGasModel(lua_State* L)
{
    pushObj!(GasModel, GasModelMT)(L, GlobalConfig.gmodel_master);
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
    lua_pushcfunction(L, &getGasModel);
    lua_setglobal(L, "getGasModel");
} // end registerGlobalConfig()
