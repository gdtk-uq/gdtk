/** readconfig.d
 * Eilmer4 compressible-flow simulation code, reading of JSON config and control files.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import core.stdc.stdlib : exit;
import std.stdio;
import std.json;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.typecons;

import util.lua;
import json_helper;
import geom;
import gas;
import kinetics;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import sblock;
import ublock;
import ssolidblock;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import grid_motion;

// Utility functions to condense the following code.
//
string update_string(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONstring(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_bool(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONbool(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_int(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONint(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_double(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONdouble(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_enum(string key, string field, string enum_from_name)
{
    return "{ // start new block scope
    auto mySaveValue = GlobalConfig."~field~";
    try {
	string name = jsonData[\""~key~"\"].str;
	GlobalConfig."~field~" = "~enum_from_name~"(name);
    } catch (Exception e) {
	GlobalConfig."~field~" = mySaveValue;
    }
    }";
}

void read_config_file()
{
    if (GlobalConfig.verbosity_level > 1) writeln("Read config file.");
    string fileName = GlobalConfig.base_file_name ~ ".config";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
	writeln("Failed to read config file: ", fileName);
	writeln("Message is: ", e.msg);
	exit(1);
    }
    JSONValue jsonData;
    try {
	jsonData = parseJSON!string(content);
    } catch (Exception e) {
	writeln("Failed to parse JSON from config file: ", fileName);
	writeln("Message is: ", e.msg);
	exit(1);
    }
    // Now that we have parsed JSON data, proceed to update those config values.
    // Note that some of the lines below are much longer than PJ would normally tolerate.
    // The trade-off for ease of reading with one line per entry won out... 
    //
    mixin(update_string("title", "title"));
    mixin(update_string("gas_model_file", "gas_model_file"));
    GlobalConfig.gmodel_master = init_gas_model(GlobalConfig.gas_model_file);
    mixin(update_bool("include_quality", "include_quality"));
    mixin(update_int("dimensions", "dimensions"));
    mixin(update_bool("axisymmetric", "axisymmetric"));
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  title: ", to!string(GlobalConfig.title));
	writeln("  gas_model_file: ", to!string(GlobalConfig.gas_model_file));
	writeln("  include_quality: ", GlobalConfig.include_quality);
	writeln("  dimensions: ", GlobalConfig.dimensions);
	writeln("  axisymmetric: ", GlobalConfig.axisymmetric);
    }

    // Parameters controlling convective update and size of storage arrays
    //
    mixin(update_enum("gasdynamic_update_scheme", "gasdynamic_update_scheme", "update_scheme_from_name"));
    GlobalConfig.n_flow_time_levels = 1 + number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
    mixin(update_enum("grid_motion", "grid_motion", "grid_motion_from_name"));
    if (GlobalConfig.grid_motion == GridMotion.none) {
	GlobalConfig.n_grid_time_levels = 1;
    } else {
	GlobalConfig.n_grid_time_levels = 1 + number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
    }
    mixin(update_bool("write_vertex_velocities", "write_vertex_velocities"));
    mixin(update_string("udf_grid_motion_file", "udf_grid_motion_file"));
    // If we have user-defined grid motion, we'll need to initialise
    // the lua_State that holds the user's function. But we can't
    // do that initialisation just yet. We have to wait until all
    // of the blocks are configured since we set that information
    // as global information available to the user. Hence, you'll
    // find that step at the very end of this function.

    // Parameters controlling convective update in detail
    //
    mixin(update_bool("separate_update_for_viscous_terms", "separate_update_for_viscous_terms"));
    mixin(update_bool("separate_update_for_k_omega_source", "separate_update_for_k_omega_source"));
    mixin(update_bool("apply_bcs_in_parallel", "apply_bcs_in_parallel"));
    mixin(update_bool("stringent_cfl", "stringent_cfl"));
    mixin(update_double("flowstate_limits_max_velocity", "flowstate_limits.max_velocity"));
    mixin(update_double("flowstate_limits_max_tke", "flowstate_limits.max_tke"));
    mixin(update_double("flowstate_limits_min_tke", "flowstate_limits.min_tke"));
    mixin(update_double("flowstate_limits_max_temp", "flowstate_limits.max_temp"));
    mixin(update_double("flowstate_limits_min_temp", "flowstate_limits.min_temp"));
    mixin(update_bool("adjust_invalid_cell_data", "adjust_invalid_cell_data"));
    mixin(update_int("max_invalid_cells", "max_invalid_cells"));
    mixin(update_int("interpolation_order", "interpolation_order"));
    mixin(update_enum("thermo_interpolator", "thermo_interpolator", "thermo_interpolator_from_name"));
    mixin(update_bool("apply_limiter", "apply_limiter"));
    mixin(update_bool("extrema_clipping", "extrema_clipping"));
    mixin(update_bool("interpolate_in_local_frame", "interpolate_in_local_frame"));
    mixin(update_bool("retain_least_squares_work_data", "retain_least_squares_work_data"));
    mixin(update_enum("flux_calculator", "flux_calculator", "flux_calculator_from_name"));
    mixin(update_double("shear_tolerance", "shear_tolerance"));
    mixin(update_double("M_inf", "M_inf"));
    mixin(update_double("compression_tolerance", "compression_tolerance"));
    mixin(update_double("shock_fitting_delay", "shock_fitting_delay"));
    mixin(update_int("shock_fitting_interpolation_order", "shock_fitting_interpolation_order"));
    mixin(update_double("shock_fitting_scale_factor", "shock_fitting_scale_factor"));
    mixin(update_bool("MHD", "MHD"));
    mixin(update_bool("divergence_cleaning", "divergence_cleaning"));
    mixin(update_double("divB_damping_length", "divB_damping_length"));

    // Checking of constraints.
    // The following checks/overrides must happen after the relevant config elements
    // have been set.
    if (GlobalConfig.grid_motion == GridMotion.shock_fitting &&
	GlobalConfig.apply_bcs_in_parallel) {
	writeln("NOTE: apply_bcs_in_parallel is set to false when shock_fitting is used.");
	GlobalConfig.apply_bcs_in_parallel = false;
    } 

    if (GlobalConfig.verbosity_level > 1) {
	writeln("  gasdynamic_update_scheme: ", gasdynamic_update_scheme_name(GlobalConfig.gasdynamic_update_scheme));
	writeln("  grid_motion: ", grid_motion_name(GlobalConfig.grid_motion));
	writeln("  shock_fitting_delay: ", GlobalConfig.shock_fitting_delay);
	writeln("  shock_fitting_interpolation_order: ", GlobalConfig.shock_fitting_interpolation_order);
	writeln("  shock_fitting_scale_factor: ", GlobalConfig.shock_fitting_scale_factor);
	writeln("  write_vertex_velocities: ", GlobalConfig.write_vertex_velocities);
	writeln("  udf_grid_motion_file: ", to!string(GlobalConfig.udf_grid_motion_file));
	writeln("  separate_update_for_viscous_terms: ", GlobalConfig.separate_update_for_viscous_terms);
	writeln("  separate_update_for_k_omega_source: ", GlobalConfig.separate_update_for_k_omega_source);
	writeln("  apply_bcs_in_parallel: ", GlobalConfig.apply_bcs_in_parallel);
	writeln("  stringent_cfl: ", GlobalConfig.stringent_cfl);
	writeln("  flowstate_limits_max_velocity: ", GlobalConfig.flowstate_limits.max_velocity);
	writeln("  flowstate_limits_max_tke: ", GlobalConfig.flowstate_limits.max_tke);
	writeln("  flowstate_limits_min_tke: ", GlobalConfig.flowstate_limits.min_tke);
	writeln("  flowstate_limits_max_temp: ", GlobalConfig.flowstate_limits.max_temp);
	writeln("  flowstate_limits_min_temp: ", GlobalConfig.flowstate_limits.min_temp);
	writeln("  adjust_invalid_cell_data: ", GlobalConfig.adjust_invalid_cell_data);
	writeln("  max_invalid_cells: ", GlobalConfig.max_invalid_cells);
	writeln("  interpolation_order: ", GlobalConfig.interpolation_order);
	writeln("  thermo_interpolator: ", thermo_interpolator_name(GlobalConfig.thermo_interpolator));
	writeln("  apply_limiter: ", GlobalConfig.apply_limiter);
	writeln("  extrema_clipping: ", GlobalConfig.extrema_clipping);
	writeln("  interpolate_in_local_frame: ", GlobalConfig.interpolate_in_local_frame);
	writeln("  retain_least_squares_work_data: ", GlobalConfig.retain_least_squares_work_data);	
	writeln("  flux_calculator: ", flux_calculator_name(GlobalConfig.flux_calculator));
	writeln("  shear_tolerance: ", GlobalConfig.shear_tolerance);
	writeln("  M_inf: ", GlobalConfig.M_inf);
	writeln("  compression_tolerance: ", GlobalConfig.compression_tolerance);
	writeln("  MHD: ", GlobalConfig.MHD);
	writeln("  divergence_cleaning: ", GlobalConfig.divergence_cleaning);
	writeln("  divB_damping_length: ", GlobalConfig.divB_damping_length);
    }
    //
    // More checking of constraints.
    //
    if (GlobalConfig.grid_motion != GridMotion.none) {
	if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_1_stage ||
	    GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_2_stage) {
	    // pass, we have a consistent selection.
	} else {
	    string msg = "We have some grid_motion but not a valid GasdynamicUpdate scheme" ~
		" for grid motion.";
	    throw new FlowSolverException(msg);
	}
    } else {
	if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_1_stage ||
	    GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_2_stage) {
	    string msg = "We have no grid_motion but have asked for a GasdynamicUpdate scheme" ~
		" for grid motion.";
	    throw new FlowSolverException(msg);
	}
    }

    // Parameters controlling viscous/molecular transport
    //
    mixin(update_bool("viscous", "viscous"));
    mixin(update_enum("spatial_deriv_calc", "spatial_deriv_calc", "spatial_deriv_calc_from_name"));
    mixin(update_enum("spatial_deriv_locn", "spatial_deriv_locn", "spatial_deriv_locn_from_name"));
    mixin(update_bool("include_ghost_cells_in_spatial_deriv_clouds", "include_ghost_cells_in_spatial_deriv_clouds"));
    mixin(update_bool("spatial_deriv_retain_lsq_work_data", "spatial_deriv_retain_lsq_work_data"));
    mixin(update_double("viscous_delay", "viscous_delay"));
    mixin(update_double("viscous_factor_increment", "viscous_factor_increment"));
    mixin(update_double("viscous_signal_factor", "viscous_signal_factor"));
    mixin(update_enum("turbulence_model", "turbulence_model", "turbulence_model_from_name"));
    mixin(update_double("turbulence_prandtl_number", "turbulence_prandtl_number"));
    mixin(update_double("turbulence_schmidt_number", "turbulence_schmidt_number"));
    mixin(update_double("max_mu_t_factor", "max_mu_t_factor"));
    mixin(update_double("transient_mu_t_factor", "transient_mu_t_factor"));
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  viscous: ", GlobalConfig.viscous);
	writeln("  spatial_deriv_calc: ", spatial_deriv_calc_name(GlobalConfig.spatial_deriv_calc));
	writeln("  spatial_deriv_locn: ", spatial_deriv_locn_name(GlobalConfig.spatial_deriv_locn));
	writeln("  include_ghost_cells_in_spatial_deriv_clouds: ", GlobalConfig.include_ghost_cells_in_spatial_deriv_clouds);
	writeln("  spatial_deriv_retain_lsq_work_data: ", GlobalConfig.spatial_deriv_retain_lsq_work_data);
	writeln("  viscous_delay: ", GlobalConfig.viscous_delay);
	writeln("  viscous_factor_increment: ", GlobalConfig.viscous_factor_increment);
	writeln("  viscous_signal_factor: ", GlobalConfig.viscous_signal_factor);
	writeln("  turbulence_model: ", turbulence_model_name(GlobalConfig.turbulence_model));
	writeln("  turbulence_prandtl_number: ", GlobalConfig.turbulence_prandtl_number);
	writeln("  turbulence_schmidt_number: ", GlobalConfig.turbulence_schmidt_number);
	writeln("  max_mu_t_factor: ", GlobalConfig.max_mu_t_factor);
	writeln("  transient_mu_t_factor: ", GlobalConfig.transient_mu_t_factor);
    }

    // User-defined source terms
    mixin(update_bool("udf_source_terms", "udf_source_terms"));
    mixin(update_string("udf_source_terms_file", "udf_source_terms_file"));
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  udf_source_terms: ", GlobalConfig.udf_source_terms);
	writeln("  udf_source_terms_file: ", to!string(GlobalConfig.udf_source_terms_file));
    }

    // Parameters controlling thermochemistry
    //
    mixin(update_bool("reacting", "reacting"));
    mixin(update_string("reactions_file", "reactions_file"));
    mixin(update_double("reaction_time_delay", "reaction_time_delay"));
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  reacting: ", GlobalConfig.reacting);
	writeln("  reactions_file: ", to!string(GlobalConfig.reactions_file));
	writeln("  reaction_time_delay: ", GlobalConfig.reaction_time_delay);
    }

    // Parameters controlling other simulation options
    //
    mixin(update_int("control_count", "control_count"));
    mixin(update_bool("block_marching", "block_marching"));
    mixin(update_int("nib", "nib"));
    mixin(update_int("njb", "njb"));
    mixin(update_int("nkb", "nkb"));
    mixin(update_bool("propagate_inflow_data", "propagate_inflow_data"));
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  control_count: ", GlobalConfig.control_count);
	writeln("  block_marching: ", GlobalConfig.block_marching);
	writeln("  nib: ", GlobalConfig.nib);
	writeln("  njb: ", GlobalConfig.njb);
	writeln("  nkb: ", GlobalConfig.nkb);
	writeln("  propagate_inflow_data: ", GlobalConfig.propagate_inflow_data);
    }

    int nhcell = getJSONint(jsonData, "nhcell", 0);
    foreach (i; 0 .. nhcell) {
	string jsonKey = format("history-cell-%d", i);
	auto hcell = getJSONintarray(jsonData, jsonKey, [0, 0]);
	GlobalConfig.hcells ~= tuple(cast(size_t) hcell[0], cast(size_t) hcell[1]);
    }
    int nsolidhcell = getJSONint(jsonData, "nsolidhcell", 0);
    foreach (i; 0 .. nsolidhcell) {
	string jsonKey = format("solid-history-cell-%d", i);
	auto hcell = getJSONintarray(jsonData, jsonKey, [0, 0]);
	GlobalConfig.solid_hcells ~= tuple(cast(size_t) hcell[0], cast(size_t) hcell[1]);
    }
    // TODO -- still have other entries such as nheatzone, nreactionzone, ...

    // Now, configure blocks that make up the flow domain.
    //
    // This is done in phases.  The blocks need valid references to LocalConfig objects
    // and the boundary conditions need valid references to Sblock objects.
    mixin(update_int("nblock", "nBlocks"));
    if (GlobalConfig.verbosity_level > 1) { writeln("  nBlocks: ", GlobalConfig.nBlocks); }
    // Set up dedicated copies of the configuration parameters for the threads.
    foreach (i; 0 .. GlobalConfig.nBlocks) {
	dedicatedConfig ~= new LocalConfig();
    }
    foreach (i; 0 .. GlobalConfig.nBlocks) {
	auto jsonDataForBlock = jsonData["block_" ~ to!string(i)];
	string blockType = getJSONstring(jsonDataForBlock, "type", "");
	switch (blockType) {
	case "SBlock": 
	    gasBlocks ~= new SBlock(i, jsonDataForBlock);
	    break;
	case "UBlock":
	    gasBlocks ~= new UBlock(i, jsonDataForBlock);
	    break;
	default:
	    throw new Error(format("Construction of block[%d], unknown type: %s",
				   i, blockType));
	} // end switch blockType
    }
    foreach (blk; gasBlocks) {
	blk.init_lua_globals();
	blk.init_boundary_conditions(jsonData["block_" ~ to!string(blk.id)]);
	if (GlobalConfig.udf_source_terms) {
	    luaL_dofile(blk.myL, GlobalConfig.udf_source_terms_file.toStringz);
	}
    } 
    // After fully constructing blocks and their boundary conditions,
    // we can optionally print their representation for checking.
    if (GlobalConfig.verbosity_level > 1) {
	foreach (i, blk; gasBlocks) { writeln("  Block[", i, "]: ", blk); }
    }
    // Read in any blocks in the solid domain.
    GlobalConfig.udfSolidSourceTerms = getJSONbool(jsonData, "udf_solid_source_terms", false);
    GlobalConfig.udfSolidSourceTermsFile = jsonData["udf_solid_source_terms_file"].str;
    GlobalConfig.nSolidBlocks = getJSONint(jsonData, "nsolidblock", 0);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  nSolidBlocks: ", GlobalConfig.nSolidBlocks);
	writeln("  udf_solid_source_terms: ", GlobalConfig.udfSolidSourceTerms);
	writeln("  udf_solid_source_terms_file: ", to!string(GlobalConfig.udfSolidSourceTermsFile));
    }
    // Set up dedicated copies of the configuration parameters for the threads.
    foreach (i; 0 .. GlobalConfig.nSolidBlocks) {
	dedicatedSolidConfig ~= new LocalConfig();
    }
    foreach (i; 0 .. GlobalConfig.nSolidBlocks) {
	solidBlocks ~= new SSolidBlock(i, jsonData["solid_block_" ~ to!string(i)]);
	if (GlobalConfig.verbosity_level > 1) {
	    writeln("  SolidBlock[", i, "]: ", solidBlocks[i]);
	}
    }
    foreach (sblk; solidBlocks) {
	sblk.initLuaGlobals();
	sblk.initBoundaryConditions(jsonData["solid_block_" ~ to!string(sblk.id)]);
	if ( GlobalConfig.udfSolidSourceTerms ) {
	    initUDFSolidSourceTerms(sblk.myL, GlobalConfig.udfSolidSourceTermsFile);
	}
    }

    // Now that the blocks are configured, we can set up the
    // lua_State for a user-defined moving grid, if needed.
    if ( GlobalConfig.grid_motion == GridMotion.user_defined ) {
	// We'll need to initialise the lua_State that holds the
	// user's function for defining grid motion.
	init_master_lua_State(GlobalConfig.udf_grid_motion_file);
    }
} // end read_config_file()

void read_control_file()
{
    if (GlobalConfig.verbosity_level > 1) writeln("read_control_file()");
    string fileName = GlobalConfig.base_file_name ~ ".control";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
	writeln("Failed to read control file: ", fileName);
	exit(1);
    }
    JSONValue jsonData;
    try {
	jsonData = parseJSON!string(content);
    } catch (Exception e) {
	writeln("Failed to parse JSON from control file: ", fileName);
	exit(1);
    }
    mixin(update_int("max_step", "max_step"));
    mixin(update_double("max_time", "max_time"));
    mixin(update_int("halt_now", "halt_now"));
    mixin(update_int("print_count", "print_count"));
    mixin(update_int("cfl_count", "cfl_count"));
    mixin(update_double("dt_init", "dt_init"));
    mixin(update_double("dt_max", "dt_max"));
    mixin(update_double("cfl_value", "cfl_value"));
    mixin(update_bool("fixed_time_step", "fixed_time_step"));
    mixin(update_double("dt_plot", "dt_plot"));
    mixin(update_double("dt_history", "dt_history"));
    //
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  max_step: ", GlobalConfig.max_step);
	writeln("  max_time: ", GlobalConfig.max_time);
	writeln("  halt_now: ", GlobalConfig.halt_now);
	writeln("  print_count: ", GlobalConfig.print_count);
	writeln("  cfl_count: ", GlobalConfig.cfl_count);
	writeln("  dt_init: ", GlobalConfig.dt_init);
	writeln("  dt_max: ", GlobalConfig.dt_max);
	writeln("  cfl_value: ", GlobalConfig.cfl_value);
	writeln("  fixed_time_step: ", GlobalConfig.fixed_time_step);
	writeln("  dt_plot: ", GlobalConfig.dt_plot);
	writeln("  dt_history: ", GlobalConfig.dt_history);
    }
} // end read_control_file()
