/** readconfig.d
 * Eilmer4 compressible-flow simulation code, reading of JSON config and control files.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

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

void read_config_file()
{
    if (GlobalConfig.verbosity_level > 1) writeln("Read config file.");
    string fileName = GlobalConfig.base_file_name ~ ".config";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
	writeln("Failed to read config file: ", fileName);
	exit(1);
    }
    JSONValue jsonData;
    try {
	jsonData = parseJSON!string(content);
    } catch (Exception e) {
	writeln("Failed to parse JSON from config file: ", fileName);
	exit(1);
    }

    // Now that we have parsed JSON data, dip into it to get config values.
    //
    GlobalConfig.title = jsonData["title"].str;
    GlobalConfig.gas_model_file = jsonData["gas_model_file"].str;
    GlobalConfig.gmodel_master = init_gas_model(GlobalConfig.gas_model_file);
    GlobalConfig.include_quality = getJSONbool(jsonData, "include_quality", false);
    GlobalConfig.dimensions = getJSONint(jsonData, "dimensions", 2);
    GlobalConfig.axisymmetric = getJSONbool(jsonData, "axisymmetric", false);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  title: ", to!string(GlobalConfig.title));
	writeln("  gas_model_file: ", to!string(GlobalConfig.gas_model_file));
	writeln("  include_quality: ", GlobalConfig.include_quality);
	writeln("  dimensions: ", GlobalConfig.dimensions);
	writeln("  axisymmetric: ", GlobalConfig.axisymmetric);
    }

    // Parameters controlling convective update and size of storage arrays
    //
    GlobalConfig.interpolation_order = getJSONint(jsonData, "interpolation_order", 2);
    try {
	string name = jsonData["gasdynamic_update_scheme"].str;
	GlobalConfig.gasdynamic_update_scheme = update_scheme_from_name(name);
    } catch (Exception e) {
	GlobalConfig.gasdynamic_update_scheme = GasdynamicUpdate.pc;
    }
    GlobalConfig.n_flow_time_levels = 1 +
	number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
    try {
	string name = jsonData["grid_motion"].str;
	GlobalConfig.grid_motion = grid_motion_from_name(name);
    } catch (Exception e) {
	GlobalConfig.grid_motion = GridMotion.none;
    }
    if (GlobalConfig.grid_motion == GridMotion.none) {
	GlobalConfig.n_grid_time_levels = 1;
    } else {
	GlobalConfig.n_grid_time_levels = 2; // [TODO] Kyle, is this correct? PJ 2016-03-04
    }
    GlobalConfig.write_vertex_velocities = 
	getJSONbool(jsonData, "write_vertex_velocities", false);
    GlobalConfig.udf_grid_motion_file = jsonData["udf_grid_motion_file"].str;
    // If we have user-defined grid motion, we'll need to initialise
    // the lua_State that holds the user's function. But we can't
    // do that initialisation just yet. We have to wait until all
    // of the blocks are configured since we set that information
    // as global information available to the user. Hence, you'll
    // find that step at the very end of this function.

    // Parameters controlling convective update in detail
    //
    GlobalConfig.separate_update_for_viscous_terms =
	getJSONbool(jsonData, "separate_update_for_viscous_terms", false);
    GlobalConfig.separate_update_for_k_omega_source =
	getJSONbool(jsonData, "separate_update_for_k_omega_source", false);
    GlobalConfig.apply_bcs_in_parallel =
	getJSONbool(jsonData, "apply_bcs_in_parallel", true);
    GlobalConfig.stringent_cfl = getJSONbool(jsonData, "stringent_cfl", false);
    GlobalConfig.adjust_invalid_cell_data =
	getJSONbool(jsonData, "adjust_invalid_cell_data", false);
    GlobalConfig.max_invalid_cells = getJSONint(jsonData, "max_invalid_cells", 0);
    try {
	string name = jsonData["thermo_interpolator"].str;
	GlobalConfig.thermo_interpolator = thermo_interpolator_from_name(name);
    } catch (Exception e) {
	GlobalConfig.thermo_interpolator = InterpolateOption.rhoe;
    }
    GlobalConfig.apply_limiter = getJSONbool(jsonData, "apply_limiter", true);
    GlobalConfig.extrema_clipping = getJSONbool(jsonData, "extrema_clipping", true);
    GlobalConfig.interpolate_in_local_frame = 
	getJSONbool(jsonData, "interpolate_in_local_frame", true);
    try {
	string name = jsonData["flux_calculator"].str;
	GlobalConfig.flux_calculator = flux_calculator_from_name(name);
    } catch (Exception e) {
	GlobalConfig.flux_calculator = FluxCalculator.adaptive;
    }
    GlobalConfig.shear_tolerance = getJSONdouble(jsonData, "shear_tolerance", 0.20);
    GlobalConfig.M_inf = getJSONdouble(jsonData, "M_inf", 0.01);
    GlobalConfig.compression_tolerance = 
	getJSONdouble(jsonData, "compression_tolerance", -0.30);

    // Parameters controlling shock fitting
    GlobalConfig.shock_fitting_delay = getJSONdouble(jsonData, "shock_fitting_delay", 0.0);
    
    GlobalConfig.MHD = getJSONbool(jsonData, "MHD", false);

    // The following checks/overrides must happen after the relevant config elements
    // have been set.
    if (GlobalConfig.grid_motion == GridMotion.shock_fitting &&
	GlobalConfig.apply_bcs_in_parallel) {
	writeln("NOTE: apply_bcs_in_parallel is set to false when shock_fitting is used.");
	GlobalConfig.apply_bcs_in_parallel = false;
    } 

    if (GlobalConfig.verbosity_level > 1) {
	writeln("  interpolation_order: ", GlobalConfig.interpolation_order);
	writeln("  gasdynamic_update_scheme: ",
		gasdynamic_update_scheme_name(GlobalConfig.gasdynamic_update_scheme));
	writeln("  grid_motion: ", grid_motion_name(GlobalConfig.grid_motion));
	writeln("  shock_fitting_delay: ", GlobalConfig.shock_fitting_delay);
	writeln("  write_vertex_velocities: ", GlobalConfig.write_vertex_velocities);
	writeln("  udf_grid_motion_file: ", GlobalConfig.udf_grid_motion_file);
	writeln("  separate_update_for_viscous_terms: ",
		GlobalConfig.separate_update_for_viscous_terms);
	writeln("  separate_update_for_k_omega_source: ",
		GlobalConfig.separate_update_for_k_omega_source);
	writeln("  apply_bcs_in_parallel: ", GlobalConfig.apply_bcs_in_parallel);
	writeln("  stringent_cfl: ", GlobalConfig.stringent_cfl);
	writeln("  adjust_invalid_cell_data: ", GlobalConfig.adjust_invalid_cell_data);
	writeln("  max_invalid_cells: ", GlobalConfig.max_invalid_cells);
	writeln("  thermo_interpolator: ",
		thermo_interpolator_name(GlobalConfig.thermo_interpolator));
	writeln("  apply_limiter: ", GlobalConfig.apply_limiter);
	writeln("  extrema_clipping: ", GlobalConfig.extrema_clipping);
	writeln("  interpolate_in_local_frame: ", GlobalConfig.interpolate_in_local_frame);
	writeln("  flux_calculator: ", flux_calculator_name(GlobalConfig.flux_calculator));
	writeln("  shear_tolerance: ", GlobalConfig.shear_tolerance);
	writeln("  M_inf: ", GlobalConfig.M_inf);
	writeln("  compression_tolerance: ", GlobalConfig.compression_tolerance);
	writeln("  MHD: ", GlobalConfig.MHD);
    }
    if (GlobalConfig.grid_motion != GridMotion.none) {
	if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_1_stage ||
	    GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_1_stage) {
	    // pass, we have a consistent selection.
	} else {
	    throw new Error("Have selected inconsistent GasdynamicUpdate scheme for grid motion.");
	}
    }

    // Parameters controlling viscous/molecular transport
    //
    GlobalConfig.viscous = getJSONbool(jsonData, "viscous", false);
    try {
	string name = jsonData["spatial_deriv_calc"].str;
	GlobalConfig.spatial_deriv_calc = spatial_deriv_calc_from_name(name);
    } catch (Exception e) {
	GlobalConfig.spatial_deriv_calc = SpatialDerivCalc.least_squares;
    }
    GlobalConfig.deriv_calc_at_vertices = getJSONbool(jsonData, "deriv_calc_at_vertices", true);
    GlobalConfig.viscous_delay = getJSONdouble(jsonData, "viscous_delay", 0.0);
    GlobalConfig.viscous_factor_increment = 
	getJSONdouble(jsonData, "viscous_factor_increment", 0.01);
    GlobalConfig.viscous_signal_factor = getJSONdouble(jsonData, "viscous_signal_factor", 1.0);
    try {
	string name = jsonData["turbulence_model"].str;
	GlobalConfig.turbulence_model = turbulence_model_from_name(name);
    } catch (Exception e) {
	GlobalConfig.turbulence_model = TurbulenceModel.none;
    }
    GlobalConfig.turbulence_prandtl_number =
	getJSONdouble(jsonData, "turbulence_prandtl_number", 0.89);
    GlobalConfig.turbulence_schmidt_number =
	getJSONdouble(jsonData, "turbulence_schmidt_number", 0.75);
    GlobalConfig.max_mu_t_factor = getJSONdouble(jsonData, "max_mu_t_factor", 300.0);
    GlobalConfig.transient_mu_t_factor = getJSONdouble(jsonData, "transient_mu_t_factor", 1.0);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  viscous: ", GlobalConfig.viscous);
	writeln("  spatial_deriv_calc: ", spatial_deriv_calc_name(GlobalConfig.spatial_deriv_calc));
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
    GlobalConfig.udf_source_terms = getJSONbool(jsonData, "udf_source_terms", false);
    GlobalConfig.udf_source_terms_file = jsonData["udf_source_terms_file"].str;
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  udf_source_terms: ", GlobalConfig.udf_source_terms);
	writeln("  udf_source_terms_file: ", to!string(GlobalConfig.udf_source_terms_file));
    }

    // Parameters controlling thermochemistry
    //
    GlobalConfig.reacting = getJSONbool(jsonData, "reacting", false);
    GlobalConfig.reactions_file = jsonData["reactions_file"].str;
    GlobalConfig.reaction_time_delay = getJSONdouble(jsonData, "reaction_time_delay", 0.0);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  reacting: ", GlobalConfig.reacting);
	writeln("  reactions_file: ", to!string(GlobalConfig.reactions_file));
	writeln("  reaction_time_delay: ", GlobalConfig.reaction_time_delay);
    }

    // Parameters controlling other simulation options
    //
    GlobalConfig.control_count = getJSONint(jsonData, "control_count", 10);
    GlobalConfig.block_marching = getJSONbool(jsonData, "block_marching", false);
    GlobalConfig.nib = getJSONint(jsonData, "nib", 1);
    GlobalConfig.njb = getJSONint(jsonData, "njb", 1);
    GlobalConfig.nkb = getJSONint(jsonData, "nkb", 1);
    GlobalConfig.propagate_inflow_data = getJSONbool(jsonData, "propagate_inflow_data", false);
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
    // TODO -- still have other entries such as nheatzone, nreactionzone, ...

    // Now, configure blocks that make up the flow domain.
    //
    // This is done in phases.  The blocks need valid references to LocalConfig objects
    // and the boundary conditions need valid references to Sblock objects.
    GlobalConfig.nBlocks = getJSONint(jsonData, "nblock", 0);
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

    GlobalConfig.max_step = getJSONint(jsonData, "max_step", 100);
    GlobalConfig.max_time = getJSONdouble(jsonData, "max_time", 1.0e-3);
    GlobalConfig.halt_now = getJSONint(jsonData, "halt_now", 0);
    GlobalConfig.print_count = getJSONint(jsonData, "print_count", 0);
    GlobalConfig.cfl_count = getJSONint(jsonData, "cfl_count", 0);
    GlobalConfig.dt_init = getJSONdouble(jsonData, "dt_init", 1.0e-6);
    GlobalConfig.dt_max = getJSONdouble(jsonData, "dt_max", 1.0-3);
    GlobalConfig.cfl_value = getJSONdouble(jsonData, "cfl_value", 0.5);
    GlobalConfig.fixed_time_step = getJSONbool(jsonData, "fixed_time_step", false);
    GlobalConfig.dt_reduction_factor = getJSONdouble(jsonData, "dt_reduction_factor", 0.2);
    GlobalConfig.write_at_step = getJSONint(jsonData, "write_at_step", 0);
    GlobalConfig.dt_plot = getJSONdouble(jsonData, "dt_plot", 1.0e-3);
    GlobalConfig.dt_history = getJSONdouble(jsonData, "dt_history", 1.0e-3);

    if (GlobalConfig.verbosity_level > 1) {
	writeln("  max_step: ", GlobalConfig.max_step);
	writeln("  max_time: ", GlobalConfig.max_time);
	writeln("  halt_now: ", GlobalConfig.halt_now);
	writeln("  print_count: ", GlobalConfig.print_count);
	writeln("  cfl_count: ", GlobalConfig.cfl_count);
	writeln("  dt_init: ", GlobalConfig.dt_init);
	writeln("  dt_max: ", GlobalConfig.dt_max);
	writeln("  cfl_value: ", GlobalConfig.cfl_value);
	writeln("  dt_reduction_factor: ", GlobalConfig.dt_reduction_factor);
	writeln("  fixed_time_step: ", GlobalConfig.fixed_time_step);
	writeln("  write_at_step: ", GlobalConfig.write_at_step);
	writeln("  dt_plot: ", GlobalConfig.dt_plot);
	writeln("  dt_history: ", GlobalConfig.dt_history);
    }
} // end read_control_file()
