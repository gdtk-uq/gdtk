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
import ssolidblock;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;

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
    GlobalConfig.dimensions = getJSONint(jsonData, "dimensions", 2);
    GlobalConfig.axisymmetric = getJSONbool(jsonData, "axisymmetric", false);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  title: ", to!string(GlobalConfig.title));
	writeln("  gas_model_file: ", to!string(GlobalConfig.gas_model_file));
	writeln("  dimensions: ", GlobalConfig.dimensions);
	writeln("  axisymmetric: ", GlobalConfig.axisymmetric);
    }

    // Parameters controlling convective update
    //
    GlobalConfig.interpolation_order = getJSONint(jsonData, "interpolation_order", 2);
    try {
	string name = jsonData["gasdynamic_update_scheme"].str;
	GlobalConfig.gasdynamic_update_scheme = update_scheme_from_name(name);
    } catch (Exception e) {
	GlobalConfig.gasdynamic_update_scheme = GasdynamicUpdate.pc;
    }
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
    GlobalConfig.moving_grid = getJSONbool(jsonData, "moving_grid", false);
    GlobalConfig.write_vertex_velocities = 
	getJSONbool(jsonData, "write_vertex_velocities", false);
	GlobalConfig.MHD = getJSONbool(jsonData, "MHD", false);

    if (GlobalConfig.verbosity_level > 1) {
	writeln("  interpolation_order: ", GlobalConfig.interpolation_order);
	writeln("  gasdynamic_update_scheme: ",
		gasdynamic_update_scheme_name(GlobalConfig.gasdynamic_update_scheme));
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
	writeln("  moving_grid: ", GlobalConfig.moving_grid);
	writeln("  write_vertex_velocities: ", GlobalConfig.write_vertex_velocities);
    }

    // Parameters controlling viscous/molecular transport
    //
    GlobalConfig.viscous = getJSONbool(jsonData, "viscous", false);
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
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  reacting: ", GlobalConfig.reacting);
	writeln("  reactions_file: ", to!string(GlobalConfig.reactions_file));
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
	gasBlocks ~= new SBlock(i, jsonData["block_" ~ to!string(i)]);
	if (GlobalConfig.verbosity_level > 1) {
	    writeln("  Block[", i, "]: ", gasBlocks[i]);
	}
    }
    foreach (blk; gasBlocks) {
	blk.init_lua_globals();
	blk.init_boundary_conditions(jsonData["block_" ~ to!string(blk.id)]);
	if (GlobalConfig.udf_source_terms) {
	    luaL_dofile(blk.myL, GlobalConfig.udf_source_terms_file.toStringz);
	}
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
