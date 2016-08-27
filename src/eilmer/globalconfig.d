/** globalconfig.d
 * A place to keep the configuration details for the simulation.
 *
 * Author: Peter J. and Rowan G.
 * First code: 2014-07-18
 * Resumed work: 2015-02-05
 */

module globalconfig;

import std.conv;
import std.stdio;
import std.string;
import std.typecons;
import core.stdc.stdlib : exit;
import std.json;
import std.file;
import std.conv;
import std.array;
import std.format;

import util.lua;
import gas;
import kinetics;
import geom;
import fvcore;
version (gpu_chem) {
    import gpu_chem;
}
import json_helper;
import globaldata;
import flowstate;
import sblock;
import ublock;
import ssolidblock;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import grid_motion;

// Symbolic names for turbulence models.
enum TurbulenceModel { none, baldwin_lomax, k_omega, spalart_allmaras }
string turbulence_model_name(TurbulenceModel i)
{
    final switch (i) {
    case TurbulenceModel.none: return "none";
    case TurbulenceModel.baldwin_lomax: return "baldwin_lomax";
    case TurbulenceModel.k_omega: return "k_omega";
    case TurbulenceModel.spalart_allmaras: return "spalart_allmaras";
    }
} // end turbulence_model_name()

TurbulenceModel turbulence_model_from_name(string name)
{
    switch (name) {
    case "none": return TurbulenceModel.none;
    case "baldwin_lomax": return TurbulenceModel.baldwin_lomax;
    case "k_omega": return TurbulenceModel.k_omega;
    case "spalart_allmaras": return TurbulenceModel.spalart_allmaras;
    default: return TurbulenceModel.none;
    }
} // end turbulence_model_from_name()

class BlockZone {
    // Note that these data are not shared across threads
    // because we will want to access them frequently at the lower levels of the code.
    Vector3 p0, p1;
    this(in Vector3 p0, in Vector3 p1) {
	this.p0 = p0;
	this.p1 = p1;
    }
    this(const BlockZone other)
    {
	p0 = other.p0;
	p1 = other.p1;
    }
    override string toString() {
	return text("BlockZone(p0=", to!string(p0), ", p1=", to!string(p1), ")");
    }
    bool is_inside(in Vector3 p, int dimensions) {
	if ( p.x >= p0.x && p.x <= p1.x &&
	     p.y >= p0.y && p.y <= p1.y ) {
	    if ( dimensions == 2 ) {
		return true;
	    } else if ( p.z >= p0.z && p.z <= p1.z ) {
		return true;
	    }
	}
	return false;
    } // end is_inside()
}

class IgnitionZone : BlockZone {
    double Tig; // temperature to apply within reaction_update ensure ignition
    this(in Vector3 p0, in Vector3 p1, double Tig) {
	super(p0, p1);
	this.Tig = Tig;
    }
    this(const IgnitionZone other)
    {
	super(other.p0, other.p1);
	Tig = other.Tig;
    }
    override string toString() {
	return text("IgnitionZone(p0=", to!string(p0), ", p1=", to!string(p1), 
		    ", Tig=", Tig, ")");
    }
}


final class GlobalConfig {
    shared static string base_file_name = "job"; // Change this to suit at run time.
    shared static string title = "Eilmer4 simulation"; // Change this to suit at run time.
    shared static string gas_model_file = "gas-model.lua";
    // Note that the following reference to the GasModel is NOT shared.
    // Do not use this one from within a thread.
    static GasModel gmodel_master; // to provide GasModel access to the Lua domain
    // Presumably, we won't be accessing this particular gas model from the 
    // individual block computations, so that parallel computations for the blocks
    // don't trip over each other.
    shared static bool include_quality = false; // if true, we include quality in the solution file  

    shared static int nBlocks = 0; // Number of blocks in the overall simulation.
    shared static int nSolidBlocks = 0; // Number of solid blocks in the overall simulation.
    shared static int dimensions = 2; // default is 2, other valid option is 3
    shared static bool axisymmetric = false;

    // Parameters controlling convective update
    shared static GasdynamicUpdate gasdynamic_update_scheme = GasdynamicUpdate.pc;
    shared static size_t n_flow_time_levels = 3;

    // Parameters related to possible motion of the grid.
    shared static grid_motion = GridMotion.none;
    shared static bool write_vertex_velocities = false;
    shared static string udf_grid_motion_file = "dummy-grid-motion-file.txt";
    static lua_State* master_lua_State;
    shared static size_t n_grid_time_levels = 1;

    // We might update some properties in with the main convective-terms
    // time-stepping function or we might choose to update them separately, 
    // like the chemistry update.
    shared static bool separate_update_for_viscous_terms = false;
    shared static bool separate_update_for_k_omega_source = false;

    // Some of the user-defined functionality depends on having access to all blocks
    // from a single thread.  For safety, in those cases, do not use parallel loops. 
    shared static bool apply_bcs_in_parallel = true;
    
    /// When decoding the array of conserved quantities, 
    /// the temperature or the density may try to go negative.  
    /// If it does and adjust_invalid_cell_data == true, the cell data
    /// is adjusted to make it reasonable.
    shared static bool adjust_invalid_cell_data = false;
    // The maximum number of bad cells (per block) 
    // which will be tolerated without complaint.
    shared static int max_invalid_cells = 0;
    shared static FlowStateLimits flowstate_limits;

    // Low order reconstruction (1) uses just the cell-centre data as left- and right-
    // flow properties in the flux calculation.
    // High-order reconstruction (2) adds a correction term to the cell-centre values
    // to approach something like a piecewise-quadratic interpolation between the
    // cell centres for structured-grids and a linear model across a cloud of cell 
    // centres for unstructured-grids.
    shared static int interpolation_order = 2; 
    // Default flow-data reconstruction includes interpolation of density 
    // and internal energy.  Other options for the thermodunamic properties
    // to be interpolated are pressure+temperature, density+temperature and
    // density+pressure.
    shared static InterpolateOption thermo_interpolator = InterpolateOption.rhoe;
    shared static bool apply_limiter = true;
    shared static bool extrema_clipping = true;
    shared static bool interpolate_in_local_frame = true; // only for structured-grid
    // The unstructured solver has a selection of limiters available
    shared static UnstructuredLimiter unstructured_limiter = UnstructuredLimiter.van_albada;

    // Default flux calculator is the adaptive mix of ausmdv and efm.
    shared static FluxCalculator flux_calculator = FluxCalculator.adaptive;

    // Set the tolerance to shear when applying the adaptive flux calculator.
    // We don't want EFM to be applied in situations of significant shear.
    // The shear value is computed as the tangential-velocity difference across an interface
    // normalised by the local sound speed.
    shared static double shear_tolerance = 0.20;

    // Reference free-stream Mach number, for use in the ausm_plus_up flux calculator.
    // Choose a value for M_inf that is good for low Mach numbers.
    // To be strictly correct, we should set this at run time
    // if an M_inf value is easily defined.
    shared static double M_inf = 0.01;

    // Set the tolerance in relative velocity change for the shock detector.
    // This value is expected to be a negative number (for compression)
    // and not too large in magnitude.
    // We have been using a value of -0.05 for years, based on some
    // early experiments with the sod and cone20 test cases, however,
    // the values may need to be tuned for other cases, especially where
    // viscous effects are important.
    shared static double compression_tolerance = -0.30;

    // With this flag on, finite-rate evolution of the vibrational energies 
    // (and in turn the total energy) is computed.
    shared static bool thermal_energy_exchange = false;

    // For including radiation energy exchange.
    //
    shared static bool radiation = false;
    shared static int radiation_update_frequency; // = 1 for every time-step
    shared static bool radiation_scaling = false;
    shared static bool halt_on_large_flow_change = false; 
    // Set to true to halt simulation when any monitor point sees a large flow change.
    shared static double tolerance_in_T;   // Temperature change for the flow change.

    shared static bool electric_field_work;

    // For Daryl Bond and Vince Wheatley's Single-fluid MHD additions.
    //
    shared static bool MHD = false;
    // Lachlan Whyborn's Divergence cleaning to go with MHD.
    shared static bool divergence_cleaning = false;
    shared static double c_h = 0.0;
    shared static double divB_damping_length = 1.0;

    // Parameters controlling viscous/molecular transport
    //
    shared static bool viscous = false; 
    // If true, viscous effects are included in the gas-dynamic update.
    shared static SpatialDerivCalc spatial_deriv_calc = SpatialDerivCalc.least_squares;
    shared static SpatialDerivLocn spatial_deriv_locn = SpatialDerivLocn.faces;
    shared static bool include_ghost_cells_in_spatial_deriv_clouds = true;
    //
    // A factor to scale the viscosity in order to achieve a soft start. 
    // The soft-start for viscous effects may be handy for impulsively-started flows.
    // A value of 1.0 means that the viscous effects are fully applied.
    shared static double viscous_factor = 1.0;
    // The amount by which to increment the viscous factor during soft-start.
    shared static double viscous_factor_increment = 0.01;
    shared static double viscous_delay = 0.0;
    // The amount of time by which to delay the shock fitting.
    shared static double shock_fitting_delay = 1.5e-3;
    // order of the special interpolation applied at the shock fitting inflow boundary
    shared static int shock_fitting_interpolation_order = 1;
    // scaling factor applied to vertices in shock fitting simulations for stability
    shared static double shock_fitting_scale_factor = 0.5;
    // When the diffusion is calculated is treated as part of the viscous calculation:
    //   false for neglecting multicomponent diffusion, 
    //   true when considering the diffusion 
    shared static bool diffusion = false; 
    // A factor to scale the diffusion in order to achieve a soft start.
    // The soft-start for diffusion effects may be handy for impulsively-started flows.
    // Note that this is separate to viscous effects.
    shared static double diffusion_factor = 1.0;
    // The amount by which to increment the diffusion factor during soft-start.
    shared static double diffusion_factor_increment = 0.01;
    shared static double diffusion_time_delay = 0.0;
    // The Lewis number when using the constant Lewis number diffusion model
    shared static double diffusion_lewis = 1.0;
    // The Schmidt number when using the constant Schmidt number diffusion model
    shared static double diffusion_schmidt = 0.7;

    shared static TurbulenceModel turbulence_model = TurbulenceModel.none;
    shared static double turbulence_prandtl_number = 0.89;
    shared static double turbulence_schmidt_number = 0.75;
    shared static double max_mu_t_factor = 300.0;
    shared static double transient_mu_t_factor = 1.0;
    static BlockZone[] turbulent_zones;

    // Indicate presence of user-defined source terms
    shared static string udf_source_terms_file = "dummy-source-terms.txt";
    shared static bool udf_source_terms = false;

    // Parameters controlling thermochemistry
    //
    // Turning on the reactions activates the chemical update function calls.
    // Chemical equilibrium simulations (via Look-Up Table) does not use this
    // chemical update function call.
    shared static bool reacting = false;
    shared static string reactions_file = "chemistry.lua";
    shared static double reaction_time_delay = 0.0;
    shared static double T_frozen; // temperature (in K) below which reactions are frozen
    shared static double T_frozen_energy; // temperature (in K) below which energy exchanges are skipped
    static BlockZone[] reaction_zones;
    shared static double ignition_time_start = 0.0;
    shared static double ignition_time_stop = 0.0;
    static IgnitionZone[] ignition_zones;
    shared static bool ignition_zone_active = false;

    // Parameters controlling other simulation options
    //
    shared static int max_step = 100;      // iteration limit
    shared static int t_level;             // time level within update
    shared static int halt_now = 0;        // flag for premature halt
    shared static int print_count = 20; // Number of steps between writing messages to console.
    shared static int control_count = 10; // Number of steps between rereading .control file.

    shared static int verbosity_level = 1;
    // Messages have a hierarchy:  // [TODO] we are not really abiding by this.
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose boundary condition messages

    shared static double max_time = 1.0e-3; // final solution time, s, set by user
    shared static double dt_init = 1.0e-6; // initial time step set by user
    shared static double dt_max = 1.0e-3; // Maximum allowable time-step, after all other considerations.
    shared static double cfl_value = 0.5; // target CFL number (worst case) set by user
    shared static bool stringent_cfl = false; 
    // If true, assume the worst with respect to cell geometry and wave speed.
    shared static double viscous_signal_factor = 1.0; // can reduce the viscous influence in CFL condition
    shared static int cfl_count = 10;  // steps between checking time step size
    shared static bool fixed_time_step = false; // set true to fix dt_allow

    shared static double dt_plot = 1.0e-3; // interval for writing soln
    shared static double dt_history = 1.0e-3; // interval for writing sample
    shared static Tuple!(size_t, size_t)[] hcells;
    shared static Tuple!(size_t, size_t)[] solid_hcells;
    
    shared static double energy_residual;      // to be monitored for steady state
    shared static Vector3 energy_residual_loc; // location of largest value
    shared static double mass_residual;
    shared static Vector3 mass_residual_loc;

    // Block-marching parameters
    shared static bool block_marching = false;
    shared static int nib = 1;
    shared static int njb = 1;
    shared static int nkb = 1;
    shared static bool propagate_inflow_data = false;

    // Parameters related to the solid domain and conjugate coupling
    shared static bool udfSolidSourceTerms = false;
    shared static string udfSolidSourceTermsFile = "dummy-solid-source-terms.txt";

    // Parameters related to the gpu chemistry mode
    version (gpu_chem) {
	static GPUChem gpuChem;
    }

    ~this()
    {
	lua_close(master_lua_State);
    }

} // end class GlobalConfig


// A place to store configuration parameters that are in memory that can be 
// dedicated to a thread.   We will want to access them frequently 
// at the lower levels of the code without having to guard them with memory barriers.
class LocalConfig {
public:
    int dimensions;
    bool axisymmetric;
    GasdynamicUpdate gasdynamic_update_scheme;
    size_t n_flow_time_levels;
    GridMotion grid_motion;
    string udf_grid_motion_file;
    size_t n_grid_time_levels;
    bool separate_update_for_viscous_terms;
    bool separate_update_for_k_omega_source;
    bool adjust_invalid_cell_data;
    FlowStateLimits flowstate_limits;
    int interpolation_order;
    InterpolateOption thermo_interpolator;
    bool apply_limiter;
    bool extrema_clipping;
    bool interpolate_in_local_frame;
    UnstructuredLimiter unstructured_limiter;
    FluxCalculator flux_calculator;
    double shear_tolerance;
    double M_inf;
    double compression_tolerance;

    bool radiation;
    bool electric_field_work;
    bool MHD;
    bool divergence_cleaning;
    double c_h;
    double divB_damping_length;

    bool viscous;
    SpatialDerivCalc spatial_deriv_calc;
    SpatialDerivLocn spatial_deriv_locn;
    bool include_ghost_cells_in_spatial_deriv_clouds;
    double viscous_factor;
    bool diffusion;
    double diffusion_factor;
    double diffusion_lewis;
    double diffusion_schmidt;

    bool stringent_cfl;
    double viscous_signal_factor;

    int shock_fitting_interpolation_order;
    double shock_fitting_scale_factor;
    
    TurbulenceModel turbulence_model;
    double turbulence_prandtl_number;
    double turbulence_schmidt_number;
    double max_mu_t_factor;
    double transient_mu_t_factor;
    BlockZone[] turbulent_zones;

    bool udf_source_terms;

    bool reacting;
    double reaction_time_delay;
    double T_frozen;
    double T_frozen_energy;
    BlockZone[] reaction_zones;

    double ignition_time_start;
    double ignition_time_stop;
    IgnitionZone[] ignition_zones;
    bool ignition_zone_active;

    GasModel gmodel;
    bool include_quality;
    ChemistryUpdate chemUpdate;

    int verbosity_level;

    this() 
    {
	dimensions = GlobalConfig.dimensions;
	axisymmetric = GlobalConfig.axisymmetric;
	gasdynamic_update_scheme = GlobalConfig.gasdynamic_update_scheme;
	n_flow_time_levels = GlobalConfig.n_flow_time_levels;
	grid_motion = GlobalConfig.grid_motion;
	udf_grid_motion_file = GlobalConfig.udf_grid_motion_file;
	n_grid_time_levels = GlobalConfig.n_grid_time_levels;
	separate_update_for_viscous_terms = GlobalConfig.separate_update_for_viscous_terms;
	separate_update_for_k_omega_source = GlobalConfig.separate_update_for_k_omega_source;
	adjust_invalid_cell_data = GlobalConfig.adjust_invalid_cell_data;
	flowstate_limits = GlobalConfig.flowstate_limits;
	interpolation_order = GlobalConfig.interpolation_order;
	thermo_interpolator = GlobalConfig.thermo_interpolator;
	apply_limiter = GlobalConfig.apply_limiter;
	extrema_clipping = GlobalConfig.extrema_clipping;
	interpolate_in_local_frame = GlobalConfig.interpolate_in_local_frame;
	unstructured_limiter = GlobalConfig.unstructured_limiter;
	flux_calculator = GlobalConfig.flux_calculator;
	shear_tolerance = GlobalConfig.shear_tolerance;
	M_inf = GlobalConfig.M_inf;
	compression_tolerance = GlobalConfig.compression_tolerance;

	radiation = GlobalConfig.radiation;
	electric_field_work = GlobalConfig.electric_field_work;
	MHD = GlobalConfig.MHD;
	divergence_cleaning = GlobalConfig.divergence_cleaning;
	c_h = GlobalConfig.c_h;
	divB_damping_length = GlobalConfig.divB_damping_length;

	viscous = GlobalConfig.viscous;
	spatial_deriv_calc = GlobalConfig.spatial_deriv_calc;
	spatial_deriv_locn = GlobalConfig.spatial_deriv_locn;
	include_ghost_cells_in_spatial_deriv_clouds = 
	    GlobalConfig.include_ghost_cells_in_spatial_deriv_clouds;
	viscous_factor = GlobalConfig.viscous_factor;
	diffusion = GlobalConfig.diffusion;
	diffusion_factor = GlobalConfig.diffusion_factor;
	diffusion_lewis = GlobalConfig.diffusion_lewis;
	diffusion_schmidt = GlobalConfig.diffusion_schmidt;

	stringent_cfl = GlobalConfig.stringent_cfl;
	viscous_signal_factor = GlobalConfig.viscous_signal_factor;

	shock_fitting_interpolation_order = GlobalConfig.shock_fitting_interpolation_order;
	shock_fitting_scale_factor = GlobalConfig.shock_fitting_scale_factor;
	
	turbulence_model = GlobalConfig.turbulence_model;
	turbulence_prandtl_number = GlobalConfig.turbulence_prandtl_number;
	turbulence_schmidt_number = GlobalConfig.turbulence_schmidt_number;
	max_mu_t_factor = GlobalConfig.max_mu_t_factor;
	transient_mu_t_factor = GlobalConfig.transient_mu_t_factor;
	foreach (bz; GlobalConfig.turbulent_zones) turbulent_zones ~= new BlockZone(bz);

	udf_source_terms = GlobalConfig.udf_source_terms;

	reacting = GlobalConfig.reacting;
	reaction_time_delay = GlobalConfig.reaction_time_delay;
	T_frozen = GlobalConfig.T_frozen;
	T_frozen_energy = GlobalConfig.T_frozen_energy;
	foreach (rz; GlobalConfig.reaction_zones) reaction_zones ~= new BlockZone(rz);

	ignition_time_start = GlobalConfig.ignition_time_start;
	ignition_time_stop = GlobalConfig.ignition_time_stop;
	ignition_zone_active = GlobalConfig.ignition_zone_active;
	foreach (iz; GlobalConfig.ignition_zones) ignition_zones ~= new IgnitionZone(iz);

	gmodel = init_gas_model(GlobalConfig.gas_model_file);
	include_quality = GlobalConfig.include_quality;
	if (GlobalConfig.reacting)
	    chemUpdate = new ChemistryUpdate(GlobalConfig.reactions_file, gmodel);

	verbosity_level = GlobalConfig.verbosity_level;
    }
} // end class LocalConfig

//-----------------------------------------------------------------------------------
// Reading from JSON file.

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
    mixin(update_enum("unstructured_limiter", "unstructured_limiter", "unstructured_limiter_from_name"));
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
    // have been set.  This is the first such check.  For details, see the function below.
    configCheckPoint1();

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
	writeln("  unstructured_limiter: ", unstructured_limiter_name(GlobalConfig.unstructured_limiter));
	writeln("  extrema_clipping: ", GlobalConfig.extrema_clipping);
	writeln("  interpolate_in_local_frame: ", GlobalConfig.interpolate_in_local_frame);
	writeln("  flux_calculator: ", flux_calculator_name(GlobalConfig.flux_calculator));
	writeln("  shear_tolerance: ", GlobalConfig.shear_tolerance);
	writeln("  M_inf: ", GlobalConfig.M_inf);
	writeln("  compression_tolerance: ", GlobalConfig.compression_tolerance);
	writeln("  MHD: ", GlobalConfig.MHD);
	writeln("  divergence_cleaning: ", GlobalConfig.divergence_cleaning);
	writeln("  divB_damping_length: ", GlobalConfig.divB_damping_length);
    }
    configCheckPoint2();

    // Parameters controlling viscous/molecular transport
    //
    mixin(update_bool("viscous", "viscous"));
    mixin(update_enum("spatial_deriv_calc", "spatial_deriv_calc", "spatial_deriv_calc_from_name"));
    mixin(update_enum("spatial_deriv_locn", "spatial_deriv_locn", "spatial_deriv_locn_from_name"));
    mixin(update_bool("include_ghost_cells_in_spatial_deriv_clouds", "include_ghost_cells_in_spatial_deriv_clouds"));
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

//
// Functions to check the compatibility of the GlobalConfig parameters.
// The individual checks are called at particular locations during the
// reading of the config parameters from the JSON .config file.
//
// Sometimes, we can make a suitable change and continue (after emitting a warning)
// but other times it is not sensible to continue, but throw an exception.
//
void configCheckPoint1()
{
    if (GlobalConfig.grid_motion == GridMotion.shock_fitting &&
	GlobalConfig.apply_bcs_in_parallel) {
	writeln("NOTE: apply_bcs_in_parallel is set to false when shock_fitting is used.");
	GlobalConfig.apply_bcs_in_parallel = false;
    }
    return;
} // end configCheckPoint1()

void configCheckPoint2()
{
    // More checking of constraints on the config parameters.
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
    return;
} // end configCheckPoint2()

void checkGlobalConfig()
// Bundle all of the checks together so that they may be conveniently applied
// at the end of processing the user's Lua input script.
{
    configCheckPoint1();
    configCheckPoint2();
    return;
}
