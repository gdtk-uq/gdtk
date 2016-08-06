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

import util.lua;
import gas;
import kinetics;
import geom;
import fvcore;
version (gpu_chem) {
    import gpu_chem;
}

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
    shared static bool spatial_deriv_retain_lsq_work_data = true;
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
    bool spatial_deriv_retain_lsq_work_data;
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
	spatial_deriv_retain_lsq_work_data = GlobalConfig.spatial_deriv_retain_lsq_work_data;
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
