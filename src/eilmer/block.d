// block.d
// Base class for blocks of cells, for use within Eilmer3.
// Peter J. 2014-07-18 first cut.

module block;

import std.conv;
import std.stdio;
import std.math;
import std.json;
import util.lua;
import geom;
import sgrid;
import gas;
import kinetics;
import globalconfig;
import globaldata;
import fvcore;
import flowstate;
import fvvertex;
import fvinterface;
import fvcell;
import viscousflux;
import bc;
import user_defined_source_terms;

enum
    nghost = 2; // Number of ghost cells surrounding the active cells.


class Block {
public:
    int id; // block identifier: assumed to be the same as the block number.
    Grid_t grid_type; // structured or unstructured
    string label;
    LocalConfig myConfig;
    lua_State* myL;

    bool active; // if true, block participates in the time integration
    // The active flag is used principally for the block-marching calculation,
    // where we want to integrate a few blocks at a time.

    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.
    double mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case
    int hncell;                 // number of sample cells
    int mncell;                 // number of monitor cells
    double[] initial_T_value; // for monitor cells to check against
    FVCell[] cells; // collection of references to be used in foreach statements.
    FVInterface[] faces; // collection of references to all in-use interfaces
    FVVertex[] vertices;
    BoundaryCondition[] bc; // collection of references to the boundary conditions

    // Although the following limits are for the structured grid, 
    // the boundary condition objects need to see them through
    // references to objects of this Block class.
    // They should be unused for unstructured-grid blocks but
    // the SBlock constructor should set appropriate values.
    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax; 
    size_t jmin, jmax;
    size_t kmin, kmax;

    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    FlowState Lft;
    FlowState Rght;

    this(int id, Grid_t grid_type, string label)
    {
	this.id = id;
	this.grid_type = grid_type;
	this.label = label;
	myConfig = dedicatedConfig[id];
	// Workspace for flux_calc method.
	Lft = new FlowState(dedicatedConfig[id].gmodel);
	Rght = new FlowState(dedicatedConfig[id].gmodel);
	// Lua interpreter for the block. 
	// It will be available for computing user-defined source terms, etc.
	myL = luaL_newstate();
	luaL_openlibs(myL);
	lua_pushinteger(myL, id); lua_setglobal(myL, "blkId");
	lua_pushinteger(myL, dedicatedConfig[id].gmodel.n_species);
	lua_setglobal(myL, "n_species");
	lua_pushinteger(myL, dedicatedConfig[id].gmodel.n_modes);
	lua_setglobal(myL, "n_modes");
    }

    ~this()
    {
	lua_close(myL);
    }

    override string toString() const { return "Block(id=" ~ to!string(id) ~ ")"; }

    abstract void init_lua_globals();
    abstract void init_boundary_conditions(JSONValue json_data);
    @nogc abstract ref FVCell get_cell(size_t i, size_t j, size_t k=0);
    @nogc abstract ref FVInterface get_ifi(size_t i, size_t j, size_t k=0);
    @nogc abstract ref FVInterface get_ifj(size_t i, size_t j, size_t k=0);
    @nogc abstract ref FVInterface get_ifk(size_t i, size_t j, size_t k=0);
    @nogc abstract ref FVVertex get_vtx(size_t i, size_t j, size_t k=0);
    abstract void init_grid_and_flow_arrays(string gridFileName);
    abstract void compute_primary_cell_geometric_data(int gtl);
    abstract void compute_distance_to_nearest_wall_for_all_cells(int gtl);
    abstract void read_grid(string filename, size_t gtl=0);
    abstract void write_grid(string filename, double sim_time, size_t gtl=0);
    abstract double read_solution(string filename, bool overwrite_geometry_data);
    abstract void write_solution(string filename, double sim_time);
    abstract void propagate_inflow_data_west_to_east();
    abstract void convective_flux();
    @nogc
    abstract void copy_into_ghost_cells(int destination_face,
					ref Block src_blk, int src_face, int src_orientation,
					int type_of_copy, bool with_encode);

    void identify_reaction_zones(int gtl)
    // Set the reactions-allowed flag for cells in this block.
    {
	size_t total_cells_in_reaction_zones = 0;
	size_t total_cells = 0;
	foreach(cell; cells) {
	    if ( myConfig.reaction_zones.length > 0 ) {
		cell.fr_reactions_allowed = false;
		foreach(rz; myConfig.reaction_zones) {
		    if ( rz.is_inside(cell.pos[gtl], myConfig.dimensions) ) {
			cell.fr_reactions_allowed = true;
		    }
		} // foreach rz
	    } else {
		cell.fr_reactions_allowed = true;
	    }
	    total_cells_in_reaction_zones += (cell.fr_reactions_allowed ? 1: 0);
	    total_cells += 1;
	} // foreach cell
	if ( myConfig.reacting && myConfig.verbosity_level >= 2 ) {
	    writeln("identify_reaction_zones(): block ", id,
		    " cells inside zones = ", total_cells_in_reaction_zones, 
		    " out of ", total_cells);
	    if ( myConfig.reaction_zones.length == 0 ) {
		writeln("Note that for no user-specified zones,",
			" the whole domain is allowed to be reacting.");
	    }
	}
    } // end identify_reaction_zones()

    void identify_turbulent_zones(int gtl)
    // Set the in-turbulent-zone flag for cells in this block.
    {
	size_t total_cells_in_turbulent_zones = 0;
	size_t total_cells = 0;
	foreach(cell; cells) {
	    if ( myConfig.turbulent_zones.length > 0 ) {
		cell.in_turbulent_zone = false;
		foreach(tz; myConfig.turbulent_zones) {
		    if ( tz.is_inside(cell.pos[gtl], myConfig.dimensions) ) {
			cell.in_turbulent_zone = true;
		    }
		} // foreach tz
	    } else {
		cell.in_turbulent_zone = true;
	    }
	    total_cells_in_turbulent_zones += (cell.in_turbulent_zone ? 1: 0);
	    total_cells += 1;
	} // foreach cell
	if ( myConfig.turbulence_model != TurbulenceModel.none && 
	     myConfig.verbosity_level >= 2 ) {
	    writeln("identify_turbulent_zones(): block ", id,
		    " cells inside zones = ", total_cells_in_turbulent_zones, 
		    " out of ", total_cells);
	    if ( myConfig.turbulent_zones.length == 0 ) {
		writeln("Note that for no user-specified zones,",
			" the whole domain is allowed to be turbulent.");
	    }
	}
    } // end identify_turbulent_zones()

    void estimate_turbulence_viscosity()
    {
	final switch (myConfig.turbulence_model) {
	case TurbulenceModel.none:
	    foreach (cell; cells) cell.turbulence_viscosity_zero();
	    return;
	case TurbulenceModel.baldwin_lomax:
	    throw new Error("need to port baldwin_lomax_turbulence_model");
	case TurbulenceModel.spalart_allmaras:
	    throw new Error("Should implement Spalart-Allmaras some day.");
	case TurbulenceModel.k_omega:
	    foreach (cell; cells) cell.turbulence_viscosity_k_omega();
	    break;

	}
	foreach (cell; cells) {
	    cell.turbulence_viscosity_factor(myConfig.transient_mu_t_factor);
	    cell.turbulence_viscosity_limit(myConfig.max_mu_t_factor);
	    cell.turbulence_viscosity_zero_if_not_in_zone();
	}
    } // end estimate_turbulence_viscosity()

    @nogc
    void set_cell_dt_chem(double dt_chem)
    {
	foreach ( cell; cells ) cell.dt_chem = dt_chem;
    }

    @nogc
    void detect_shock_points()
    // Detects shocks by looking for compression between adjacent cells.
    //
    // The velocity component normal to the cell interfaces
    // is used as the indicating variable.
    {
	// Change in normalised velocity to indicate a shock.
	// A value of -0.05 has been found suitable to detect the levels of
	// shock compression observed in the "sod" and "cone20" test cases.
	// It may need to be tuned for other situations, especially when
	// viscous effects are important.
	double tol = myConfig.compression_tolerance;
	// First, work across interfaces and locate shocks using the (local) normal velocity.
	foreach (iface; faces) {
	    auto cL = iface.left_cells[0];
	    auto cR = iface.right_cells[0];
	    double uL = cL.fs.vel.x * iface.n.x + cL.fs.vel.y * iface.n.y + cL.fs.vel.z * iface.n.z;
	    double uR = cR.fs.vel.x * iface.n.x + cR.fs.vel.y * iface.n.y + cR.fs.vel.z * iface.n.z;
	    double aL = cL.fs.gas.a;
	    double aR = cR.fs.gas.a;
	    double a_min = fmin(aL, aR);
	    iface.fs.S = ((uR - uL) / a_min < tol);
	}
	// Finally, mark cells as shock points if any of their interfaces are shock points.
	foreach (cell; cells) {
	    cell.fs.S = false;
	    foreach (face; cell.iface) {
		if (face.fs.S) cell.fs.S = true;
	    }
	}
    } // end detect_shock_points()

    int count_invalid_cells(int gtl)
    // Returns the number of cells that contain invalid data,
    // optionally patching bad cell data as it goes.
    //
    // This data can be identified by the density of internal energy 
    // being on the minimum limit or the velocity being very large.
    {
	int number_of_invalid_cells = 0;
	foreach(cell; cells) {
	    if ( cell.fs.check_data(cell.pos[0]) == false ) {
		++number_of_invalid_cells;
		writefln("count_invalid_cells: block_id = %d, cell at %g,%g,%g\n",
			 id, cell.pos[0].x, cell.pos[0].y, cell.pos[0].z);
		writeln(cell);
		if ( myConfig.adjust_invalid_cell_data ) {
		    // We shall set the cell data to something that
		    // is valid (and self consistent).
		    FlowState[] neighbour_flows;
		    printf( "Adjusting cell data to a local average.\n" );
		    foreach (i; 0 .. cell.iface.length) {
			auto face = cell.iface[i];
			auto other_cell = (cell.outsign[i] == 1) ? face.right_cells[0] : face.left_cells[0];
			if (other_cell.fs.check_data(other_cell.pos[0])) neighbour_flows ~= other_cell.fs;
		    }
		    if (neighbour_flows.length == 0) {
			throw new Error(text("Block::count_invalid_cells(): " ~
					     "There were no valid neighbours " ~
					     "to replace flow data in cell."));
		    }
		    cell.fs.copy_average_values_from(neighbour_flows, myConfig.gmodel);
		    cell.encode_conserved(gtl, 0, omegaz);
		    cell.decode_conserved(gtl, 0, omegaz);
		    writefln("after flow-data replacement: block_id = %d, cell @ %g,%g,%g\n",
			     id, cell.pos[0].x, cell.pos[0].y, cell.pos[0].z);
		    writeln(cell);
		} // end adjust_invalid_cell_data 
	    } // end of if invalid data...
	} // foreach cell
	return number_of_invalid_cells;
    } // end count_invalid_cells()

    @nogc
    void flow_property_derivatives(int gtl)
    {
	if (myConfig.dimensions == 2) {
	    final switch ( myConfig.spatial_deriv_calc ) {
	    case SpatialDerivCalc.least_squares:
		foreach(vtx; vertices) { gradients_xy_leastsq(vtx, myConfig.diffusion); }
		break;
	    case SpatialDerivCalc.divergence:
		foreach(vtx; vertices) { gradients_xy_div(vtx, myConfig.diffusion); }
	    } // end switch
	} else {
	    foreach(vtx; vertices) { gradients_xyz_leastsq(vtx, myConfig.diffusion); }
	} // end if (myConfig.dimensions
    } // end flow_property_derivatives()

    @nogc
    void clear_fluxes_of_conserved_quantities()
    {
	foreach (iface; faces) iface.F.clear_values();
    }

    void viscous_flux()
    {
	auto vfwork = new ViscousFluxData(myConfig);
	foreach (iface; faces) {
	    vfwork.average_vertex_values(iface);
	    vfwork.viscous_flux_calc(iface);
	}
    }

    @nogc
    void init_residuals()
    // Initialization of data for later computing residuals.
    {
	mass_residual = 0.0;
	mass_residual_loc = Vector3(0.0, 0.0, 0.0);
	energy_residual = 0.0;
	energy_residual_loc = Vector3(0.0, 0.0, 0.0);
	foreach(FVCell cell; cells) {
	    cell.rho_at_start_of_step = cell.fs.gas.rho;
	    cell.rE_at_start_of_step = cell.U[0].total_energy;
	}
    } // end init_residuals()

    @nogc
    void compute_residuals(int gtl)
    // Compute the residuals using previously stored data.
    //
    // The largest residual of density for all cells was the traditional way
    // mbcns/Elmer estimated the approach to steady state.
    // However, with the splitting up of the increments for different physical
    // processes, this residual calculation code needed a bit of an update.
    // Noting that the viscous-stress, chemical and radiation increments
    // do not affect the mass within a cell, we now compute the residuals 
    // for both mass and (total) energy for all cells, the record the largest
    // with their location. 
    {
	mass_residual = 0.0;
	mass_residual_loc = Vector3(0.0, 0.0, 0.0);
	energy_residual = 0.0;
	energy_residual_loc = Vector3(0.0, 0.0, 0.0);
	foreach(FVCell cell; cells) {
	    double local_residual = (cell.fs.gas.rho - cell.rho_at_start_of_step) 
		/ cell.fs.gas.rho;
	    local_residual = fabs(local_residual);
	    if ( local_residual > mass_residual ) {
		mass_residual = local_residual;
		mass_residual_loc.refx = cell.pos[gtl].x;
		mass_residual_loc.refy = cell.pos[gtl].y;
		mass_residual_loc.refz = cell.pos[gtl].z;
	    }
	    // In the following line, the zero index is used because,
	    // at the end of the gas-dynamic update, that index holds
	    // the updated data.
	    local_residual = (cell.U[0].total_energy - cell.rE_at_start_of_step) 
		/ cell.U[0].total_energy;
	    local_residual = fabs(local_residual);
	    if ( local_residual > energy_residual ) {
		energy_residual = local_residual;
		energy_residual_loc.refx = cell.pos[gtl].x;
		energy_residual_loc.refy = cell.pos[gtl].y;
		energy_residual_loc.refz = cell.pos[gtl].z;
	    }
	} // for cell
    } // end compute_residuals()

    double determine_time_step_size(double dt_current)
    // Compute the local time step limit for all cells in the block.
    // The overall time step is limited by the worst-case cell.
    {
	double cfl_value = GlobalConfig.cfl_value;
	double dt_local;
	double cfl_local;
	double signal;
	double cfl_allow; // allowable CFL number, t_order dependent
	double dt_allow;
	double cfl_min, cfl_max;
	// The following limits allow the simulation of the sod shock tube
	// to get just a little wobbly around the shock.
	// Lower values of cfl should be used for a smooth solution.
	switch (number_of_stages_for_update_scheme(myConfig.gasdynamic_update_scheme)) {
	case 1: cfl_allow = 0.9; break;
	case 2: cfl_allow = 1.2; break;
	case 3: cfl_allow = 1.6; break;
	default: cfl_allow = 0.9;
	}
	bool first = true;
	foreach(FVCell cell; cells) {
	    signal = cell.signal_frequency();
	    cfl_local = dt_current * signal; // Current (Local) CFL number
	    dt_local = cfl_value / signal; // Recommend a time step size.
	    if ( first ) {
		cfl_min = cfl_local;
		cfl_max = cfl_local;
		dt_allow = dt_local;
		first = false;
	    } else {
		cfl_min = fmin(cfl_min, cfl_local);
		cfl_max = fmax(cfl_max, cfl_local);
		dt_allow = fmin(dt_allow, dt_local);
	    }
	} // foreach cell
	if ( cfl_max < 0.0 || cfl_max > cfl_allow ) {
	    writeln( "determine_time_step_size(): bad CFL number was encountered");
	    writeln( "    cfl_max=", cfl_max, " for Block ", id);
	    writeln( "    If this cfl_max value is not much larger than 1.0,");
	    writeln( "    your simulation could probably be restarted successfully");
	    writeln( "    with some minor tweaking.");
	    writeln( "    That tweaking should probably include a reduction");
	    writeln( "    in the size of the initial time-step, dt");
	    writeln( "    If this job is a restart/continuation of an old job, look in");
	    writeln( "    the old-job.finish file for the value of dt at termination.");
	    throw new Error(text("Bad cfl number encountered cfl_max=", cfl_max));
	}
	return dt_allow;
    } // end determine_time_step_size()

    void write_history(string filename, double sim_time, bool write_header=false)
    {
	throw new Error("[TODO] Block.write_history() not implemented yet.");
    }

    void applyPreReconAction(double t, int gtl, int ftl)
    {
	foreach(boundary; bc) { boundary.applyPreReconAction(t, gtl, ftl); }
    }

    void applyPostConvFluxAction(double t, int gtl, int ftl)
    {
	foreach(boundary; bc) { boundary.applyPostConvFluxAction(t, gtl, ftl); }
    }

    void applyPreSpatialDerivAction(double t, int gtl, int ftl)
    {
	foreach(boundary; bc) { boundary.applyPreSpatialDerivAction(t, gtl, ftl); }
    }

    void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
	foreach(boundary; bc) { boundary.applyPostDiffFluxAction(t, gtl, ftl); }
    }
} // end class Block
