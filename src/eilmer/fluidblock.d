// fluidblock.d
// Base class for blocks of cells, for use within the Eilmer flow solver.
// Peter J. 2014-07-18 first cut.

module fluidblock;

import std.conv;
import std.stdio;
import std.math;
import std.json;
import util.lua;
import nm.complex;
import nm.number;
import nm.bbla;
import nm.smla;
import util.lua_service;
import gas.luagas_model;
import geom;
import gas;
import kinetics;
import globalconfig;
import globaldata;
import fvcore;
import flowstate;
import fvvertex;
import fvinterface;
import fvcell;
import flowgradients;
import bc;
import user_defined_source_terms;
import conservedquantities;
import lua_helper;
import grid_motion;
import grid_deform;
import sfluidblock; // needed for some special-case processing, below

// To distinguish ghost cells from active cells, we start their id values at
// an arbitrarily high value.  It seem high to me (PJ) but feel free to adjust it
// if you start using grids larger I expect.
enum ghost_cell_start_id = 1_000_000_000;


// The flow solver handles structured- and unstructured-grid blocks via this base class.
// Mostly, we view the block as an unstructured bag of cells because that requires least
// knowledge in the calling code.
class FluidBlock {
public:
    int id; // block identifier: assumed to be the same as the block number.
    Grid_t grid_type; // structured or unstructured
    string label;
    LocalConfig myConfig;
    lua_State* myL;
    //
    bool active; // if true, block participates in the time integration
    // The active flag is used principally for the block-marching calculation,
    // where we want to integrate a few blocks at a time.
    //
    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.
    number mass_balance; // domain mass balance used to monitor for steady state
    number L2_residual; // L2 norm of the global residual
    number mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case
    ConservedQuantities Linf_residuals;
    number c_h, divB_damping_length; //divergence cleaning parameters for MHD
    int mncell;                 // number of monitor cells
    number[] initial_T_value; // for monitor cells to check against
    //
    // Collections of cells, vertices and faces are held as arrays of references.
    // These allow us to conveniently work through the items via foreach statements.
    FVCell[] cells;
    FVInterface[] faces;
    FVVertex[] vertices;
    BoundaryCondition[] bc; // collection of references to the boundary conditions
    //
    // We need to know the number of cells even if the grid is not read
    // for this block in the local process.
    size_t ncells_expected;
    size_t globalCellIdStart = 0; // needed to compute globalCellId
    //
    // Sometimes we need to look up cells and faces that are attached to a vertex.
    size_t[][] cellIndexListPerVertex;
    size_t[][] faceIndexListPerVertex;
    //
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    FlowState Lft;
    FlowState Rght;

    // list of vertex id's that makeup the fluidblock boundary
    // (used in the grid deformation methods in conjunction with
    // the shape sensitivity calculator).
    size_t[] boundaryVtxIndexList;

    // Shape sensitivity calculator workspace.
    version(shape_sensitivity) {
        immutable size_t MAX_PERTURBED_INTERFACES = 800;
        FVCell cellSave;
        FVInterface[MAX_PERTURBED_INTERFACES] ifaceP;

	size_t[] local_pcell_global_coord_list;
	size_t[][] local_ecell_global_coord_list;
	number[][] local_entry_list;
    
        // local objective function evaluation
        number locObjFcn;
        // arrays used to temporarily store data during construction of the flow Jacobian transpose 
        number[][] aa;
        size_t[][] ja;
        // local effects matrix for flow Jacobian transpose.
        // dimensions: [# local cells x # primitive vars] X [# local cells x # primitive vars]
        SMatrix!number JlocT;
        // external effects matrix for flow Jacobian transpose.
        // dimensions: [# local boundary cells x # primitive vars] X [# global cells x # primitive vars]
        SMatrix!number JextT;
        // Matrix used in preconditioning (low order, local, flow Jacobian).
        SMatrix!number P;
        SMatrix!number A; // Jacobian (w.r.t conserved variables)
        SMatrix!number Aext; // Jacobian (w.r.t conserved variables)
        // objective function senstivity w.r.t primitive variables
        number[] f;
        number[] b;
        // adjoint variables
        number[] psi;
        number[] delpsi;
        // residual sensitivity w.r.t. design variables (transposed)
        Matrix!number rT;           
        // local dot product of the residual sensitivity w.r.t. design variables (transposed) with the adjoint variables
        number[] rTdotPsi;
        // These arrays and matrices are directly tied to using the
        // GMRES iterative solver (use some directly from steady-state solver).
        number[] Z, z, wext, zext;
        //
        // Make a block-local copy of conserved quantities info
        size_t nConserved;
        size_t MASS;
        size_t X_MOM;
        size_t Y_MOM;
        size_t Z_MOM;
        size_t TOT_ENERGY;
        size_t TKE;
        size_t OMEGA;
    }
    
    version(steady_state)
    {
    // Work-space for steady-state solver
    // These arrays and matrices are directly tied to using the
    // GMRES iterative solver.
    SMatrix!number JcT; // transposed Jacobian (w.r.t conserved variables)
    ConservedQuantities maxRate, residuals;
    double normAcc, dotAcc;
    size_t nvars;
    Matrix!number Minv;
    number[] FU, dU, Dinv, r0, x0;
    number[] v, w, zed;
    number[] g0, g1;
    Matrix!number Q1;
    Matrix!number V;
    }

    this(int id, Grid_t grid_type, size_t ncells, string label)
    {
        this.id = id;
        this.grid_type = grid_type;
        this.ncells_expected = ncells;
        this.label = label;
        myConfig = dedicatedConfig[id];
        Linf_residuals = new ConservedQuantities(dedicatedConfig[id].gmodel.n_species,
                                                 dedicatedConfig[id].gmodel.n_modes);
        // Workspace for flux_calc method.
        Lft = new FlowState(dedicatedConfig[id].gmodel);
        Rght = new FlowState(dedicatedConfig[id].gmodel);
        // Lua interpreter for the block. 
        // It will be available for computing user-defined source terms.
        myL = luaL_newstate();
        luaL_openlibs(myL);
        registerGasModel(myL, LUA_GLOBALSINDEX);
        lua_pushinteger(myL, id);
        lua_setglobal(myL, "blkId");
        pushObj!(GasModel, GasModelMT)(myL, dedicatedConfig[id].gmodel);
        lua_setglobal(myL, "gmodel");
        lua_pushinteger(myL, dedicatedConfig[id].gmodel.n_species);
        lua_setglobal(myL, "n_species");
        lua_pushinteger(myL, dedicatedConfig[id].gmodel.n_modes);
        lua_setglobal(myL, "n_modes");
        // Although we make the helper functions available within 
        // the block-specific Lua interpreter, we should use 
        // those functions only in the context of the master thread.
        setSampleHelperFunctions(myL);
        setGridMotionHelperFunctions(myL);
    }

    ~this()
    {
        lua_close(myL);
    }

    override string toString() const { return "Block(id=" ~ to!string(id) ~ ")"; }
    @nogc size_t globalCellId(size_t localCellId) { return globalCellIdStart + localCellId; }

    @nogc abstract int get_interpolation_order();
    @nogc abstract void set_interpolation_order(int order);
    abstract void init_lua_globals();
    abstract void init_boundary_conditions(JSONValue json_data);
    @nogc abstract void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found);
    abstract void init_grid_and_flow_arrays(string gridFileName);
    @nogc abstract void compute_primary_cell_geometric_data(size_t gtl);
    @nogc abstract void compute_least_squares_setup(size_t gtl);
    @nogc abstract void sync_vertices_from_underlying_grid(size_t gtl=0);
    @nogc abstract void sync_vertices_to_underlying_grid(size_t gtl=0);
    abstract void read_new_underlying_grid(string fileName);
    abstract void write_underlying_grid(string fileName);
    abstract double read_solution(string filename, bool overwrite_geometry_data);
    abstract void write_solution(string filename, double sim_time);
    @nogc abstract void propagate_inflow_data_west_to_east();
    @nogc abstract void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0);
    @nogc abstract void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0);

    @nogc
    void identify_reaction_zones(int gtl)
    // Adjust the reactions-allowed flag for cells in this block.
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
        debug {
            if (myConfig.reacting && myConfig.verbosity_level >= 2) {
                writeln("identify_reaction_zones(): block ", id,
                        " cells inside zones = ", total_cells_in_reaction_zones, 
                        " out of ", total_cells);
                if (myConfig.reaction_zones.length == 0) {
                    writeln("Note that for no user-specified zones,",
                            " the whole domain is allowed to be reacting.");
                }
            }
        }
    } // end identify_reaction_zones()

    @nogc
    void identify_suppress_reconstruction_zones()
    // Adjust the in_suppress_reconstruction-zone flag for faces in this block.
    {
        foreach(f; faces) {
            f.in_suppress_reconstruction_zone = false;
            if (myConfig.suppress_reconstruction_zones.length > 0 ) {
                foreach(srz; myConfig.suppress_reconstruction_zones) {
                    if (srz.is_inside(f.pos, myConfig.dimensions)) {
                        f.in_suppress_reconstruction_zone = true;
                    }
                } // foreach srz
            }
            if (myConfig.suppress_reconstruction_at_boundaries && f.is_on_boundary) {
                f.in_suppress_reconstruction_zone = true;
            }
        } // foreach face
        //
        auto sfb = cast(SFluidBlock) this;
        if (sfb) {
            // Some special processing for Structured-grid blocks.
            if (myConfig.suppress_radial_reconstruction_at_xaxis &&
                myConfig.dimensions == 2 && myConfig.axisymmetric) {
                immutable double tol = 1.0e-9;
                // Work along each boundary and suppress reconstruction for faces that are
                // on the axis of symmetry or just one off the axis (but parallel to the axis).
                for (size_t i = sfb.imin; i <= sfb.imax; ++i) {
                    auto f = sfb.get_ifj!()(i,sfb.jmin);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifj!()(i,sfb.jmin+1);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
                for (size_t i = sfb.imin; i <= sfb.imax; ++i) {
                    auto f = sfb.get_ifj!()(i,sfb.jmax+1);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifj!()(i,sfb.jmax);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
                for (size_t j = sfb.jmin; j <= sfb.jmax; ++j) {
                    auto f = sfb.get_ifi!()(sfb.imin,j);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifi!()(sfb.imin+1,j);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
                for (size_t j = sfb.jmin; j <= sfb.jmax; ++j) {
                    auto f = sfb.get_ifi!()(sfb.imax+1,j);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifi!()(sfb.imax,j);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
            } // end if (myConfig...)
        } // end if(sfb)
        //
        debug {
            size_t total_faces_in_suppress_reconstruction_zones = 0;
            size_t total_faces = 0;
            foreach(f; faces) {
                total_faces_in_suppress_reconstruction_zones += (f.in_suppress_reconstruction_zone ? 1: 0);
                total_faces += 1;
            } // foreach face
            if (myConfig.verbosity_level >= 2) {
                writeln("identify_suppress_reconstruction_zones(): block ", id,
                        " faces inside zones = ", total_faces_in_suppress_reconstruction_zones, 
                        " out of ", total_faces);
                if (myConfig.suppress_reconstruction_zones.length == 0) {
                    writeln("Note that for no user-specified zones,",
                            " nominally reconstruction is allowed for the whole domain",
                            " but may be suppressed at the boundaries and/or near the axis of symmetry.");
                }
            }
        }
    } // end identify_suppress_reconstruction_zones()

    @nogc
    void identify_turbulent_zones(int gtl)
    // Adjust the in-turbulent-zone flag for cells in this block.
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
        debug {
            if (myConfig.turbulence_model != TurbulenceModel.none && 
                myConfig.verbosity_level >= 2) {
                writeln("identify_turbulent_zones(): block ", id,
                        " cells inside zones = ", total_cells_in_turbulent_zones, 
                        " out of ", total_cells);
                if (myConfig.turbulent_zones.length == 0) {
                    writeln("Note that for no user-specified zones,",
                            " the whole domain is allowed to be turbulent.");
                }
            }
        }
    } // end identify_turbulent_zones()

    @nogc
    void estimate_turbulence_viscosity()
    {
        final switch (myConfig.turbulence_model) {
        case TurbulenceModel.none:
            foreach (cell; cells) { cell.turbulence_viscosity_zero(); }
            return;
        case TurbulenceModel.baldwin_lomax:
            throw new FlowSolverException("need to port baldwin_lomax_turbulence_model");
        case TurbulenceModel.spalart_allmaras:
            throw new FlowSolverException("Should implement Spalart-Allmaras some day.");
        case TurbulenceModel.k_omega:
            version(komega) {
                foreach (cell; cells) { cell.turbulence_viscosity_k_omega(); }
            }
            break;
        } // end switch
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
            auto cL = iface.left_cell;
            auto cR = iface.right_cell;
            if (cL && cR) {
                number uL = geom.dot(cL.fs.vel, iface.n);
                number uR = geom.dot(cR.fs.vel, iface.n);
                number aL = cL.fs.gas.a;
                number aR = cR.fs.gas.a;
                number a_min = (aL < aR) ? aL : aR;
                iface.fs.S = ((uR - uL)/a_min) < tol;
            } else {
                iface.fs.S = 0;
            }
        }
        // Finally, mark cells as shock points if any of their interfaces are shock points.
        foreach (cell; cells) {
            cell.fs.S = false;
            foreach (face; cell.iface) {
                if (face.fs.S) cell.fs.S = true;
            }
        }
    } // end detect_shock_points()

    int count_invalid_cells(int gtl, int ftl)
    // Returns the number of cells that contain invalid data,
    // optionally patching bad cell data as it goes.
    //
    // Since this function may be called at the end of each stage of the gasdynamic update,
    // we must patch the conserved quantities at the appropriate flow-level, as well as
    // patching the flow data.
    //
    // Bad flow data can be identified by the density of internal energy 
    // being on the minimum limit or the velocity being very large.
    // There is also a flag to indicate that the thermo data is dodgy from an earlier
    // call to one of the thermochemical update functions.
    {
        int number_of_invalid_cells = 0;
        foreach(cell; cells) {
            if (cell.data_is_bad || cell.fs.check_data(cell.pos[0], myConfig.flowstate_limits) == false) {
                ++number_of_invalid_cells;
                if (myConfig.report_invalid_cells) {
                    writefln("count_invalid_cells: block_id=%d, cell_id=%d at pos %s\n",
                             id, cell.id, to!string(cell.pos[gtl]));
                    writeln(cell);
                }
                if (myConfig.adjust_invalid_cell_data) {
                    // We shall set the cell data to something that
                    // is valid (and self consistent).
                    FlowState[] neighbour_flows;
                    if (myConfig.report_invalid_cells) {
                        writeln("Adjusting cell data to a local average.");
                    }
                    foreach (i; 0 .. cell.iface.length) {
                        auto face = cell.iface[i];
                        auto other_cell = (cell.outsign[i] == 1) ? face.right_cell : face.left_cell;
                        if (other_cell && other_cell.contains_flow_data &&
                            other_cell.fs.check_data(other_cell.pos[gtl], myConfig.flowstate_limits)) {
                            neighbour_flows ~= other_cell.fs;
                        }
                    }
                    if (neighbour_flows.length == 0) {
                        string msg = "Block::count_invalid_cells(): There were no valid neighbours " ~
                            "to replace flow data in cell.";
                        if (!myConfig.report_invalid_cells) {
                            msg ~= "\nTo get more information, rerun with config.report_invalid_cells=true";
                        }
                        throw new FlowSolverException(msg);
                    }
                    cell.fs.copy_average_values_from(neighbour_flows, myConfig.gmodel);
                    scale_mass_fractions(cell.fs.gas.massf, 0.0, 0.9); // big assertion-error-tolerance
                    cell.data_is_bad = false; // assume that we've fixed it at this point.
                    cell.encode_conserved(gtl, ftl, omegaz);
                    if (0 != cell.decode_conserved(gtl, ftl, omegaz)) {
                        string msg = "Block::count_invalid_cells(): " ~
                            "Tried to replace flow data in cell but it's still bad.";
                        if (!myConfig.report_invalid_cells) {
                            msg ~= "\nTo get more information, rerun with config.report_invalid_cells=true";
                        }
                        throw new FlowSolverException(msg);
                    }
                    if (myConfig.report_invalid_cells) {
                        writefln("after flow-data replacement: block_id = %d, cell pos %.18e,%.18e,%.18e\n",
                                 id, cell.pos[gtl].x, cell.pos[gtl].y, cell.pos[gtl].z);
                        writeln(cell);
                    }
                } // end adjust_invalid_cell_data 
            } // end of if invalid data...
        } // foreach cell
        return number_of_invalid_cells;
    } // end count_invalid_cells()

    @nogc
    void flow_property_spatial_derivatives(int gtl)
    {
        final switch (myConfig.spatial_deriv_locn) {
        case SpatialDerivLocn.vertices:
            if (myConfig.dimensions == 2) {
                final switch (myConfig.spatial_deriv_calc) {
                case SpatialDerivCalc.least_squares:
                    foreach(vtx; vertices) { 
                        vtx.grad.gradients_leastsq(vtx.cloud_fs, vtx.cloud_pos, vtx.ws_grad);
                    }
                    break;
                case SpatialDerivCalc.divergence:
                    foreach(vtx; vertices) {
                        vtx.grad.gradients_xy_div(vtx.cloud_fs, vtx.cloud_pos);
                    }
                } // end switch
            } else {
                // Have only least-squares in 3D.
                foreach(vtx; vertices) {
                    vtx.grad.gradients_leastsq(vtx.cloud_fs, vtx.cloud_pos, vtx.ws_grad);
                }
            }
            foreach (iface; faces) {
                iface.average_vertex_deriv_values();
            }
            break;
        case SpatialDerivLocn.faces:
            if (myConfig.dimensions == 2) {
                final switch (myConfig.spatial_deriv_calc) {
                case SpatialDerivCalc.least_squares:
                    foreach(iface; faces) { 
                        iface.grad.gradients_leastsq(iface.cloud_fs, iface.cloud_pos, iface.ws_grad);
                    }
                    break;
                case SpatialDerivCalc.divergence:
                    foreach(iface; faces) {
                        iface.grad.gradients_xy_div(iface.cloud_fs, iface.cloud_pos);
                    }
                } // end switch
            } else { //3D
                final switch (myConfig.spatial_deriv_calc) {
                case SpatialDerivCalc.least_squares:
                    foreach(iface; faces) {
                        iface.grad.gradients_leastsq(iface.cloud_fs, iface.cloud_pos, iface.ws_grad);
                    }
                    break;
                case SpatialDerivCalc.divergence:
                    foreach(iface; faces) {
                        assert(0, "divergence thereom not implemented for 3D");
                    }
                } // end switch
            } // end if (myConfig.dimensions)
            break;
        case SpatialDerivLocn.cells:
            foreach(cell; cells) {
                cell.grad.gradients_leastsq(cell.cloud_fs, cell.cloud_pos, cell.ws_grad);
            }
        } // end switch (myConfig.spatial_deriv_locn)
    } // end flow_property_spatial_derivatives()

    @nogc
    void exchange_boundary_spatial_deriv_data(int gtl)
    {
        final switch (myConfig.spatial_deriv_locn) {
        case SpatialDerivLocn.vertices:
            break;
        case SpatialDerivLocn.faces:
            break;
        case SpatialDerivLocn.cells:
            if (grid_type == Grid_t.structured_grid) {
            } else {
                foreach (bcond; bc) {
                    if (bcond.ghost_cell_data_available == false) { continue; }
                    // Proceed to do some work only if we have ghost cells in which to insert gradients.
                    bool found_mapped_cell_bc = false;
                    foreach (gce; bcond.preReconAction) {
                        auto mygce = cast(GhostCellMappedCellCopy)gce;
                        if (mygce && !myConfig.in_mpi_context) {
                            found_mapped_cell_bc = true;
                            // There is a mapped-cell backing the ghost cell, so we can copy its gradients.
                            foreach (i, f; bcond.faces) {
                                // Only FVCell objects in an unstructured-grid are expected to have
                                // precomputed gradients.  There will be an initialized reference
                                // in the FVCell object of a structured-grid block, so we need to
                                // test and avoid copying from such a reference.
                                auto mapped_cell_grad = mygce.get_mapped_cell(i).grad;
                                if (bcond.outsigns[i] == 1) {
                                    if (mapped_cell_grad) {
                                        f.right_cell.grad.copy_values_from(mapped_cell_grad);
                                    } else {
                                        // Fall back to looking over the face for suitable gradient data.
                                        f.right_cell.grad.copy_values_from(f.left_cell.grad);
                                    }
                                } else {
                                    if (mapped_cell_grad) {
                                        f.left_cell.grad.copy_values_from(mapped_cell_grad);
                                    } else {
                                        f.left_cell.grad.copy_values_from(f.right_cell.grad);
                                    }
                                }
                            } // end foreach f
                        } // end if (mygce)
                    } // end foreach gce
                    if (!found_mapped_cell_bc) {
                        // There are no other cells backing the ghost cells on this boundary.
                        // Fill in ghost-cell gradients from the other side of the face.
                        foreach (i, f; bcond.faces) {
                            if (bcond.outsigns[i] == 1) {
                                f.right_cell.grad.copy_values_from(f.left_cell.grad);
                            } else {
                                f.left_cell.grad.copy_values_from(f.right_cell.grad);
                            }
                        } // end foreach f
                    } // end if !found_mapped_cell_bc
                } // end foreach bcond
            }
        } // end switch (myConfig.spatial_deriv_locn)
    } // end flow_property_spatial_derivatives()
    
    @nogc
    void clear_fluxes_of_conserved_quantities()
    {
        foreach (iface; faces) { iface.F.clear(); }
    }

    @nogc
    void viscous_flux()
    {
        foreach (iface; faces) { iface.viscous_flux_calc(); } 
    }

    @nogc
    void init_residuals()
    // Initialization of data for later computing residuals.
    {
        mass_residual = 0.0;
        mass_residual_loc.clear();
        energy_residual = 0.0;
        energy_residual_loc.clear();
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
        mass_residual_loc.clear();
        energy_residual = 0.0;
        energy_residual_loc.clear();
        foreach(FVCell cell; cells) {
            number local_residual = (cell.fs.gas.rho - cell.rho_at_start_of_step) / cell.fs.gas.rho;
            local_residual = fabs(local_residual);
            if ( local_residual > mass_residual ) {
                mass_residual = local_residual;
                mass_residual_loc.set(cell.pos[gtl]);
            }
            // In the following line, the zero index is used because,
            // at the end of the gas-dynamic update, that index holds
            // the updated data.
            local_residual = (cell.U[0].total_energy - cell.rE_at_start_of_step) / cell.U[0].total_energy;
            local_residual = fabs(local_residual);
            if ( local_residual > energy_residual ) {
                energy_residual = local_residual;
                energy_residual_loc.set(cell.pos[gtl]);
            }
        } // for cell
    } // end compute_residuals()

    @nogc
    void compute_Linf_residuals()
    // Compute Linf residuals for conserved quantities.
    // This is similar to the calculation above of
    // residual, but this differs by a factor of the timestep size
    // because here the residual is taken as R(U) = dU/dt.
    // We will assume that dUdt[0] is up-to-date.
    {
        Linf_residuals.copy_values_from(cells[0].dUdt[0]);
        Linf_residuals.mass = fabs(Linf_residuals.mass);
        Linf_residuals.momentum.set(fabs(Linf_residuals.momentum.x),
                                    fabs(Linf_residuals.momentum.y),
                                    fabs(Linf_residuals.momentum.z));
        Linf_residuals.total_energy = fabs(Linf_residuals.total_energy);
        foreach (cell; cells) {
            Linf_residuals.mass = fmax(Linf_residuals.mass, fabs(cell.dUdt[0].mass));
            Linf_residuals.momentum.set(fmax(Linf_residuals.momentum.x, fabs(cell.dUdt[0].momentum.x)),
                                        fmax(Linf_residuals.momentum.y, fabs(cell.dUdt[0].momentum.y)),
                                        fmax(Linf_residuals.momentum.z, fabs(cell.dUdt[0].momentum.z)));
            Linf_residuals.total_energy = fmax(Linf_residuals.total_energy, fabs(cell.dUdt[0].total_energy));
        }
    } // end compute_Linf_residuals()

    @nogc
    void compute_L2_residual()
    {
        L2_residual = 0.0;
        foreach (cell; cells) {
            L2_residual += fabs(cell.dUdt[0].mass)^^2;
	    L2_residual += fabs(cell.dUdt[0].momentum.x)^^2;
	    L2_residual += fabs(cell.dUdt[0].momentum.y)^^2;
	    L2_residual += fabs(cell.dUdt[0].momentum.z)^^2;
	    L2_residual += fabs(cell.dUdt[0].total_energy)^^2;
	    if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
		L2_residual += fabs(cell.dUdt[0].tke)^^2;
		L2_residual += fabs(cell.dUdt[0].omega)^^2;
	    } 
	    L2_residual = sqrt(L2_residual);
        }
    } // end compute_Linf_residuals()

    @nogc
    void compute_mass_balance()
    {
        mass_balance = 0.0;
        foreach(boundary; bc) {
            if (boundary.type != "exchange_over_full_face" && boundary.type != "exchange_using_mapped_cells") {
                foreach(i, face; boundary.faces) {
                    mass_balance += boundary.outsigns[i] * face.F.mass * face.area[0];
                }
            }
        }
    } // end compute_Linf_residuals()
    
    @nogc
    void residual_smoothing_dUdt(size_t ftl)
    {
        assert(ftl < cells[0].dUdt.length, "inconsistent flow time level and allocated dUdt");
        foreach (c; cells) {
            c.dUdt_copy[0].copy_values_from(c.dUdt[ftl]);
            c.dUdt_copy[1].copy_values_from(c.dUdt[ftl]);
        }
        double eps = GlobalConfig.residual_smoothing_weight;
        if (GlobalConfig.residual_smoothing_type == ResidualSmoothingType.explicit) { 
            foreach (c; cells) {
                double total = 1.0;
                foreach (i, f; c.iface) {
                    total += eps;
                    auto other_cell = (c.outsign[i] > 0.0) ? f.right_cell : f.left_cell;
                    if (other_cell && other_cell.contains_flow_data) {
                        c.dUdt[ftl].add(other_cell.dUdt_copy[0], eps);
                    }
                }
                c.dUdt[ftl].scale(1.0/total);
            }
        } else { // ResidualSmoothingType.implicit
            int n_iters = GlobalConfig.residual_smoothing_iterations;
            // perform Jacobi iterations
            foreach (ni; 0..n_iters) {
                foreach (c; cells) {
                    double n = 0;
                    c.dUdt_copy[1].copy_values_from(c.dUdt[ftl]);
                    foreach (i, f; c.iface) {
                        auto other_cell = (c.outsign[i] > 0.0) ? f.right_cell : f.left_cell;
                        if (other_cell && other_cell.contains_flow_data) {
                            n += 1;
                            c.dUdt_copy[1].add(other_cell.dUdt_copy[0], eps);
                        }
                    }
                    double scale = 1.0+n*eps;
                    c.dUdt_copy[1].scale(1.0/scale);
                }
                foreach (c; cells) {
                    c.dUdt_copy[0].copy_values_from(c.dUdt_copy[1]);
                }
            }
            // replace residual with smoothed residual
            foreach (c; cells) {
                c.dUdt[ftl].copy_values_from(c.dUdt_copy[1]);
            }
        }
    } // end residual_smoothing_dUdt()
    
    @nogc
    double update_c_h(double dt_current)
    // Update the c_h value for the divergence cleaning mechanism.
    {
        double min_L_for_block, cfl_local, cfl_max;
        bool first = true;
        foreach(FVCell cell; cells) {
            // Search for the minimum length scale and the maximum CFL value in the block.
            if (first) {
                min_L_for_block = cell.L_min.re;
                cfl_local = cell.signal_frequency() * dt_current;
                cfl_max = cfl_local;
                first = false;
            } else {
                min_L_for_block = fmin(cell.L_min.re, min_L_for_block);
                cfl_local = cell.signal_frequency() * dt_current;
                cfl_max = fmax(cfl_local, cfl_max);
            }
        }
        return cfl_max * min_L_for_block / dt_current;
    } // end update_c_h()

    @nogc
    double[2] determine_time_step_size(double dt_current, bool check_cfl)
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
        double cfl_adjust = 0.5; // Adjust cfl_max (and hence dt_allow) if cfl_max > cfl_allow
        // The following limits allow the simulation of the sod shock tube
        // to get just a little wobbly around the shock.
        // Lower values of cfl should be used for a smooth solution.
        switch (number_of_stages_for_update_scheme(myConfig.gasdynamic_update_scheme)) {
        case 1: cfl_allow = 0.9; break;
        case 2: cfl_allow = 1.2; break;
        case 3: cfl_allow = 1.6; break;
        default: cfl_allow = 0.9;
        }
        // when using residual smoothing (implicit) we should be able to achieve a higher stable CFL
        // so let's relax the cfl_allow
        if (myConfig.residual_smoothing &&
            myConfig.with_local_time_stepping &&
            GlobalConfig.residual_smoothing_type == ResidualSmoothingType.implicit) cfl_allow *= 10.0;
        int local_time_stepping_limit_factor = myConfig.local_time_stepping_limit_factor; // for local time-stepping we limit the larger time-steps by a factor of the smallest timestep
        bool first = true;
        foreach(FVCell cell; cells) {
            signal = cell.signal_frequency();
            cfl_local = dt_current * signal; // Current (Local) CFL number
            dt_local = cfl_value / signal; // Recommend a time step size.
            cell.dt_local = fmin(dt_local, dt_current*local_time_stepping_limit_factor); // set local time-step in cell
            cell.dt_local = fmin(cell.dt_local, GlobalConfig.dt_max); // Limit the largest local time-step to a set input value
	    if (first) {
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
        if (check_cfl && (cfl_max < 0.0 || cfl_max > cfl_allow)) {
            string msg = "Bad cfl number encountered";
            debug { msg ~= text(" cfl_max=", cfl_max, " for FluidBlock ", id); }
            debug { writeln(msg); } // Write out warning message when running in debug mode
            cfl_max = cfl_adjust*cfl_allow;// If cfl_max exceeds cfl_allow, simply reduce the
                                    // cfl_max to cfl_adjust*cfl_allow. A value of 0.5 seems 
                                    // to work robustly. Values to 0.7 also work, beyond this 
                                    // and code begins to crash due to numerical instability.
                                    // Results in auto-limitation of the time step/cfl
            dt_allow = cfl_max/signal; // Reduce dt_allow according to new rescaled cfl_max
            //throw new FlowSolverException(msg); // Previous code threw an error and halted
        }
        return [dt_allow, cfl_max];
    } // end determine_time_step_size()

    void applyPreReconAction(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPreReconAction(t, gtl, ftl); }
    }

    void applyPostConvFluxAction(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPostConvFluxAction(t, gtl, ftl); }
    }

    void applyPreSpatialDerivActionAtBndryFaces(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPreSpatialDerivActionAtBndryFaces(t, gtl, ftl); }
    }

    void applyPreSpatialDerivActionAtBndryCells(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPreSpatialDerivActionAtBndryCells(t, gtl, ftl); }
    }

    void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPostDiffFluxAction(t, gtl, ftl); }
    }

    version(steady_state) {
    void allocate_GMRES_workspace()
    {
        size_t nConserved = nConservedQuantities;
        int n_species = GlobalConfig.gmodel_master.n_species();
        int n_modes = GlobalConfig.gmodel_master.n_modes();
        maxRate = new ConservedQuantities(n_species, n_modes);
        residuals = new ConservedQuantities(n_species, n_modes);

        size_t mOuter = to!size_t(GlobalConfig.sssOptions.maxOuterIterations);
        size_t mInner = to!size_t(GlobalConfig.sssOptions.nInnerIterations);
        size_t n = nConserved*cells.length;
        nvars = n;
        // Now allocate arrays and matrices
        FU.length = n;
        dU.length = n; dU[] = to!number(0.0);
        r0.length = n;
        x0.length = n;
        Dinv.length = n;
        v.length = n;
        w.length = n;
        zed.length = n;
        g0.length = mOuter+1;
        g1.length = mOuter+1;
        //h_outer.length = mOuter+1;
        //hR_outer.length = mOuter+1;
        V = new Matrix!number(n, mOuter+1);
        //H0_outer = new Matrix!number(mOuter+1, mOuter);
        //H1_outer = new Matrix!number(mOuter+1, mOuter);
        //Gamma_outer = new Matrix!number(mOuter+1, mOuter+1);
        //Q0_outer = new Matrix!number(mOuter+1, mOuter+1);
        Q1 = new Matrix!number(mOuter+1, mOuter+1);
    }
    }
} // end class FluidBlock
