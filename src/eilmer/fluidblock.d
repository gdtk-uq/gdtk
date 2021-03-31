// fluidblock.d
// Base class for blocks of cells, for use within the Eilmer flow solver.
// Peter J. 2014-07-18 first cut.
// Now a derived class from the Block base class
// Kyle A. Damm 2020-02-11

module fluidblock;

import std.algorithm;
import std.conv;
import std.stdio;
import std.string;
import std.math;
import std.json;
import std.format;
import std.array;
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
import json_helper;
import grid_motion;
import grid_deform;
import sfluidblock; // needed for some special-case processing, below
import shockdetectors;
import block;
import jacobian;
version(mpi_parallel) {
    import mpi;
}

// version(matplotlib) {
// import plt = matplotlibd.pyplot;
// import std.format;
// }



// To distinguish ghost cells from active cells, we start their id values at
// an arbitrarily high value.  It seem high to me (PJ) but feel free to adjust it
// if you start using grids larger I expect.
enum ghost_cell_start_id = 1_000_000_000;


// The flow solver handles structured- and unstructured-grid blocks via this base class.
// Mostly, we view the block as an unstructured bag of cells because that requires least
// knowledge in the calling code.
class FluidBlock : Block {
public:
    Grid_t grid_type; // structured or unstructured
    bool may_be_turbulent; // if true, the selected turbulence model is active
                           // within this block.
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
    size_t n_ghost_cell_layers;
    //
    // Sometimes we need to look up cells and faces that are attached to a vertex.
    size_t[][] cellIndexListPerVertex;
    size_t[][] faceIndexListPerVertex;
    //
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    FlowState Lft;
    FlowState Rght;
    //
    // Super-time-stepping parameters used when applying flexible stages per block
    int s_RKL; // number of super-step
    double dt_parab; // allowable parabolic time-step
    // list of vertex id's that makeup the fluidblock boundary
    // (used in the grid deformation methods in conjunction with
    // the shape sensitivity calculator).
    size_t[] boundaryVtxIndexList;

    // Shape sensitivity calculator workspace.
    version(shape_sensitivity) {
        FlowJacobianT flowJacobianT;
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
        size_t SPECIES;
    }

    version(nk_accelerator)
    {
    // Work-space for Newton-Krylov accelerator
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

    this(int id, Grid_t grid_type, size_t ncells, size_t n_ghost_cell_layers, string label)
    {
        super(id, label);
        this.grid_type = grid_type;
        this.ncells_expected = ncells;
        this.n_ghost_cell_layers = n_ghost_cell_layers;
        myConfig = dedicatedConfig[id];
    }

    void init_workspace()
    {
        Linf_residuals = new ConservedQuantities(dedicatedConfig[id].n_species,
                                                 dedicatedConfig[id].n_modes);
        // Workspace for flux_calc method.
        Lft = new FlowState(dedicatedConfig[id].gmodel);
        Rght = new FlowState(dedicatedConfig[id].gmodel);
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
    @nogc abstract void propagate_inflow_data_west_to_east();
    @nogc abstract void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0, FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = []);
    @nogc abstract void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0, FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = []);

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
                foreach (i; 0 .. sfb.nic) {
                    auto f = sfb.get_ifj(i,0);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifj(i,1);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
                foreach (i; 0 .. sfb.nic) {
                    auto f = sfb.get_ifj(i,sfb.njc);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifj(i,sfb.njc-1);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
                foreach (j; 0 .. sfb.njc) {
                    auto f = sfb.get_ifi(0,j);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifi(1,j);
                        f1.in_suppress_reconstruction_zone = true;
                    }
                }
                foreach (j; 0 .. sfb.njc) {
                    auto f = sfb.get_ifi(sfb.nic,j);
                    if (f.pos.y < tol) {
                        f.in_suppress_reconstruction_zone = true;
                        auto f1 = sfb.get_ifi(sfb.nic-1,j);
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
            if (may_be_turbulent) {
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
            } else {
                cell.in_turbulent_zone = false;
            }
            total_cells_in_turbulent_zones += (cell.in_turbulent_zone ? 1: 0);
            total_cells += 1;
        } // foreach cell
        debug {
            if (myConfig.turb_model.isTurbulent &&
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

    // Increment the local DFT in each cell.
    // Do we need garbage collection here? I don't think so?
    @nogc
    void increment_DFT(size_t DFT_step) {
        foreach(cell; cells) {
            cell.increment_local_DFT(DFT_step);
        }
    }

    @nogc
    void estimate_turbulence_viscosity(FVCell[] cell_list = [])
    {
        version(turbulence) { // Exit instantly if turbulence capability disabled
        if (cell_list.length == 0) { cell_list = cells; }
        foreach (cell; cell_list) {
            cell.turbulence_viscosity();
            cell.turbulence_viscosity_factor(myConfig.transient_mu_t_factor);
            cell.turbulence_viscosity_limit(myConfig.max_mu_t_factor);
            cell.turbulence_viscosity_zero_if_not_in_zone();
        }
        }
    } // end estimate_turbulence_viscosity()

    @nogc
    void set_cell_dt_chem(double dt_chem)
    {
        foreach (cell; cells) { cell.dt_chem = dt_chem; }
    }

    @nogc
    void detect_shock_points()
    // Detects shocks by looking for compression between adjacent cells.
    //
    {
        number comp_tol = myConfig.compression_tolerance;
        number shear_tol = myConfig.shear_tolerance;
        final switch (myConfig.shock_detector) {
        case ShockDetector.PJ:
            foreach (iface; faces) {
                iface.fs.S = PJ_ShockDetector(iface, comp_tol, shear_tol);
            }
            break;
        }
    } // end detect_shock_points()

    @nogc
    void diffuse_shock_marker()
    {
        // Spatially diffuse the shock marker.
        // To avoid order of operation effects, we load the cell based shock
        // detector value into faces so we have static data to work from.
        foreach (face; faces) {
            if (face.fs.S == 1.0) continue;
            face.fs.S = fmax(face.left_cell.fs.S, face.right_cell.fs.S);
        }
        // now use the face shock detector value to update the cell values
        foreach (cell; cells) {
            if (cell.fs.S == 1.0) continue;
            number sum = 0.0;
            int count = 0;
            foreach (face; cell.iface) {
                sum += face.fs.S;
                count++;
            }
            cell.fs.S = sum / count;
        }
    } // end diffuse_shock_marker()

    @nogc
    void shock_cells_to_faces()
    // Mark faces as shocked if either attached cell is shocked
    //
    {
        foreach(face; faces) {
            face.fs.S = fmax(face.left_cell.fs.S, face.right_cell.fs.S);
        }
    } // end shock_cells_to_faces()

    @nogc
    void shock_faces_to_cells()
    // Mark cells as shocked if any of their interfaces are shocked
    //
    {
        foreach (cell; cells) {
            cell.fs.S = 0.0;
            foreach (face; cell.iface) {
                cell.fs.S = fmax(cell.fs.S, face.fs.S);
            }
        }
    } // end shock_faces_to_cells()

    @nogc
    void enforce_strict_shock_detector()
    // mark faces as shocked if either left or right cell is shocked
    {
        foreach(face; faces) {
            if (face.fs.S > 0.0) {
                face.fs.S = 1.0;
                continue;
            }
            // need to check for existence of cell before checking its shock state
            if (face.left_cell) {
                if (face.left_cell.fs.S > 0.0) {
                    face.fs.S = 1.0;
                    continue;
                }
            }
            if (face.right_cell) {
                if (face.right_cell.fs.S > 0.0) {
                    face.fs.S = 1.0;
                    continue;
                }
            }
        }
    } // end enforce_strict_shock_detector

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
            if (cell.data_is_bad || cell.fs.check_data(cell.pos[0], myConfig) == false) {
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
                            other_cell.fs.check_data(other_cell.pos[gtl], myConfig)) {
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
                    version(multi_species_gas){
                        scale_mass_fractions(cell.fs.gas.massf, 0.0, 0.9); // big assertion-error-tolerance
                    }
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

    // Per discussion with PJ/KD on 191129 most of this routine is on the chopping block (NNG)
    // Possibly move to lsq Cell centered gradients only, followed to average_cell_deriv_values
    // on faces in simcore AFTER calling exchange_ghost_cell_boundary_viscous_gradient_data
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

            // We've finished computing gradients at vertices, now copy them around if needed
            // Interfaces need them for reconstruction/viscous fluxes
            foreach (iface; faces) {
                iface.average_vertex_deriv_values();
            }
            // Turbulence models will need cell centered gradients (also for axisymmetric!)
            if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent)){
                foreach (cell; cells) {
                    cell.average_vertex_deriv_values();
                }
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

            // Finished computing gradients at interfaces, now copy them around if needed
            // Turbulence models and axisymmetric source terms need cell centered gradients
            if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent)){
                foreach (cell; cells) {
                    cell.average_interface_deriv_values();
                }
            }
            break;

        case SpatialDerivLocn.cells:
            foreach(cell; cells) {
                cell.grad.gradients_leastsq(cell.cloud_fs, cell.cloud_pos, cell.ws_grad);
            }
            // Cell centered gradients are transformed to interfaces in simcore.

        } // end switch (myConfig.spatial_deriv_locn)
    } // end flow_property_spatial_derivatives()

    @nogc
    void clear_fluxes_of_conserved_quantities()
    {
        foreach (iface; faces) { iface.F.clear(); }
    }

    @nogc
    void viscous_flux(FVInterface[] face_list = [])
    {
        if (face_list.length == 0) { face_list = faces; }
        foreach (iface; face_list) { iface.viscous_flux_calc(); }
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
	    version(turbulence) {
            foreach(it; 0 .. myConfig.turb_model.nturb){
                L2_residual += fabs(cell.dUdt[0].rhoturb[it])^^2;
            }
	    }
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
    double[3] determine_time_step_size(double dt_current, bool check_cfl)
    // Compute the local time step limit for all cells in the block.
    // The overall time step is limited by the worst-case cell.
    {
        // for STS (hyp = hyperbolic/convective, parab = parabolic/viscous)
        double signal_hyp;
        double signal_parab;
	double dt_allow_hyp;
	double dt_allow_parab;
        //
        double cfl_value = GlobalConfig.cfl_value;
        double dt_local;
        double cfl_local;
        double signal;
       	//
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
        // when using implicit residual smoothing we should be able to achieve a higher stable CFL
        // so let's relax the cfl_allow
        if (myConfig.residual_smoothing &&
            myConfig.with_local_time_stepping &&
            GlobalConfig.residual_smoothing_type == ResidualSmoothingType.implicit) cfl_allow *= 10.0;
        // for local time-stepping we limit the larger time-steps by a factor of the smallest timestep
        int local_time_stepping_limit_factor = myConfig.local_time_stepping_limit_factor;
        bool first = true;
        foreach(FVCell cell; cells) {
            signal = cell.signal_frequency();
	    if (myConfig.with_super_time_stepping) {
		signal_hyp = cell.signal_hyp.re;
		signal_parab = cell.signal_parab.re;
		if (first) {
		    dt_allow_hyp = cfl_value / signal_hyp;
		    dt_allow_parab = cfl_value / signal_parab;
		    first = false;
		} else {
                    dt_allow_hyp = fmin(dt_allow_hyp, cfl_value / signal_hyp);
		    dt_allow_parab = fmin(dt_allow_parab, cfl_value / signal_parab);
		}
		dt_allow = fmin(dt_allow_hyp, GlobalConfig.dt_max); // set the allowable time-step based on hyperbolic time-step
            } else {
		// no STS
		dt_allow_hyp = 0;
		dt_allow_parab = 0;
		//
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
	    }
        } // foreach cell
        if (myConfig.with_super_time_stepping == false && check_cfl && (cfl_max < 0.0 || cfl_max > cfl_allow)) {
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
        return [dt_allow, cfl_max, dt_allow_parab];
    } // end determine_time_step_size()

    void applyPreReconAction(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPreReconAction(t, gtl, ftl); }
    }

    void applyPreReconAction(double t, int gtl, int ftl, FVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.which_boundary == f.bc_id) { boundary.applyPreReconAction(t, gtl, ftl, f); }
        }
    }

    void applyPostConvFluxAction(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPostConvFluxAction(t, gtl, ftl); }
    }

    void applyPostConvFluxAction(double t, int gtl, int ftl, FVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.which_boundary == f.bc_id) { boundary.applyPostConvFluxAction(t, gtl, ftl, f); }
        }
    }

    void applyPreSpatialDerivActionAtBndryFaces(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPreSpatialDerivActionAtBndryFaces(t, gtl, ftl); }
    }

    void applyPreSpatialDerivActionAtBndryFaces(double t, int gtl, int ftl, FVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.which_boundary == f.bc_id) { boundary.applyPreSpatialDerivActionAtBndryFaces(t, gtl, ftl, f); }
        }
    }

    void applyPreSpatialDerivActionAtBndryCells(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPreSpatialDerivActionAtBndryCells(t, gtl, ftl); }
    }

    void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
        foreach(boundary; bc) { boundary.applyPostDiffFluxAction(t, gtl, ftl); }
    }

    void applyPostDiffFluxAction(double t, int gtl, int ftl, FVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.which_boundary == f.bc_id) { boundary.applyPostDiffFluxAction(t, gtl, ftl, f); }
        }
    }

    version(shape_sensitivity) {

    void initialize_transpose_jacobian(size_t spatial_order_of_jacobian)
    {
        /*
          This method initializes the transpose flow Jacobian matrix attached the FluidBlock object.
          We gather the cell/interface residual stencils for the interior and ghost cells at this point,
          since we need to know the number of expected entries in the Jacobian matrix to pre-size the
          sparse matrix arrays.
         */

        size_t nentry = 0;
        // gather the expected number of non-zero entries in the flow Jacobian
        foreach (cell; cells) {
            cell.gather_residual_stencil_lists(spatial_order_of_jacobian);
            nentry += cell.cell_list.length;
        }
        shared size_t my_nConserved = GlobalConfig.cqi.nConservedQuantities;
        size_t ncells = cells.length;
        flowJacobianT = new FlowJacobianT(myConfig.dimensions, my_nConserved, spatial_order_of_jacobian, nentry, ncells);

        // we will gather the ghost cell residual stencil lists
        // at this point as well for convenience, we need them
        // when we apply the boundary condition corrections later
        foreach ( bndary; bc ) {
            if ( bndary.type == "exchange_using_mapped_cells" || bndary.type == "exchange_over_full_face") { continue; }
            foreach ( iface, face; bndary.faces) {
                FVCell ghost_cell; FVCell cell;
                if (bndary.outsigns[iface] == 1) {
                    cell = face.left_cell;
                    ghost_cell = face.right_cell;
                } else {
                    cell = face.right_cell;
                    ghost_cell = face.left_cell;
                }
                ghost_cell.gather_residual_stencil_lists_for_ghost_cells(spatial_order_of_jacobian, cell.cell_cloud);
            }
        }
    } // end initialize_transpose_jacobian()

    void evaluate_transpose_jacobian()
    {
        /*
          Higher level method used to evaluate the transpose flow Jacobian attached to the FluidBlock object.
         */

        // temporarily change interpolation order
        shared int interpolation_order_save = GlobalConfig.interpolation_order;
        myConfig.interpolation_order = to!int(flowJacobianT.spatial_order);

        // fill out the rows of the Jacobian for a cell
        flowJacobianT.prepare_crs_indexes();
        foreach(cell; cells) { evaluate_block_row_of_jacobian(cell); }

        // add boundary condition corrections to boundary cells
        apply_transpose_jacobian_bcs();

        // return the interpolation order to its original state
        myConfig.interpolation_order = interpolation_order_save;

    } // end evaluate_transpose_jacobian()

    void evaluate_block_row_of_jacobian(FVCell pcell)
    {
        // Make a stack-local copy of conserved quantities info
        size_t nConserved = GlobalConfig.cqi.nConservedQuantities;
        size_t MASS = GlobalConfig.cqi.mass;
        size_t X_MOM = GlobalConfig.cqi.xMom;
        size_t Y_MOM = GlobalConfig.cqi.yMom;
        size_t Z_MOM = GlobalConfig.cqi.zMom;
        size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
        size_t TKE = GlobalConfig.cqi.tke;
        size_t SPECIES = GlobalConfig.cqi.species;
        number EPS = flowJacobianT.eps;
        int gtl = 0; int ftl = 1;
        pcell.U[ftl].copy_values_from(pcell.U[0]);

        // perturb either the cell flow state or conserved quantities and
        // evaluate the perturbed residuals for the cells in the residual stencil
        if ( flowJacobianT.wrtConserved ) {
            mixin(computeResidualSensitivity("pcell", "mass", "MASS", true, true));
            mixin(computeResidualSensitivity("pcell", "momentum.refx", "X_MOM", false, true));
            mixin(computeResidualSensitivity("pcell", "momentum.refy", "Y_MOM", false, true));
            if ( myConfig.dimensions == 3 ) { mixin(computeResidualSensitivity("pcell", "momentum.refz", "Z_MOM", false, true)); }
            mixin(computeResidualSensitivity("pcell", "total_energy", "TOT_ENERGY", true, true));
            foreach(ii; 0 .. myConfig.turb_model.nturb) { mixin(computeResidualSensitivity("pcell", "rhoturb[ii]", "TKE+ii", false, true)); }
            version(multi_species_gas){
            if (myConfig.n_species > 1) {
                foreach(ii; 0 .. myConfig.gmodel.n_species) { mixin(computeResidualSensitivity("pcell", "massf[ii]", "SPECIES+ii", false, true)); }
            }
            }
        } else { // wrtPrimitive
            mixin(computeResidualSensitivity("pcell", "gas.rho", "MASS", true, false));
            mixin(computeResidualSensitivity("pcell", "vel.refx", "X_MOM", false, false));
            mixin(computeResidualSensitivity("pcell", "vel.refy", "Y_MOM", false, false));
            if ( myConfig.dimensions == 3 ) { mixin(computeResidualSensitivity("pcell", "vel.refz", "Z_MOM", false, false)); }
            mixin(computeResidualSensitivity("pcell", "gas.p", "TOT_ENERGY", true, false));
            foreach(ii; 0 .. myConfig.turb_model.nturb) { mixin(computeResidualSensitivity("pcell", "turb[ii]", "TKE+ii", false, false)); }
            version(multi_species_gas){
            if (myConfig.n_species > 1) {
                foreach(ii; 0 .. myConfig.gmodel.n_species) { mixin(computeResidualSensitivity("pcell", "gas.massf[ii]", "SPECIES+ii", false, false)); }
            }
            }
        }

        // we now populate the pre-sized sparse matrix representation of the flow Jacobian
        // note that nConserved and nPrimitive are interchangeable
        size_t jidx; // column index into the matrix
        // loop through nConserved rows
        for (size_t ip = 0; ip < nConserved; ++ip) {
            // loop through cells that will have non-zero entries
            foreach(cell; pcell.cell_list) {
                // loop through nConserved columns for each effected cell
                for ( size_t jp = 0; jp < nConserved; ++jp ) {
                    assert(cell.id < ghost_cell_start_id, "Oops, we expect to not find a ghost cell at this point.");
                    jidx = cell.id*nConserved + jp; // column index
                    // note we are inherently performing a transpose operation on the next line,
                    // we do this since it appeared more natural to build up the
                    // global transpose flow Jacobian in compressed row storage format
                    flowJacobianT.local.aa[flowJacobianT.aa_idx] = cell.dQdU[jp][ip];
                    // fill out the sparse matrix idexes ready for the next entry in the row
                    flowJacobianT.aa_idx += 1;
                    flowJacobianT.local.ja[flowJacobianT.ja_idx] = jidx;
                    flowJacobianT.ja_idx += 1;
                }
            }
            // prepare the sparse matrix for a new row
            flowJacobianT.local.ia[flowJacobianT.ia_idx] = flowJacobianT.aa_idx;
            flowJacobianT.ia_idx += 1;
        }
    } // end evaluate_block_row_of_jacobian()

    static string computeResidualSensitivity(string cellName, string varName, string posInArray, bool includeThermoUpdate, bool wrtConserved)
    {
        string codeStr ="{";

        if ( wrtConserved ) {
            codeStr ~= ""~cellName~".U[ftl]."~varName~" += EPS;
                        "~cellName~".decode_conserved(gtl, ftl, 0.0); ";
        } else {
            codeStr ~= " "~cellName~".fs."~varName~" += EPS; ";
        }

        if ( includeThermoUpdate ) {
            codeStr ~="
            myConfig.gmodel.update_thermo_from_rhop("~cellName~".fs.gas);
            myConfig.gmodel.update_trans_coeffs("~cellName~".fs.gas);
            myConfig.gmodel.update_sound_speed("~cellName~".fs.gas);
            ";
        }

        codeStr ~="
        evalRHS(gtl, ftl, "~cellName~".cell_list, "~cellName~".face_list);

        foreach (i, cell; "~cellName~".cell_list) {
            cell.dQdU[MASS][" ~ posInArray ~ "] = cell.dUdt[ftl].mass.im/EPS.im;
            cell.dQdU[X_MOM][" ~ posInArray ~ "] = cell.dUdt[ftl].momentum.x.im/EPS.im;
            cell.dQdU[Y_MOM][" ~ posInArray ~ "] = cell.dUdt[ftl].momentum.y.im/EPS.im;
            if (myConfig.dimensions == 3) { cell.dQdU[Z_MOM][" ~ posInArray ~ "] = cell.dUdt[ftl].momentum.z.im/EPS.im; }
            cell.dQdU[TOT_ENERGY][" ~ posInArray ~ "] = cell.dUdt[ftl].total_energy.im/EPS.im;
            foreach(it; 0 .. myConfig.turb_model.nturb) { cell.dQdU[TKE+it][" ~ posInArray ~ "] = cell.dUdt[ftl].rhoturb[it].im/EPS.im; }
            version(multi_species_gas){
            if (myConfig.n_species > 1) {
                foreach(sp; 0 .. myConfig.gmodel.n_species) { cell.dQdU[SPECIES+sp][" ~ posInArray ~ "] = cell.dUdt[ftl].massf[sp].im/EPS.im; }
            }
            }
         }

         foreach(cell; "~cellName~".cell_list) {
             cell.fs.clear_imaginary_components();
             cell.U[ftl].clear_imaginary_components();
             cell.dUdt[ftl].clear_imaginary_components();
         }
         "~cellName~".fs.clear_imaginary_components();
         "~cellName~".U[ftl].clear_imaginary_components();
         "~cellName~".dUdt[ftl].clear_imaginary_components();
         foreach(face; "~cellName~".face_list) { face.fs.clear_imaginary_components(); } ";
         if ( wrtConserved ) {
             codeStr ~= ""~cellName~".decode_conserved(gtl, ftl, 0.0); ";
         }
         codeStr ~= "evalRHS(gtl, ftl, "~cellName~".cell_list, "~cellName~".face_list);
         }";
         return codeStr;
    } // end computeResidualSensitivity()

    void evalRHS(int gtl, int ftl, ref FVCell[] cell_list, FVInterface[] iface_list)
    /*
     *  This method evaluates the RHS residual on a subset of cells for a given FluidBlock.
     *  It is used when constructing the numerical Jacobian.
     *  Its effect should replicate evalRHS() in steadystatecore.d for a subset of cells.
     */
    {

        foreach(iface; iface_list) iface.F.clear();
        foreach(cell; cell_list) cell.clear_source_vector();

        bool do_reconstruction = ( flowJacobianT.spatial_order > 1 );

        convective_flux_phase0(do_reconstruction, gtl, cell_list, iface_list);
        convective_flux_phase1(do_reconstruction, gtl, cell_list, iface_list);

        foreach(f; iface_list) {
            if (f.is_on_boundary) { applyPostConvFluxAction(0.0, gtl, ftl, f); }
        }

        // Viscous flux update
        if (myConfig.viscous) {

            foreach(f; iface_list) {
                if (f.is_on_boundary) { applyPreSpatialDerivActionAtBndryFaces(0.0, gtl, ftl, f); }
            }

            // currently only for least-squares at faces
            // TODO: generalise for all spatial gradient methods
            foreach(c; cell_list) {
                c.grad.gradients_leastsq(c.cloud_fs, c.cloud_pos, c.ws_grad); // flow_property_spatial_derivatives(0);
            }

            // we need to average cell-centered spatial (/viscous) gradients to get approximations of the gradients
            // at the cell interfaces before the viscous flux calculation.
            if (myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                foreach(f; iface_list) {
                    f.average_cell_deriv_values(gtl);
                }
            }
            estimate_turbulence_viscosity(cell_list);
            viscous_flux(iface_list);
            foreach(f; iface_list) {
                if (f.is_on_boundary) { applyPostDiffFluxAction(0.0, gtl, ftl, f); }
            }
        }

        foreach (i, cell; cell_list) {
            cell.add_inviscid_source_vector(gtl, 0.0);
            if (myConfig.viscous) {
                cell.add_viscous_source_vector();
            }
            if (myConfig.reacting) {
                cell.add_chemistry_source_vector();
            }
            if (myConfig.udf_source_terms) {
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                // can't call a blk/this inside the FluidBlock class
                /*
                if (grid_type == Grid_t.structured_g31rid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                    i_cell = ijk_indices[0];
                    j_cell = ijk_indices[1];
                    k_cell = ijk_indices[2];
                }
                */
                addUDFSourceTermsToCell(myL, cell, 0,
                                        0.0, myConfig,
                                        id, i_cell, j_cell, k_cell);
            }
            cell.time_derivatives(gtl, ftl);
        }

    } // end evalRHS()

    void apply_transpose_jacobian_bcs() {
        /*
          This method accounts for the boundary conditions for the boundary cell entries in the transpose flow Jacobian.

          To calculate a boundary correction dRdU (or dRdQ), we calculate the sensitivity of a
          ghost cell flow state (q) with respect to a perturbation of the interior cells
          conserved quantity (U) or flow state (Q) depending on whether the Jacobian is
          calculated with respect to conserved or primitive flow variables. This sensitivity
          is denoted as dQdq or dUdQ. We also calculate the residual sensitivity of interior
          cells with respect to perturbing the ghost cell flow state, dRdq. Then we multiply
          these results to obtain the local sensitivity matrix to be added onto the current entries
          in the Jacobian, e.g. dRdU = dRdq * dqdU or dRdQ = dRdq * dqdQ.
          Refer to Kyle's thesis for more details.

          TODO: this function is a bit busy, think about breaking it up into smaller pieces.
         */

        // Make a stack-local copy of conserved quantities info
        size_t nConserved = GlobalConfig.cqi.nConservedQuantities;
        size_t MASS = GlobalConfig.cqi.mass;
        size_t X_MOM = GlobalConfig.cqi.xMom;
        size_t Y_MOM = GlobalConfig.cqi.yMom;
        size_t Z_MOM = GlobalConfig.cqi.zMom;
        size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
        size_t TKE = GlobalConfig.cqi.tke;
        size_t SPECIES = GlobalConfig.cqi.species;
        number EPS = flowJacobianT.eps;
        int gtl = 0; int ftl = 1;

        foreach ( bndary; bc ) {
            if ( bndary.type == "exchange_using_mapped_cells" || bndary.type == "exchange_over_full_face") { continue; }
            foreach ( bi, bface; bndary.faces) {
                FVCell ghost_cell; FVCell pcell;
                if (bndary.outsigns[bi] == 1) {
                    pcell = bface.left_cell;
                    ghost_cell = bface.right_cell;
                } else {
                    pcell = bface.right_cell;
                    ghost_cell = bface.left_cell;
                }
                pcell.U[ftl].copy_values_from(pcell.U[0]);

                // calculate the sensitivity of the ghost cell with respect to
                // a perturbation of the interior cell it shares an interface with
                if ( flowJacobianT.wrtConserved ) {
                    mixin(computeGhostCellSensitivity("mass", "MASS", true, true));
                    mixin(computeGhostCellSensitivity("momentum.refx", "X_MOM", false, true));
                    mixin(computeGhostCellSensitivity("momentum.refy", "Y_MOM", false, true));
                    if ( myConfig.dimensions == 3 ) { mixin(computeGhostCellSensitivity("momentum.refz", "Z_MOM", false, true)); }
                    mixin(computeGhostCellSensitivity("total_energy", "TOT_ENERGY", true, true));
                    foreach(ii; 0 .. myConfig.turb_model.nturb) { mixin(computeGhostCellSensitivity("rhoturb[ii]", "TKE+ii", false, true)); }
                    version(multi_species_gas){
                    if (myConfig.n_species > 1) {
                        foreach(ii; 0 .. myConfig.gmodel.n_species) { mixin(computeGhostCellSensitivity("massf[ii]", "SPECIES+ii", false, true)); }
                    }
                    }
                } else {
                    mixin(computeGhostCellSensitivity("gas.rho", "MASS", true, false));
                    mixin(computeGhostCellSensitivity("vel.refx", "X_MOM", false, false));
                    mixin(computeGhostCellSensitivity("vel.refy", "Y_MOM", false, false));
                    if ( myConfig.dimensions == 3 ) { mixin(computeGhostCellSensitivity("vel.refz", "Z_MOM", false, false)); }
                    mixin(computeGhostCellSensitivity("gas.p", "TOT_ENERGY", true, false));
                    foreach(ii; 0 .. myConfig.turb_model.nturb) { mixin(computeGhostCellSensitivity("turb[ii]", "TKE+ii", false, false)); }
                    version(multi_species_gas){
                    if (myConfig.n_species > 1) {
                        foreach(ii; 0 .. myConfig.gmodel.n_species) { mixin(computeGhostCellSensitivity("gas.massf[ii]", "SPECIES+ii", false, false)); }
                    }
                    }
                }

                // calculate the sensivitiy of the interior cells residuals with respect to
                // a perturbation of the ghost cell flow state
                mixin(computeResidualSensitivity("ghost_cell", "gas.rho", "MASS", true, false));
                mixin(computeResidualSensitivity("ghost_cell", "vel.refx", "X_MOM", false, false));
                mixin(computeResidualSensitivity("ghost_cell", "vel.refy", "Y_MOM", false, false));
                if ( myConfig.dimensions == 3 ) { mixin(computeResidualSensitivity("ghost_cell", "vel.refz", "Z_MOM", false, false)); }
                mixin(computeResidualSensitivity("ghost_cell", "gas.p", "TOT_ENERGY", true, false));
                foreach(ii; 0 .. myConfig.turb_model.nturb) { mixin(computeResidualSensitivity("ghost_cell", "turb[ii]", "TKE+ii", false, false)); }
                version(multi_species_gas){
                if (myConfig.n_species > 1) {
                    foreach(ii; 0 .. myConfig.gmodel.n_species) { mixin(computeResidualSensitivity("ghost_cell", "gas.massf[ii]", "SPECIES+ii", false, false)); }
                }
                }
                // multiply the sensitivity matrices and add the corrections to the relevant flow Jacobian entries
                foreach(bcell; ghost_cell.cell_list) { //
                    assert(bcell.id < ghost_cell_start_id, "Oops, we expect to not find a ghost cell at this point.");
                    for ( size_t ip = 0; ip < nConserved; ++ip ) {
                        for ( size_t jp = 0; jp < nConserved; ++jp ) {
                            number entry = bcell.dQdU[ip][jp];
                            flowJacobianT.dRdq[ip][jp] = entry;
                        }
                    }
                    // perform matrix-matrix multiplication
                    for (size_t i = 0; i < nConserved; i++) {
                        for (size_t j = 0; j < nConserved; j++) {
                            flowJacobianT.dRdQ[i][j] = 0;
                            for (size_t k = 0; k < nConserved; k++) {
                                flowJacobianT.dRdQ[i][j] += flowJacobianT.dRdq[i][k]*flowJacobianT.dqdQ[k][j];
                            }
                        }
                    }

                    // add correction to boundary entry in Jacobian
                    size_t I, J;
                    for ( size_t ip = 0; ip < nConserved; ++ip ) {
                        I = bcell.id*nConserved + ip; // column index
                        for ( size_t jp = 0; jp < nConserved; ++jp ) {
                            J = pcell.id*nConserved + jp; // row index
                            //writeln(I, ", ", J, ", ", flowJacobianT.local[J,I], ", ", flowJacobianT.dRdQ[ip][jp]);
                            flowJacobianT.local[J,I] = flowJacobianT.local[J,I] + flowJacobianT.dRdQ[ip][jp];
                        }
                    }
                }
            } // foreach ( bi, bface; bndary.faces)
        } // foreach ( bndary; bc )
    } // end apply_transpose_jacobian_bcs()

    static string computeGhostCellSensitivity(string varName, string posInArray, bool includeThermoUpdate, bool wrtConserved)
    {
        string codeStr = "{ ";
        if ( wrtConserved ) {
            codeStr ~= "pcell.U[ftl]."~varName~" += EPS;
                        pcell.decode_conserved(gtl, ftl, 0.0); ";
        } else {
            codeStr ~= "pcell.fs."~varName~" += EPS; ";
        }
        if ( includeThermoUpdate ) {
            codeStr ~= "myConfig.gmodel.update_thermo_from_rhop(pcell.fs.gas);";
        }

        codeStr ~= "
        if (bc[bface.bc_id].preReconAction.length > 0) { bc[bface.bc_id].applyPreReconAction(0.0, 0, 0, bface); }

        flowJacobianT.dqdQ[MASS][" ~ posInArray ~ "] = ghost_cell.fs.gas.rho.im/(EPS.im);
        flowJacobianT.dqdQ[X_MOM][" ~ posInArray ~ "] = ghost_cell.fs.vel.x.im/(EPS.im);
        flowJacobianT.dqdQ[Y_MOM][" ~ posInArray ~ "] = ghost_cell.fs.vel.y.im/(EPS.im);
        if (myConfig.dimensions == 3) { flowJacobianT.dqdQ[Z_MOM][" ~ posInArray ~ "] = ghost_cell.fs.vel.z.im/(EPS.im); }
        flowJacobianT.dqdQ[TOT_ENERGY][" ~ posInArray ~ "] = ghost_cell.fs.gas.p.im/(EPS.im);
        foreach(it; 0 .. myConfig.turb_model.nturb){ flowJacobianT.dqdQ[TKE+it][" ~ posInArray ~ "] = ghost_cell.fs.turb[it].im/(EPS.im); }
        version(multi_species_gas){
        if (myConfig.n_species > 1) {
            foreach(sp; 0 .. myConfig.gmodel.n_species){ flowJacobianT.dqdQ[SPECIES+sp][" ~ posInArray ~ "] = ghost_cell.fs.gas.massf[sp].im/(EPS.im); }
        }
        }
        foreach ( ref c; [pcell, ghost_cell] ) {
            c.fs.clear_imaginary_components();
            c.U[ftl].clear_imaginary_components();
            c.dUdt[ftl].clear_imaginary_components();
        }
        }";
        return codeStr;
    }

    // The following two methods are used to verify the numerical Jacobian implementation.
    void verify_transpose_jacobian()
    {

        // we perform a residual evaluation to ensure the ghost cells are filled with good data
        import steadystate_core;
        steadystate_core.evalRHS(0.0, 0);

        // calculate the numerical Jaacobian
        initialize_transpose_jacobian(1);
        evaluate_transpose_jacobian();
        assert(flowJacobianT !is null, "Oops, we expect a flowJacobianT object to be attached to the fluidblock.");
        size_t nConserved = GlobalConfig.cqi.nConservedQuantities;

        // temporarily change interpolation order
        shared int interpolation_order_save = GlobalConfig.interpolation_order;
        myConfig.interpolation_order = to!int(flowJacobianT.spatial_order);

        // create an arbitrary unit vector
        number[] vec;
        vec.length = cells.length*nConserved;
        foreach ( i, ref val; vec) { val = i+1; }

        // normalise the vector
        number norm = 0.0;
        foreach( i; 0..vec.length) { norm += vec[i]*vec[i]; }
        norm = sqrt(norm);
        foreach( ref val; vec) { val = val/norm; }

        // result vectors
        number[] sol1;
        sol1.length = vec.length;
        number[] sol2;
        sol2.length = vec.length;

        // explicit multiplication of J*vec
        foreach ( i; 0..vec.length) {
            sol1[i] = 0.0;
            foreach ( j; 0..vec.length) {
                sol1[i] += flowJacobianT.local[j,i]*vec[j];
            }
        }

        // Frechet derivative of J*vec
        steadystate_core.evalRHS(0.0, 0);
        evalConservativeJacobianVecProd(vec, sol2);

        // write out results
        string fileName = "jacobian_test.output";
        auto outFile = File(fileName, "w");
        foreach( i; 0..v.length ) {
            size_t id = i/nConserved;
            outFile.writef("%d    %d    %.16e    %.16e    %.16f    %.16f \n", i, id, fabs((sol1[i]-sol2[i])/sol1[i]).re, fabs(sol1[i]-sol2[i]).re, sol1[i].re, sol2[i].re);
        }

        // return the interpolation order to its original state
        myConfig.interpolation_order = interpolation_order_save;

        // stop the program at this point
        import core.runtime;
        Runtime.terminate();
    } // end verify_transpose_jacobian

    void evalConservativeJacobianVecProd(number[] vec, ref number[] sol) {
        size_t nConserved = GlobalConfig.cqi.nConservedQuantities;
        size_t MASS = GlobalConfig.cqi.mass;
        size_t X_MOM = GlobalConfig.cqi.xMom;
        size_t Y_MOM = GlobalConfig.cqi.yMom;
        size_t Z_MOM = GlobalConfig.cqi.zMom;
        size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
        size_t TKE = GlobalConfig.cqi.tke;
        size_t SPECIES = GlobalConfig.cqi.species;
        double EPS = flowJacobianT.eps.im;

        // We perform a Frechet derivative to evaluate J*D^(-1)v
        size_t nturb = myConfig.turb_model.nturb;
        size_t nsp = myConfig.gmodel.n_species;
        clear_fluxes_of_conserved_quantities();
        foreach (cell; cells) cell.clear_source_vector();
        int cellCount = 0;
        foreach (cell; cells) {
            cell.U[1].copy_values_from(cell.U[0]);
            cell.U[1].mass += complex(0.0, EPS*vec[cellCount+MASS].re);
            cell.U[1].momentum.refx += complex(0.0, EPS*vec[cellCount+X_MOM].re);
            cell.U[1].momentum.refy += complex(0.0, EPS*vec[cellCount+Y_MOM].re);
            if ( myConfig.dimensions == 3 ) { cell.U[1].momentum.refz += complex(0.0, EPS*vec[cellCount+Z_MOM].re); }
            cell.U[1].total_energy += complex(0.0, EPS*vec[cellCount+TOT_ENERGY].re);
            foreach(it; 0 .. nturb) { cell.U[1].rhoturb[it] += complex(0.0, EPS*vec[cellCount+TKE+it].re); }
            version(multi_species_gas){
            if (myConfig.n_species > 1) {
                foreach(sp; 0 .. nsp) { cell.U[1].massf[sp] += complex(0.0, EPS*vec[cellCount+SPECIES+sp].re); }
            }
            }
            cell.decode_conserved(0, 1, 0.0);
            cellCount += nConserved;
        }
        import steadystate_core;
        steadystate_core.evalRHS(0.0, 1);
        //evalRHS(cells, ifaces);
        cellCount = 0;
        foreach (cell; cells) {
            sol[cellCount+MASS] = cell.dUdt[1].mass.im/EPS;
            sol[cellCount+X_MOM] = cell.dUdt[1].momentum.x.im/EPS;
            sol[cellCount+Y_MOM] = cell.dUdt[1].momentum.y.im/EPS;
            if ( myConfig.dimensions == 3 ) { sol[cellCount+Z_MOM] = cell.dUdt[1].momentum.z.im/EPS; }
            sol[cellCount+TOT_ENERGY] = cell.dUdt[1].total_energy.im/EPS;
            foreach(it; 0 .. nturb) { sol[cellCount+TKE+it] = cell.dUdt[1].rhoturb[it].im/EPS; }
            version(multi_species_gas){
            if (myConfig.n_species > 1) {
                foreach(sp; 0 .. nsp) { sol[cellCount+SPECIES+sp] = cell.dUdt[1].massf[sp].im/EPS; }
            }
            }
            cellCount += nConserved;
        }
    }
    } // end version(shape_sensitivity)

    version(nk_accelerator) {
    void allocate_GMRES_workspace()
    {
        size_t nConserved = GlobalConfig.cqi.nConservedQuantities;
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
