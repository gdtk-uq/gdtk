// fluidblock.d
// Base class for blocks of cells, for use within the Eilmer flow solver.
// Peter J. 2014-07-18 first cut.
// Now a derived class from the Block base class
// Kyle A. Damm 2020-02-11

module lmr.fluidblock;

import std.algorithm;
import std.array;
import std.conv;
import std.format;
import std.json;
import std.math;
import std.stdio;
import std.string;
version(mpi_parallel) {
    import mpi;
}

import lmr.grid_motion;
import lmr.grid_motion_udf;
import lmr.grid_motion_shock_fitting;

import gas.luagas_model;
import gas;
import geom;
import kinetics;
import nm.bbla;
import nm.number;
import nm.smla;
import ntypes.complex;
import util.json_helper;
import util.lua;
import util.lua_service;

import lmr.bc;
import lmr.block;
import lmr.conservedquantities;
import lmr.coredata;
import lmr.flowgradients;
import lmr.flowstate;
import lmr.fluidfvcell;
import lmr.fvinterface;
import lmr.fvvertex;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.grid_motion;
import lmr.jacobian;
import lmr.mass_diffusion : MassDiffusionModel;
import lmr.lua_helper;
import lmr.sfluidblock; // needed for some special-case processing, below
import lmr.shockdetectors;
import lmr.user_defined_source_terms;

// version(diagnostics) {
// import plt = matplotlibd.pyplot;
// import std.format;
// }

// The flow solver handles structured- and unstructured-grid blocks via this base class.
// Mostly, we view the block as an unstructured bag of cells because that requires least
// knowledge in the calling code.
class FluidBlock : Block {
public:
    size_t ncells;
    size_t nfaces;
    Grid_t grid_type; // structured or unstructured
    bool may_be_turbulent; // if true, the selected turbulence model is active
                           // within this block.
    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.
    double mass_balance; // domain mass balance used to assess convergence of steady-state simulations
    number L2_residual; // L2 norm of the global residual
    number mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case
    ConservedQuantities Linf_residuals;
    number c_h, divB_damping_length; //divergence cleaning parameters for MHD
    int mncell;                 // number of monitor cells
    number[] initial_T_value; // for monitor cells to check against
    number[] jx,jy,jz,hs;          // temporary storage for computing diffusion fluxes
    //
    // Collections of cells, vertices and faces are held as arrays of references.
    // These allow us to conveniently work through the items via foreach statements.
    FluidFVCell[] cells;
    FVInterface[] faces;
    FVVertex[] vertices;
    BoundaryCondition[] bc; // collection of references to the boundary conditions
    //
    // Dense core datastructures
    FluidCellData celldata;
    FVInterfaceData facedata;

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
    FlowState* Lft;
    FlowState* Rght;
    //
    // Super-time-stepping parameters used when applying flexible stages per block
    int s_RKL; // number of super-step
    double dt_parab; // allowable parabolic time-step
    // list of vertex id's that makeup the fluidblock boundary
    // (used in the grid deformation methods in conjunction with
    // the shape sensitivity calculator).
    size_t[] boundaryVtxIndexList;
    //
    // Workspace for transient-solver implicit update.
    Matrix!double crhs;
    double[] dRUdU;
    ConservedQuantities U0save, RU0;
    //
    // Workspace for transient-solver implicit update and NK accelerator.
    // source terms for finite-rate chemistry
    number[] thermochem_source;
    //
    version(newton_krylov)
    {
    double dtMin;
    double omegaLocal;
    FlowState* fs_save;

    // storage for a precondition matrix
    FlowJacobian flowJacobian;

    // Work-space for Newton-Krylov accelerator
    // These arrays and matrices are directly tied to using the
    // GMRES iterative solver.
    SMatrix!double JcT; // transposed Jacobian (w.r.t conserved variables)
    ConservedQuantities maxRate, residuals;
    number[] maxVel; // used for computing scale factors for steady-state shock-fitting
    double normAcc, dotAcc;
    size_t nvars;
    Matrix!number Minv;
    double[] R, R0, dU, DinvR, r0, x0, rhs;
    double[] v, w, zed;
    double[] g0, g1;
    Matrix!double Q1;
    double[] VT;

    // memory for fGMRES
    double[] ZT; // the preconditioned Krylov vectors used by fGMRES
    double[] VT_inner; // memory for the Krylov vectors used in the inner iterations of fGMRES
    double[] inner_x0;
    double[] inner_g1;
    double[] inner_R;
    }

    this(int id, string label)
    {
        super(id, label);
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
        auto cqi = dedicatedConfig[id].cqi;
        Linf_residuals = new_ConservedQuantities(cqi.n);
        // Workspace for flux_calc method.
        Lft = new FlowState(dedicatedConfig[id].gmodel, dedicatedConfig[id].turb_model.nturb);
        Rght = new FlowState(dedicatedConfig[id].gmodel, dedicatedConfig[id].turb_model.nturb);
        // Workspace for implicit updates of the thermochemistry.
        version(multi_species_gas) {
            if (myConfig.reacting) {
                thermochem_source.length = cqi.n_species;
            }
            // Workspace for viscous diffusion
            if (dedicatedConfig[id].turb_model.isTurbulent ||
               (dedicatedConfig[id].mass_diffusion_model != MassDiffusionModel.none)) {
                jx.length = cqi.n_species;
                jy.length = cqi.n_species;
                jz.length = cqi.n_species;
                hs.length = cqi.n_species;
            }
        }
        version(multi_T_gas) {
            if (cqi.n_modes > 0) {
                thermochem_source.length += cqi.n_modes;
            }
        }
        version(newton_krylov) {
            fs_save = new FlowState(dedicatedConfig[id].gmodel, dedicatedConfig[id].turb_model.nturb);
        }
    }

    override string toString() const { return "Block(id=" ~ to!string(id) ~ ")"; }
    @nogc size_t globalCellId(size_t localCellId) { return globalCellIdStart + localCellId; }
    abstract JSONValue get_header();

    @nogc abstract int get_interpolation_order();
    @nogc abstract void set_interpolation_order(int order);
    abstract void init_lua_globals();
    abstract void init_boundary_conditions(JSONValue json_data);
    @nogc abstract void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found);
    abstract void init_grid_and_flow_arrays(string gridFileName);
    @nogc abstract void compute_primary_cell_geometric_data(size_t gtl);
    @nogc abstract void precompute_stencil_data(size_t gtl);
    @nogc abstract void compute_least_squares_setup(size_t gtl);
    @nogc abstract void sync_vertices_from_underlying_grid(size_t gtl=0);
    @nogc abstract void sync_vertices_to_underlying_grid(size_t gtl=0);
    abstract void read_new_underlying_grid(string fileName);
    abstract void write_underlying_grid(string fileName);
    @nogc abstract void propagate_inflow_data_west_to_east();
    @nogc abstract void set_face_flowstates_to_averages_from_cells();
    @nogc abstract void convective_flux_phase0_legacy(
                                               bool allow_high_order_interpolation, size_t gtl=0,
                                               FluidFVCell[] cell_list = [], FVInterface[] iface_list = [],
                                               FVVertex[] vertex_list = []);
    @nogc abstract void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0);
    @nogc abstract void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0,
                                               FluidFVCell[] cell_list = [], FVInterface[] iface_list = [],
                                               FVVertex[] vertex_list = []);
    @nogc abstract void convective_flux_phase2(bool allow_high_order_interpolation, size_t gtl=0,
                                               FluidFVCell[] cell_list = [], FVInterface[] iface_list = [],
                                               FVVertex[] vertex_list = []);
    abstract size_t[] get_cell_write_indices();
    abstract void eval_udf_source_vectors(double simTime, size_t gtl, size_t[] cell_idxs=[]);

    void allocate_dense_celldata(size_t ncells, size_t nghost, size_t neq, size_t nftl)
    {
    /*
        Both kinds of blocks now share a structure of densely packed core flow data. This routine
        allocates the storage for these structures, attempting to keep related bits of data together
        on the heap.

        In an older draft of this code, I allocated all of the normal cells in one loop,
        then added the ghost cells onto the end later. This caused all of the pointers
        in the FVCell classes to no longer be pointing to the right objects, which was
        very bad. The solution to this was to use .reserve on the arrays before
        filling them, which prevented the pointers from being copied. Now however,
        we do the cells and the ghost cells all at once, but we still do .reserve
        because it should be faster.

        @author: Nick Gibbons
    */
        auto gmodel = myConfig.gmodel;
        size_t nturb = myConfig.turb_model.nturb;
        size_t ngtl = myConfig.n_grid_time_levels;

        celldata.all_cell_idxs.length = ncells;
        celldata.halo_cell_ids.length = (ncells + nghost);
        celldata.halo_face_ids.length = (ncells + nghost);
        foreach(i; 0 .. ncells) celldata.all_cell_idxs[i] = i;

        celldata.nfaces.length = ncells;
        celldata.dt_local.length = ncells;
        celldata.wall_distances.length = ncells;
        celldata.data_is_bad.length = ncells;
        celldata.in_turbulent_zone.length = ncells;
        celldata.lengths.length = ncells + nghost;

        celldata.areas.length     = ngtl*(ncells + nghost);
        celldata.volumes.length   = ngtl*(ncells + nghost);
        celldata.positions.length = ngtl*(ncells + nghost);

        celldata.U0.length = (ncells + nghost)*neq*nftl;
        if (nftl>1) celldata.U1.length = (ncells + nghost)*neq*nftl;
        if (nftl>2) celldata.U2.length = (ncells + nghost)*neq*nftl;
        if (nftl>3) celldata.U3.length = (ncells + nghost)*neq*nftl;
        if (nftl>4) celldata.U4.length = (ncells + nghost)*neq*nftl;
        celldata.dUdt0.length = (ncells + nghost)*neq*nftl;
        if (nftl>1) celldata.dUdt1.length = (ncells + nghost)*neq*nftl;
        if (nftl>2) celldata.dUdt2.length = (ncells + nghost)*neq*nftl;
        if (nftl>3) celldata.dUdt3.length = (ncells + nghost)*neq*nftl;
        if (nftl>4) celldata.dUdt4.length = (ncells + nghost)*neq*nftl;
        celldata.source_terms.length = (ncells + nghost)*neq;

        celldata.flowstates.reserve(ncells + nghost);
        celldata.gradients.reserve(ncells + nghost);
        celldata.workspaces.reserve(ncells + nghost);
        foreach (n; 0 .. ncells+nghost) celldata.flowstates ~= FlowState(gmodel, nturb);
        foreach (n; 0 .. ncells+nghost) celldata.gradients ~= FlowGradients(myConfig);
        foreach (n; 0 .. ncells+nghost) celldata.workspaces ~= WLSQGradWorkspace();

        version(newton_krylov) {
            celldata.cell_jacobians.length = (ncells+nghost)*neq*neq;
            celldata.saved_gradients.reserve(ncells + nghost);
            foreach (n; 0 .. ncells+nghost) celldata.saved_gradients ~= FlowGradients(myConfig);
            celldata.saved_source_terms.length = (ncells + nghost)*neq;
            celldata.doNotPerturb.length = (ncells+nghost);
        }
    }

    void allocate_dense_facedata(size_t nfaces, size_t nbfaces, size_t neq, size_t nftl)
    {
    /*
        Both kinds of blocks now share a structure of densely packed core flow
        data. This routine allocates the storage for these structures,
        attempting to keep related bits of data together on the heap.

        @author: Nick Gibbons
    */
        auto gmodel = myConfig.gmodel;
        size_t nturb = myConfig.turb_model.nturb;
        size_t ngtl = myConfig.n_grid_time_levels;

        facedata.all_face_idxs.length = nfaces;
        foreach(i; 0 .. nfaces) facedata.all_face_idxs[i] = i;

        facedata.positions.length = nfaces;
        facedata.areas.length = nfaces*ngtl;
        facedata.grid_velocities.length = nfaces*ngtl;
        facedata.normals.length = nfaces;
        facedata.tangents1.length = nfaces;
        facedata.tangents2.length = nfaces;
        facedata.left_interior_only.length = nfaces;
        facedata.right_interior_only.length = nfaces;
        facedata.stencil_idxs.length = nfaces;
        facedata.fluxes.length = nfaces*neq;
        facedata.f2c.length = nfaces;

        if (myConfig.interpolation_order>=2) facedata.l2r2_interp_data.length = nfaces;
        if (myConfig.interpolation_order==3) facedata.l3r3_interp_data.length = nfaces;

        facedata.flowstates.reserve(nfaces);
        foreach (n; 0 .. nfaces) facedata.flowstates ~= FlowState(gmodel, nturb);
        facedata.gradients.reserve(nfaces);
        foreach (n; 0 .. nfaces) facedata.gradients ~= FlowGradients(myConfig);
        facedata.workspaces.reserve(nfaces);
        foreach (n; 0 .. nfaces) facedata.workspaces ~= WLSQGradWorkspace();
    }

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
            if (myConfig.suppress_reconstruction_zones.length > 0 ) {
                foreach(srz; myConfig.suppress_reconstruction_zones) {
                    if (srz.is_inside(f.pos, myConfig.dimensions)) {
                        f.in_suppress_reconstruction_zone = true;
                    }
                } // foreach srz
            }
            if (myConfig.suppress_reconstruction_at_boundaries && f.is_on_boundary) {
                if (!startsWith(bc[f.bc_id].type, "exchange_")) {
                    f.in_suppress_reconstruction_zone = true;
                }
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
    void identify_suppress_viscous_stresses_zones()
    // Adjust the in_suppress_reconstruction-zone flag for faces in this block.
    {
        foreach(f; faces) {
            f.in_suppress_viscous_stresses_zone = false;
            if (myConfig.suppress_viscous_stresses_zones.length > 0 ) {
                foreach(svgz; myConfig.suppress_viscous_stresses_zones) {
                    if (svgz.is_inside(f.pos, myConfig.dimensions)) {
                        f.in_suppress_viscous_stresses_zone = true;
                    }
                } // foreach srz
            }
        } // foreach face
        //
        debug {
            size_t total_faces_in_suppress_viscous_stresses_zones = 0;
            size_t total_faces = 0;
            foreach(f; faces) {
                total_faces_in_suppress_viscous_stresses_zones += (f.in_suppress_viscous_stresses_zone ? 1: 0);
                total_faces += 1;
            } // foreach face
            if (myConfig.verbosity_level >= 2) {
                writeln("identify_suppress_viscous_stresses_zones(): block ", id,
                        " faces inside zones = ", total_faces_in_suppress_viscous_stresses_zones,
                        " out of ", total_faces);
                if (myConfig.suppress_viscous_stresses_zones.length == 0) {
                    writeln("Note that for no user-specified zones,",
                            " nominally viscous stresses are allowed for the whole domain",
                            " but may be suppressed at the boundaries and/or near the axis of symmetry.");
                }
            }
        }
    } // end identify_suppress_viscous_stresses_zones()

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
    void estimate_turbulence_viscosity(FluidFVCell[] cell_list = [])
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
    void evaluate_electrical_conductivity(FluidFVCell[] cell_list = [])
    {
        if (cell_list.length == 0) { cell_list = cells; }
        foreach (cell; cell_list) {
            cell.evaluate_electrical_conductivity();
        }
    }

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
        if (comp_tol>0.0) throw new Error("compression_tolerance should be negative!");

        final switch (myConfig.shock_detector) {
        case ShockDetector.PJ:
            foreach (iface; faces) {
                iface.fs.S = PJ_ShockDetector(iface, comp_tol, shear_tol);
            }
            break;
        case ShockDetector.NNG:
            // Re-use comp_tol as our tunable parameter. It needs to be
            // positive however. From the users point of view, larger (more
            // negative) values make the shock detector more sensitive, and
            // hence on more often. This is consistent with the PJ behaviour.
            auto gm = myConfig.gmodel;
            double Mx = -1.0*comp_tol.re;
            foreach (iface; faces) {
                iface.fs.S = NNG_ShockDetector(gm, *(iface.left_cell.fs),
                                                   *(iface.right_cell.fs),
                                                   iface.n, Mx);
            }
            break;
        case ShockDetector.KAD:
            double weight = comp_tol.re;
            foreach (iface; faces) {
                iface.fs.S = KAD_ShockDetector(*(iface.left_cell.fs),
                                               *(iface.right_cell.fs),
                                               weight);
            }
            break;
        }

        // Set the outflow interfaces to be shocked for high-order simulations with boundary layers
        // Drops back to the more robust flux calculators at the outflow to help prevent numerical noise
        // propagating back upstream
        if (myConfig.damped_outflow) {
            foreach (bndry; bc) {
                if (bndry.type == "outflow_simple_extrapolate") {
                    foreach (face; bndry.faces) {
                        face.fs.S = 1;
                    }
                }
            }
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
            if (!face.left_cell || !face.left_cell.is_interior_to_domain) {
                face.fs.S = face.right_cell.fs.S;
            } else if (!face.right_cell || !face.right_cell.is_interior_to_domain) {
                face.fs.S = face.left_cell.fs.S;
            } else {
                face.fs.S = fmax(face.left_cell.fs.S, face.right_cell.fs.S);
            }
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
    // Note: Changed 26/07/22 by NNG to add guards for missing ghost cells.
    //       Consider revising for efficiency, by preallocating a list of left onlies
    //       and a list of right onlies, and a list of faces with both.
    {
        foreach (face; faces) {
            if (!face.left_cell || !face.left_cell.is_interior_to_domain) {
                face.fs.S = face.right_cell.fs.S;
            } else if (!face.right_cell || !face.right_cell.is_interior_to_domain) {
                face.fs.S = face.left_cell.fs.S;
            } else {
                face.fs.S = fmax(face.left_cell.fs.S, face.right_cell.fs.S);
            }
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
                    FlowState*[] neighbour_flows;
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
                        vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, *(vtx.ws_grad));
                    }
                    break;
                case SpatialDerivCalc.divergence:
                    foreach(vtx; vertices) {
                        vtx.grad.gradients_xy_div(myConfig, vtx.cloud_fs, vtx.cloud_pos);
                    }
                } // end switch
            } else {
                // Have only least-squares in 3D.
                foreach(vtx; vertices) {
                    vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, *(vtx.ws_grad));
                }
            }

            // We've finished computing gradients at vertices, now copy them around if needed
            // Interfaces need them for reconstruction/viscous fluxes
            foreach (iface; faces) {
                iface.average_vertex_deriv_values();
            }
            // Turbulence models will need cell centered gradients (also for axisymmetric!)
            if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent || myConfig.save_viscous_gradients)){
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
                        iface.grad.gradients_leastsq(myConfig, iface.cloud_fs, *(iface.ws_grad));
                    }
                    break;
                case SpatialDerivCalc.divergence:
                    foreach(iface; faces) {
                        iface.grad.gradients_xy_div(myConfig, iface.cloud_fs, iface.cloud_pos);
                    }
                } // end switch
            } else { //3D
                final switch (myConfig.spatial_deriv_calc) {
                case SpatialDerivCalc.least_squares:
                    foreach(iface; faces) {
                        iface.grad.gradients_leastsq(myConfig, iface.cloud_fs, *(iface.ws_grad));
                    }
                    break;
                case SpatialDerivCalc.divergence:
                    foreach(iface; faces) {
                        assert(0, "divergence thereom not implemented when computing gradients at faces in 3D");
                    }
                } // end switch
            } // end if (myConfig.dimensions)

            // Finished computing gradients at interfaces, now copy them around if needed
            // Turbulence models and axisymmetric source terms need cell centered gradients
            if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent || myConfig.save_viscous_gradients)){
                foreach (cell; cells) {
                    cell.average_interface_deriv_values();
                }
            }
            break;

        case SpatialDerivLocn.cells:
            final switch (myConfig.spatial_deriv_calc) {
            case SpatialDerivCalc.least_squares:
                compute_gradients_at_cells_leastsq();
                break;
            case SpatialDerivCalc.divergence:
                foreach(iface; faces) {
                    assert(0, "divergence thereom not implemented when computing gradients at cell centres.");
                }
            } // end switch
            //
            // Cell centered gradients are transformed to interfaces in simcore.
            //
        } // end switch (myConfig.spatial_deriv_locn)
    } // end flow_property_spatial_derivatives()

    @nogc
    void compute_gradients_at_cells_leastsq(size_t[] cell_idxs=[])
    {
        immutable bool is3D = (myConfig.dimensions == 3);
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool doSpecies = myConfig.turb_model.isTurbulent || myConfig.mass_diffusion_model != MassDiffusionModel.none;

        if (cell_idxs.length==0) cell_idxs = celldata.all_cell_idxs;
        foreach(idx; cell_idxs){
            celldata.gradients[idx].gradients_at_cells_leastsq(
                celldata.flowstates[idx], facedata.flowstates, celldata.c2f[idx],
                celldata.workspaces[idx], celldata.nfaces[idx],
                is3D, nsp, nmodes, nturb, doSpecies);
        }
    }

    @nogc
    void clear_fluxes_of_conserved_quantities()
    {
        foreach (iface; faces) { iface.F.clear(); }
    }

    @nogc
    void average_turbulent_transprops_to_faces(FVInterface[] face_list)
    {
        if (myConfig.turb_model.isTurbulent){
            foreach (iface; face_list) { iface.average_turbulent_transprops(); }
        }
    }

    @nogc
    void average_turbulent_transprops_to_faces()
    {
        if (myConfig.turb_model.isTurbulent){
            foreach (iface; faces) { iface.average_turbulent_transprops(); }
        }
    }

    @nogc
    void viscous_flux(FVInterface[] face_list)
    {
        foreach (iface; face_list) { iface.viscous_flux_calc(); }
    }

   @nogc
   void viscous_flux(size_t[] face_idxs=[])
   {
       immutable size_t n_species      = myConfig.n_species;
       immutable size_t n_modes        = myConfig.n_modes;
       immutable size_t nturb          = myConfig.turb_model.nturb;
       immutable bool is3d             = myConfig.dimensions == 3;
       immutable bool isTurbulent      = myConfig.turb_model.isTurbulent;
       immutable bool axisymmetric     = myConfig.axisymmetric;
       immutable double viscous_factor = myConfig.viscous_factor;
       immutable size_t neq            = myConfig.cqi.n;
       immutable bool laminarDiffusion = myConfig.mass_diffusion_model != MassDiffusionModel.none;
       immutable double Sc_t = myConfig.turbulence_schmidt_number;

       if (face_idxs.length==0) face_idxs = facedata.all_face_idxs;
       foreach(fid; face_idxs){
           navier_stokes_viscous_fluxes(facedata.flowstates[fid], facedata.gradients[fid],
                             myConfig, n_modes, nturb,
                             isTurbulent, axisymmetric, is3d,
                             facedata.normals[fid], facedata.positions[fid].y.re,
                             facedata.fluxes[fid*neq .. (fid+1)*neq]);
       }
       if ((isTurbulent || laminarDiffusion)&&(n_species>1)) {
           foreach(fid; face_idxs){
               diffusion_viscous_fluxes(facedata.flowstates[fid], facedata.gradients[fid],
                                 myConfig, n_species, n_modes,
                                 Sc_t, laminarDiffusion, isTurbulent,
                                 facedata.normals[fid], jx, jy, jz, hs,
                                 facedata.fluxes[fid*neq .. (fid+1)*neq]);
           }
       }
   }

    @nogc
    void init_residuals()
    // Initialization of data for later computing residuals.
    {
        auto cqi = myConfig.cqi;
        mass_residual = 0.0;
        mass_residual_loc.clear();
        energy_residual = 0.0;
        energy_residual_loc.clear();
        foreach(FluidFVCell cell; cells) {
            cell.rho_at_start_of_step = cell.fs.gas.rho;
            cell.rE_at_start_of_step = cell.U[0][cqi.totEnergy];
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
        auto cqi = myConfig.cqi;
        mass_residual = 0.0;
        mass_residual_loc.clear();
        energy_residual = 0.0;
        energy_residual_loc.clear();
        foreach(FluidFVCell cell; cells) {
            number local_residual = (cell.fs.gas.rho - cell.rho_at_start_of_step) / cell.fs.gas.rho;
            local_residual = fabs(local_residual);
            if ( local_residual > mass_residual ) {
                mass_residual = local_residual;
                mass_residual_loc.set(cell.pos[gtl]);
            }
            // In the following line, the zero index is used because,
            // at the end of the gas-dynamic update, that index holds
            // the updated data.
            number currentEnergy = cell.U[0][cqi.totEnergy];
            local_residual = (currentEnergy - cell.rE_at_start_of_step) / currentEnergy;
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
        auto cqi = myConfig.cqi;
        Linf_residuals.copy_values_from(cells[0].dUdt[0]);
        foreach (i; 0 .. cqi.n) {
            Linf_residuals[i] = fabs(Linf_residuals[i]);
            foreach (cell; cells) {
                Linf_residuals[i] = fmax(Linf_residuals[i], fabs(cell.dUdt[0][i]));
            }
        }
    } // end compute_Linf_residuals()

    @nogc
    void compute_L2_residual()
    {
        auto cqi = myConfig.cqi;
        L2_residual = 0.0;
        foreach (cell; cells) {
            if (cqi.mass==0) L2_residual += fabs(cell.dUdt[0][cqi.mass])^^2;
	    L2_residual += fabs(cell.dUdt[0][cqi.xMom])^^2;
	    L2_residual += fabs(cell.dUdt[0][cqi.yMom])^^2;
	    if (cqi.threeD) { L2_residual += fabs(cell.dUdt[0][cqi.zMom])^^2; }
	    L2_residual += fabs(cell.dUdt[0][cqi.totEnergy])^^2;
	    version(turbulence) {
                foreach(it; 0 .. myConfig.turb_model.nturb){
                    L2_residual += fabs(cell.dUdt[0][cqi.rhoturb+it])^^2;
                }
	    }
            version(multi_species_gas) {
                foreach(isp; 0 .. myConfig.gmodel.n_species){
                    L2_residual += fabs(cell.dUdt[0][cqi.species+isp])^^2;
                }
	    }
            version(multi_T_gas) {
                foreach(imode; 0 .. myConfig.gmodel.n_modes){
                    L2_residual += fabs(cell.dUdt[0][cqi.modes+imode])^^2;
                }
	    }
        }
    } // end compute_L2_residual()

    @nogc
    void compute_mass_balance()
    {
        auto cqi = myConfig.cqi;
        mass_balance = 0.0;
        foreach(boundary; bc) {
            if (boundary.type != "exchange_over_full_face" && boundary.type != "exchange_using_mapped_cells") {
                foreach(i, face; boundary.faces) {
                    double massflux=0.0;
                    if (cqi.mass==0) {
                        massflux = face.F[cqi.mass].re;
                    } else {
                        foreach(isp; 0 .. cqi.n_species) massflux += face.F[cqi.species+isp].re;
                    }
                    mass_balance += boundary.outsigns[i] * massflux * face.area[0].re;
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
                    auto other_cell = (c.outsign[i] == 1) ? f.right_cell : f.left_cell;
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
                        auto other_cell = (c.outsign[i] == 1) ? f.right_cell : f.left_cell;
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
        bool gridFlag = grid_type == Grid_t.structured_grid;
        foreach(FluidFVCell cell; cells) {
            // Search for the minimum length scale and the maximum CFL value in the block.
            if (first) {
                min_L_for_block = cell.L_min.re;
                cfl_local = cell.signal_frequency(gridFlag) * dt_current;
                cfl_max = cfl_local;
                first = false;
            } else {
                min_L_for_block = fmin(cell.L_min.re, min_L_for_block);
                cfl_local = cell.signal_frequency(gridFlag) * dt_current;
                cfl_max = fmax(cfl_local, cfl_max);
            }
        }
        return cfl_max * min_L_for_block / dt_current;
    } // end update_c_h()

    @nogc
    double[3] determine_time_step_size(double dt_current, double cfl_value, bool check_cfl)
    // Compute the local time step limit for all cells in the block.
    // The overall time step is limited by the worst-case cell.
    {
        // for STS (hyp = hyperbolic/convective, parab = parabolic/viscous)
        double signal_hyp;
        double signal_parab;
	double dt_allow_hyp;
	double dt_allow_parab;
        //
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
        if (is_explicit_update_scheme(myConfig.gasdynamic_update_scheme)) {
            switch (number_of_stages_for_update_scheme(myConfig.gasdynamic_update_scheme)) {
            case 1: cfl_allow = 0.9; break;
            case 2: cfl_allow = 1.2; break;
            case 3: cfl_allow = 1.6; break;
            default: cfl_allow = 0.9;
            }
            // When using implicit residual smoothing we should be able to achieve a higher stable CFL
            // so let's relax the cfl_allow value..
            if (myConfig.residual_smoothing && myConfig.with_local_time_stepping &&
                GlobalConfig.residual_smoothing_type == ResidualSmoothingType.implicit) {
                cfl_allow *= 10.0;
            }
        } else {
            // [TODO] PJ 2021-05-17 Implicit update schemes should run with cfl >> 1
            // but we don't have much experience to know how far the cfl can be pushed.
            // 100.0 seems a little too far for the LaRC flat plate with turbulent Mach 5 flow.
            // 2021-09-17 It has been reported that 500.0 is ok for a nozzle flow.
            cfl_allow = 600.0;
        }
        // for local time-stepping we limit the larger time-steps by a factor of the smallest timestep
        int local_time_stepping_limit_factor = myConfig.local_time_stepping_limit_factor;
        bool first = true;
        bool gridFlag = grid_type == Grid_t.structured_grid;
        foreach(FluidFVCell cell; cells) {
            signal = cell.signal_frequency(gridFlag);
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
                // Set the allowable time-step based on hyperbolic time-step.
                dt_allow = fmin(dt_allow_hyp, GlobalConfig.dt_max);
                cell.dt_hyper = cfl_value / signal_hyp;
                cell.dt_parab = cfl_value / signal_parab;
            } else {
                // no STS
                dt_allow_hyp = 0;
                dt_allow_parab = 0;
                // Compute current (Local) CFL number and recommend a time step.
                cfl_local = dt_current * signal;
                dt_local = cfl_value / signal;
                // Set local time-step in cell.
                cell.dt_local = fmin(dt_local, dt_current*local_time_stepping_limit_factor);
                // Limit the largest local time-step to a set input value.
                cell.dt_local = fmin(cell.dt_local, GlobalConfig.dt_max);
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
            debug { writeln("Bad cfl number encountered cfl_max=", cfl_max, " for FluidBlock ", id); }
            // If cfl_max exceeds cfl_allow, simply reduce the cfl_max to cfl_adjust*cfl_allow.
            // For the low-order explicit updating schemes, a value of 0.5 seems
            // to work robustly. Values to 0.7 also work, beyond this
            // and code begins to crash due to numerical instability.
            // Results in auto-limitation of the time step/cfl
            cfl_max = cfl_adjust*cfl_allow;
            // Reduce dt_allow according to new rescaled cfl_max
            dt_allow = cfl_max/signal;
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

    version(newton_krylov) {

    void initialize_jacobian(int spatial_order_of_jacobian, double sigma, int iluFill=-1)
    {
        /*
          This method initializes the flow Jacobian matrix attached the FluidBlock object.
          We gather the cell/interface residual stencils for the interior and ghost cells at this point,
          since we need to know the number of expected entries in the Jacobian matrix to pre-size the
          sparse matrix arrays.
         */

        size_t nentry = 0;
        GasModel gmodel = cast(GasModel) myConfig.gmodel;
        if (gmodel is null) { gmodel = GlobalConfig.gmodel_master; }
        fs_save = new FlowState(gmodel, myConfig.turb_model.nturb);

        // gather the expected number of non-zero entries in the flow Jacobian
        foreach (cell; cells) {
            if (cell.doNotPerturb) {
                cell.cell_list ~= cell; 
            }
            else {
                cell.gather_residual_stencil_lists(spatial_order_of_jacobian);
            }
            nentry += cell.cell_list.length;
        }

        // we will gather the ghost cell residual stencil lists
        // at this point as well for convenience, we need them
        // when we apply the boundary condition corrections later
        foreach ( bndary; bc ) {
            if (!bndary.ghost_cell_data_available) { continue; }
            if (bndary.type == "exchange_using_mapped_cells" ||
                bndary.type == "exchange_over_full_face" ||
                bndary.type == "do_nothing") {
                continue;
            }
            foreach ( iface, face; bndary.faces) {
                FluidFVCell ghost_cell, cell;
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

        jacobian_nonzero_pattern(spatial_order_of_jacobian, nentry, cells.length, sigma, iluFill);
    } // end initialize_jacobian()

    void jacobian_nonzero_pattern(int spatial_order_of_jacobian, size_t nentry, size_t ncells, double sigma, int iluFill=-1) {
        /*
          This routine will fill out the sparse matrix with the appropriate non-zero structure
        */

        // We first construct a temporary sparse matrix to store the point structure of the system
        // i.e. each nConserved x nConserved block is represented by a single entry
        FlowJacobian ptJac = new FlowJacobian(0.0, myConfig.dimensions, 1, spatial_order_of_jacobian, nentry, ncells);
        ptJac.local.ia[ptJac.ia_idx] = 0; // the first entry will always be filled
        ptJac.ia_idx += 1;
        foreach (pcell; cells) {
            // loop through cells that will have non-zero entries
            foreach(cell; pcell.cell_list) {
                assert(!cell.is_ghost, "Oops, we expect to not find a ghost cell at this point.");
                size_t jidx = cell.id; // column index
                // populate entry with a place holder value
                ptJac.local.aa[ptJac.aa_idx] = 1.0;
                // fill out the sparse matrix ready for the next entry in the row
                ptJac.aa_idx += 1;
                ptJac.local.ja[ptJac.ja_idx] = jidx;
                ptJac.ja_idx += 1;
            }
            // prepare the sparse matrix for a new row
            ptJac.local.ia[ptJac.ia_idx] = ptJac.aa_idx;
            ptJac.ia_idx += 1;
        }

        // For ILU(p>0) we need to modify the sparsity pattern
        if (iluFill > 0) {

            // [TODO] think about making the lev matrix sparse as well KAD 2022-03-31
            // construct a level matrix
            int p = iluFill;
            int n = to!int(ptJac.local.ia.length)-1;
            int[][] lev; // fill levels
            lev.length = n;
            foreach ( i; 0..n) lev[i].length = n;
            // assign initial fill levels
            foreach ( i; 0 .. n ) {
                foreach ( j; 0 .. n ) {
                    if (ptJac.local[i,j] < 1.0) lev[i][j] = n-1;
                    else lev[i][j] = 0;
                }
            }

            // apply symbolic phase algorithm
            foreach ( i; 1 .. n ) { // Begin from 2nd row
                foreach ( k; 0 .. i ) {
                    if (lev[i][k] <= p) {
                        foreach ( j ; k+1..n) {
                            lev[i][j] = min(lev[i][j], lev[i][k]+lev[k][j]+1);
                        }
                    }
                }
            }

            // modify the sparse matrix non-zero pattern based on the level matrix
            foreach ( i; 0..n) {
                foreach ( j; 0..n) {
                    // add new entry
                    if (lev[i][j] <= p) { ptJac.local[i,j] = 1.0; }
                }
            }
        }

        // we now construct the full sparse matrix by replacing each entry with an nConserved x nConserved block
        auto cqi = myConfig.cqi;
        auto nConserved = cqi.n;

        size_t n_vtx = (myConfig.grid_motion != GridMotion.none) ? vertices.length : 0;
        flowJacobian = new FlowJacobian(sigma, myConfig.dimensions, nConserved, spatial_order_of_jacobian, ptJac.local.aa.length, ncells, n_vtx);
        flowJacobian.local.ia[flowJacobian.ia_idx] = 0;
        flowJacobian.ia_idx += 1;
        size_t n = ptJac.local.ia.length-1;
        foreach ( rowi; 0 .. n ) { // foreach row in the point sparse matrix
            // loop through nConserved rows
            for (size_t ip = 0; ip < nConserved; ++ip) {
                foreach (colj; ptJac.local.ja[ptJac.local.ia[rowi] .. ptJac.local.ia[rowi+1]] ) { // foreach non-zero col in the point sparse matrix
                    for (size_t jp = 0; jp < nConserved; ++jp) {
                        size_t jidx = colj*nConserved + jp; // column index
                        // populate entry with a place holder value
                        flowJacobian.local.aa[flowJacobian.aa_idx] = 1.0;
                        // fill out the sparse matrix ready for the next entry in the row
                        flowJacobian.aa_idx += 1;
                        flowJacobian.local.ja[flowJacobian.ja_idx] = jidx;
                        flowJacobian.ja_idx += 1;
                    }
                }
                // prepare the sparse matrix for a new row
                flowJacobian.local.ia[flowJacobian.ia_idx] = flowJacobian.aa_idx;
                flowJacobian.ia_idx += 1;
            }
        }

        /*
        // Output some handy debugging information about the sparse matrices
        writef("blk[%d] has %d entries in the sparse matrix \n", id, ptJac.local.aa.length);
        writef("blk[%d] has %d entries in the block sparse matrix \n", id, flowJacobian.local.aa.length);

        // point sparse matrix
        string filename = "b" ~ to!string(id) ~ "_ilu.dat";
        File outFile;
        outFile = File(filename, "w");
        foreach(val; ptJac.local.aa) {
            outFile.writef("%.16f    ", val.re);
        }
        outFile.writef("\n");
        foreach(val; ptJac.local.ja) {
            outFile.writef("%d    ", val);
        }
        outFile.writef("\n");
        foreach(val; ptJac.local.ia) {
            outFile.writef("%d    ", val);
        }

        // full sparse matrix
        filename = "b" ~ to!string(id) ~ "_bfilu.dat";
        File outFile2;
        outFile2 = File(filename, "w");
        foreach(val; flowJacobian.local.aa) {
            outFile2.writef("%.16f    ", val.re);
        }
        outFile2.writef("\n");
        foreach(val; flowJacobian.local.ja) {
            outFile2.writef("%d    ", val);
        }
        outFile2.writef("\n");
        foreach(val; flowJacobian.local.ia) {
            outFile2.writef("%d    ", val);
        }
        */
    } // end jacobian_nonzero_pattern()

    void evaluate_jacobian()
    {
        /*
          Higher level method used to evaluate the flow Jacobian attached to the FluidBlock object.
         */

        // zero entries
        flowJacobian.local.aa[] = 0.0;

        // temporarily change flux calculator
        auto flux_calculator_save = myConfig.flux_calculator;
        myConfig.flux_calculator = myConfig.flux_calculator; // TODO: add preconditionMatrixFluxCalculator back

        // temporarily change interpolation order
        shared int interpolation_order_save = GlobalConfig.interpolation_order;
        myConfig.interpolation_order = to!int(flowJacobian.spatial_order);

        // copy some data for later use
        if (myConfig.viscous) {
            foreach(cell; cells) { cell.grad_save.copy_values_from(*(cell.grad)); }
        }
        if (flowJacobian.spatial_order >= 2) {
            foreach(cell; cells) { cell.gradients_save.copy_values_from(*(cell.gradients)); }
        }
        if (myConfig.reacting) {
            foreach(cell; cells) {
                cell.clear_source_vector();
                cell.add_thermochemical_source_vector(thermochem_source, 1.0);
                cell.Q_save.copy_values_from(cell.Q);
            }
        }

        // the real-valued finite difference needs a base residual (R0)
        version(complex_numbers) { } // do nothing
        else {
            foreach(cell; cells) {
                if (cell.doNotPerturb) continue;
                evalRHS(0, 0, cell.cell_list, cell.face_list, cell);
            }
        }

        // fill out the rows of the Jacobian for a cell
        foreach(cell; cells) {
            if (cell.doNotPerturb) continue;
                /*
                debug {
                    writefln("cell.id= %d, pos= %s ", cell.id, cell.pos[0]);
                    foreach (cl; cell.cell_list) writef("[%d: %s] ", cl.id, cl.pos[0]);
                    writeln();
                }
                */
            evaluate_cell_contribution_to_jacobian(cell);
        }

        // add boundary condition corrections to boundary cells
        apply_jacobian_bcs();

        // return flux calculator to its original state
        myConfig.flux_calculator = flux_calculator_save;

        // return the interpolation order to its original state
        myConfig.interpolation_order = interpolation_order_save;
    } // end evaluate_jacobian()

    void evaluate_cell_contribution_to_jacobian(FluidFVCell pcell)
    {
        auto cqi = myConfig.cqi; // was GlobalConfig.cqi;
        auto nConserved = cqi.n;
        auto eps0 = flowJacobian.eps;
        number eps;
        int ftl = 1; int gtl = 0;

        // save a copy of the flowstate and copy conserved quantities
        fs_save.copy_values_from(*(pcell.fs));
        pcell.U[ftl].copy_values_from(pcell.U[0]);

        // perturb the current cell's conserved quantities
        // and then evaluate the residuals for each cell in
        // the local domain of influence
        foreach(j; 0..nConserved) {

            // peturb conserved quantity
            version(complex_numbers) { eps = complex(0.0, eps0.re); }
            else { eps = eps0*fabs(pcell.U[ftl][j]) + eps0; }
            pcell.U[ftl][j] += eps;
            pcell.decode_conserved(gtl, ftl, 0.0);

            // evaluate perturbed residuals in local stencil
            evalRHS(gtl, ftl, pcell.cell_list, pcell.face_list, pcell);

            // fill local Jacobians
            foreach (cell; pcell.cell_list) {
                foreach(i; 0..nConserved) {
                    version(complex_numbers) { cell.dRdU[i][j] = cell.dUdt[ftl][i].im/eps.im; }
                    else { cell.dRdU[i][j] = (cell.dUdt[ftl][i]-cell.dUdt[0][i])/eps; }
                }
            }

            // return cell to original state
            version(complex_numbers) {
                pcell.U[ftl].copy_values_from(pcell.U[0]);
                pcell.fs.copy_values_from(*fs_save);
                if (myConfig.viscous) {
                    foreach (cell; pcell.cell_list) { cell.grad.copy_values_from(*(cell.grad_save)); }
                }
                if (flowJacobian.spatial_order >= 2) {
                    foreach(cell; pcell.cell_list) { cell.gradients.copy_values_from(*(cell.gradients_save)); }
                }
            } else {
                // TODO: for structured grids employing a real-valued low-order numerical Jacobian as
                // the precondition matrix the above code doesn't completely clean up the effect of the
                // perturbuation. Note that this does not occur for unstructured grids.
                // This is currently under investigation, in the interim we will apply the more costly
                // method of re-evaluating evalRHS with the unperturbed conserved state. KAD 2023-08-31
                pcell.U[ftl].copy_values_from(pcell.U[0]);
                pcell.decode_conserved(gtl, 0, 0.0);
                evalRHS(gtl, 0, pcell.cell_list, pcell.face_list, pcell);
            }
        }

        // we now populate the pre-sized sparse matrix representation of the flow Jacobian
        size_t jidx; // column index into the matrix
        // loop through nConserved rows
        for (size_t ip = 0; ip < nConserved; ++ip) {
            // loop through cells that will have non-zero entries
            foreach(cell; pcell.cell_list) {
                // loop through nConserved columns for each effected cell
                for ( size_t jp = 0; jp < nConserved; ++jp ) {
                    assert(!cell.is_ghost, "Oops, we expect to not find a ghost cell at this point.");
                    size_t I = cell.id*nConserved + ip;
                    size_t J = pcell.id*nConserved + jp;
                    flowJacobian.local[I,J] = cell.dRdU[ip][jp].re;
                }
            }
        }
    } // end evaluate_cell_contribution_to_jacobian()
    void apply_jacobian_bcs() {
        /*
          This method accounts for the boundary conditions for the boundary cell entries in the flow Jacobian.

          To calculate the boundary correction terms dR/dU:

          Step 1. calculate du/dU
          The sensitivity of the ghost cell conserved quantity (u) with respect to a perturbation of the interior cells conserved quantity (U).
          Step 2. calculate dR/du
          The sensitivity of the interior cells residuals (R) with respect to a perturbation of the ghost cells conserved quantity (U).
          Step 3. calculate dR/dU
          dRdU = dRdu * dudU

          Refer to Kyle's thesis for more details.
          TODO: this function is a bit busy, think about breaking it up into smaller pieces.
         */
        auto cqi = myConfig.cqi; // The block object has its own.
        auto nConserved = cqi.n;
        auto eps0 = flowJacobian.eps;
        number eps;
        int gtl = 0; int ftl = 1;

        foreach ( bndary; bc ) {
            if (!bndary.ghost_cell_data_available) { continue; }
            if (bndary.type == "exchange_using_mapped_cells" ||
                bndary.type == "exchange_over_full_face" ||
                bndary.type == "do_nothing") {
                continue;
            }
            foreach ( bi, bface; bndary.faces) {
                FluidFVCell ghost_cell, pcell;
                if (bndary.outsigns[bi] == 1) {
                    pcell = bface.left_cell;
                    ghost_cell = bface.right_cell;
                } else {
                    pcell = bface.right_cell;
                    ghost_cell = bface.left_cell;
                }

                // Step 1. Calculate du/dU
                // u: ghost cell conserved quantities
                // U: interior cell conserved quantities

                // save a copy of the flowstate and copy conserved quantities
                fs_save.copy_values_from(*(pcell.fs));
                pcell.U[ftl].copy_values_from(pcell.U[0]);
                ghost_cell.encode_conserved(gtl, 0, 0.0);

                foreach(idx; 0..nConserved) {

                    // peturb conserved quantity
                    version(complex_numbers) { eps = complex(0.0, eps0.re); }
                    else { eps = eps0*fabs(pcell.U[ftl][idx]) + eps0; }
                    pcell.U[ftl][idx] += eps;
                    pcell.decode_conserved(gtl, ftl, 0.0);

                    // update (ghost cell) boundary conditions
                    if (bc[bface.bc_id].preReconAction.length > 0) { bc[bface.bc_id].applyPreReconAction(0.0, 0, 0, bface); }
                    ghost_cell.encode_conserved(gtl, ftl, 0.0);

                    // fill local Jacobian
                    foreach(jdx; 0..nConserved) {
                        version(complex_numbers) { flowJacobian.dudU[jdx][idx] = ghost_cell.U[ftl][jdx].im/(eps.im); }
                        else { flowJacobian.dudU[jdx][idx] = (ghost_cell.U[ftl][jdx]-ghost_cell.U[0][jdx])/eps; }
                    }

                    // return cells to original state
                    pcell.U[ftl].copy_values_from(pcell.U[0]);
                    pcell.fs.copy_values_from(*fs_save);

                    // update (ghost cell) boundary conditions
                    if (bc[bface.bc_id].preReconAction.length > 0) { bc[bface.bc_id].applyPreReconAction(0.0, 0, 0, bface); }
                    ghost_cell.encode_conserved(gtl, ftl, 0.0);
                }

                // Step 2. Calculate dR/du
                // R: residual of interior cells conserved quantities
                // u: ghost cell conserved quantities

                // save a copy of the flowstate and copy conserved quantities
                fs_save.copy_values_from(*(ghost_cell.fs));
                ghost_cell.U[ftl].copy_values_from(ghost_cell.U[0]);

                foreach(idx; 0..nConserved) {

                    // peturb conserved quantity
                    version(complex_numbers) { eps = complex(0.0, eps0.re); }
                    else { eps = eps0*fabs(ghost_cell.U[ftl][idx]) + eps0; }
                    ghost_cell.U[ftl][idx] += eps;
                    ghost_cell.decode_conserved(gtl, ftl, 0.0);

                    // evaluate perturbed residuals in local stencil
                    evalRHS(gtl, ftl, ghost_cell.cell_list, ghost_cell.face_list, ghost_cell);

                    // fill local Jacobians
                    foreach (cell; ghost_cell.cell_list) {
                        foreach(jdx; 0..nConserved) {
                            version(complex_numbers) { cell.dRdU[jdx][idx] = cell.dUdt[ftl][jdx].im/eps.im; }
                            else { cell.dRdU[jdx][idx] = (cell.dUdt[ftl][jdx]-cell.dUdt[0][jdx])/eps; }
                        }
                    }

                    // return cell to original state
                    version(complex_numbers) {
                        ghost_cell.U[ftl].copy_values_from(ghost_cell.U[0]);
                        ghost_cell.fs.copy_values_from(*fs_save);
                        if (myConfig.viscous) {
                            foreach (cell; ghost_cell.cell_list) { cell.grad.copy_values_from(*(cell.grad_save)); }
                        }
                        if (flowJacobian.spatial_order >= 2) {
                            foreach(cell; ghost_cell.cell_list) { cell.gradients.copy_values_from(*(cell.gradients_save)); }
                        }
                    } else {
                        // TODO: for structured grids employing a real-valued low-order numerical Jacobian as
                        // the precondition matrix the above code doesn't completely clean up the effect of the
                        // perturbuation. Note that this does not occur for unstructured grids.
                        // This is currently under investigation, in the interim we will apply the more costly
                        // method of re-evaluating evalRHS with the unperturbed conserved state. KAD 2023-08-31
                        ghost_cell.U[ftl].copy_values_from(ghost_cell.U[0]);
                        ghost_cell.decode_conserved(gtl, 0, 0.0);
                        evalRHS(gtl, 0, ghost_cell.cell_list, ghost_cell.face_list, ghost_cell);
                    }
                }

                // Step 3. Calculate dR/dU and add corrections to Jacobian
                // dR/dU = dR/du * du/dU
                foreach(bcell; ghost_cell.cell_list) { //
                    assert(!bcell.is_ghost, "Oops, we expect to not find a ghost cell at this point.");
                    size_t I, J;
                    for ( size_t i = 0; i < nConserved; ++i ) {
                        I = bcell.id*nConserved + i; // column index
                        for ( size_t j = 0; j < nConserved; ++j ) {
                            J = pcell.id*nConserved + j; // row index
                            for (size_t k = 0; k < nConserved; k++) {
                                flowJacobian.local[I,J] = flowJacobian.local[I,J] + bcell.dRdU[i][k].re*flowJacobian.dudU[k][j];
                            }
                        }
                    }
                }
            } // foreach ( bi, bface; bndary.faces)
        } // foreach ( bndary; bc )
    } // end apply_jacobian_bcs()

    void evalRHS(int gtl, int ftl, ref FluidFVCell[] cell_list, FVInterface[] iface_list, FluidFVCell pcell)
    /*
     *  This method evaluates the RHS residual on a subset of cells for a given FluidBlock.
     *  It is used when constructing the numerical Jacobian.
     *  Its effect should replicate evalRHS() in steadystatecore.d for a subset of cells.
     */
    {
        // assume static grid for the preconditioner jacobian for the moment.
        // if (GlobalConfig.grid_motion == GridMotion.shock_fitting) {
        //     foreach (i, fba; fluidBlockArrays) {
        //         if (fba.shock_fitting) { compute_vtx_velocities_for_sf(fba, ftl); }
        //     }

        //     foreach (blk; localFluidBlocksBySize) {
        //         compute_avg_face_vel(blk, ftl);        
        //     }
        // }


        foreach(iface; iface_list) iface.F.clear();
        foreach(cell; cell_list) cell.clear_source_vector();

        bool do_reconstruction = ( flowJacobian.spatial_order > 1 );

        // convective flux update
        convective_flux_phase0_legacy(do_reconstruction, gtl, cell_list, iface_list);
        convective_flux_phase1(do_reconstruction, gtl, cell_list, iface_list);
        convective_flux_phase2(do_reconstruction, gtl, cell_list, iface_list);

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
                c.grad.gradients_leastsq(myConfig, c.cloud_fs, *(c.ws_grad)); // flow_property_spatial_derivatives(0);
            }

            // we need to average cell-centered spatial (/viscous) gradients to get approximations of the gradients
            // at the cell interfaces before the viscous flux calculation.
            if (myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                foreach(f; iface_list) {
                    f.average_cell_deriv_values(gtl);
                }
            }
            estimate_turbulence_viscosity(cell_list);
            average_turbulent_transprops_to_faces(iface_list);
            viscous_flux(iface_list);
            foreach(f; iface_list) {
                if (f.is_on_boundary) { applyPostDiffFluxAction(0.0, gtl, ftl, f); }
            }
        }

        // source terms and flux integration
        // the limit_factor is used to slowly increase the magnitude of the
        // thermochemical source terms from 0 to 1 for problematic reacting flows
        double limit_factor = 1.0;
        if (myConfig.nsteps_of_chemistry_ramp > 0) {
            double S = SimState.step/to!double(myConfig.nsteps_of_chemistry_ramp);
            limit_factor = min(1.0, S);
        }
        if (myConfig.conductivity_model) {
            foreach (cell; cell_list) cell.evaluate_electrical_conductivity();
        }
        foreach (i, cell; cell_list) {
            cell.add_inviscid_source_vector(gtl, 0.0);
            if (myConfig.viscous) {
                cell.add_viscous_source_vector();
            }
            if (myConfig.reacting) {
                // NOTE: we only need to evaluate the chemical source terms for the perturb cell
                //       this saves us a lot of unnecessary computations
                if (cell.id == pcell.id) {
                    cell.add_thermochemical_source_vector(thermochem_source, limit_factor);
                } else {
                    foreach (j; 0 .. myConfig.cqi.n) { cell.Q[j] += cell.Q_save[j]; }
                }
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
                getUDFSourceTermsForCell(myL, cell, 0, 0.0, myConfig, id, i_cell, j_cell, k_cell);
                cell.add_udf_source_vector();
            }
            cell.time_derivatives(gtl, ftl);
        }

    } // end evalRHS()

    void allocate_GMRES_workspace(int maxLinearSolverIterations, bool useRealValuedFrechetDerivative)
    {
        size_t nConserved = GlobalConfig.cqi.n;
        int n_species = GlobalConfig.gmodel_master.n_species();
        int n_modes = GlobalConfig.gmodel_master.n_modes();
        maxRate = new ConservedQuantities(nConserved);
        residuals = new ConservedQuantities(nConserved);

        size_t m = to!size_t(maxLinearSolverIterations);
        size_t n = nConserved*cells.length;

        // add the number of grid degrees of freedom to the total
        // number of variables
        if (myConfig.grid_motion != GridMotion.none) {
            foreach (ref vtx; vertices) {
                if (vtx.solve_position) { n += myConfig.dimensions; } 
            }

            maxVel.length = myConfig.dimensions;
        }

        nvars = n;
        // Now allocate arrays and matrices
        R.length = n;
        if (useRealValuedFrechetDerivative) { R0.length = n; }
        dU.length = n; dU[] = 0.0;
        r0.length = n;
        x0.length = n;
        DinvR.length = n;
        rhs.length = n;
        v.length = n;
        w.length = n;
        zed.length = n;
        g0.length = m+1;
        g1.length = m+1;
        VT.length = (m+1)*n;
        Q1 = new Matrix!double(m+1, m+1);
    }

    void allocate_FGMRES_workspace(int maxLinearSolverIterations,
                                   int maxPreconditioningIterations,
                                   bool useRealValuedFrechetDerivative)
    {
        size_t nConserved = GlobalConfig.cqi.n;
        int n_species = GlobalConfig.gmodel_master.n_species();
        int n_modes = GlobalConfig.gmodel_master.n_modes();
        maxRate = new ConservedQuantities(nConserved);
        residuals = new ConservedQuantities(nConserved);

        size_t m = to!size_t(maxLinearSolverIterations);
        size_t m_inner = to!size_t(maxPreconditioningIterations);
        size_t n = nConserved*cells.length;

        // add the number of grid degrees of freedom to the total
        // number of variables
        if (myConfig.grid_motion != GridMotion.none) {
            foreach (ref vtx; vertices) {
                if (vtx.solve_position) { n += myConfig.dimensions; } 
            }

            maxVel.length = myConfig.dimensions;
        }

        nvars = n;
        // Now allocate arrays and matrices
        R.length = n;
        if (useRealValuedFrechetDerivative) { R0.length = n; }
        dU.length = n; dU[] = 0.0;
        r0.length = n;
        x0.length = n;
        DinvR.length = n;
        rhs.length = n;
        v.length = n;
        w.length = n;
        zed.length = n;
        g0.length = m+1;
        g1.length = m+1;
        VT.length = (m+1)*n;
        Q1 = new Matrix!double(m+1, m+1);

        ZT.length = (m+1)*n;
        VT_inner.length = (m_inner+1)*n;
        inner_x0.length = n;
        inner_g1.length = m_inner+1;
    }

    } // end version(newton_krylov)

    void evalRU(double t, int gtl, int ftl, FluidFVCell c, bool do_reconstruction, double reaction_fraction)
    {
        // This method evaluates the R(U) for a single cell.
        // It is used when constructing the numerical Jacobian.
        // Adapted from Kyle's evalRHS().
        //
        // [TODO] PJ 2021-05-15 Might be good to move this code to the FluidFVCell class
        // but we will need to make the flux-calculation functions more cell centric.
        //
        foreach(f; c.iface) { f.F.clear(); }
        c.clear_source_vector();
        //
        FluidFVCell[1] cell_list = [c];
        convective_flux_phase0_legacy(do_reconstruction, gtl, cell_list, c.iface);
        convective_flux_phase1(do_reconstruction, gtl, cell_list, c.iface);
        convective_flux_phase2(do_reconstruction, gtl, cell_list, c.iface);

        foreach(f; c.iface) {
            if (f.is_on_boundary) { applyPostConvFluxAction(t, gtl, ftl, f); }
        }
        c.add_inviscid_source_vector(gtl, omegaz);
        //
        if (myConfig.viscous) {
            // Prepare spatial derivatives for viscous flux.
            // Note that with the cell-by-cell calculation process
            // that is associated with the point-implicit update,
            // the "natural" place for evaluating the derivatives will be at the faces.
            // If another location is selected, the values will still be needed at the faces
            // in order to compute viscous fluxes.
            // In such a case, nearby gradient values will be averaged to get an estimate
            // at the face and some of the values being fed into the average will be old.
            // For the initial step, they might even be zero.
            foreach(f; c.iface) {
                if (f.is_on_boundary) { applyPreSpatialDerivActionAtBndryFaces(t, gtl, ftl, f); }
            }
            final switch (myConfig.spatial_deriv_locn) {
            case SpatialDerivLocn.vertices:
                foreach (vtx; c.vtx) {
                    if (myConfig.dimensions == 2) {
                        if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                            vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, *(vtx.ws_grad));
                        } else {
                            vtx.grad.gradients_xy_div(myConfig, vtx.cloud_fs, vtx.cloud_pos);
                        }
                    } else {
                        // Have only least-squares in 3D.
                        vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, *(vtx.ws_grad));
                    }
                } // end foreach vtx
                // We've finished computing gradients at vertices, now copy them around if needed.
                // Interfaces need them for reconstruction/viscous fluxes.
                foreach (f; c.iface) { f.average_vertex_deriv_values(); }
                // Turbulence models will need cell centered gradients (also for axisymmetric!)
                if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent || myConfig.save_viscous_gradients)){
                    c.average_vertex_deriv_values();
                }
                break;
            case SpatialDerivLocn.faces:
                foreach(f; c.iface) {
                    if (myConfig.dimensions == 2) {
                        if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                            f.grad.gradients_leastsq(myConfig, f.cloud_fs, *(f.ws_grad));
                        } else {
                            f.grad.gradients_xy_div(myConfig, f.cloud_fs, f.cloud_pos);
                        }
                    } else {
                        // 3D has only least-squares available.
                        f.grad.gradients_leastsq(myConfig, f.cloud_fs, *(f.ws_grad));
                    } // end if (myConfig.dimensions)
                } // end foreach f
                // Finished computing gradients at faces, now copy them around if needed.
                // Turbulence models and axisymmetric source terms need cell centered gradients
                if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent || myConfig.save_viscous_gradients)){
                    c.average_interface_deriv_values();
                }
                break;
            case SpatialDerivLocn.cells:
                c.grad.gradients_leastsq(myConfig, c.cloud_fs, *(c.ws_grad));
                // We need to average cell-centered spatial (/viscous) gradients
                // to get approximations of the gradients at the faces for the viscous flux.
                foreach(f; c.iface) { f.average_cell_deriv_values(gtl); }
            } // end switch spatial_deriv_locn
            estimate_turbulence_viscosity(cell_list);
            average_turbulent_transprops_to_faces(c.iface);
            viscous_flux(c.iface);
            foreach(f; c.iface) {
                if (f.is_on_boundary) { applyPostDiffFluxAction(t, gtl, ftl, f); }
            }
            c.add_viscous_source_vector();
        } // end if viscous
        if (myConfig.reacting && myConfig.chemistry_update == ChemistryUpdateMode.integral) {
            c.add_thermochemical_source_vector(thermochem_source, reaction_fraction);
        }
        if (myConfig.udf_source_terms) { c.add_udf_source_vector(); }
        c.time_derivatives(gtl, ftl);
    } // end evalRU()

    @nogc
    void time_derivatives(int gtl, int ftl, size_t[] cell_idxs=[])
    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    // gtl : (grid-time-level) flow derivatives are evaluated at this grid level
    // ftl : (flow-time-level) specifies where computed derivatives are to be stored.
    //       0: Start of stage-1 update.
    //       1: End of stage-1.
    //       2: End of stage-2.
    {

        immutable size_t neq = myConfig.cqi.n;
        immutable size_t ngtl = myConfig.n_grid_time_levels;
        if (cell_idxs.length == 0) cell_idxs = celldata.all_cell_idxs;

        // Doesn't allocate new memory and seems to be very fast (NNG)
        number[] dUdt;
        if (ftl==0) { dUdt = celldata.dUdt0; }
        if (ftl==1) { dUdt = celldata.dUdt1; }
        if (ftl==2) { dUdt = celldata.dUdt2; }
        if (ftl==3) { dUdt = celldata.dUdt3; }
        if (ftl==4) { dUdt = celldata.dUdt4; }

        foreach(cidx; cell_idxs){
            // Note this is the number of faces that the cell cidx has, not the total number
            size_t nfaces = celldata.nfaces[cidx];
            number vol_inv = 1.0 / celldata.volumes[cidx*ngtl + gtl]; // Cell volume (inverted).

            foreach(j; 0 .. neq){
                number surface_integral = to!number(0.0);
                size_t idx = cidx*neq + j;

                foreach(i; 0 .. nfaces){
                    size_t fidx = celldata.c2f[cidx][i];
                    number area = celldata.outsigns[cidx][i]*facedata.areas[fidx*ngtl + gtl];
                    surface_integral -= facedata.fluxes[fidx*neq + j] * area;
                }
                dUdt[idx] = vol_inv*surface_integral + celldata.source_terms[idx];
            }
        }
        // GCL for steady-state
        if ((myConfig.solverMode == SolverMode.steady) && (myConfig.grid_motion != GridMotion.none)) {
            moving_grid_correction(gtl, ftl, cell_idxs);
        }

        if (myConfig.cqi.MHD && myConfig.MHD_static_field) {
            zero_magnetic_field_derivatives(gtl, ftl, cell_idxs);
        }
    } // end time_derivatives()

    @nogc
    void moving_grid_correction(int gtl, int ftl, size_t[] cell_idxs)
    /*
        For steady-state moving grid calculations, we need to correct the
        residual for the change in conserved quantities due to grid movement.
        This correction isn't needed for the transient moving grid sims, which
        handle this in a different, and probably better, way.

        @author: Nick Gibbons and Rob Watt (May 2025)
    */
    {

        immutable size_t neq = myConfig.cqi.n;
        immutable size_t ngtl = myConfig.n_grid_time_levels;

        number[] U; number[] dUdt;
        if (ftl==0) { dUdt = celldata.dUdt0; U = celldata.U0; }
        if (ftl==1) { dUdt = celldata.dUdt1; U = celldata.U1; }
        if (ftl==2) { dUdt = celldata.dUdt2; U = celldata.U2; }
        if (ftl==3) { dUdt = celldata.dUdt3; U = celldata.U3; }
        if (ftl==4) { dUdt = celldata.dUdt4; U = celldata.U4; }


        foreach(cidx; cell_idxs){
            // Note this is the number of faces that the cell cidx has, not the total number
            size_t nfaces = celldata.nfaces[cidx];
            number vol_inv = 1.0 / celldata.volumes[cidx*ngtl + gtl]; // Cell volume (inverted).

            foreach(j; 0 .. neq){
                number surface_integral = to!number(0.0);
                size_t idx = cidx*neq + j;
                foreach(i; 0 .. nfaces){
                    size_t fidx = celldata.c2f[cidx][i];
                    number area = celldata.outsigns[cidx][i]*facedata.areas[fidx*ngtl + gtl];
                    number vn = dot(facedata.normals[fidx], facedata.grid_velocities[fidx]);
                    surface_integral -= U[idx] * vn * area;
                }
                dUdt[idx] += vol_inv*surface_integral;
            }
        }
    }

    @nogc
    void zero_magnetic_field_derivatives(int gtl, int ftl, size_t[] cell_idxs)
    /*
        Implement MHD_static_field by zeroing the bits of the conservation
        equations that are involved in MHD>

        @author: Nick Gibbons and Peter Jacobs (May 2025)
    */
    {

        immutable size_t neq = myConfig.cqi.n;
        immutable size_t ngtl = myConfig.n_grid_time_levels;
        auto cqi = myConfig.cqi;

        number[] dUdt;
        if (ftl==0) { dUdt = celldata.dUdt0; }
        if (ftl==1) { dUdt = celldata.dUdt1; }
        if (ftl==2) { dUdt = celldata.dUdt2; }
        if (ftl==3) { dUdt = celldata.dUdt3; }
        if (ftl==4) { dUdt = celldata.dUdt4; }


        foreach(cidx; cell_idxs){
            size_t idx = cidx*neq;
            dUdt[idx + cqi.xB] = 0.0;
            dUdt[idx + cqi.yB] = 0.0;
            dUdt[idx + cqi.zB] = 0.0;
            dUdt[idx + cqi.psi] = 0.0;
            dUdt[idx + cqi.divB] = 0.0;
        }
    }

    @nogc
    void add_udf_source_vectors(size_t[] cell_list=[]){
    /*
        To save some computation, it is possible to compute the
        user defined source terms once per timestep, store them
        in Qudf (in the cells) and then apply them to the RHS
        for each stage of the explicit update here.

        However, for Kyle's temporal MMS and some other specific
        applications, we need to sometimes recalculate the terms'
        each stage. This neccessitated a separation of this
        routine from eval_udf_source_vectors, which are defined
        separately in sfluidblock and ufluidblock.

        @author: Nick Gibbons (May 2025)
    */
        if (myConfig.udf_source_terms) {
            if (cell_list.length==0) cell_list = celldata.all_cell_idxs;
            foreach (i; cell_list) {
                auto cell = cells[i];
                cell.add_udf_source_vector();
            }
        }
    }
} // end class FluidBlock
