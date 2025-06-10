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
import flowstate;
import fvvertex;
import fvinterface;
import fvcell;
import flowgradients;
import bc;
import fluxcalc;
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
import fluidblockio;
import fluidblockio_new;
import mass_diffusion;
import onedinterp;
version(mpi_parallel) {
    import mpi;
}

// version(diagnostics) {
// import plt = matplotlibd.pyplot;
// import std.format;
// }



// We no longer use a large number for the cell IDs to indicate ghost cells.
// Instead, consult the is_ghost_cell member. (NNG, Feb 23) 
//enum ghost_cell_start_id = 1_000_000_000;


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
    number mass_balance; // domain mass balance used to monitor for steady state
    number L2_residual; // L2 norm of the global residual
    number mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case
    ConservedQuantities Linf_residuals;
    number c_h, divB_damping_length; //divergence cleaning parameters for MHD
    int mncell;                 // number of monitor cells
    number[] initial_T_value; // for monitor cells to check against
    number[] jx,jy,jz,hs;
    //
    // Collections of cells, vertices and faces are held as arrays of references.
    // These allow us to conveniently work through the items via foreach statements.
    FVCell[] cells;
    FVInterface[] faces;
    FVVertex[] vertices;
    BoundaryCondition[] bc; // collection of references to the boundary conditions
    // Densified core datastructures
    FVCellData celldata;
    FVInterfaceData facedata;
    FVVertexData vertexdata;

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
    ConservedQuantities flux;
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
    // Shape sensitivity calculator workspace.
    version(shape_sensitivity) {
    FlowJacobianT flowJacobian;
    immutable size_t MAX_PERTURBED_INTERFACES = 800;
    FlowState fsSave;
    FlowGradients gradSave;
    FVCell cellSave;
    FVInterface[MAX_PERTURBED_INTERFACES] ifaceP;
    FlowState[MAX_PERTURBED_INTERFACES] fsP;
    FlowGradients[MAX_PERTURBED_INTERFACES] gradP;
    //
    size_t[] local_pcell_global_coord_list;
    size_t[][] local_ecell_global_coord_list;
    number[][] local_entry_list;
    //
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

    FlowState* fs_save;

    // storage for a precondition matrix
    FlowJacobian flowJacobian;

    // Work-space for Newton-Krylov accelerator
    // These arrays and matrices are directly tied to using the
    // GMRES iterative solver.
    SMatrix!double JcT; // transposed Jacobian (w.r.t conserved variables)
    ConservedQuantities maxRate, residuals;
    double normAcc, dotAcc;
    size_t nvars;
    Matrix!double Minv;
    double[] FU, dU, DinvR, r0, x0, rhs;
    double[] v, w, zed;
    double[] g0, g1;
    Matrix!double Q1;
    double[] VT;
    size_t Vstride;
    }

    FluidBlockIO[] block_io; // io handlers

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
        version(nk_accelerator) {
            fs_save = new FlowState(dedicatedConfig[id].gmodel, dedicatedConfig[id].turb_model.nturb);
        }
    }

    void add_IO()
    {
        if (!is_legacy_format(GlobalConfig.flow_format))
            block_io = get_fluid_block_io(this);
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
    @nogc abstract void update_nonshared_ghost_cell_positions(size_t gtl);
    @nogc abstract void precompute_stencil_data(size_t gtl);
    @nogc abstract void compute_least_squares_setup(size_t gtl);
    @nogc abstract void sync_vertices_from_underlying_grid(size_t gtl=0);
    @nogc abstract void sync_vertices_to_underlying_grid(size_t gtl=0);
    abstract void read_new_underlying_grid(string fileName);
    abstract void write_underlying_grid(string fileName);
    @nogc abstract void propagate_inflow_data_west_to_east();
    @nogc abstract void convective_flux_phase0new(bool allow_high_order_interpolation, size_t[] cell_idxs=[], size_t[] face_idxs=[]);
    @nogc abstract void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl,
                                               FVCell[] cell_list = [], FVInterface[] iface_list = [],
                                               FVVertex[] vertex_list = []);
    @nogc abstract void convective_flux_phase1new(bool allow_high_order_interpolation, size_t[] cell_idxs=[], size_t[] face_idxs=[]);
    @nogc abstract void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0,
                                               FVCell[] cell_list = [], FVInterface[] iface_list = [],
                                               FVVertex[] vertex_list = []);
    @nogc abstract void convective_flux_phase2new(bool allow_high_order_interpolation, size_t[] cell_idxs=[], size_t[] face_idxs=[]);
    @nogc abstract void convective_flux_phase2(bool allow_high_order_interpolation, size_t gtl=0,
                                               FVCell[] cell_list = [], FVInterface[] iface_list = [],
                                               FVVertex[] vertex_list = []);
    abstract size_t[] get_cell_write_indices();
    abstract void eval_udf_source_vectors(double simTime, size_t[] cell_idxs=[]);

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

        TODO: Some of this stuff doesn't need to be allocated in the ghost cells
        @author: Nick Gibbons
    */
        auto gmodel = myConfig.gmodel;
        size_t nturb = myConfig.turb_model.nturb;

        // This, apparently stupid, array is an experiment
        celldata.all_cell_idxs.length = ncells;
        celldata.halo_cell_ids.length = (ncells + nghost);
        celldata.halo_face_ids.length = (ncells + nghost);
        foreach(i; 0 .. ncells) celldata.all_cell_idxs[i] = i;

        celldata.c2f.length = ncells;
        celldata.outsigns.length=ncells;
        celldata.nfaces.length = ncells;
        celldata.dt_local.length = ncells;
        celldata.areas.length = ncells + nghost;
        celldata.wall_distances.length = ncells;
        celldata.data_is_bad.length = ncells;
        celldata.in_turbulent_zone.length = ncells;
        celldata.volumes.length = ncells + nghost;
        celldata.lengths.length = ncells + nghost;
        celldata.positions.length = ncells + nghost;
        celldata.U0.length = (ncells + nghost)*neq*nftl;
        celldata.cell_cloud_indices.length = ncells;
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

        version(nk_accelerator) {
            celldata.saved_gradients.reserve(ncells + nghost);
            foreach (n; 0 .. ncells+nghost) celldata.saved_gradients ~= FlowGradients(myConfig);
            celldata.saved_source_terms.length = (ncells + nghost)*neq;
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

        // This, apparently stupid, array is an experiment
        facedata.all_face_idxs.length = nfaces;
        foreach(i; 0 .. nfaces) facedata.all_face_idxs[i] = i;

        facedata.positions.length = nfaces;
        facedata.areas.length = nfaces;
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
        facedata.gradients.reserve(nfaces);
        facedata.workspaces.reserve(nfaces);
        foreach (n; 0 .. nfaces) facedata.flowstates ~= FlowState(gmodel, nturb);
        foreach (n; 0 .. nfaces) facedata.gradients ~= FlowGradients(myConfig);
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

        foreach(i, cell; cells){
            celldata.in_turbulent_zone[i] = cell.in_turbulent_zone;
        }
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

    @nogc void set_face_flowstates_to_averages_from_cells(size_t[] face_idxs=[])
    {
        if (face_idxs.length==0) face_idxs = facedata.all_face_idxs;

        foreach(idx; face_idxs){
            size_t l = facedata.f2c[idx].left;
            size_t r = facedata.f2c[idx].right;

            if (facedata.left_interior_only[idx]) {
                facedata.flowstates[idx] = celldata.flowstates[l];

            } else if (facedata.right_interior_only[idx]) {
                facedata.flowstates[idx] = celldata.flowstates[r];

            } else {
                // Assume we have both cells, for a shared or internal boundary interface
                facedata.flowstates[idx].copy_average_values_from(celldata.flowstates[l], celldata.flowstates[r]);
            }
        } // end main face loop
    }

    // Increment the local DFT in each cell.
    // Do we need garbage collection here? I don't think so?
    @nogc
    void increment_DFT(size_t DFT_step) {
        foreach(cell; cells) {
            cell.increment_local_DFT(DFT_step);
        }
    }

    @nogc
    void estimate_turbulence_viscosity(size_t[] cell_idxs = [])
    {
    version(turbulence) { // Exit instantly if turbulence capability disabled
        auto gmodel = myConfig.gmodel;
        if (cell_idxs.length == 0) cell_idxs = celldata.all_cell_idxs;

        foreach (i; cell_idxs) {
            if ( celldata.in_turbulent_zone[i] ) {
                celldata.flowstates[i].mu_t = myConfig.turb_model.turbulent_viscosity(celldata.flowstates[i],
                                                                                      celldata.gradients[i],
                                                                                      celldata.positions[i].y,
                                                                                      celldata.wall_distances[i]);
                celldata.flowstates[i].k_t = myConfig.turb_model.turbulent_conductivity(celldata.flowstates[i], gmodel);
            } else {
                /* Presume this part of the flow is laminar; clear turbulence quantities. */
                celldata.flowstates[i].mu_t= 0.0;
                celldata.flowstates[i].k_t = 0.0;
            }
        }
    }
    } // end estimate_turbulence_viscosity()

    @nogc
    void estimate_turbulence_viscosity(FVCell[] cell_list)
    {
        version(turbulence) { // Exit instantly if turbulence capability disabled
        if (cell_list.length == 0) { cell_list = cells; }
        foreach (cell; cell_list) {
            cell.turbulence_viscosity();
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
        case ShockDetector.NNG:
            auto gm = myConfig.gmodel;
            if (!facedata.stencil_idxs) throw new Error("NNG Shock detector needs a 4 cell stencil");
            foreach (idx; 0 .. nfaces) {
                size_t l = facedata.f2c[idx].left;
                size_t r = facedata.f2c[idx].right;
                number S0 = NNG_ShockDetector(gm, celldata.flowstates[l],  celldata.flowstates[r], facedata.normals[idx], comp_tol.re);

                size_t LL = facedata.stencil_idxs[idx].L1;
                size_t  L = facedata.stencil_idxs[idx].L0;
                size_t  R = facedata.stencil_idxs[idx].R0;
                size_t RR = facedata.stencil_idxs[idx].R1;

                number S1 = NNG_DiscoDectector1(gm, celldata.flowstates[LL], celldata.flowstates[L], celldata.flowstates[R], celldata.flowstates[RR]);
                number S2 = NNG_DiscoDectector2(gm, celldata.flowstates[LL], celldata.flowstates[L], celldata.flowstates[R], celldata.flowstates[RR]);
                facedata.flowstates[idx].S = fmax(fmax(S0, S1), S2);
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

    int count_invalid_cells_safe(int gtl, int ftl)
    {
        try {
            return count_invalid_cells(gtl, ftl);
        } catch (FlowSolverException e) {
            writefln("Adjust invalid cells failed: %s", e.msg);
            return to!int(ncells);
        }
    }

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
            if (cell.data_is_bad || celldata.data_is_bad[cell.id] || cell.fs.check_data(cell.pos[0], myConfig) == false) {
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
                    celldata.data_is_bad[cell.id] = false;
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
                        vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, vtx.cloud_pos, *(vtx.ws_grad));
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
                    vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, vtx.cloud_pos, *(vtx.ws_grad));
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
                        iface.grad.gradients_leastsq(myConfig, iface.cloud_fs, iface.cloud_pos, *(iface.ws_grad));
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
                        iface.grad.gradients_leastsq(myConfig, iface.cloud_fs, iface.cloud_pos, *(iface.ws_grad));
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
                throw new Error("divergence thereom not implemented when computing gradients at cell centres.");
            } // end switch
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
        facedata.fluxes.clear();
    }

    @nogc
    void clear_cell_source_vectors()
    {
        celldata.source_terms.clear();
    }

    @nogc
    void average_turbulent_transprops_to_faces(FVInterface[] face_list)
    {
        if (myConfig.turb_model.isTurbulent){
            foreach (iface; face_list) { iface.average_turbulent_transprops(); }
        }
    }

    @nogc
    void average_turbulent_transprops_to_faces(size_t[] face_idxs=[])
    {
        if (!myConfig.turb_model.isTurbulent) return;
        if (face_idxs.length==0) face_idxs = facedata.all_face_idxs;

        foreach(idx; face_idxs){
            size_t l = facedata.f2c[idx].left;
            size_t r = facedata.f2c[idx].right;

            if (facedata.left_interior_only[idx]) {
                facedata.flowstates[idx].mu_t = celldata.flowstates[l].mu_t;
                facedata.flowstates[idx].k_t =  celldata.flowstates[l].k_t;

            } else if (facedata.right_interior_only[idx]) {
                facedata.flowstates[idx].mu_t = celldata.flowstates[r].mu_t;
                facedata.flowstates[idx].k_t =  celldata.flowstates[r].k_t;

            } else {
                // Assume we have both cells, for a shared or internal boundary interface
                facedata.flowstates[idx].mu_t = 0.5*(celldata.flowstates[l].mu_t + celldata.flowstates[r].mu_t);
                facedata.flowstates[idx].k_t =  0.5*(celldata.flowstates[l].k_t  + celldata.flowstates[r].k_t);
            }
        } // end main face loop
    }

    @nogc
    void viscous_flux(FVInterface[] face_list)
    {
        if (face_list.length == 0) { face_list = faces; }
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
        foreach(FVCell cell; cells) {
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
                    number massflux=0.0;
                    if (cqi.mass==0) {
                        massflux = face.F[cqi.mass];
                    } else {
                        foreach(isp; 0 .. cqi.n_species) massflux += face.F[cqi.species+isp];
                    }
                    mass_balance += boundary.outsigns[i] * massflux * face.area[0];
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
        foreach(i, cell; cells) {
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
                    // Set the allowable time-step based on hyperbolic time-step.
                    dt_allow = fmin(dt_allow_hyp, GlobalConfig.dt_max);
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
            celldata.dt_local[i] = cell.dt_local;
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

    version(nk_accelerator) {

    void initialize_jacobian(int spatial_order_of_jacobian, double sigma)
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
            cell.gather_residual_stencil_lists(spatial_order_of_jacobian);
            nentry += cell.cell_list.length;
        }

        // we will gather the ghost cell residual stencil lists
        // at this point as well for convenience, we need them
        // when we apply the boundary condition corrections later
        foreach ( bndary; bc ) {
            if (!bndary.ghost_cell_data_available) { continue; }
            if (bndary.type == "exchange_using_mapped_cells" || bndary.type == "exchange_over_full_face") { continue; }
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

        jacobian_nonzero_pattern(spatial_order_of_jacobian, nentry, cells.length, sigma);

        // Setup a helper array for use with the smla.solve later
        flowJacobian.local_diags.length = flowJacobian.local.ia.length-1;
        get_smatrix_diag_indices(flowJacobian.local, flowJacobian.local_diags);
    } // end initialize_jacobian()

    void jacobian_nonzero_pattern(int spatial_order_of_jacobian, size_t nentry, size_t ncells, double sigma) {
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
                assert(!cell.is_ghost_cell, "Oops, we expect to not find a ghost cell at this point.");
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
        if (myConfig.sssOptions.iluFill > 0) {

            // [TODO] think about making the lev matrix sparse as well KAD 2022-03-31
            // construct a level matrix
            int p = GlobalConfig.sssOptions.iluFill;
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
        // remove the conserved mass variable for multi-species gas
        //if (cqi.n_species > 1) { nConserved -= 1; }

        flowJacobian = new FlowJacobian(sigma, myConfig.dimensions, nConserved, spatial_order_of_jacobian, ptJac.local.aa.length, ncells);
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
        myConfig.flux_calculator = myConfig.sssOptions.preconditionMatrixFluxCalculator;

        // temporarily change interpolation order
        shared int interpolation_order_save = GlobalConfig.interpolation_order;
        myConfig.interpolation_order = to!int(flowJacobian.spatial_order);

        // copy some data for later use
        if (myConfig.viscous) {
            foreach(id; 0 .. ncells) celldata.saved_gradients[id].copy_values_from(celldata.gradients[id]);
        }
        if (myConfig.reacting) {
            clear_cell_source_vectors();
            eval_thermochem_source_vector(0);
            size_t cncq = ncells*myConfig.cqi.n;
            foreach(i; 0 .. cncq) celldata.saved_source_terms[i] = celldata.source_terms[i];
        }

        // the real-valued finite difference needs a base residual (R0)
        version(complex_numbers) { } // do nothing
        else { // TODO: Shouldn't this be done upstairs using the real evalRHS?
            foreach(i; 0 .. ncells) { evalRHS2(0, 0, celldata.halo_cell_ids[i], celldata.halo_face_ids[i], i); }
        }

        // fill out the rows of the Jacobian for a cell
        foreach(i; celldata.all_cell_idxs) { evaluate_cell_contribution_to_jacobian(i); }

        // add boundary condition corrections to boundary cells
        apply_jacobian_bcs();

        // return flux calculator to its original state
        myConfig.flux_calculator = flux_calculator_save;

        // return the interpolation order to its original state
        myConfig.interpolation_order = interpolation_order_save;
    } // end evaluate_jacobian()

    void evaluate_cell_contribution_to_jacobian(size_t pid)
    {
        auto cqi = myConfig.cqi; // was GlobalConfig.cqi;
        auto nConserved = cqi.n;
        size_t idx = nConserved*pid;

        auto eps0 = flowJacobian.eps;
        number eps;
        int ftl = 1; int gtl = 0;

        // save a copy of the flowstate and copy conserved quantities
        fs_save.copy_values_from(celldata.flowstates[pid]);
        foreach(j; 0..nConserved) celldata.U1[idx + j] = celldata.U0[idx + j];

        // perturb the current cell's conserved quantities
        // and then evaluate the residuals for each cell in
        // the local domain of influence
        foreach(j; 0 .. nConserved) {

            // peturb conserved quantity
            version(complex_numbers) { eps = complex(0.0, eps0.re); }
            else { eps = eps0*fabs(celldata.U1[idx + j]) + eps0; }
            celldata.U1[idx + j] += eps;
            decode_conserved(celldata.positions[pid], celldata.U1[idx .. idx+nConserved], celldata.flowstates[pid], omegaz, pid, myConfig);

            // evaluate perturbed residuals in local stencil
            evalRHS2(gtl, ftl, celldata.halo_cell_ids[pid], celldata.halo_face_ids[pid], pid);

            // fill local Jacobians
            foreach (cid; celldata.halo_cell_ids[pid]) {
                size_t cidx = cid*nConserved;
                foreach(i; 0 .. nConserved) {

                    double cell_dRdU_ij;
                    version(complex_numbers) { cell_dRdU_ij = celldata.dUdt1[cidx + i].im/eps.im; }
                    else { cell_dRdU_ij = (celldata.dUdt1[cidx + i] - celldata.dUdt0[cidx + i])/eps; }

                    size_t I = cidx + i;
                    size_t J = idx + j;
                    flowJacobian.local[I,J] = cell_dRdU_ij;
                }
            }

            // return cell to original state
            foreach(k; 0..nConserved) celldata.U1[idx + k] = celldata.U0[idx + k]; // Needed?
            celldata.flowstates[pid].copy_values_from(*fs_save);
            if (myConfig.viscous) {
                foreach (cid; celldata.halo_cell_ids[pid]) {
                    celldata.gradients[cid].copy_values_from(celldata.saved_gradients[cid]);
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
            if ( bndary.type == "exchange_using_mapped_cells" || bndary.type == "exchange_over_full_face") { continue; }
            foreach ( bi, bface; bndary.faces) {
                FVCell ghost_cell; FVCell pcell;
                size_t fid = bface.id;
                size_t ghost_cell_id; size_t pcell_id;
                if (bndary.outsigns[bi] == 1) {
                    pcell_id = facedata.f2c[fid].left;
                    ghost_cell_id = facedata.f2c[fid].right;
                } else {
                    pcell_id = facedata.f2c[fid].right;
                    ghost_cell_id = facedata.f2c[fid].left;
                }
                size_t pidx = nConserved*pcell_id;
                size_t gidx = nConserved*ghost_cell_id;

                // Step 1. Calculate du/dU
                // u: ghost cell conserved quantities
                // U: interior cell conserved quantities

                // save a copy of the flowstate and copy conserved quantities
                fs_save.copy_values_from(celldata.flowstates[pcell_id]);
                foreach(j; 0..nConserved) celldata.U1[pidx + j] = celldata.U0[pidx + j];
                encode_conserved(celldata.positions[ghost_cell_id], celldata.flowstates[ghost_cell_id], omegaz, myConfig, celldata.U0[gidx .. gidx+nConserved]);

                foreach(idx; 0..nConserved) {

                    // peturb conserved quantity
                    version(complex_numbers) { eps = complex(0.0, eps0.re); }
                    else { eps = eps0*fabs(celldata.U1[pidx + idx]) + eps0; }
                    celldata.U1[pidx + idx] += eps;
                    decode_conserved(celldata.positions[pcell_id], celldata.U1[pidx .. pidx+nConserved], celldata.flowstates[pcell_id], omegaz, pcell_id, myConfig);

                    // update (ghost cell) boundary conditions
                    if (bndary.preReconAction.length > 0) { bndary.applyPreReconAction(0.0, 0, 0, bface); }
                    encode_conserved(celldata.positions[ghost_cell_id], celldata.flowstates[ghost_cell_id], omegaz, myConfig, celldata.U1[gidx .. gidx+nConserved]);

                    // fill local Jacobian
                    foreach(jdx; 0..nConserved) {
                        version(complex_numbers) { flowJacobian.dudU[jdx][idx] = celldata.U1[gidx + jdx].im/eps.im; }
                        else { flowJacobian.dudU[jdx][idx] = (celldata.U1[gidx + jdx]-celldata.U0[gidx + jdx])/eps; }
                    }

                    // return cells to original state
                    foreach(k; 0..nConserved) celldata.U1[pidx + k] = celldata.U0[pidx + k]; // Needed?
                    celldata.flowstates[pcell_id].copy_values_from(*fs_save);

                    // update (ghost cell) boundary conditions
                    if (bndary.preReconAction.length > 0) { bndary.applyPreReconAction(0.0, 0, 0, bface); }
                    encode_conserved(celldata.positions[ghost_cell_id], celldata.flowstates[ghost_cell_id], omegaz, myConfig, celldata.U1[gidx .. gidx+nConserved]);
                }

                // Step 2. Calculate dR/du
                // R: residual of interior cells conserved quantities
                // u: ghost cell conserved quantities

                // save a copy of the flowstate and copy conserved quantities
                fs_save.copy_values_from(celldata.flowstates[ghost_cell_id]);
                foreach(j; 0..nConserved) celldata.U1[gidx + j] = celldata.U0[gidx + j];

                foreach(idx; 0..nConserved) {

                    // peturb conserved quantity
                    version(complex_numbers) { eps = complex(0.0, eps0.re); }
                    else { eps = eps0*fabs(celldata.U1[gidx + idx]) + eps0; }
                    celldata.U1[gidx + idx] += eps;
                    decode_conserved(celldata.positions[ghost_cell_id], celldata.U1[gidx .. gidx+nConserved], celldata.flowstates[ghost_cell_id], omegaz, pidx, myConfig);

                    // evaluate perturbed residuals in local stencil
                    evalRHS2(gtl, ftl, celldata.halo_cell_ids[ghost_cell_id], celldata.halo_face_ids[ghost_cell_id], ghost_cell_id);

                    // fill local Jacobians
                    foreach (cid; celldata.halo_cell_ids[ghost_cell_id]) {
                        foreach(jdx; 0..nConserved) {
                            double dRdu_ji;
                            version(complex_numbers) { dRdu_ji = celldata.dUdt1[cid*nConserved + jdx].im/eps.im; }
                            else { dRdu_ji = (celldata.dUdt1[cid*nConserved + jdx]-celldata.dUdt0[cid*nConserved + jdx])/eps; }

                            // Step 3. Calculate dR/dU and add corrections to Jacobian
                            // dR/dU = dR/du * du/dU
                            size_t J = cid*nConserved + jdx; // column index
                            foreach(kdx; 0 .. nConserved) {
                                size_t K = pcell_id*nConserved + kdx; // row index
                                flowJacobian.local[J,K] = flowJacobian.local[J,K] + dRdu_ji*flowJacobian.dudU[idx][kdx];
                            }
                        }
                    }

                    // return cell to original state
                    foreach(j; 0..nConserved) celldata.U1[gidx + j] = celldata.U0[gidx + j];
                    celldata.flowstates[ghost_cell_id].copy_values_from(*fs_save);
                    if (myConfig.viscous) {
                        foreach (cid; celldata.halo_cell_ids[ghost_cell_id]) {
                            celldata.gradients[cid].copy_values_from(celldata.saved_gradients[cid]);
                        }
                    }
                }
            } // foreach ( bi, bface; bndary.faces)
        } // foreach ( bndary; bc )
    } // end apply_jacobian_bcs()

    void evalRHS2(int gtl, int ftl, size_t[] halo_cell_ids, size_t[] halo_face_ids, size_t pid)
    /*
     *  This method evaluates the RHS residual on a subset of cells for a given FluidBlock.
     *  It is used when constructing the numerical Jacobian.
     *  Its effect should replicate evalRHS() in steadystatecore.d for a subset of cells.
     */
    {
        auto cqi = myConfig.cqi;
        immutable size_t ncq = myConfig.cqi.n;

        foreach(id; halo_face_ids) { foreach(i; 0 .. ncq) facedata.fluxes[id*ncq + i] = 0.0; }
        foreach(id; halo_cell_ids) { foreach(i; 0 .. ncq) celldata.source_terms[id*ncq + i] = 0.0; }

        bool do_reconstruction = ( flowJacobian.spatial_order > 1 );

        // convective flux update
        convective_flux_phase0new(do_reconstruction, halo_cell_ids, halo_face_ids);
        convective_flux_phase1new(do_reconstruction, halo_cell_ids, halo_face_ids);
        convective_flux_phase2new(do_reconstruction, halo_cell_ids, halo_face_ids);

        foreach(id; halo_face_ids) {
            if (facedata.left_interior_only[id] || facedata.right_interior_only[id]) {
                applyPostConvFluxAction(0.0, gtl, ftl, faces[id]);
            }
        }

        // Viscous flux update
        if (myConfig.viscous) {

            foreach(id; halo_face_ids) {
                if (facedata.left_interior_only[id] || facedata.right_interior_only[id]) {
                    applyPreSpatialDerivActionAtBndryFaces(0.0, gtl, ftl, faces[id]); 
                }
            }

            // currently only for least-squares at faces
            compute_gradients_at_cells_leastsq(halo_cell_ids);
            average_lsq_cell_derivs_to_faces(halo_face_ids);
            estimate_turbulence_viscosity(halo_cell_ids);
            average_turbulent_transprops_to_faces(halo_face_ids);
            viscous_flux(halo_face_ids);

            foreach(id; halo_face_ids) {
                if (facedata.left_interior_only[id] || facedata.right_interior_only[id]) {
                    applyPostDiffFluxAction(0.0, gtl, ftl, faces[id]);
                }
            }
        }

        // Compute source terms. Note we only compute the thermochem source for the pcell,
        // which should always be the first entry in halo_cell_ids. This is because
        // the thermochem sources have no spatial dependance, so they do not change
        // when pcell is perturbed.
        size_t[1] this_cell_id = [pid];

        eval_fluid_source_vectors(omegaz, halo_cell_ids);
        eval_thermochem_source_vector(SimState.step, this_cell_id[0 .. 1]);
        eval_udf_source_vectors(0.0, halo_cell_ids);

        // Compute time derivatives
        time_derivatives(gtl, ftl, halo_cell_ids);

    } // end evalRHS()


    void evalRHS(int gtl, int ftl, ref FVCell[] cell_list, FVInterface[] iface_list, FVCell pcell)
    /*
     *  This method evaluates the RHS residual on a subset of cells for a given FluidBlock.
     *  It is used when constructing the numerical Jacobian.
     *  Its effect should replicate evalRHS() in steadystatecore.d for a subset of cells.
     */
    {
        foreach(iface; iface_list) iface.F.clear();
        foreach(cell; cell_list) cell.clear_source_vector();

        bool do_reconstruction = ( flowJacobian.spatial_order > 1 );

        // convective flux update
        convective_flux_phase0(do_reconstruction, gtl, cell_list, iface_list);
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
                c.grad.gradients_leastsq(myConfig, c.cloud_fs, c.cloud_pos, *(c.ws_grad)); // flow_property_spatial_derivatives(0);
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
        foreach (i, cell; cell_list) {
            cell.add_inviscid_source_vector(gtl, omegaz);
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


    // The following two methods are used to verify the numerical Jacobian implementation.
    void verify_jacobian(double sigma)
    {
        // we perform a residual evaluation to ensure the ghost cells are filled with good data
        import steadystate_core;
        steadystate_core.evalRHS(0.0, 0);

        // calculate the numerical Jacobian
        initialize_jacobian(GlobalConfig.interpolation_order, sigma);
        evaluate_jacobian();
        assert(flowJacobian !is null, "Oops, we expect a flowJacobian object to be attached to the fluidblock.");
        size_t nConserved = GlobalConfig.cqi.n;
        // remove the conserved mass variable for multi-species gas
        //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }

        // create an arbitrary unit vector
        double[] vec;
        vec.length = cells.length*nConserved;
        foreach ( i, ref val; vec) { val = i+1; }

        // normalise the vector
        double norm = 0.0;
        foreach( i; 0..vec.length) { norm += vec[i]*vec[i]; }
        norm = sqrt(norm);
        foreach( ref val; vec) { val = val/norm; }

        // result vectors
        double[] sol1;
        sol1.length = vec.length;
        double[] sol2;
        sol2.length = vec.length;

        // explicit multiplication of J*vec
        nm.smla.multiply(flowJacobian.local, vec, sol1);

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

        // stop the program at this point
        import core.stdc.stdlib : exit;
        exit(0);
    } // end verify_jacobian

    void evalConservativeJacobianVecProd(double[] vec, ref double[] sol) {
        size_t nConserved = GlobalConfig.cqi.n;
        auto EPS = GlobalConfig.sssOptions.sigma0;

        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach(i; 0 .. nvars) celldata.U1[i] = celldata.U0[i];

        foreach(i; 0 .. nvars) {
            version(complex_numbers) {
                celldata.U1[i] += complex(0.0, EPS.re*vec[i].re);
            } else {
                celldata.U1[i] += EPS*vec[i];
            }
        }
        foreach(id; 0 .. ncells) {
            size_t idx = id*nConserved;
            decode_conserved(celldata.positions[id], celldata.U1[idx .. idx+nConserved], celldata.flowstates[id], omegaz, id, myConfig);
        }

        import steadystate_core;
        clear_fluxes_of_conserved_quantities();
        clear_cell_source_vectors();
        steadystate_core.evalRHS(0.0, 1);

        foreach(i; 0 .. nvars) {
            version(complex_numbers) {
                sol[i] = celldata.dUdt1[i].im/EPS;
            } else {
                sol[i] = (celldata.dUdt1[i]-celldata.dUdt0[i])/EPS;
            }
        }
    }
    } // end version(nk_accelerator)

    version(nk_accelerator) {
    void allocate_GMRES_workspace()
    {
        size_t nConserved = GlobalConfig.cqi.n;
        // remove the conserved mass variable for multi-species gas
        //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
        int n_species = GlobalConfig.gmodel_master.n_species();
        int n_modes = GlobalConfig.gmodel_master.n_modes();
        maxRate = new_ConservedQuantities(nConserved);
        residuals = new_ConservedQuantities(nConserved);

        size_t mOuter = to!size_t(GlobalConfig.sssOptions.maxOuterIterations);
        size_t mInner = to!size_t(GlobalConfig.sssOptions.nInnerIterations);
        size_t n = nConserved*cells.length;
        nvars = n;
        // Now allocate arrays and matrices
        FU.length = n;
        dU.length = n; dU[] = 0.0;
        r0.length = n;
        x0.length = n;
        DinvR.length = n;
        rhs.length = n;
        v.length = n;
        w.length = n;
        zed.length = n;
        g0.length = mOuter+1;
        g1.length = mOuter+1;
        //h_outer.length = mOuter+1;
        //hR_outer.length = mOuter+1;
        VT.length = (mOuter+1)*n;
        //H0_outer = new Matrix!number(mOuter+1, mOuter);
        //H1_outer = new Matrix!number(mOuter+1, mOuter);
        //Gamma_outer = new Matrix!number(mOuter+1, mOuter+1);
        //Q0_outer = new Matrix!number(mOuter+1, mOuter+1);
        Q1 = new Matrix!double(mOuter+1, mOuter+1);
    }
    }

    void update_aux(double dt, double time, size_t step)
    {
        foreach (cell ; cells) {
            foreach(aux ; cell.aux_cell_data) {
                aux.update(cell, dt, time, step);
            }
        }
    }


    void evalRU(double t, int gtl, int ftl, FVCell c, bool do_reconstruction, double reaction_fraction)
    {
        // This method evaluates the R(U) for a single cell.
        // It is used when constructing the numerical Jacobian.
        // Adapted from Kyle's evalRHS().
        //
        // [TODO] PJ 2021-05-15 Might be good to move this code to the FVCell class
        // but we will need to make the flux-calculation functions more cell centric.
        //
        foreach(f; c.iface) { f.F.clear(); }
        c.clear_source_vector();
        //
        FVCell[1] cell_list = [c];
        convective_flux_phase0(do_reconstruction, gtl, cell_list, c.iface);
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
                            vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, vtx.cloud_pos, *(vtx.ws_grad));
                        } else {
                            vtx.grad.gradients_xy_div(myConfig, vtx.cloud_fs, vtx.cloud_pos);
                        }
                    } else {
                        // Have only least-squares in 3D.
                        vtx.grad.gradients_leastsq(myConfig, vtx.cloud_fs, vtx.cloud_pos, *(vtx.ws_grad));
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
                            f.grad.gradients_leastsq(myConfig, f.cloud_fs, f.cloud_pos, *(f.ws_grad));
                        } else {
                            f.grad.gradients_xy_div(myConfig, f.cloud_fs, f.cloud_pos);
                        }
                    } else {
                        // 3D has only least-squares available.
                        f.grad.gradients_leastsq(myConfig, f.cloud_fs, f.cloud_pos, *(f.ws_grad));
                    } // end if (myConfig.dimensions)
                } // end foreach f
                // Finished computing gradients at faces, now copy them around if needed.
                // Turbulence models and axisymmetric source terms need cell centered gradients
                if (myConfig.axisymmetric || (myConfig.turb_model.isTurbulent || myConfig.save_viscous_gradients)){
                    c.average_interface_deriv_values();
                }
                break;
            case SpatialDerivLocn.cells:
                c.grad.gradients_leastsq(myConfig, c.cloud_fs, c.cloud_pos, *(c.ws_grad));
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
        if (cell_idxs.length == 0) cell_idxs = celldata.all_cell_idxs;

        foreach(cidx; cell_idxs){
            // Note this is the number of faces that the cell cidx has, not the total number
            size_t nfaces = celldata.nfaces[cidx];
            number vol_inv = 1.0 / celldata.volumes[cidx]; // Cell volume (inverted).

            foreach(j; 0 .. neq){
                number surface_integral = to!number(0.0);
                foreach(i; 0 .. nfaces){
                    size_t fidx = celldata.c2f[cidx][i];
                    number area = celldata.outsigns[cidx][i]*facedata.areas[fidx];
                    surface_integral -= facedata.fluxes[fidx*neq + j] * area;
                }

                size_t idx = cidx*neq + j;
                if (ftl==0) celldata.dUdt0[idx] = vol_inv*surface_integral + celldata.source_terms[idx];
                if (ftl==1) celldata.dUdt1[idx] = vol_inv*surface_integral + celldata.source_terms[idx];
                if (ftl==2) celldata.dUdt2[idx] = vol_inv*surface_integral + celldata.source_terms[idx];
                if (ftl==3) celldata.dUdt3[idx] = vol_inv*surface_integral + celldata.source_terms[idx];
            }
        }
    } // end time_derivatives()

    @nogc void eval_fluid_source_vectors(double omegaz=0.0, size_t[] cell_idxs=[])
    {
        auto cqi = myConfig.cqi;
        size_t ncq = myConfig.cqi.n;
        immutable bool rotating_frame = omegaz!=0.0;
        immutable bool axisymmetric = myConfig.axisymmetric;
        immutable bool gravity_non_zero = myConfig.gravity_non_zero;
        immutable bool axiviscous = myConfig.viscous && myConfig.axisymmetric;
        immutable bool turbulent = myConfig.viscous && myConfig.turb_model.isTurbulent;
        immutable size_t cqiyMom = cqi.yMom;
        immutable Vector3 gravity = myConfig.gravity;
        immutable size_t cqirhoturb = cqi.rhoturb;
        immutable size_t n_turb = cqi.n_turb;

        if (cell_idxs.length==0) cell_idxs = celldata.all_cell_idxs;
        foreach(i; cell_idxs) {
            size_t idx = i*ncq;
            number[] Q = celldata.source_terms[idx .. idx+ncq]; // view or reference to source_terms data structure

            if (rotating_frame)
                add_rotating_frame_source_vector(myConfig, celldata.positions[i], celldata.flowstates[i], celldata.areas[i], celldata.volumes[i], omegaz, Q);

            if (axisymmetric)
                add_axisymmetric_source_vector(cqiyMom, celldata.flowstates[i], celldata.areas[i], celldata.volumes[i], Q);

            if (gravity_non_zero)
                add_gravitational_source_vector(gravity, *cqi, celldata.flowstates[i], Q);

            if (axiviscous)
                add_axisymmetric_viscous_source_vector(celldata.positions[i].y, cqiyMom, celldata.areas[i], celldata.volumes[i],
                                                   celldata.flowstates[i], celldata.gradients[i], Q);

            if (turbulent){
                if (celldata.in_turbulent_zone[i]) {
                    number[] rhoturb = celldata.source_terms[idx+cqirhoturb .. idx+cqirhoturb+n_turb];
                    myConfig.turb_model.source_terms(celldata.flowstates[i], celldata.gradients[i], celldata.positions[i].y,
                                                     celldata.wall_distances[i], celldata.lengths[i], rhoturb);
                }
            }
        }
        return;
    }

    @nogc void eval_thermochem_source_vector(int step, size_t[] cell_idxs=[])
    {
        if (!myConfig.reacting) return;
        auto cqi = myConfig.cqi;
        size_t ncq = myConfig.cqi.n;

        double limit_factor = 1.0;
        // the limit_factor is used to slowly increase the magnitude of the
        // thermochemical source terms from 0 to 1 for problematic reacting flows
        if (myConfig.nsteps_of_chemistry_ramp > 0) {
            double S = step/to!double(myConfig.nsteps_of_chemistry_ramp);
            limit_factor = min(1.0, S);
        }

        if (cell_idxs.length==0) cell_idxs = celldata.all_cell_idxs;
        foreach(i; cell_idxs) {
            size_t idx = i*ncq;
            number[] Q = celldata.source_terms[idx .. idx+ncq]; // view or reference to source_terms data structure
            add_thermochemical_source_vector(myConfig, thermochem_source, limit_factor, celldata.flowstates[i], Q);
        }
        return;
    }

    @nogc void average_lsq_cell_derivs_to_faces(size_t[] face_idxs=[]){
    /*
        Fast averageing of the cell-based derivatives, using the expression from
        A. Haselbacher, J. Blazek, Accurate and efficient discretization
        of Navier-Stokes equations on mixed grids, AIAA Journal 38
        (2000) 2094–2102. doi:10.2514/2.871.

        NNG, April 2023
    */
        size_t nsp = myConfig.n_species;
        size_t nmodes = myConfig.n_modes;
        size_t nturb = myConfig.turb_model.nturb;
        bool is3D = (myConfig.dimensions == 3);
        if (face_idxs.length==0) face_idxs = facedata.all_face_idxs;

        foreach(idx; face_idxs){
            size_t l = facedata.f2c[idx].left;
            size_t r = facedata.f2c[idx].right;

            if (facedata.left_interior_only[idx]) {
                facedata.gradients[idx].copy_values_from(celldata.gradients[l]);

            } else if (facedata.right_interior_only[idx]) {
                facedata.gradients[idx].copy_values_from(celldata.gradients[r]);

            } else {
                // Assume we have both cells, for a shared or internal boundary interface
                apply_haschelbacher_averaging(celldata.positions[l], celldata.positions[r], facedata.normals[idx],
                                              celldata.gradients[l], celldata.gradients[r],
                                              celldata.flowstates[l],celldata.flowstates[r],
                                              facedata.gradients[idx], nsp, nmodes, nturb, is3D);
            }
        } // end main face loop
    } // end average_lsq_cell_derivs_to_faces routine

    @nogc
    void first_order_flux_calc(size_t gtl, size_t[] face_idxs)
    {
        immutable size_t neq = myConfig.cqi.n;
        Vector3 gvel;
        gvel.clear();
        foreach(idx; face_idxs){
            size_t l = facedata.stencil_idxs[idx].L0;
            size_t r = facedata.stencil_idxs[idx].R0;

            Lft.copy_values_from(celldata.flowstates[l]);
            Rght.copy_values_from(celldata.flowstates[r]);

            facedata.flowstates[idx].copy_average_values_from(*Lft, *Rght);

            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
    }

    @nogc
    void second_order_flux_calc(size_t gtl, size_t[] face_idxs)
    {
    /*
        Eilmer's classic second-order piece-wise parabolic reconstruction, using a
        4 cell symmetric stencil.
        Notes: This routine currently ignores the supress_reconstruction face parameter.
    */
        immutable bool hpl = myConfig.apply_heuristic_pressure_based_limiting;
        immutable size_t neq = myConfig.cqi.n;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool is3D = (myConfig.dimensions == 3);
        immutable bool MHD = myConfig.MHD;
        immutable bool apply_limiter = myConfig.apply_limiter;
        immutable bool extrema_clipping = myConfig.extrema_clipping;
        immutable InterpolateOption ti = myConfig.thermo_interpolator;

        number beta = 1.0;
        Vector3 gvel;
        gvel.clear();

        foreach(idx; face_idxs){
            size_t L1 = facedata.stencil_idxs[idx].L1;
            size_t L0 = facedata.stencil_idxs[idx].L0;
            size_t R0 = facedata.stencil_idxs[idx].R0;
            size_t R1 = facedata.stencil_idxs[idx].R1;
            if (hpl) beta =
                compute_heuristic_pressure_limiter(celldata.flowstates[L1].gas.p,
                                                   celldata.flowstates[L0].gas.p,
                                                   celldata.flowstates[R0].gas.p,
                                                   celldata.flowstates[R1].gas.p);

            Lft.copy_values_from(celldata.flowstates[L0]);
            Rght.copy_values_from(celldata.flowstates[R0]);

            interp_l2r2(celldata.flowstates[L1], celldata.flowstates[L0],
                        celldata.flowstates[R0], celldata.flowstates[R1],
                        facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                        facedata.l2r2_interp_data[idx], nsp, nmodes, nturb,
                        ti, MHD, apply_limiter, extrema_clipping,
                        myConfig, *Lft, *Rght, beta);
            facedata.flowstates[idx].copy_average_values_from(*Lft, *Rght);

            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
        return;
    }

    @nogc
    void third_order_flux_calc(size_t gtl, size_t[] face_idxs)
    {
    /*
        Experimental high-order parabolic reconstruction, using a 6 cell
        symmetric stencil, and possibly the unusual ASF flux calculator.
    */
        immutable bool hpl = myConfig.apply_heuristic_pressure_based_limiting;
        immutable size_t neq = myConfig.cqi.n;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool is3D = (myConfig.dimensions == 3);
        immutable bool MHD = myConfig.MHD;
        immutable bool apply_limiter = myConfig.apply_limiter;
        immutable bool extrema_clipping = myConfig.extrema_clipping;
        immutable InterpolateOption ti = myConfig.thermo_interpolator;

        number beta = 1.0;
        Vector3 gvel;
        gvel.clear();

        foreach(idx; face_idxs){
            size_t L2 = facedata.stencil_idxs[idx].L2;
            size_t L1 = facedata.stencil_idxs[idx].L1;
            size_t L0 = facedata.stencil_idxs[idx].L0;
            size_t R0 = facedata.stencil_idxs[idx].R0;
            size_t R1 = facedata.stencil_idxs[idx].R1;
            size_t R2 = facedata.stencil_idxs[idx].R2;
            if (hpl) beta =
                compute_heuristic_pressure_limiter(celldata.flowstates[L2].gas.p,
                                                   celldata.flowstates[L1].gas.p,
                                                   celldata.flowstates[L0].gas.p,
                                                   celldata.flowstates[R0].gas.p,
                                                   celldata.flowstates[R1].gas.p,
                                                   celldata.flowstates[R2].gas.p);

            Lft.copy_values_from(celldata.flowstates[L0]);
            Rght.copy_values_from(celldata.flowstates[R0]);

            interp_l3r3(celldata.flowstates[L2], celldata.flowstates[L1],
                        celldata.flowstates[L0], celldata.flowstates[R0],
                        celldata.flowstates[R1], celldata.flowstates[R2],
                        facedata.l3r3_interp_data[idx], nsp, nmodes, nturb,
                        ti, MHD, apply_limiter, extrema_clipping,
                        myConfig, *Lft, *Rght, beta);

            facedata.flowstates[idx].copy_average_values_from(*Lft, *Rght);

            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
        return;
    }

    @nogc
    void asf_flux_calc(size_t gtl, size_t[] face_idxs, bool adaptive=false)
    {
    /*
        Low dissipation flux calculator, using the unusual ASF method.
    */
        immutable size_t neq = myConfig.cqi.n;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool is3D = (myConfig.dimensions == 3);

        number factor = 1.0;
        foreach(idx; face_idxs){
            size_t L1 = facedata.stencil_idxs[idx].L1;
            size_t L0 = facedata.stencil_idxs[idx].L0;
            size_t R0 = facedata.stencil_idxs[idx].R0;
            size_t R1 = facedata.stencil_idxs[idx].R1;

            facedata.flowstates[idx].copy_average_values_from(celldata.flowstates[L0], celldata.flowstates[R0]);

            if (adaptive) factor = 1.0 - facedata.flowstates[idx].S;

            // This routine uses and reuses a "flux" vector that belongs to the SFluidBlock
            // This is because we need to transform it to the global face, without messing
            // with the components of facedata.fluxes
            flux.clear();
            ASF_242(celldata.flowstates[L1], celldata.flowstates[L0],
                    celldata.flowstates[R0], celldata.flowstates[R1],
                    myConfig, nsp, nmodes, nturb, 
                    facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                    flux, factor);
            foreach(i; 0 .. neq) {
                facedata.fluxes[idx*neq +i] += flux[i];
            }
        }
        return;
    }

    @nogc void block_decode_conserved(int gtl, int ftl)
    {
        number[] U;
        if (ftl==0) { U = celldata.U0; }
        if (ftl==1) { U = celldata.U1; }
        if (ftl==2) { U = celldata.U2; }
        if (ftl==3) { U = celldata.U3; }
        if (ftl==4) { U = celldata.U4; }

        size_t ncq = myConfig.cqi.n;
        foreach (icell; 0 .. ncells) {
            size_t idx = icell*ncq;
            decode_conserved(celldata.positions[icell], U[idx .. idx+ncq],
                             celldata.flowstates[icell], omegaz, icell, myConfig);
        }
    }

    @nogc void block_decode_conserved_safe(int gtl, int ftl)
    {
        number[] U;
        if (ftl==0) { U = celldata.U0; }
        if (ftl==1) { U = celldata.U1; }
        if (ftl==2) { U = celldata.U2; }
        if (ftl==3) { U = celldata.U3; }
        if (ftl==4) { U = celldata.U4; }

        size_t ncq = myConfig.cqi.n;
        foreach (icell; 0 .. ncells) {
            size_t idx = icell*ncq;
            celldata.data_is_bad[icell] = false;
            try {
                decode_conserved(celldata.positions[icell], U[idx .. idx+ncq],
                                 celldata.flowstates[icell], omegaz, icell, myConfig);
            }  catch (FlowSolverException e) {
                celldata.data_is_bad[icell] = true;
            }
        }
    }

} // end class FluidBlock

@nogc void apply_haschelbacher_averaging(Vector3 lpos, Vector3 rpos, Vector3 n,
                                         ref const FlowGradients lgrad, ref const FlowGradients rgrad,
                                         ref const FlowState lfs, ref const FlowState rfs, ref FlowGradients grad,
                                         size_t nsp, size_t nmodes, size_t nturb, bool is3D){

    // vector from left-cell-centre to right-cell-centre
    number ex = rpos.x - lpos.x;
    number ey = rpos.y - lpos.y;
    number ez = rpos.z - lpos.z;
    // ehat
    number emag = sqrt(ex*ex + ey*ey + ez*ez);
    number ehatx = ex/emag;
    number ehaty = ey/emag;
    number ehatz = ez/emag;
    // ndotehat
    number nx = n.x;
    number ny = n.y;
    number nz = n.z;
    number ndotehat = nx*ehatx + ny*ehaty + nz*ehatz;
    nx /= ndotehat;
    ny /= ndotehat;
    nz /= ndotehat;

    grad.vel[0] = haselbacher(lgrad.vel[0],
                              rgrad.vel[0],
                              lfs.vel.x,
                              rfs.vel.x,
                              nx, ny, nz, ehatx, ehaty, ehatz, emag);

    grad.vel[1] = haselbacher(lgrad.vel[1],
                              rgrad.vel[1],
                              lfs.vel.y,
                              rfs.vel.y,
                              nx, ny, nz, ehatx, ehaty, ehatz, emag);

    if (is3D) {
        grad.vel[2] = haselbacher(lgrad.vel[2],
                                  rgrad.vel[2],
                                  lfs.vel.z,
                                  rfs.vel.z,
                                  nx, ny, nz, ehatx, ehaty, ehatz, emag);
    }
    version(multi_species_gas) {
        foreach (isp; 0 .. nsp) {
            grad.massf[isp] = haselbacher(lgrad.massf[isp],
                                          rgrad.massf[isp],
                                          lfs.gas.massf[isp],
                                          rfs.gas.massf[isp],
                                          nx, ny, nz, ehatx, ehaty, ehatz, emag);
        }
    }
    grad.T = haselbacher(lgrad.T,
                         rgrad.T,
                         lfs.gas.T,
                         rfs.gas.T,
                         nx, ny, nz, ehatx, ehaty, ehatz, emag);
    version(multi_T_gas) {
        foreach (imode; 0 .. nmodes) {
            grad.T_modes[imode] = haselbacher(lgrad.T_modes[imode],
                                              rgrad.T_modes[imode],
                                              lfs.gas.T_modes[imode],
                                              rfs.gas.T_modes[imode],
                                              nx, ny, nz, ehatx, ehaty, ehatz, emag);
        }
    }
    version(turbulence) {
        foreach(i; 0 .. nturb) {
            grad.turb[i] = haselbacher(lgrad.turb[i],
                                       rgrad.turb[i],
                                       lfs.turb[i],
                                       rfs.turb[i],
                                       nx, ny, nz, ehatx, ehaty, ehatz, emag);
        }
    } // end turbulence
}


@nogc pure number[3] haselbacher(number[3] Lgrad, number[3] Rgrad, number L, number R, number nx, number ny, number nz, number ehatx, number ehaty, number ehatz, number emag){
    number[3] grad;
    number avgdotehat = 0.5*(Lgrad[0]+Rgrad[0])*ehatx +
                        0.5*(Lgrad[1]+Rgrad[1])*ehaty +
                        0.5*(Lgrad[2]+Rgrad[2])*ehatz;
    number jump = avgdotehat - (R - L)/emag;
    grad[0] = 0.5*(Lgrad[0]+Rgrad[0]) - jump*(nx);
    grad[1] = 0.5*(Lgrad[1]+Rgrad[1]) - jump*(ny);
    grad[2] = 0.5*(Lgrad[2]+Rgrad[2]) - jump*(nz);
    return grad;
}
