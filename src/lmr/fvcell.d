/**
 * fvcell.d
 * Finite-volume cell class for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcell;

import std.conv;
import std.string;
import std.array;
import std.format;
import std.stdio;
import std.math;
import std.algorithm;
import nm.complex;
import nm.number;
import nm.bbla;
import geom;
import gas;
import kinetics;
import flowstate;
import flowgradients;
import conservedquantities;
import fvvertex;
import fvinterface;
import globalconfig;
import lsqinterp;
import gas.fuel_air_mix;
import globaldata : SimState;
import turbulence;
import celldata;

import kinetics.chemistry_update;
import kinetics.reaction_mechanism;

version(debug_chem) {
    GasState savedGasState;
}


// The following functions are used at compile time.
// Look for mixin statements further down in the file.
string avg_over_vtx_list(string quantity, string result)
{
    string code = result ~ " = 0.0; ";
    code ~= "foreach(v; vtx) { " ~ result ~ " += v." ~ quantity ~ "; } ";
    code ~= result ~ " /= vtx.length;";
    return code;
}
string avg_over_iface_list(string quantity, string result)
{
    string code = result ~ " = 0.0; ";
    code ~= "foreach(face; iface) { " ~ result ~ " += face." ~ quantity ~ "; } ";
    code ~= result ~ " /= iface.length;";
    return code;
}

class FVCell {
public:
    int id;  // allows us to work out where, in the block, the cell is
    bool data_is_bad; // Set to false at the start of an update.
    // Reset and checked at points through the update so that we don't stagger on
    // with bad data poisoning the simulation.
    //
    // local time-stepping
    double dt_local;
    double t_local;
    // super time-stepping
    number signal_hyp;
    number signal_parab;
    number dt_hyper;
    number dt_parab;
    //
    bool fr_reactions_allowed; // if true, will call thermochemical_increment
    double dt_chem; // acceptable time step for finite-rate chemistry
    double dt_therm; // acceptable time step for thermal relaxation
    bool in_turbulent_zone; // if true, we will keep the turbulence viscosity
    number base_qdot; // base-level of heat addition to cell, W/m**3
    // Geometry
    Vector3[] pos; // Centre x,y,z-coordinates for time-levels, m,m,m
    number[3] lengths; // length in i,j,k index direction
    number L_min;   // minimum length scale for cell
    number L_max;   // maximum length scale for cell
    // Connections
    FVInterface[] iface;  // references to defining interfaces of cell
    int[] outsign; // +1 if iface is outward-facing; -1 for an inward-facing iface
    FVVertex[] vtx;  // references to vertices for quad (2D) and hexahedral (3D) cells
    FVCell[] cell_cloud; // references to neighbouring cells
    // More geometry
    number[] volume; // Cell volume for time-levels (per unit depth or radian in 2D), m**3
    number[] areaxy; // (x,y)-plane area for time-levels, m**2
    // Flow
    // Although most do, some boundary conditions will not fill in
    // valid flow state data for the ghost cell. The following flag
    // is used for the unstructured-grid code to determine if we
    // should add the cell to the list of points in the cloud about
    // an interface location.
    // [TODO] PJ 2016-04-23, Consider if we should use this flag in
    // the context of structured grids also.
    // [TODO] PJ 2016-10-18, We should really only have cells where we truly
    // expect to have gas.  If we eliminate the use of ghost cells for boundaries
    // where there is no gas on the otherside of the "wall",
    // this flag might not be needed.
    // Such a change will need a big rework of boundary condition code.
    bool contains_flow_data;
    bool is_interior_to_domain; // true if the cell is interior to the flow domain
    bool allow_k_omega_update = true; // turbulent wall functions may turn this off
    FlowState fs; // Flow properties
    ConservedQuantities[] U;  // Conserved flow quantities for the update stages.
    ConservedQuantities[] dUdt; // Time derivatives for the update stages.
    ConservedQuantities Q; // source (or production) terms
    ConservedQuantities Qudf; // source terms from user-defined function (temporary storage, each update)
    ConservedQuantities[2] dUdt_copy; // for residual smoothing
    // for unstructured grids, we may be doing high-order reconstruction
    LSQInterpWorkspace* ws;
    LSQInterpGradients* gradients; // we only need these workspaces for the unstructured
                                  // solver, they are instantiated in ufluidblock.d
    // Viscous-flux-related quantities.
    FlowGradients* grad;
    WLSQGradWorkspace* ws_grad;
    Vector3*[] cloud_pos; // Positions of flow points for gradients calculation.
    FlowState*[] cloud_fs; // References to flow states at those points.

    // Heat flux used in the loosely coupled CHT solver (we currently store this here for convenience).
    number heat_transfer_into_solid; // TODO: We should think about whether there is a more appropriate place to store this. KAD 2022-11-08

    // Terms for loose-coupling of radiation.
    number Q_rad_org;
    number f_rad_org;
    number Q_rE_rad; // Rate of energy addition to cell via radiation.
    number Q_rE_rad_save; // Presently, the radiation source term is calculated
                          // at the first update stage.  We need to retain that
                          // value for all of the update stages.
    // Data for computing residuals.
    number rho_at_start_of_step, rE_at_start_of_step;
    // distance to nearest viscous wall (only computed if turb_model.needs_dwall)
    number dwall;

    // For use with LU-SGS solver/preconditioner (note: we don't need complex numbers here)
    number[] LU;
    number[] dUk;
    number[] dF;
    number[] scalar_diag_inv;
    Matrix!number dFdU;
    Matrix!number dFdU_rotated;

    // Arrays to store the local DFT values
    // Lengths are known at run time (globalconfig.n_DFT_modes) but not at compile time, handle lengths later
    number[] DFT_local_real;
    number[] DFT_local_imag;

    // array of auxiliary data
    AuxCellData[] aux_cell_data;

    // Electromagnetic Field Variables
    double electric_potential;
    double[2] electric_field;

    // Shape sensitivity calculator workspace
    FVCell[] cell_list;            // list of cells in the residual stencil
    FVInterface[] face_list;       // list of faces in the residual stencil
    version(newton_krylov) {
        number[][] dRdU;
        ConservedQuantities Q_save;
        FlowGradients* grad_save;
        LSQInterpGradients* gradients_save;
        Matrix!number dConservative;
    }

    // 2021-03-12: Note that we have moved the IO functions for the cell into fluidblockio.d
    // in preparation for the new block-level IO code.  We thus need to be able to access the
    // cell's local reference to the config data from over there.
public:
    LocalConfig myConfig;

public:
    @disable this();

    this(LocalConfig myConfig, bool allocate_spatial_deriv_lsq_workspace=false, int id_init=-1)
    {
        this.myConfig = myConfig;
        id = id_init;
        contains_flow_data = false; // initial presumption to be adjusted later
        is_interior_to_domain = false;
        pos.length = myConfig.n_grid_time_levels;
        volume.length = myConfig.n_grid_time_levels;
        areaxy.length = myConfig.n_grid_time_levels;

        GasModel gmodel = cast(GasModel) myConfig.gmodel;
        if (gmodel is null) { gmodel = GlobalConfig.gmodel_master; }

        int n_species = myConfig.n_species;
        int n_modes = myConfig.n_modes;
        double T = 300.0;
        double[] T_modes; foreach(i; 0 .. n_modes) { T_modes ~= 300.0; }
        double[] turb_init;
        foreach(i; 0 .. myConfig.turb_model.nturb)
            turb_init ~= myConfig.turb_model.turb_limits(i).re;
        fs = FlowState(gmodel, 100.0e3, T, T_modes, Vector3(0.0,0.0,0.0), turb_init);
        size_t ncq = myConfig.cqi.n; // number of conserved quantities
        foreach(i; 0 .. myConfig.n_flow_time_levels) {
            U ~= new_ConservedQuantities(ncq);
            U[i].clear();
            dUdt ~= new_ConservedQuantities(ncq);
        }
        Q = new_ConservedQuantities(ncq);
        Q.clear();
        Qudf = new_ConservedQuantities(ncq);
        Qudf.clear();
        if (myConfig.residual_smoothing) {
            dUdt_copy[0] = new_ConservedQuantities(ncq);
            dUdt_copy[1] = new_ConservedQuantities(ncq);
        }
        grad = new FlowGradients(myConfig);
        if (allocate_spatial_deriv_lsq_workspace) {
            ws_grad = new WLSQGradWorkspace();
        }
        //
        version(newton_krylov) {
            grad_save = new FlowGradients(myConfig);
            gradients_save = new LSQInterpGradients(myConfig.n_species, myConfig.n_modes, myConfig.turb_model.nturb);
            Q_save = new_ConservedQuantities(ncq);
            dRdU.length = ncq; // number of conserved variables
            foreach (ref a; dRdU) a.length = ncq;
            foreach (i; 0..dRdU.length) {
                foreach (j; 0..dRdU[i].length) {
                    dRdU[i][j] = 0.0;
                }
            }
        }
        version(debug_chem) {
            // The savedGasState is a module-level variable.
            // It only needs to be initialised when debug_chem mode
            // is on AND it only required initialisation once.
            savedGasState = GasState(gmodel);
        }
        //
        // some data structures used in the LU-SGS solver
        version(steady_state) {
            size_t nConserved = myConfig.cqi.n;
            scalar_diag_inv.length = nConserved;
            dFdU = new Matrix!number(nConserved,nConserved);
            dFdU.zeros;
            dFdU_rotated = new Matrix!number(nConserved,nConserved);
            dFdU_rotated.zeros;
            dF.length = nConserved;
            dUk.length = nConserved;
            dUk[] = to!number(0.0);
            LU.length = nConserved;
        }
        //
        DFT_local_real.length = myConfig.DFT_n_modes;
        DFT_local_imag.length = myConfig.DFT_n_modes;
        //
        // generate auxiliary data items
        aux_cell_data = AuxCellData.get_aux_cell_data_items(myConfig);
    }

    this(LocalConfig myConfig, in Vector3 pos, FlowState fs, in number volume, int id_init=-1)
    // stripped down initialisation
    {
        id = id_init;
        this.myConfig = myConfig;
        this.pos.length = 1;
        this.pos[0] = pos;
        this.volume.length = 1;
        this.volume[0] = volume;
        this.fs = fs;

        // generate auxiliary data items
        aux_cell_data = AuxCellData.get_aux_cell_data_items(myConfig);
    }

    // length in the i-index direction
    @property @nogc number iLength() const { return lengths[0]; }
    @property @nogc number iLength(number l) { return lengths[0] = l; }

    // length in the j-index direction
    @property @nogc number jLength() const { return lengths[1]; }
    @property @nogc number jLength(number l) { return lengths[1] = l; }

    // length in the k-index direction
    @property @nogc number kLength() const { return lengths[2]; }
    @property @nogc number kLength(number l) { return lengths[2] = l; }

    @nogc
    void copy_values_from(FVCell other, int type_of_copy)
    {
        switch ( type_of_copy ) {
        case CopyDataOption.minimal_flow:
            fs.copy_values_from(other.fs);
            break;
        case CopyDataOption.all_flow:
            fs.copy_values_from(other.fs);
            Q.copy_values_from(other.Q);
            foreach(i; 0 .. other.myConfig.n_flow_time_levels) {
                U[i].copy_values_from(other.U[i]);
                dUdt[i].copy_values_from(other.dUdt[i]);
            }
            break;
        case CopyDataOption.grid:
            foreach(i; 0 .. other.myConfig.n_grid_time_levels) {
                pos[i].set(other.pos[i]);
                volume[i] = other.volume[i];
                areaxy[i] = other.areaxy[i];
            }
            iLength = other.iLength;
            jLength = other.jLength;
            kLength = other.kLength;
            L_min = other.L_min;
            L_max = other.L_max;
            break;
        case CopyDataOption.cell_lengths_only:
            iLength = other.iLength;
            jLength = other.jLength;
            kLength = other.kLength;
            L_min = other.L_min;
            L_max = other.L_max;
            break;
        case CopyDataOption.all:
        default:
            // [TODO] really need to think about what needs to be copied...
            id = other.id;
            is_interior_to_domain = other.is_interior_to_domain;
            myConfig = other.myConfig;
            foreach(i; 0 .. other.myConfig.n_grid_time_levels) {
                pos[i].set(other.pos[i]);
                volume[i] = other.volume[i];
                areaxy[i] = other.areaxy[i];
            }
            iLength = other.iLength;
            jLength = other.jLength;
            kLength = other.kLength;
            L_min = other.L_min;
            L_max = other.L_max;
            fs.copy_values_from(other.fs);
            Q.copy_values_from(other.Q);
            grad.copy_values_from(*(other.grad));
            foreach(i; 0 .. other.myConfig.n_flow_time_levels) {
                U[i].copy_values_from(other.U[i]);
                dUdt[i].copy_values_from(other.dUdt[i]);
            }
        } // end switch
    }

    @nogc
    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
        pos[to_level] = pos[from_level];
        volume[to_level] = volume[from_level];
        areaxy[to_level] = areaxy[from_level];
        // When working over all cells in a block, the following copies
        // will no doubt do some doubled-up work, but it should be otherwise benign.
        foreach(ref face; iface) {
            if (face) face.copy_grid_level_to_level(from_level, to_level);
        }
        foreach(ref v; vtx) {
            if (v) v.copy_grid_level_to_level(from_level, to_level);
        }
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "FVCell(";
        repr ~= "id=" ~ to!string(id);
        repr ~= ", universe_blk_id=" ~ to!string(myConfig.universe_blk_id);
        repr ~= ", pos=" ~ to!string(pos);
        repr ~= ", iface_ids=["; foreach (f; iface) { repr ~= format("%d,", f.id); } repr ~= "]";
        repr ~= ", outsigns=["; foreach (osgn; outsign) { repr ~= format("%d,", osgn); } repr ~= "]";
        repr ~= ", vtx_ids=["; foreach (v; vtx) { repr ~= format("%d,", v.id); } repr ~= "]";
        repr ~= ",\n... volume=" ~ to!string(volume);
        repr ~= ", areaxy=" ~ to!string(areaxy);
        repr ~= ", iLength=" ~ to!string(iLength);
        repr ~= ", jLength=" ~ to!string(jLength);
        repr ~= ", kLength=" ~ to!string(kLength);
        repr ~= ", L_min=" ~ to!string(L_min);
        repr ~= ", dt_chem=" ~ to!string(dt_chem);
        repr ~= ", dt_therm=" ~ to!string(dt_therm);
        repr ~= ", in_turbulent_zone=" ~ to!string(in_turbulent_zone);
        repr ~= ", fr_reactions_allowed=" ~ to!string(fr_reactions_allowed);
        repr ~= ", contains_flow_data=" ~ to!string(contains_flow_data);
        repr ~= ", allow_k_omega_update=" ~ to!string(allow_k_omega_update);
        repr ~= ",\n... fs=" ~ to!string(fs);
        repr ~= ",\n... U=" ~ to!string(U);
        repr ~= ",\n... dUdt=" ~ to!string(dUdt);
        repr ~= ")";
        return to!string(repr);
    }

    @nogc
    int universe_blk_id()
    {
        return myConfig.universe_blk_id;
    }

    @nogc
    void update_2D_geometric_data(size_t gtl, bool axisymmetric)
    {
        string msg = "FVCell.update_2D_geometric_data(): ";
        number vol, xyplane_area, iL, jL;
        switch (vtx.length) {
        case 3:
            xyplane_triangle_cell_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl],
                                             pos[gtl], xyplane_area, iL, jL, L_min);
            iLength = iL; jLength = jL;
            break;
        case 4:
            xyplane_quad_cell_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                                         vtx[2].pos[gtl], vtx[3].pos[gtl],
                                         pos[gtl], xyplane_area, iL, jL, L_min);
            iLength = iL; jLength = jL;
            break;
        default:
            debug { msg ~= format("Unhandled number of vertices: %d", vtx.length); }
            throw new FlowSolverException(msg);
        } // end switch
        // Cell Volume.
        if (axisymmetric) {
            // Volume per radian = centroid y-ordinate * cell area
            vol = xyplane_area * pos[gtl].y;
        } else {
            // Assume unit depth in the z-direction.
            vol = xyplane_area;
        }
        if (vol < 0.0) {
            msg = "Negative cell volume";
            debug {
                msg ~= format(" for cell[%d]= %g\n", id, vol);
                foreach (i; 0 .. vtx.length) {
                    msg ~= format("   vtx[%d].pos[%d]=%s\n", i, gtl, vtx[i].pos[gtl].toString);
                }
            }
            throw new FlowSolverException(msg);
        }
        volume[gtl] = vol;
        areaxy[gtl] = xyplane_area;
        kLength = to!number(0.0);
        L_max = fmax(iLength, jLength);
    } // end update_2D_geometric_data()

    @nogc
    void update_3D_geometric_data(size_t gtl, bool true_centroid)
    {
        string msg = "FVCell.update_3D_geometric_data(): ";
        number iL, jL, kL;
        switch (vtx.length) {
        case 4:
            tetrahedron_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                                   vtx[2].pos[gtl], vtx[3].pos[gtl],
                                   pos[gtl], volume[gtl], L_min);
            iLength = L_min; jLength = L_min; kLength = L_min;
            break;
        case 8:
            hex_cell_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl], vtx[3].pos[gtl],
                                vtx[4].pos[gtl], vtx[5].pos[gtl], vtx[6].pos[gtl], vtx[7].pos[gtl],
                                true_centroid, pos[gtl], volume[gtl], iL, jL, kL);
            iLength = iL; jLength = jL; kLength = kL;
            L_min = min(iLength, jLength, kLength);
            break;
        case 5:
            pyramid_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl], vtx[3].pos[gtl],
                               vtx[4].pos[gtl], true_centroid, pos[gtl], volume[gtl], L_min);
            iLength = L_min; jLength = L_min; kLength = L_min;
            break;
        case 6:
            wedge_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl],
                             vtx[3].pos[gtl], vtx[4].pos[gtl], vtx[5].pos[gtl],
                             true_centroid, pos[gtl], volume[gtl], L_min);
            iLength = L_min; jLength = L_min; kLength = L_min;
            break;
        default:
            debug { msg ~= format("Unhandled number of vertices: %d", vtx.length); }
            throw new FlowSolverException(msg);
        } // end switch
        if (volume[gtl] <= 0.0) {
            debug {
                msg ~= format("Invalid volume %g for cell %d in block %d at pos %s",
                              volume[gtl], id, myConfig.universe_blk_id, pos[gtl]);
                msg ~= format(" Lmin=%g vtx.length=%d", L_min, vtx.length);
                foreach (i; 0 .. vtx.length) {
                    msg ~= format(" vtx[%d].pos=%s", i, vtx[i].pos[gtl]);
                }
            }
            throw new FlowSolverException(msg);
        }
        L_max = fmax(fmax(iLength, jLength), kLength);
    } // end update_3D_geometric_data()

    void replace_flow_data_with_average(in FVCell[] others)
    {
        auto gmodel = myConfig.gmodel;
        size_t n = others.length;
        if (n == 0) throw new FlowSolverException("Need to average from a nonempty array.");
        FlowState*[] fsList;
        // We need to be honest and not to fiddle with the other gas states.
        foreach(other; others) {
            if ( this is other ) {
                throw new FlowSolverException("Must not include destination in source list.");
            }
            fsList ~= cast(FlowState*)&(other.fs);
        }
        fs.copy_average_values_from(fsList, gmodel);
        // Accumulate from a clean slate and then divide.
        Q_rE_rad = 0.0;
        foreach(other; others) {
            Q_rE_rad += other.Q_rE_rad;
        }
        Q_rE_rad /= n;
    } // end replace_flow_data_with_average()

    @nogc
    void encode_conserved(int gtl, int ftl, double omegaz)
    // gtl = grid time level
    // ftl = flow time level
    {
        auto cqi = myConfig.cqi;
        ConservedQuantities myU = U[ftl];
        number myrho = fs.gas.rho;
        // Mass per unit volume.
        myU[cqi.mass] = myrho;
        // Momentum per unit volume.
        myU[cqi.xMom] = fs.gas.rho*fs.vel.x;
        myU[cqi.yMom] = fs.gas.rho*fs.vel.y;
        if (cqi.threeD) { myU[cqi.zMom] = fs.gas.rho*fs.vel.z; }
        version(MHD) {
            // Magnetic field
            if (cqi.MHD) {
                myU[cqi.xB] = fs.B.x;
                myU[cqi.yB] = fs.B.y;
                myU[cqi.zB] = fs.B.z;
                myU[cqi.psi] = fs.psi;
                myU[cqi.divB] = fs.divB;
            }
        }
        // Total Energy / unit volume
        number u = myConfig.gmodel.internal_energy(fs.gas);
        number ke = 0.5*(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y+fs.vel.z*fs.vel.z);
        myU[cqi.totEnergy] = fs.gas.rho*(u + ke);
        version(turbulence) {
            if (cqi.turb) {
                foreach(i; 0 .. myConfig.turb_model.nturb){
                    myU[cqi.rhoturb+i] = fs.gas.rho * fs.turb[i];
                }
                myU[cqi.totEnergy] += fs.gas.rho * myConfig.turb_model.turbulent_kinetic_energy(fs);
            }
        }
        version(MHD) {
            if (myConfig.MHD) {
                number me = 0.5*(fs.B.x*fs.B.x + fs.B.y*fs.B.y + fs.B.z*fs.B.z);
                myU[cqi.totEnergy] += me;
            }
        }
        version(multi_T_gas) {
            // Other internal energies: energy in mode per unit volume.
            foreach(imode; 0 .. cqi.n_modes) {
                myU[cqi.modes+imode] = fs.gas.rho*fs.gas.u_modes[imode];
            }
        }
        if (omegaz != 0.0) {
            // Rotating frame.
            // Finally, we adjust the total energy to make rothalpy.
            // We do this last because the gas models don't know anything
            // about rotating frames and we don't want to mess their
            // energy calculations around.
            number rho = fs.gas.rho;
            number x = pos[gtl].x;
            number y = pos[gtl].y;
            number rsq = x*x + y*y;
            // The conserved quantity is rothalpy. I = E - (u**2)/2
            // where rotating frame velocity  u = omegaz * r.
            myU[cqi.totEnergy] -= rho*0.5*omegaz*omegaz*rsq;
        }
        version(multi_species_gas) {
            // Species densities: mass of species is per unit volume.
            if (cqi.n_species > 1) {
                foreach(isp; 0 .. cqi.n_species) {
                    myU[cqi.species+isp] = fs.gas.rho*fs.gas.massf[isp];
                }
            }
        }
        assert(U[ftl][cqi.mass] > 0.0, "invalid density in conserved quantities vector" ~
               " at end of FVCell.encode_conserved().");
        return;
    } // end encode_conserved()

    @nogc
    int decode_conserved(int gtl, int ftl, double omegaz,
                         double tolerance=0.0, double assert_error_tolerance=0.1)
    {
        auto cqi = myConfig.cqi;
        auto gmodel = myConfig.gmodel;
        ConservedQuantities myU = U[ftl];
        // The conserved quantities are carried as quantity per unit volume.
        // mass / unit volume = density
        number rho = 0.0;
        if (cqi.mass == 0) {
            rho = myU[cqi.mass];
        } else {
            // if cqi.mass \= 0, then we assume that we have dropped the mass continuity equation from the solver,
            // in that case, the total density is a sum of the species densities
            foreach(isp; 0 .. cqi.n_species) {
                // if (myU[cqi.species+isp] < 0.0) myU[cqi.species+isp] = to!number(0.0);
                // We were originally imposing positivity on the species densities here to prevent slightly negative mass fractions,
                // this was causing noisy convergence for some problems using the steady-state accelerator due to the artificial nature
                // of suppressing the negative mass fractions after an update (i.e. the numerical solution really should have negative
                // mass fractions due to the numerical methods employed).
                // We are trialing the removal of this enforcement (note that this will only take effect for the steady-state solver here).
                // KAD 2022-08-22.
                rho += myU[cqi.species+isp];
            }
        }
        if (!(rho.re > 0.0)) {
            if (myConfig.adjust_invalid_cell_data) {
                data_is_bad = true;
                // We can do nothing more with the present data but the caller may
                // be able to replace the data with other nearby-cell data.
                return -1;
            } else {
                debug {
                    writeln("FVCell.decode_conserved(): Density invalid in conserved quantities.");
                    writeln("  universe-blk-id= ", myConfig.universe_blk_id, " cell-id= ", id);
                    writeln("  x= ", pos[gtl].x, " y= ", pos[gtl].y, " z= ", pos[gtl].z);
                    writeln("  gas= ", fs.gas);
                    writeln("  ftl= ", ftl, " (flow-time-level)");
                    writeln("  U[ftl]= ", myU);
                    writeln("  U[0]= ", U[0]);
                    writeln("  interfaces:", iface.length);
                    foreach(i, f; iface) { writeln("    iface[", i, "]= ", f); }
                }
                throw new FlowSolverException("Bad cell with negative mass.");
            }
        } // end if mass is not positive
        fs.gas.rho = rho; // This is limited to nonnegative and finite values.
        number dinv = 1.0 / rho;

        // Velocities from momenta.
        number zMom = (cqi.threeD) ? myU[cqi.zMom] : to!number(0.0);
        fs.vel.set(myU[cqi.xMom]*dinv, myU[cqi.yMom]*dinv, zMom*dinv);
        version(MHD) {
            // Magnetic field.
            if (cqi.MHD) {
                fs.B.set(myU[cqi.xB], myU[cqi.yB], myU[cqi.zB]);
                fs.psi = myU[cqi.psi];
                fs.divB = myU[cqi.divB];
            }
        }
        // Divide up the total energy per unit volume.
        number rE;
        if (omegaz != 0.0) {
            // Rotating frame.
            // The conserved quantity is rothalpy so we need to convert
            // back to enthalpy to do the rest of the decode.
            number x = pos[gtl].x;
            number y = pos[gtl].y;
            number rsq = x*x + y*y;
            rE = myU[cqi.totEnergy] + rho*0.5*omegaz*omegaz*rsq;
        } else {
            // Non-rotating frame.
            rE = myU[cqi.totEnergy];
        }
        version(MHD) {
            number me = 0.0;
            if (cqi.MHD) {
                me = 0.5*(fs.B.x*fs.B.x + fs.B.y*fs.B.y + fs.B.z*fs.B.z);
                rE -= me;
            }
        }
        // Start with the total energy, then take out the other components.
        // Internal energy is what remains.
        number u = rE * dinv;
        version(turbulence) {
            if (cqi.turb) {
                if (allow_k_omega_update) {
                    foreach(i; 0 .. myConfig.turb_model.nturb) {
                        if (isNaN(myU[cqi.rhoturb+i]))
                            throw new FlowSolverException("Turbulent quantity is Not A Number.");
                        // enforce positivity of turbulence model conserved quantities
                        if (myU[cqi.rhoturb+i] < 0.0) {
                            foreach (j; 0 .. myConfig.turb_model.nturb) {
                                // taking the absolute value has been observed to be more
                                // stable than clipping them to a small value (e.g. 1.0e-10)
                                myU[cqi.rhoturb+j] = fabs(myU[cqi.rhoturb+j]);
                            }
                        }
                        fs.turb[i] = myU[cqi.rhoturb+i] * dinv;
                    }
                }
                u -= myConfig.turb_model.turbulent_kinetic_energy(fs);
            }
        }
        // Remove kinetic energy for bulk flow.
        number ke = 0.5*(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y + fs.vel.z*fs.vel.z);
        u -= ke;
        // Other energies, if any.
        version(multi_T_gas) {
            number u_other = 0.0;
            foreach(imode; 0 .. gmodel.n_modes) { fs.gas.u_modes[imode] = myU[cqi.modes+imode] * dinv; }
            foreach(ei; fs.gas.u_modes) { u_other += ei; }
            fs.gas.u = u - u_other;
        } else {
            fs.gas.u = u;
        }
        // Thermochemical species, if appropriate.
        version(multi_species_gas) {
            // scale the species densities such that they sum to the total density (this ensures massfs sum to 1)
            // if cqi.mass \= 0, then we have dropped the mass continuity equation, so this step would be redundant
            if (myConfig.n_species > 1 && cqi.mass == 0) {
                number rhos_sum = 0.0;
                foreach(isp; 0 .. cqi.n_species) {
                    // impose positivity on the species densities
                    if (myU[cqi.species+isp] < 0.0) myU[cqi.species+isp] = to!number(0.0);
                    rhos_sum += myU[cqi.species+isp];
                }
                if (fabs(rhos_sum - rho) > assert_error_tolerance) {
                    string msg = "Sum of species densities far from total density";
                    debug {
                        msg ~= format(": fabs(rhos_sum - rho) = %s \n", fabs(rhos_sum - rho));
                        msg ~= format("  assert_error_tolerance = %s \n", assert_error_tolerance);
                        msg ~= format("  tolerance = %s \n", tolerance);
                        msg ~= "]\n";
                    }
                    throw new FlowSolverException(msg);
                }
                if ( fabs(rhos_sum - rho) > tolerance ) {
                    number scale_factor = rho/rhos_sum;
                    foreach(isp; 0 .. cqi.n_species) { myU[cqi.species+isp] *= scale_factor; }
                }
            } // end if (myConfig.n_species > 1 && cqi.mass == 0)
            try {
                if (cqi.n_species > 1) {
                    foreach(isp; 0 .. cqi.n_species) { fs.gas.massf[isp] = myU[cqi.species+isp] * dinv; }
                } else {
                    fs.gas.massf[0] = 1.0;
                }
		if (myConfig.sticky_electrons) { gmodel.balance_charge(fs.gas); }
	    } catch (GasModelException err) {
		if (myConfig.adjust_invalid_cell_data) {
		    data_is_bad = true;
		    return -2;
		} else {
		    string msg = "Bad cell with mass fractions that do not add correctly.";
		    debug {
			msg ~= format(" scale_mass_fractions exception with message:\n  %s", err.msg);
			msg ~= format("The decode_conserved() failed for cell: %d\n", id);
			msg ~= format("  This cell is located at: %s\n", pos[0]);
			msg ~= format("  This cell is located in block: %d\n", myConfig.universe_blk_id);
			msg ~= format("  The gas state before thermo update is:\n   fs.gas %s", fs.gas);
		    }
		    throw new FlowSolverException(msg);
		} // end if
	    } // end catch
        }
        //
        // Fill out the other variables: P, T, a, and viscous transport coefficients.
        try {
            try {
                gmodel.update_thermo_from_rhou(fs.gas);
            } catch (GasModelException err) {
                // Oops, it seems that the thermo update has failed to work
                // using the internal energy and density that have been
                // decoded from the current conserved quantities.
                if (myConfig.ignore_low_T_thermo_update_failure && (rho > 0.0)) {
                    // This small-energy, hopefully-transient error may get
                    // washed out of the flow field, so let's try to keep going.
                    // We reset the thermo data to an acceptable low-T state
                    // and make the current conserved quantities consistent.
                    fs.gas.T = myConfig.suggested_low_T_value;
                    version(multi_T_gas) {
                        foreach(i; 0 .. gmodel.n_modes) {
                            fs.gas.T_modes[i] = myConfig.suggested_low_T_value;
                        }
                    }
                    gmodel.update_thermo_from_rhoT(fs.gas);
                    encode_conserved(gtl, ftl, omegaz);
                } else {
                    // We do not ignore the thermo update error at this point.
                    throw err;
                }
            }
            if (fs.gas.T<=0.0) throw new FlowSolverException("update_thermo returned negative temperature.");
            version(multi_T_gas) {
                foreach(i; 0 .. gmodel.n_modes) {
                    if (fs.gas.T_modes[i]<=0.0) throw new FlowSolverException("update_thermo returned negative T_modes.");
                }
            }
            gmodel.update_sound_speed(fs.gas);
            if (myConfig.viscous) gmodel.update_trans_coeffs(fs.gas);
        } catch (GasModelException err) {
            if (myConfig.adjust_invalid_cell_data) {
                data_is_bad = true;
                return -2;
            } else {
                string msg = "Bad cell with failed thermodynamic update.";
                debug {
                    msg ~= format(" thermodynamic update exception with message:\n  %s", err.msg);
                    msg ~= format("The decode_conserved() failed for cell: %d\n", id);
                    msg ~= format("  This cell is located at: %s\n", pos[0]);
                    msg ~= format("  This cell is located in block: %d\n", myConfig.universe_blk_id);
                    msg ~= format("  The gas state after the failed update is:\n   fs.gas %s", fs.gas);
                }
                throw new FlowSolverException(msg);
            } // end if
        } // end catch
        //
        if (myConfig.radiation_energy_dump_allowed &&
            fs.gas.T > myConfig.radiation_energy_dump_temperature_limit) {
            // Dump excess energy and blame radiation.
            fs.gas.T = myConfig.radiation_energy_dump_temperature_limit;
            gmodel.update_thermo_from_rhoT(fs.gas);
            encode_conserved(gtl, ftl, omegaz);
            gmodel.update_sound_speed(fs.gas);
            if (myConfig.viscous) { gmodel.update_trans_coeffs(fs.gas); }
        }
        return 0; // success
    } // end decode_conserved()

    @nogc
    void time_derivatives(int gtl, int ftl)
    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    // gtl : (grid-time-level) flow derivatives are evaluated at this grid level
    // ftl : (flow-time-level) specifies where computed derivatives are to be stored.
    //       0: Start of stage-1 update.
    //       1: End of stage-1.
    //       2: End of stage-2.
    {
        auto my_dUdt = dUdt[ftl];
        auto cqi = myConfig.cqi;
        number vol_inv = 1.0 / volume[gtl]; // Cell volume (inverted).
        foreach (j; 0 .. cqi.n) {
            number surface_integral = to!number(0.0);
            // Integrate the fluxes across the interfaces that bound the cell.
            foreach(i; 0 .. iface.length) {
                number area = outsign[i] * iface[i].area[gtl];
                surface_integral -= iface[i].F[j] * area;
            }
            // Then evaluate the derivatives of conserved quantities.
            // Conserved quantities are stored per-unit-volume.
            my_dUdt[j] = vol_inv * surface_integral + Q[j];
        }
    } // end time_derivatives()


    @nogc
    void thermochemical_increment(double dt)
    // Use the finite-rate chemistry module to update the species fractions
    // and the other thermochemical properties.
    {
        if (!fr_reactions_allowed || fs.gas.T <= myConfig.T_frozen) return;
        number T_save = fs.gas.T;
        if (myConfig.ignition_zone_active) {
            // When active, replace gas temperature with an effective ignition temperature
            foreach(zone; myConfig.ignition_zones) {
                if (zone.is_inside(pos[0], myConfig.dimensions)) { fs.gas.T = zone.Tig; }
            }
        }

        number[maxParams] params;
        if ((cast(FuelAirMix) myConfig.gmodel) !is null) {
            // for this gas model thermochemical reactor we need turbulence info
            if (params.length < 1) { throw new Error("params vector too short."); }
            version(turbulence) {
                params[0]=myConfig.turb_model.turbulent_signal_frequency(fs);
            } else {
                throw new Error("FuelAirMix needs komega capability.");
            }
        }

        // Take a copy of dt_chem since it will be modified during the update.
        // The copy is useful to print to the screen if there's a failure of the
        // chemistry update.
        double dt_chem_save = dt_chem;

        if (myConfig.sticky_electrons) { myConfig.gmodel.balance_charge(fs.gas); }

        version(debug_chem) {
            savedGasState.copy_values_from(fs.gas);
        }

        try {
            myConfig.thermochemUpdate(fs.gas, dt, dt_chem, params);
            if (myConfig.ignition_zone_active) {
                // Restore actual gas temperature
                fs.gas.T = T_save;
            }
        } catch(ThermochemicalReactorUpdateException err) {
            // It's probably worth one more try but setting dt_chem = -1.0 to give
            // the ODE solver a fresh chance to find a good timestep.
            dt_chem = -1.0;
            try {
                 myConfig.thermochemUpdate(fs.gas, dt, dt_chem, params);
                 if (myConfig.ignition_zone_active) {
                     // Restore actual gas temperature
                     fs.gas.T = T_save;
                 }
            } catch(ThermochemicalReactorUpdateException err) {
                string msg = "The thermochemical_increment() failed.";
                debug {
                    msg ~= format("\nFOR CELL: %d\n", id);
                    msg ~= format("CAUGHT: %s\n", err.msg);
                    msg ~= format("This cell is located at: %s\n", pos[0]);
                    msg ~= format("This cell is located in block: %d\n", myConfig.universe_blk_id);
                    msg ~= format("The cell's id is: %d\n", id);
                    msg ~= format("The flow timestep is: %12.6e\n", dt);
                    msg ~= format("The initial attempted dt_chem is: %12.6e\n", dt_chem_save);
                    version(debug_chem) {
                        msg ~= format("The gas state BEFORE thermochemUpdate was:\n %s\n", savedGasState);
                    }
                    msg ~= format("The gas state AFTER the failed update is:\n   fs.gas %s", fs.gas);
                }
                throw new FlowSolverException(msg);
            }
        }

        // The update only changes mass fractions; we need to impose
        // a thermodynamic constraint based on a call to the equation of state.
        try {
            myConfig.gmodel.update_thermo_from_rhou(fs.gas);
        }
        catch (Exception err) {
            string msg = "The thermochemical_increment() failed update_thermo_from_rhou";
            debug {
                msg ~= format("\nfor cell: %d\n", id);
                msg ~= format("caught %s", err.msg);
                msg ~= format("This cell is located at: %s\n", pos[0]);
                msg ~= format("This cell is located in block: %d\n", myConfig.universe_blk_id);
                msg ~= "This failure occurred when trying to update the thermo state after\n";
                msg ~= "computing the species change due to chemical reactions.\n";
                version(debug_chem) {
                    msg ~= format("The gas state BEFORE thermochemUpdate was:\n %s\n", savedGasState);
                }
                msg ~= format("The present gas state is:\n   fs.gas %s", fs.gas);
            }
            throw new FlowSolverException(msg);
        }

        // If we are doing a viscous sim, we'll need to ensure
        // viscous properties are up-to-date
        if (myConfig.viscous) myConfig.gmodel.update_trans_coeffs(fs.gas);
        // [TODO] if ( myConfig.diffusion ) myConfig.gmodel.update_diffusion_coeffs(fs.gas);

        // Finally, we have to manually update the conservation quantities
        // for the gas-dynamics time integration.
        auto cqi = myConfig.cqi;
        version(multi_species_gas) {
            // Species densities: mass of species isp per unit volume.
            if (cqi.n_species > 1) {
                foreach(isp; 0 .. fs.gas.massf.length) {
                    U[0][cqi.species+isp] = fs.gas.rho * fs.gas.massf[isp];
                }
            }
        }
        version(multi_T_gas) {
            // Independent energies energy: Joules per unit volume.
            foreach(imode; 0 .. fs.gas.u_modes.length) {
                U[0][cqi.modes+imode] = fs.gas.rho * fs.gas.u_modes[imode];
            }
        }
    } // end thermochemical_increment()

    @nogc
    double signal_frequency()
    // Remember to use stringent_cfl=true for unstructured-grid.
    {
        number signal = 0; // Signal speed is something like a frequency, with units 1/s.
        //
        // Check the convective/wave-driven time step limit first,
        // then add a component to ensure viscous stability.
        //
        // Look at gas-dynamic signal speeds along each face.
        // This works for gas-dynamics only (not MHD), on a structured grid.
        //
        // Get the local normal velocities by rotating the local frame of reference.
        // Also, compute the velocity magnitude and recall the minimum length.
        number un_N = fabs(fs.vel.dot(iface[Face.north].n));
        number un_E = fabs(fs.vel.dot(iface[Face.east].n));
        // just in case we are given a non-hex cell
        size_t third_face = min(Face.top, iface.length-1);
        number un_T = (myConfig.dimensions == 3) ? fabs(fs.vel.dot(iface[third_face].n)) : to!number(0.0);
        if (myConfig.stringent_cfl) {
            // Compute the speed with the shortest length and the highest speed.
            number un_max = fmax(un_N, un_E);
            number minLength = fmin(iLength, jLength);
            if (myConfig.dimensions == 3) {
                un_max = fmax(un_max, un_T);
                minLength = fmin(minLength, kLength);
            }
            signal = (un_max + fs.gas.a) / minLength;
        } else {
            // Compute the signal speed in each index direction.
            number signalN = (un_N + fs.gas.a) / jLength;
            signal = fmax(signal, signalN);
            number signalE = (un_E + fs.gas.a) / iLength;
            signal = fmax(signal, signalE);
            if (myConfig.dimensions == 3) {
                number signalT = (un_T + fs.gas.a) / kLength;
                signal = fmax(signal, signalT);
            }
        }
        this.signal_hyp = signal; // store hyperbolic signal for STS
        // Factor for the viscous time limit.
        // See Swanson, Turkel and White (1991)
        // This factor is not included if viscosity is zero.
        if (myConfig.viscous && (fs.gas.mu > 10.0e-23)) {
            auto gmodel = myConfig.gmodel;
            number gam_eff = gmodel.gamma(fs.gas);
            // Need to sum conductivities for thermal nonequilibrium.
            number k_total = fs.gas.k;
            version(multi_T_gas) {
                foreach(k_value; fs.gas.k_modes) { k_total += k_value; }
            }
            number Prandtl = fs.gas.mu * gmodel.Cp(fs.gas) / k_total;
            signal += 4.0 * myConfig.viscous_factor * (fs.gas.mu + fs.mu_t)
                * gam_eff / (Prandtl * fs.gas.rho)
                * 1.0/(L_min^^2) * myConfig.viscous_signal_factor;
        }
        this.signal_parab = signal - this.signal_hyp; // store parabolic signal for STS
        version(turbulence) {
            number turbulent_signal = myConfig.turb_model.turbulent_signal_frequency(fs);
            turbulent_signal *= myConfig.turbulent_signal_factor;
            signal = fmax(signal, turbulent_signal);
            this.signal_parab = fmax(signal_parab, turbulent_signal);
        }
        version(MHD) {
            if (myConfig.MHD) {
                assert(myConfig.stringent_cfl, "MHD seems to only works if stringent_cfl is used.");
                // Gas dynamics speed
                // Ignoring flow and index directions, make the worst case assumptions.
                number u_mag_sq = (fs.vel.x)^^2 + (fs.vel.y)^^2;
                if (myConfig.dimensions == 3) { u_mag_sq += (fs.vel.z)^^2; }
                number u_mag = sqrt(u_mag_sq);
                // MHD signal speed
                number B_mag_sq = (fs.B.x)^^2 + (fs.B.y)^^2 + (fs.B.z)^^2;
                number ca2 = B_mag_sq / fs.gas.rho;
                number cfast = sqrt(ca2 + (fs.gas.a)^^2);
                signal = fmax(signal, (u_mag+cfast)/L_min);
            }
        }
        return signal.re;
    } // end signal_frequency()

    @nogc
    void turbulence_viscosity_zero()
    {
        fs.mu_t = 0.0;
        fs.k_t = 0.0;
    }

    @nogc
    double divergence_damping_factor(double dt, double c_h, double divB_damping_length)
    //Divergence factor factor used to scale the cleaning factor psi after each timestep.
    {
        double c_h2 = c_h * c_h;
        double c_p2 = 0.18 * divB_damping_length * c_h;
        return exp(-(c_h2 / c_p2) * dt);
    }

    @nogc
    void turbulence_viscosity_zero_if_not_in_zone()
    {
        if ( in_turbulent_zone ) {
            /* Do nothing, leaving the turbulence quantities as set. */
        } else {
            /* Presume this part of the flow is laminar; clear turbulence quantities. */
            fs.mu_t = 0.0;
            fs.k_t = 0.0;
        }
    }

    @nogc
    void turbulence_viscosity_limit(double factor)
    // Limit the turbulent viscosity to reasonable values relative to
    // the local molecular viscosity.
    // In shock started flows, we seem to get crazy values on the
    // starting shock structure and the simulations do not progress.
    {
        fs.mu_t = fmin(fs.mu_t, factor * fs.gas.mu);
        fs.k_t = fmin(fs.k_t, factor * fs.gas.k); // ASSUMPTION re k
    }

    @nogc
    void turbulence_viscosity_factor(double factor)
    // Scale the turbulent viscosity to model effects
    // such as not-fully-developed turbulence that might be expected
    // in short-duration transient flows.
    {
        fs.mu_t *= factor;
        fs.k_t *= factor;
    }

    @nogc
    void turbulence_viscosity()
    {
        auto gmodel = myConfig.gmodel;
        fs.mu_t = myConfig.turb_model.turbulent_viscosity(fs, *grad, pos[0].y, dwall);
        fs.k_t = myConfig.turb_model.turbulent_conductivity(fs, gmodel);
    }

    /*
    Old k-omega stuff moved (see NNG 11/02/20)
    */

    @nogc
    void clear_source_vector()
    // When doing the gasdynamic update stages, the source vector values
    // are accumulated for the inviscid and then viscous terms, so we
    // have to start with a clean slate, so to speak.
    {
        Q.clear();
    }

    @nogc
    void add_inviscid_source_vector(int gtl, double omegaz=0.0)
    // Add the components of the source vector, Q, for inviscid flow.
    //
    // Currently, the axisymmetric equations include the
    // pressure contribution to the y-momentum equation
    // here rather than in the boundary fluxes.
    // By default, assume 2D-planar, or 3D-Cartesian flow.
    {
        auto cqi = myConfig.cqi;
        //
        if (omegaz != 0.0) {
            // Rotating frame.
            number rho = fs.gas.rho;
            number x = pos[gtl].x;
            number y = pos[gtl].y;
            number wx = fs.vel.x;
            number wy = fs.vel.y;
            // Coriolis and centrifugal forces contribute to momenta.
            Q[cqi.xMom] += rho * (omegaz*omegaz*x + 2.0*omegaz*wy);
            Q[cqi.yMom] += rho * (omegaz*omegaz*y - 2.0*omegaz*wx);
            // There is no contribution to the energy equation in the rotating frame
            // because it is implicit in the use of rothalpy as the conserved quantity.
        }
        if (myConfig.axisymmetric) {
            // For axisymmetric flow:
            // pressure contribution from the Front and Back (radial) interfaces.
            Q[cqi.yMom] += fs.gas.p * areaxy[gtl] / volume[gtl];
        }
        if (myConfig.gravity_non_zero) {
            number rho = fs.gas.rho;
            // Force per unit volume.
            Q[cqi.xMom] += myConfig.gravity.x * rho;
            Q[cqi.yMom] += myConfig.gravity.y * rho;
            if (cqi.threeD) { Q[cqi.zMom] += myConfig.gravity.z * rho; }
            // Work done per unit volume.
            Q[cqi.totEnergy] += rho * dot(myConfig.gravity, fs.vel);
        }
        // Species production (other than chemistry).
        // For the chemistry and other-internal energy exchange,
        // see thermochemical_increment().
        // Individual energies (other than energy exchange)
        // Radiation can potentially be removed from both the electronic and
        // total energy source terms.
        if (myConfig.radiation) {
            // Radiative source term should be already calculated
            // Add value to total energy
            // FIX-ME: - assuming electronic mode is the last in the vector of energies
            //         - what about Q_renergies[0]?
            Q[cqi.totEnergy] += Q_rE_rad;
            version(multi_T_gas) {
                // Q[cqi.modes+cqi.n_modes-1] += Q_rE_rad; // FIX-ME old C++ code
            }
        }
        return;
    } // end add_inviscid_source_vector()

    @nogc
    void add_viscous_source_vector()
    {
        auto cqi = myConfig.cqi;
        //
        if (myConfig.axisymmetric) {
            // For viscous, axisymmetric flow:
            number v_over_y = fs.vel.y / pos[0].y;
            number dudx=grad.vel[0][0];
            number dvdy=grad.vel[1][1];

            number mu  = fs.gas.mu + fs.mu_t;
            mu *= myConfig.viscous_factor;
            number lmbda = -2.0/3.0 * mu;
            number tau_00 = 2.0 * mu * v_over_y + lmbda * (dudx + dvdy + v_over_y);
            // Y-Momentum; viscous stress contribution from the front and Back interfaces.
            // Note that these quantities are approximated at the
            // mid-point of the cell face and so should never be
            // singular -- at least I hope that this is so.
            Q[cqi.yMom] -= tau_00 * areaxy[0] / volume[0];
        } // end if ( myConfig.axisymmetric )

        version(turbulence) {
            if (in_turbulent_zone) {
                number[] rhoturb = Q[cqi.rhoturb .. cqi.rhoturb+cqi.n_turb];
                myConfig.turb_model.source_terms(fs, *grad, pos[0].y, dwall, L_min, L_max, rhoturb);
            }
        }

        if (myConfig.electric_field_work) {
            // Work done on electrons due to electric field induced by charge separation
            // on scales less than the Debye length
            // FIXME: Only consistent with ambipolar diffusion. Currently this is up to
            //        the user to enforce.
            // Estimate electron pressure gradient as average of all vertices then
            // use approximation for work done on electrons: u dot div(pe)
            // number udivpe, dpedx, dpedy, dpedz;
            // if ( myConfig.dimensions == 2 ) {
            //  mixin(avg_over_vtx_list("grad.pe.x", "dpedx"));
            //  mixin(avg_over_vtx_list("grad.pe.y", "dpedy"));
            //  udivpe = fs.vel.x * dpedx + fs.vel.y * dpedy;
            // } else {
            //  mixin(avg_over_vtx_list("grad.pe.x", "dpedx"));
            //  mixin(avg_over_vtx_list("grad.pe.y", "dpedy"));
            //  mixin(avg_over_vtx_list("grad.pe.z", "dpedz"));
            //  udivpe = fs.vel.x * dpedx + fs.vel.y * dpedy + fs.vel.z * dpedz;
            // }
            // // [TODO] FIXME: Assuming the free electron energy is included in the last mode
            // Q.energies.back() += udivpe * myConfig.diffusion_factor;
        } // end if ( myConfig.electric_field_work )
        return;
    } // end add_viscous_source_vector()


    @nogc
    void add_udf_source_vector()
    {
        Q.add(Qudf);
    }

    @nogc
    void add_thermochemical_source_vector(number[] thermochem_source, double reaction_factor)
    {
        // Bug fixes to this function by NNG on 17/10/22. Previously we skipped
        // calling eval_source_terms if n_species==1, which caused issues with
        // multi-temperature single species calculations. The new code now
        // includes a dummy entry in thermochem_source[0] if nspecies==1, but does
        // not update Q with the results in that case.

        auto cqi = myConfig.cqi;
        if (fs.gas.T <= myConfig.T_frozen) { return; } // TODO: What about T_frozen_energy?
        myConfig.thermochemUpdate.eval_source_terms(myConfig.gmodel, fs.gas, thermochem_source);

        version(multi_species_gas) {
            if (cqi.n_species > 1) {
                foreach(sp; 0 .. cqi.n_species) { Q[cqi.species+sp] += reaction_factor*thermochem_source[sp]; }
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. cqi.n_modes) {
                Q[cqi.modes+imode] += reaction_factor*thermochem_source[$-cqi.n_modes+imode];
            }
        }
    }

    @nogc
    number calculate_wall_Reynolds_number(int which_boundary, GasModel gmodel)
    // [TODO] unstructured-grid adaption to be done, however,
    // this function is not presently used because we have not ported the
    // writing of boundary flow and heat-transfer conditions.
    {
        FVInterface IFace = iface[which_boundary];
        gmodel.update_thermo_from_rhoT(IFace.fs.gas); // Note that we adjust IFace here.
        number a_wall = IFace.fs.gas.a;
        number cell_width = 0.0;
        if ( which_boundary == Face.east || which_boundary == Face.west )
            cell_width = iLength;
        else if ( which_boundary == Face.north || which_boundary == Face.south )
            cell_width = jLength;
        else if ( which_boundary == Face.top || which_boundary == Face.bottom )
            cell_width = kLength;
        number Re_wall = IFace.fs.gas.rho * a_wall * cell_width / IFace.fs.gas.mu;
        return Re_wall;
    } // end calculate_wall_Reynolds_number()

    @nogc
    void store_rad_scaling_params()
    // Store parameters for (re-)scaling of radiative source term.
    // Simple rho x T**4 scaling seems to be adequate.
    {
        // 1. Store the freshly computed radiative flux as the 'original'
        Q_rad_org = Q_rE_rad;
        // 2. Compute the scaling factor based on local gas properties
        // NOTE: - The idea is that f_rad_org is proportional to actual value
        number T = fs.gas.T;
        if ( Q_rad_org <= 0.0 ) {
            // This cell is a net emitter
            f_rad_org = fs.gas.rho * pow(T, 4);
        } else if ( Q_rad_org > 0.0 ) {
            // This cell is a net absorber
            f_rad_org = fs.gas.rho / pow(T, 4);
        }
    } // end store_rad_scaling_params()

    @nogc
    void rescale_Q_rE_rad()
    {
        // 1. Compute the current scaling factor based on local gas properties
        number T = fs.gas.T;
        number f_rad_new = 1.0;
        if ( Q_rad_org <= 0.0 ) {
            // This cell is a net emitter
            f_rad_new = fs.gas.rho * pow(T, 4);
        }
        else if ( Q_rad_org > 0.0 ) {
            // This cell is a net absorber
            f_rad_new = fs.gas.rho / pow(T, 4);
        }
        // 2. (Re-)scale the original source term
        Q_rE_rad = ( f_rad_new / f_rad_org ) * Q_rad_org;
    } // end rescale_Q_rE_rad()

    @nogc
    void reset_Q_rad_to_zero()
    {
        Q_rE_rad = 0.0;
    } // end reset_Q_rad_to_zero()

    @nogc
    number rad_scaling_ratio()
    {
        // 1. Compute the current scaling factor based on local gas properties
        number T = fs.gas.T;
        number f_rad = 1.0;
        if ( Q_rE_rad <= 0.0 ) {
            // This cell is a net emitter
            f_rad = fs.gas.rho * pow(T, 4);
        }
        else if ( Q_rE_rad > 0.0 ) {
            // This cell is a net absorber
            f_rad = fs.gas.rho / pow(T, 4);
        }
        return fabs( f_rad - f_rad_org ) / f_rad_org;
    } // end rad_scaling_ratio()

    @nogc
    void average_vertex_deriv_values()
    {
        grad.copy_values_from(*(vtx[0].grad));
        foreach (i; 1 .. vtx.length) grad.accumulate_values_from(*(vtx[i].grad));
        grad.scale_values_by(to!number(1.0/vtx.length));
    } // end average_vertex_deriv_values()

    @nogc
    void rotate_gradients(const(double[]) Rmatrix)
    {
        // Rotate velocity gradients
        // This is a gradient of a vector, so the result is a tensor. Apply T*A*transpose(T)
        // where T is the rotation matrix and A is the velocity gradient tensor.
        // Some literature reports it as T*A*inv(T), fortunately the rotation matrices have
        // the property that the inverse and the tranpose are equivalent i.e. its an
        // orthogonal matrix

        // Perform the T*A operation
        number old_x, old_y, old_z;
        foreach (i; 0 .. 3) {
            old_x = grad.vel[i][0]; old_y = grad.vel[i][1]; old_z = grad.vel[i][2];

            grad.vel[i][0] = Rmatrix[0] * old_x + Rmatrix[1] * old_y + Rmatrix[2] * old_z;
            grad.vel[i][1] = Rmatrix[3] * old_x + Rmatrix[4] * old_y + Rmatrix[5] * old_z;
            grad.vel[i][2] = Rmatrix[6] * old_x + Rmatrix[7] * old_y + Rmatrix[8] * old_z;
        }

        // Now apply the *inv(T)
        foreach (i; 0 .. 3) {
            old_x = grad.vel[0][i]; old_y = grad.vel[1][i]; old_z = grad.vel[2][i];

            grad.vel[0][i] = Rmatrix[0] * old_x + Rmatrix[1] * old_y + Rmatrix[2] * old_z;
            grad.vel[1][i] = Rmatrix[3] * old_x + Rmatrix[4] * old_y + Rmatrix[5] * old_z;
            grad.vel[2][i] = Rmatrix[6] * old_x + Rmatrix[7] * old_y + Rmatrix[8] * old_z;
        }

        // Rotate temperature gradients- just a vector rotation
        old_x = grad.T[0]; old_y = grad.T[1]; old_z = grad.T[2];

        grad.T[0] = Rmatrix[0] * old_x + Rmatrix[1] * old_y + Rmatrix[2] * old_z;
        grad.T[1] = Rmatrix[3] * old_x + Rmatrix[4] * old_y + Rmatrix[5] * old_z;
        grad.T[2] = Rmatrix[6] * old_x + Rmatrix[7] * old_y + Rmatrix[8] * old_z;

        version(multi_species_gas) {
            foreach (i; 0 .. myConfig.n_species) {
                old_x = grad.massf[i][0]; old_y = grad.massf[i][1]; old_z = grad.massf[i][2];

                grad.massf[i][0] = Rmatrix[0] * old_x + Rmatrix[1] * old_y + Rmatrix[2] * old_z;
                grad.massf[i][1] = Rmatrix[3] * old_x + Rmatrix[4] * old_y + Rmatrix[5] * old_z;
                grad.massf[i][2] = Rmatrix[6] * old_x + Rmatrix[7] * old_y + Rmatrix[8] * old_z;
            }
        }

        version(multi_T_gas) {
            foreach (i; 0 .. myConfig.n_modes) {
                old_x = grad.T_modes[i][0]; old_y = grad.T_modes[i][1]; old_z = grad.T_modes[i][2];

                grad.T_modes[i][0] = Rmatrix[0] * old_x + Rmatrix[1] * old_y + Rmatrix[2] * old_z;
                grad.T_modes[i][1] = Rmatrix[3] * old_x + Rmatrix[4] * old_y + Rmatrix[5] * old_z;
                grad.T_modes[i][2] = Rmatrix[6] * old_x + Rmatrix[7] * old_y + Rmatrix[8] * old_z;
            }
        }

        version(turbulence) {
            foreach (i; 0 .. myConfig.turb_model.nturb) {
                old_x = grad.turb[i][0]; old_y = grad.turb[i][1]; old_z = grad.turb[i][2];

                grad.turb[i][0] = Rmatrix[0] * old_x + Rmatrix[1] * old_y + Rmatrix[2] * old_z;
                grad.turb[i][1] = Rmatrix[3] * old_x + Rmatrix[4] * old_y + Rmatrix[5] * old_z;
                grad.turb[i][2] = Rmatrix[6] * old_x + Rmatrix[7] * old_y + Rmatrix[8] * old_z;
            }
        }
    }
    // End rotate_gradients

    @nogc
    void average_interface_deriv_values()
    {
        grad.copy_values_from(*(iface[0].grad));
        foreach (i; 1 .. iface.length) grad.accumulate_values_from(*(iface[i].grad));
        grad.scale_values_by(to!number(1.0/iface.length));
    } // end average_interface_deriv_values()

    // Think this should be fine as nogc? Taking transform of pressure in this example
    @nogc
    void increment_local_DFT(size_t DFT_step) {
        // If it's the first step, we should set the values rather than incrementing
        if (DFT_step == 0) {
            foreach (i; 0..myConfig.DFT_n_modes) {
                DFT_local_real[i] = cos(2 * std.math.PI * i * DFT_step / myConfig.DFT_n_modes) * fs.gas.p;
                DFT_local_imag[i] = sin(2 * std.math.PI * i * DFT_step / myConfig.DFT_n_modes) * fs.gas.p;
            }
        } else {
            foreach (i; 0..myConfig.DFT_n_modes) {
                DFT_local_real[i] += cos(2 * std.math.PI * i * DFT_step / myConfig.DFT_n_modes) * fs.gas.p;
                DFT_local_imag[i] -= sin(2 * std.math.PI * i * DFT_step / myConfig.DFT_n_modes) * fs.gas.p;
            }
        }
    }

    void gather_residual_stencil_lists(int spatial_order_of_jacobian)
    {
        /*
          This function gathers references to the interfaces and cells
          that make up the residual stencil for a cell needed for the
          flow Jacobian construction. These stencils can be thought of
          in terms of what neighbouring cells will have perturbed residuals
          in the event this cells flow state or conserved quantities are
          perturbed.

          Note that we need the cells to be in numerical order
          according to their local id for entry into the flow
          Jacobian later

          TODO: extend to handle structured grid solver
         */

        FVCell[] unordered_cell_list;  // TODO: think about possibly pre-sizing this array
        size_t[size_t] cell_pos_array; // this is used to retrieve a cell from the unordered list

        bool include_viscous_effects = myConfig.viscous;
        // when using the flow Jacobian as a precondition matrix for methods such as GMRES it has
        // been observed that only filling entries of the nearest-neighbours provides more
        // robust preconditioning even if viscous effects would suggest a larger stencil is
        // required. Note that we will still apply the viscous effects later when forming
        // the flow Jacobian, we are in effect just dropping some of the Jacobian entries
        // by reducing the stencil footprint.

        // add this cell
        size_t[] cell_ids;
        unordered_cell_list ~= cell_cloud[0];
        cell_pos_array[cell_cloud[0].id] = unordered_cell_list.length-1;
        cell_ids ~= cell_cloud[0].id;
        bool nearest_neighbours = false;
        if ( spatial_order_of_jacobian >= 0) { nearest_neighbours = true; }

        if (nearest_neighbours) {
            // this first order stencil adds the nearest neighbours
            // Note that if the use_extended_stencil is set to true, then there will be more than just
            // nearest neighbour cells in the cell cloud, so care must be taken to gather the correct cells
            size_t np = iface.length + 1;
            foreach (i; 1 .. np) {
                FVCell cell = cell_cloud[i];
                bool cell_exists = cell_ids.canFind(cell.id);
                if (!cell_exists && cell.id < 1_000_000_000 && is_interior_to_domain) {
                    unordered_cell_list ~= cell;
                    cell_pos_array[cell.id] = unordered_cell_list.length-1;
                    cell_ids ~= cell.id;
                }
            } // finished gathering cells

        }

        bool extended_neighbours = false;
        if ( nearest_neighbours &&
             ( (spatial_order_of_jacobian >= 2) ||
               (spatial_order_of_jacobian >= 1 && include_viscous_effects) ) ) { extended_neighbours = true; }

        if (extended_neighbours) {
            // second order (&/or viscous, or first order with viscous effects)
            // stencil adds the nearest neighbours of the nearest neighbours

            // gather additional cells
            foreach (icell; 1 .. cell_cloud.length) {
                foreach (cell; cell_cloud[icell].cell_cloud) {
                    bool cell_exists = cell_ids.canFind(cell.id);
                    if (!cell_exists && cell.id < 1_000_000_000 && is_interior_to_domain) {
                        unordered_cell_list ~= cell;
                        cell_pos_array[cell.id] = unordered_cell_list.length-1;
                        cell_ids ~= cell.id;
                    }
                }
            } // finished gathering cells
        }

        // now sort the cells
        cell_ids.sort();
        foreach (id; cell_ids) { cell_list ~= unordered_cell_list[cell_pos_array[id]]; }

        // gather the interfaces of those cells
        size_t[] face_ids;
        foreach (cell; cell_list) {
            foreach (face; cell.iface) {
                bool face_exists = face_ids.canFind(face.id);
                if (!face_exists) {
                    face_list ~= face;
                    face_ids ~= face.id;
                }
            }
        } // finished gathering faces

    } // end gather_residual_stencil_lists()

    void gather_residual_stencil_lists_for_ghost_cells(int spatial_order_of_jacobian, FVCell[] neighbour_cell_cloud)
    {
        /*
          This function gathers references to the interfaces and cells
          that make up the residual stencil for a ghost cell along domain
          boundaries (e.g. all boundary conditions except FullFaceCopy and
          MappedCellCopy BCs). These stencils can be thought of
          in terms of what interior neighbouring cells will have perturbed residuals
          in the event this ghost cells flow state or conserved quantities are
          perturbed.

          Note that we DO NOT need the cells to be in numerical order
          for this stencil

          TODO: extend to handle structured grid solver
         */

        bool include_viscous_effects = myConfig.viscous;
        bool nearest_neighbours_only = false;
        if ( (spatial_order_of_jacobian == 1 && include_viscous_effects == false) ||
             spatial_order_of_jacobian == 0) { nearest_neighbours_only = true; }

        // this first order stencil includes the ghost cells nearest neighbours
        cell_list ~= neighbour_cell_cloud[0]; // this is the interior cell that shares an interface

        bool extended_neighbours = false;
        if ( (spatial_order_of_jacobian >= 2) ||
             (spatial_order_of_jacobian >= 1 && include_viscous_effects) ) { extended_neighbours = true; }

        if (extended_neighbours) {
            // second order (&/or viscous, or first order with viscous effects) stencil includes the cell
            // and its nearest neighbours as well as the nearest neighbours of the nearest neighbours

            // gather cells
            foreach (c; neighbour_cell_cloud) {
                if ( c.id != id && c.id != neighbour_cell_cloud[0].id && c.id < 1_000_000_000 && c.is_interior_to_domain) {
                    cell_list ~= c;
                }
            } // finished gathering cells
        }

        // gather faces
        size_t[] face_ids;
        foreach (c; cell_list) {
            foreach (face; c.iface) {
                bool face_exists = face_ids.canFind(face.id);
                if (!face_exists) {
                    face_list ~= face;
                    face_ids ~= face.id;
                }
            }
        } // finished gathering faces

    } // end gather_residual_stencil_lists_for_ghost_cells()

} // end class FVCell
