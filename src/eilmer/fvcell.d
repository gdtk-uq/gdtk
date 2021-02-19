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
import fvcore;
import flowstate;
import flowgradients;
import conservedquantities;
import fvvertex;
import fvinterface;
import globalconfig;
import lsqinterp;
import gas.fuel_air_mix;
import simcore : SimState;
import turbulence;

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
    //
    bool fr_reactions_allowed; // if true, will call thermochemical_increment
    double dt_chem; // acceptable time step for finite-rate chemistry
    double dt_therm; // acceptable time step for thermal relaxation
    bool in_turbulent_zone; // if true, we will keep the turbulence viscosity
    number base_qdot; // base-level of heat addition to cell, W/m**3
    // Geometry
    Vector3[] pos; // Centre x,y,z-coordinates for time-levels, m,m,m
    number iLength; // length in the i-index direction
    number jLength; // length in the j-index direction
    number kLength; // length in the k-index direction
    number L_min;   // minimum length scale for cell
    number L_max;   // maximum length scale for cell
    // Connections
    FVInterface[] iface;  // references to defining interfaces of cell
    double[] outsign; // +1.0 if iface is outward-facing; -1.0 for an inward-facing iface
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
    ConservedQuantities[2] dUdt_copy; // for residual smoothing
    // for unstructured grids, we may be doing high-order reconstruction
    LSQInterpWorkspace ws;
    LSQInterpGradients gradients; // we only need these workspaces for the unstructured
                                  // solver, they are instantiated in ufluidblock.d
    // Viscous-flux-related quantities.
    FlowGradients grad;
    WLSQGradWorkspace ws_grad;
    Vector3*[] cloud_pos; // Positions of flow points for gradients calculation.
    FlowState[] cloud_fs; // References to flow states at those points.
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

    // source terms for finite-rate chemistr
    number[] chem_conc;
    number[] chem_rates;
    number[] chem_source;
    ReactionMechanism rmech;
    
    // For use with LU-SGS solver/preconditioner (note: we don't need complex numbers here)
    number[] LU;
    number[] dUk;
    number[] dF;
    number[] scalar_diag_inv;
    Matrix!number dFdU;
    Matrix!number dFdU_rotated;

    // Shape sensitivity calculator workspace
    FVCell[] cell_list;            // list of cells in the residual stencil
    FVInterface[] face_list;       // list of faces in the residual stencil
    version(shape_sensitivity) {
       	size_t[] pcell_global_coord_list;
	size_t[][] ecell_global_coord_list;
	number[][] entry_list;

	size_t global_id;
	number[][] dqdQ;
        number[][] dQdU;
        // stencil of effected cells & faces used in forming the flow Jacobian
        FVCell[] jacobian_cell_stencil;
        FVInterface[] jacobian_face_stencil;
        // arrays used to temporarily store data intended for the neighbouring block
        // during construction of the external portion of the flow Jacobian.
        size_t[] idList;
        number[] aa;
        // block-diagonal contribution to Jacobian used in steady-state solver pre-conditioner
        Matrix!number dPrimitive;
        Matrix!number dConservative;
        int[] pivot;
    }

private:
    LocalConfig myConfig;

public:
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
        fs = new FlowState(gmodel, 100.0e3, T, T_modes, Vector3(0.0,0.0,0.0));
        foreach(i; 0 .. myConfig.n_flow_time_levels) {
            U ~= new ConservedQuantities(n_species, n_modes);
            U[i].clear();
            dUdt ~= new ConservedQuantities(n_species, n_modes);
        }
        Q = new ConservedQuantities(n_species, n_modes);
        Q.clear();
        if (myConfig.residual_smoothing) {
            dUdt_copy[0] = new ConservedQuantities(n_species, n_modes);
            dUdt_copy[1] = new ConservedQuantities(n_species, n_modes);
        }
        grad = new FlowGradients(myConfig);
        if (allocate_spatial_deriv_lsq_workspace) {
            ws_grad = new WLSQGradWorkspace();
        }
        version(shape_sensitivity) {
            dQdU.length = GlobalConfig.cqi.nConservedQuantities; // number of conserved variables
            dqdQ.length = GlobalConfig.cqi.nConservedQuantities;
	    foreach (ref a; dQdU) a.length = GlobalConfig.cqi.nConservedQuantities;
	    foreach (ref a; dqdQ) a.length = GlobalConfig.cqi.nConservedQuantities;
            foreach (i; 0..dQdU.length) {
                foreach (j; 0..dQdU[i].length) {
                    dQdU[i][j] = 0.0;
		    dqdQ[i][j] = 0.0;
                }
            }
        }
        version(debug_chem) {
            // The savedGasState is a module-level variable.
            // It only needs to be initialised when debug_chem mode
            // is on AND it only required initialisation once.
            savedGasState = new GasState(gmodel);
        }

        // some data structures used in the LU-SGS solver
        version(steady_state) {

            if (myConfig.reacting) {
                chem_source.length = n_species;
                chem_conc.length = n_species;
                chem_rates.length = n_species;
                
                auto myChemUpdate = cast(ChemistryUpdate) myConfig.thermochemUpdate;
                if (myChemUpdate !is null) { 
                    rmech = myChemUpdate.rmech.dup();
                } else {
                    throw new Exception("Opps, incorrect ThermochemicalReactor.");
                }
            }
            
            size_t nConserved = myConfig.cqi.nConservedQuantities;
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

    }

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
            grad.copy_values_from(other.grad);
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
        repr ~= ", outsigns=["; foreach (osgn; outsign) { repr ~= format("%.18e,", osgn); } repr ~= "]";
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
        number vol, xyplane_area;
        switch (vtx.length) {
        case 3:
            xyplane_triangle_cell_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl],
                                             pos[gtl], xyplane_area, iLength, jLength, L_min);
            break;
        case 4:
            xyplane_quad_cell_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                                         vtx[2].pos[gtl], vtx[3].pos[gtl],
                                         pos[gtl], xyplane_area, iLength, jLength, L_min);
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
        kLength = 0.0;
        L_max = fmax(iLength, jLength);
    } // end update_2D_geometric_data()

    @nogc
    void update_3D_geometric_data(size_t gtl)
    {
        string msg = "FVCell.update_3D_geometric_data(): ";
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
                                pos[gtl], volume[gtl], iLength, jLength, kLength);
            L_min = min(iLength, jLength, kLength);
            break;
        case 5:
            pyramid_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl], vtx[3].pos[gtl],
                               vtx[4].pos[gtl], pos[gtl], volume[gtl], L_min);
            iLength = L_min; jLength = L_min; kLength = L_min;
            break;
        case 6:
            wedge_properties(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[2].pos[gtl],
                             vtx[3].pos[gtl], vtx[4].pos[gtl], vtx[5].pos[gtl],
                             pos[gtl], volume[gtl], L_min);
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
        FlowState[] fsList;
        // We need to be honest and not to fiddle with the other gas states.
        foreach(other; others) {
            if ( this is other ) {
                throw new FlowSolverException("Must not include destination in source list.");
            }
            fsList ~= cast(FlowState)other.fs;
        }
        fs.copy_average_values_from(fsList, gmodel);
        // Accumulate from a clean slate and then divide.
        Q_rE_rad = 0.0;
        foreach(other; others) {
            Q_rE_rad += other.Q_rE_rad;
        }
        Q_rE_rad /= n;
    } // end replace_flow_data_with_average()

    void scan_values_from_string(string buffer, ref string[] varNameList, bool fixedOrder,
                                 GasModel gmodel, bool overwrite_geometry_data)
    // Note that the position data is read into grid_time_level 0.
    {
        Vector3 new_pos;
        number new_volume;
        if (fixedOrder) {
            scan_cell_data_from_fixed_order_string
                (buffer,
                 new_pos, new_volume, fs,
                 Q_rad_org, f_rad_org, Q_rE_rad,
                 myConfig.with_local_time_stepping, dt_local, dt_chem, dt_therm,
                 myConfig.include_quality, myConfig.MHD,
                 myConfig.divergence_cleaning, myConfig.radiation,
                 myConfig.turb_model.nturb);
        } else {
            scan_cell_data_from_variable_order_string
                (buffer, varNameList, gmodel, myConfig.turb_model,
                 new_pos, new_volume, fs,
                 Q_rad_org, f_rad_org, Q_rE_rad,
                 myConfig.with_local_time_stepping, dt_local, dt_chem, dt_therm,
                 myConfig.include_quality, myConfig.MHD,
                 myConfig.divergence_cleaning, myConfig.radiation);
        }
        if (overwrite_geometry_data) {
            pos[0].set(new_pos);
            volume[0] = new_volume;
        }
    } // end scan_values_from_string()

    void read_values_from_raw_binary(ref File fin, bool overwrite_geometry_data)
    // Note that the position data is read into grid_time_level 0.
    {
        Vector3 new_pos;
        number new_volume;
        raw_binary_to_cell_data(fin, new_pos, new_volume, fs,
                                Q_rad_org, f_rad_org, Q_rE_rad,
                                myConfig.with_local_time_stepping, dt_local, dt_chem, dt_therm,
                                myConfig.include_quality, myConfig.MHD,
                                myConfig.divergence_cleaning, myConfig.radiation,
                                myConfig.turb_model.nturb);
        if (overwrite_geometry_data) {
            pos[0].set(new_pos);
            volume[0] = new_volume;
        }
    } // end read_values_from_raw_binary()

    string write_values_to_string() const
    {
        return cell_data_as_string(pos[0], volume[0], fs,
                                   Q_rad_org, f_rad_org, Q_rE_rad,
                                   myConfig.with_local_time_stepping, dt_local, dt_chem, dt_therm,
                                   myConfig.include_quality, myConfig.MHD,
                                   myConfig.divergence_cleaning, myConfig.radiation,
                                   myConfig.turb_model.nturb);
    } // end write_values_to_string()

    void write_values_to_raw_binary(ref File fout) const
    {
        cell_data_to_raw_binary(fout, pos[0], volume[0], fs,
                                Q_rad_org, f_rad_org, Q_rE_rad,
                                myConfig.with_local_time_stepping, dt_local, dt_chem, dt_therm,
                                myConfig.include_quality, myConfig.MHD,
                                myConfig.divergence_cleaning, myConfig.radiation,
                                myConfig.turb_model.nturb);
    } // end write_values_to_raw_binary()

    string write_residuals_to_string() const
    {
        auto writer = appender!string();
        version(complex_numbers) {
            formattedWrite(writer, "%.18e %.18e %.18e %.18e",
                           -dUdt[0].mass.re, -dUdt[0].momentum.x.re,
                           -dUdt[0].momentum.y.re, -dUdt[0].total_energy.re);
        } else {
            formattedWrite(writer, "%.18e %.18e %.18e %.18e",
                           -dUdt[0].mass, -dUdt[0].momentum.x,
                           -dUdt[0].momentum.y, -dUdt[0].total_energy);
        }
        return writer.data;
    }

    @nogc
    void encode_conserved(int gtl, int ftl, double omegaz)
    // gtl = grid time level
    // ftl = flow time level
    {
        ConservedQuantities myU = U[ftl];
        number myrho = fs.gas.rho;
        // Mass per unit volume.
        myU.mass = myrho;
        // Momentum per unit volume.
        myU.momentum.set(fs.gas.rho*fs.vel.x, fs.gas.rho*fs.vel.y, fs.gas.rho*fs.vel.z);
        version(MHD) {
            // Magnetic field
            myU.B.set(fs.B);
            myU.psi = fs.psi;
            myU.divB = fs.divB;
        }
        // Total Energy / unit volume
        number u = myConfig.gmodel.internal_energy(fs.gas);
        number ke = 0.5*(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y+fs.vel.z*fs.vel.z);
        myU.total_energy = fs.gas.rho*(u + ke);
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                myU.rhoturb[i] = fs.gas.rho * fs.turb[i];
            }
            myU.total_energy += fs.gas.rho * myConfig.turb_model.turbulent_kinetic_energy(fs);
        }
        version(MHD) {
            if (myConfig.MHD) {
                number me = 0.5*(fs.B.x*fs.B.x + fs.B.y*fs.B.y + fs.B.z*fs.B.z);
                myU.total_energy += me;
            }
        }
        version(multi_T_gas) {
            // Other internal energies: energy in mode per unit volume.
            foreach(imode; 0 .. myU.energies.length) {
                myU.energies[imode] = fs.gas.rho*fs.gas.u_modes[imode];
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
            myU.total_energy -= rho*0.5*omegaz*omegaz*rsq;
        }
        version(multi_species_gas) {
            // Species densities: mass of species is per unit volume.
            foreach(isp; 0 .. myConfig.n_species) {
                myU.massf[isp] = fs.gas.rho*fs.gas.massf[isp];
            }
        }
        assert(U[ftl].mass > 0.0, "invalid density in conserved quantities vector" ~
               " at end of FVCell.encode_conserved().");
        return;
    } // end encode_conserved()

    @nogc
    int decode_conserved(int gtl, int ftl, double omegaz)
    {
        auto gmodel = myConfig.gmodel;
        ConservedQuantities myU = U[ftl];
        // The conserved quantities are carried as quantity per unit volume.
        // mass / unit volume = density
        if (!(myU.mass > 0.0)) {
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
        number rho = myU.mass;
        fs.gas.rho = rho; // This is limited to nonnegative and finite values.
        number dinv = 1.0 / rho;
        // Velocities from momenta.
        fs.vel.set(myU.momentum.x*dinv, myU.momentum.y*dinv, myU.momentum.z*dinv);
        version(MHD) {
            // Magnetic field.
            fs.B.set(myU.B);
            fs.psi = myU.psi;
            fs.divB = myU.divB;
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
            rE = myU.total_energy + rho*0.5*omegaz*omegaz*rsq;
        } else {
            // Non-rotating frame.
            rE = myU.total_energy;
        }
        version(MHD) {
            number me = 0.0;
            if ( myConfig.MHD ) { me = 0.5*(fs.B.x*fs.B.x + fs.B.y*fs.B.y + fs.B.z*fs.B.z); }
            rE -= me;
        }
        // Start with the total energy, then take out the other components.
        // Internal energy is what remains.
        number u = rE * dinv;
        version(turbulence) {
            if (allow_k_omega_update) {
                foreach(i; 0 .. myConfig.turb_model.nturb) {
                    // for stability, we enforce tke and omega to be positive.
                    // This approach is referred to as clipping in Chisholm's (2007) thesis:
                    // A fully coupled Newton-Krylov solver with a one-equation turbulence model.
                    // to prevent division by 0.0 set variables to a very small positive value.
                    fs.turb[i] = myU.rhoturb[i] * dinv;
                    if (fs.turb[i] < 0.0) fs.turb[i] = 1.0e-10;
                }
            }
            u -= myConfig.turb_model.turbulent_kinetic_energy(fs);
        }
        // Remove kinetic energy for bulk flow.
        number ke = 0.5*(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y + fs.vel.z*fs.vel.z);
        u -= ke;
        // Other energies, if any.
        version(multi_T_gas) {
            number u_other = 0.0;
            foreach(imode; 0 .. gmodel.n_modes) { fs.gas.u_modes[imode] = myU.energies[imode] * dinv; }
            foreach(ei; fs.gas.u_modes) { u_other += ei; }
            fs.gas.u = u - u_other;
        } else {
            fs.gas.u = u;
        }
        // Thermochemical species, if appropriate.
        version(multi_species_gas) {
            try {
		foreach(isp; 0 .. myConfig.n_species) { fs.gas.massf[isp] = myU.massf[isp] * dinv; }
		if (myConfig.sticky_electrons) { gmodel.balance_charge(fs.gas); }
		if (myConfig.n_species > 1) { scale_mass_fractions(fs.gas.massf); }
	    } catch (GasModelException err) {
		if (myConfig.adjust_invalid_cell_data) {
		    data_is_bad = true;
		    return -2;
		} else {
		    string msg = "Bad cell with mass fractions that do not add correctly.";
		    debug {
			msg ~= format("scale_mass_fractions exception with message:\n  %s", err.msg);
			msg ~= format("The decode_conserved() failed for cell: %d\n", id);
			msg ~= format("This cell is located at: %s\n", pos[0]);
			msg ~= format("This cell is located in block: %d\n", myConfig.universe_blk_id);
			msg ~= format("The gas state before thermo update is:\n   fs.gas %s", fs.gas);
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
            gmodel.update_sound_speed(fs.gas);
            if (myConfig.viscous) gmodel.update_trans_coeffs(fs.gas);
        } catch (GasModelException err) {
            if (myConfig.adjust_invalid_cell_data) {
                data_is_bad = true;
                return -2;
            } else {
                string msg = "Bad cell with failed thermodynamic update.";
                debug {
                    msg ~= format("Thermodynamic update exception with message:\n  %s", err.msg);
                    msg ~= format("The decode_conserved() failed for cell: %d\n", id);
                    msg ~= format("This cell is located at: %s\n", pos[0]);
                    msg ~= format("This cell is located in block: %d\n", myConfig.universe_blk_id);
                    msg ~= format("The gas state after the failed update is:\n   fs.gas %s", fs.gas);
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
            if (myConfig.viscous) gmodel.update_trans_coeffs(fs.gas);
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
        //
        // Mass.
        number integral_mass = 0.0;
        // Momentum.
        number integral_momx = 0.0;
        number integral_momy = 0.0;
        number integral_momz = 0.0;
        version(MHD) {
            // Magnetic Field.
            number integral_Bx = 0.0;
            number integral_By = 0.0;
            number integral_Bz = 0.0;
            // Divergence of the magnetic field; it is not actually a time-derivative
            // but seems to be the best way to calculate it. (Lachlan W.)
            my_dUdt.divB = 0.0;
            number integral_psi = 0.0;
        }
        // Total Energy.
        number integral_E = 0.0;
        version(turbulence) {
            size_t nturb = myConfig.turb_model.nturb;
            number[2] integral_rhoturb;
            foreach(ref t; integral_rhoturb) t = 0.0;
        }
        version(multi_species_gas) {
            // Time-derivative for individual species.
            // The conserved quantity is the mass per unit
            // volume of species isp and
            // the fluxes are mass/unit-time/unit-area.
            // Units of DmassfDt are 1/sec.
            immutable uint max_species = 32;
            number[max_species] integral_species;
            if (myConfig.n_species > max_species) { throw new Error("oops too many chemical species for work array"); }
            foreach(j; 0 .. myConfig.n_species) { integral_species[j] = 0.0; }
        }
        version(multi_T_gas) {
            // Individual energies.
            immutable uint max_modes = 5;
            number[max_modes] integral_modes;
            uint nmodes = myConfig.n_modes;
            if (nmodes > max_modes) { throw new Error("oops too many energy modes for work array"); }
            foreach(j; 0 .. nmodes) { integral_modes[j] = 0.0; }
        }
        //
        // Integrate the fluxes across the interfaces that bound the cell.
        foreach(i; 0 .. iface.length) {
            ConservedQuantities* myF = &(iface[i].F);
            number area = outsign[i]*iface[i].area[gtl];
            //
            integral_mass -= myF.mass*area;
            integral_momx -= myF.momentum.x*area;
            integral_momy -= myF.momentum.y*area;
            if ((myConfig.dimensions == 3) || ( myConfig.MHD )) {
                // require z-momentum for MHD even in 2D
                integral_momz -= myF.momentum.z*area;
            }
            version(MHD) {
                if (myConfig.MHD) {
                    integral_Bx -= myF.B.x*area;
                    integral_By -= myF.B.y*area;
                    integral_Bz -= myF.B.z*area;
                    my_dUdt.divB += myF.divB*area;
                    if (myConfig.divergence_cleaning) {
                        integral_psi -= myF.psi*area;
                    }
                }
            }
            integral_E -= myF.total_energy*area;
            version(turbulence) {
                foreach(j; 0 .. nturb) {
                    integral_rhoturb[j] -= myF.rhoturb[j]*area;
                }
            }
            version(multi_species_gas) {
                foreach(j; 0 .. myConfig.n_species) { integral_species[j] -= myF.massf[j]*area; }
            }
            version(multi_T_gas) {
                foreach(j; 0 .. nmodes) { integral_modes[j] -= myF.energies[j]*area; }
            }
        } // end foreach iface
        //
        // Finally, evaluate the derivatives of conserved quantities.
        // These are quantity-per-unit-volume.
        number vol_inv = 1.0 / volume[gtl]; // Cell volume (inverted).
        my_dUdt.mass = vol_inv*integral_mass + Q.mass;
        my_dUdt.momentum.set(vol_inv*integral_momx + Q.momentum.x,
                             vol_inv*integral_momy + Q.momentum.y,
                             vol_inv*integral_momz + Q.momentum.z);
        version(MHD) {
            if (myConfig.MHD) {
                if (myConfig.MHD_static_field) {
                    // then the magnetic field won't change...
                    my_dUdt.B.set(Q.B.x, Q.B.y, Q.B.z);
                } else {
                    my_dUdt.B.set(vol_inv*integral_Bx + Q.B.x,
                                  vol_inv*integral_By + Q.B.y,
                                  vol_inv*integral_Bz + Q.B.z);
                }
                if (myConfig.divergence_cleaning) {
                    my_dUdt.psi = vol_inv*integral_psi + Q.psi;
                }
            } else {
                my_dUdt.B.clear();
                my_dUdt.psi = 0.0;
                my_dUdt.divB = 0.0;
            }
        }
        my_dUdt.total_energy = vol_inv*integral_E + Q.total_energy;
        version(turbulence) {
            foreach(i; 0 .. nturb) {
                my_dUdt.rhoturb[i] = vol_inv*integral_rhoturb[i] + Q.rhoturb[i];
            }
        }
        version(multi_species_gas) {
            foreach(j; 0 .. myConfig.n_species) {
                my_dUdt.massf[j] = vol_inv*integral_species[j] + Q.massf[j];
            }
        }
        version(multi_T_gas) {
            foreach(j; 0 .. nmodes) {
                my_dUdt.energies[j] = vol_inv*integral_modes[j] + Q.energies[j];
            }
        }
    } // end time_derivatives()


    void rkl1_stage_update_for_flow_on_fixed_grid1(double dt, int j, int s, bool with_local_time_stepping) 
    {
        ConservedQuantities dUdt0;
        ConservedQuantities U0;
        ConservedQuantities U1;
        ConservedQuantities U2;
        U0 = U[0];
        U1 = U[1];
        dUdt0 = dUdt[0];

        // coefficients
        double muj; double vuj; double muj_tilde;
        muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);
        muj = 1.0;
        vuj = 0.0;

        U1.mass = U0.mass + muj_tilde*dt*dUdt0.mass;
        U1.momentum.set(U0.momentum.x + muj_tilde*dt*dUdt0.momentum.x,
                        U0.momentum.y + muj_tilde*dt*dUdt0.momentum.y,
                        U0.momentum.z + muj_tilde*dt*dUdt0.momentum.z);
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U1.B.set(U0.B.x + muj_tilde*dt*dUdt0.B.x,
                         U0.B.y + muj_tilde*dt*dUdt0.B.y,
                         U0.B.z + muj_tilde*dt*dUdt0.B.z);
                if (myConfig.divergence_cleaning) {
                    U1.psi = U0.psi + muj_tilde*dt*dUdt0.psi;
                    U1.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U1.B.clear();
                U1.psi = 0.0;
            }
        }
        U1.total_energy = U0.total_energy + muj_tilde*dt*dUdt0.total_energy;
        version(turbulence) {
        foreach(i; 0 .. myConfig.turb_model.nturb){
                U1.rhoturb[i] = U0.rhoturb[i] + muj_tilde*dt*dUdt0.rhoturb[i];
                U1.rhoturb[i] = fmax(U1.rhoturb[i],  U0.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. U0.massf.length) {
                U1.massf[isp] = U0.massf[isp] + muj_tilde*dt*dUdt0.massf[isp];
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U0.energies.length) {
                U1.energies[imode] = U0.energies[imode] + muj_tilde*dt*dUdt0.energies[imode];
            }
        }
        // shuffle time-levels
        U[0] = U0;
        U[1] = U1;
        return;
    } // end rkl1_stage_update_for_flow_on_fixed_grid1()

    void rkl1_stage_update_for_flow_on_fixed_grid2(double dt, int j, int s, bool with_local_time_stepping) 
    {
        ConservedQuantities dUdt0;
        ConservedQuantities U0;
        ConservedQuantities U1;
        ConservedQuantities U2;
        U0 = U[0];
        U1 = U[1];
        U2 = U[2];
        dUdt0 = dUdt[1];

        // coefficients
        double muj; double vuj; double muj_tilde;
        muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);
        muj = (2.0*j-1)/j;
        vuj = (1.0-j)/j;

        U2.mass = muj*U1.mass + vuj*U0.mass + muj_tilde*dt*dUdt0.mass;
        U2.momentum.set(muj*U1.momentum.x + vuj*U0.momentum.x + muj_tilde*dt*dUdt0.momentum.x,
                        muj*U1.momentum.y + vuj*U0.momentum.y + muj_tilde*dt*dUdt0.momentum.y,
                        muj*U1.momentum.z + vuj*U0.momentum.z + muj_tilde*dt*dUdt0.momentum.z);
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U2.B.set(muj*U1.B.x + vuj*U0.B.x + muj_tilde*dt*dUdt0.B.x,
                         muj*U1.B.y + vuj*U0.B.y + muj_tilde*dt*dUdt0.B.y,
                         muj*U1.B.z + vuj*U0.B.z + muj_tilde*dt*dUdt0.B.z);
                if (myConfig.divergence_cleaning) {
                    U2.psi = muj*U1.psi + vuj*U0.psi + muj_tilde*dt*dUdt0.psi;
                    U2.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U2.B.clear();
                U2.psi = 0.0;
            }
        }
        U2.total_energy = muj*U1.total_energy + vuj*U0.total_energy + muj_tilde*dt*dUdt0.total_energy;
        version(turbulence) {
        foreach(i; 0 .. myConfig.turb_model.nturb){
                U2.rhoturb[i] = muj*U1.rhoturb[i] + vuj*U0.rhoturb[i] + muj_tilde*dt*dUdt0.rhoturb[i];
                U2.rhoturb[i] = fmax(U2.rhoturb[i],  U0.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. U2.massf.length) {
                U2.massf[isp] = muj*U1.massf[isp] + vuj*U0.massf[isp] + muj_tilde*dt*dUdt0.massf[isp];
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U2.energies.length) {
                U2.energies[imode] = muj*U1.energies[imode] + vuj*U0.energies[imode] + muj_tilde*dt*dUdt0.energies[imode];
            }
        }
        // shuffle time-levels
        U[0] = U0;
        U[1] = U1;
        U[2] = U2;
        return;
    } // end rkl1_stage_update_for_flow_on_fixed_grid2()

    void rkl2_stage_update_for_flow_on_fixed_grid1(double dt, int j, int s, bool with_local_time_stepping)
    {
        ConservedQuantities dUdt0;
        ConservedQuantities U0;
        ConservedQuantities U1;
        ConservedQuantities U2;
        U0 = U[0];
        U1 = U[1];
        dUdt0 = dUdt[0];

        // coefficients
        double muj; double vuj; double muj_tilde;
        muj_tilde = 4.0/(3.0*(s*s+s-2.0));

        U1.mass = U0.mass + muj_tilde*dt*dUdt0.mass;
        U1.momentum.set(U0.momentum.x + muj_tilde*dt*dUdt0.momentum.x,
                        U0.momentum.y + muj_tilde*dt*dUdt0.momentum.y,
                        U0.momentum.z + muj_tilde*dt*dUdt0.momentum.z);
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U1.B.set(U0.B.x + muj_tilde*dt*dUdt0.B.x,
                         U0.B.y + muj_tilde*dt*dUdt0.B.y,
                         U0.B.z + muj_tilde*dt*dUdt0.B.z);
                if (myConfig.divergence_cleaning) {
                    U1.psi = U0.psi + muj_tilde*dt*dUdt0.psi;
                    U1.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U1.B.clear();
                U1.psi = 0.0;
            }
        }
        U1.total_energy = U0.total_energy + muj_tilde*dt*dUdt0.total_energy;
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U1.rhoturb[i] = U0.rhoturb[i] + muj_tilde*dt*dUdt0.rhoturb[i];
                U1.rhoturb[i] = fmax(U1.rhoturb[i], U0.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. U0.massf.length) {
                U1.massf[isp] = U0.massf[isp] + muj_tilde*dt*dUdt0.massf[isp];
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U0.energies.length) {
                U1.energies[imode] = U0.energies[imode] + muj_tilde*dt*dUdt0.energies[imode];
            }
        }

        // make a copy of the initial conserved quantities
        U[3].copy_values_from(U[0]);
        return;
    } // end rkl2_stage_update_for_flow_on_fixed_grid1()

    void rkl2_stage_update_for_flow_on_fixed_grid2(double dt, int j, int s, bool with_local_time_stepping)
    {
        ConservedQuantities dUdt0;
        ConservedQuantities U0;
        ConservedQuantities U1;
        ConservedQuantities U2;
        ConservedQuantities U3;
        ConservedQuantities dUdtO;
        U0 = U[0];
        U1 = U[1];
        U2 = U[2];
        U3 = U[3];
        dUdt0 = dUdt[1];
        dUdtO = dUdt[0];

        // coefficients
        double ajm1; double bj; double bjm1, bjm2; double muj; double vuj; double muj_tilde; double gam_tilde;

        if (j == 2) {
            bj = 1.0/3.0;
            bjm1 = 1.0/3.0;
            bjm2 = 1.0/3.0;
        } else if (j == 3) {
            bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
            bjm1 = 1.0/3.0;
            bjm2 = 1.0/3.0;
        } else if (j == 4) {
            bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
            bjm1 = ((j-1.0)*(j-1.0)+(j-1.0)-2.0)/(2.0*(j-1.0)*((j-1.0)+1.0));
            bjm2 = 1.0/3.0;
        } else {
            bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
            bjm1 = ((j-1.0)*(j-1.0)+(j-1.0)-2.0)/(2.0*(j-1.0)*((j-1.0)+1.0));
            bjm2 = ((j-2.0)*(j-2.0)+(j-2.0)-2.0)/(2.0*(j-2.0)*((j-2.0)+1.0));
        }
        ajm1 = 1.0-bjm1;
        muj_tilde = (4.0*(2.0*j-1.0))/(j*(s*s+s-2.0)) * (bj/bjm1);
        gam_tilde = -ajm1*muj_tilde;
        muj = (2*j-1.0)/(j) * (bj/bjm1);
        vuj = -(j-1.0)/(j) * (bj/bjm2);

        U2.mass = muj*U1.mass + vuj*U0.mass + (1.0-muj-vuj)*U3.mass + muj_tilde*dt*dUdt0.mass + gam_tilde*dt*dUdtO.mass;
        U2.momentum.set(muj*U1.momentum.x + vuj*U0.momentum.x + (1.0-muj-vuj)*U3.momentum.x + muj_tilde*dt*dUdt0.momentum.x + gam_tilde*dt*dUdtO.momentum.x,
                        muj*U1.momentum.y + vuj*U0.momentum.y + (1.0-muj-vuj)*U3.momentum.y + muj_tilde*dt*dUdt0.momentum.y + gam_tilde*dt*dUdtO.momentum.y,
                        muj*U1.momentum.z + vuj*U0.momentum.z + (1.0-muj-vuj)*U3.momentum.z + muj_tilde*dt*dUdt0.momentum.z + gam_tilde*dt*dUdtO.momentum.z);
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U2.B.set(muj*U1.B.x + vuj*U0.B.x + (1.0-muj-vuj)*U3.B.x + muj_tilde*dt*dUdt0.B.x + gam_tilde*dt*dUdtO.B.x,
                         muj*U1.B.y + vuj*U0.B.y + (1.0-muj-vuj)*U3.B.y + muj_tilde*dt*dUdt0.B.y + gam_tilde*dt*dUdtO.B.y,
                         muj*U1.B.z + vuj*U0.B.z + (1.0-muj-vuj)*U3.B.z + muj_tilde*dt*dUdt0.B.z + gam_tilde*dt*dUdtO.B.z);
                if (myConfig.divergence_cleaning) {
                    U2.psi = muj*U1.psi + vuj*U0.psi + (1.0-muj-vuj)*U3.psi + muj_tilde*dt*dUdt0.psi + gam_tilde*dt*dUdtO.psi;
                    U2.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U2.B.clear();
                U2.psi = 0.0;
            }
        }
        U2.total_energy = muj*U1.total_energy + vuj*U0.total_energy + (1.0-muj-vuj)*U3.total_energy + muj_tilde*dt*dUdt0.total_energy + gam_tilde*dt*dUdtO.total_energy;
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U2.rhoturb[i] = muj*U1.rhoturb[i] + vuj*U0.rhoturb[i] + (1.0-muj-vuj)*U3.rhoturb[i] + muj_tilde*dt*dUdt0.rhoturb[i] + gam_tilde*dt*dUdtO.rhoturb[i];
                U2.rhoturb[i] = fmax(U2.rhoturb[i], U0.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. U2.massf.length) {
                U2.massf[isp] = muj*U1.massf[isp] + vuj*U0.massf[isp] + (1.0-muj-vuj)*U3.massf[isp] + muj_tilde*dt*dUdt0.massf[isp] + gam_tilde*dt*dUdtO.massf[isp];
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U2.energies.length) {
                U2.energies[imode] = muj*U1.energies[imode] + vuj*U0.energies[imode] + muj_tilde*dt*dUdt0.energies[imode] + muj_tilde*dt*dUdt0.energies[imode] + gam_tilde*dt*dUdtO.energies[imode];
            }
        }
	return;
    } // end rkl2_stage_update_for_flow_on_fixed_grid2()

    @nogc
    void stage_1_update_for_flow_on_fixed_grid(double dt, bool force_euler, bool with_local_time_stepping)
    {
        // use the local-time step
        if (with_local_time_stepping) dt = this.dt_local;

        ConservedQuantities dUdt0 = dUdt[0];
        ConservedQuantities U0 = U[0];
        ConservedQuantities U1 = U[1];
        double gamma_1 = 1.0; // for normal Predictor-Corrector or Euler update.
        // In some parts of the code (viscous updates, k-omega updates)
        // we use this function as an Euler update even when the main
        // gasdynamic_update_scheme is of higher order.
        // force_euler is set true for these situations.
        if (!force_euler) {
            final switch (myConfig.gasdynamic_update_scheme) {
            case GasdynamicUpdate.euler:
            case GasdynamicUpdate.moving_grid_1_stage:
            case GasdynamicUpdate.moving_grid_2_stage:
            case GasdynamicUpdate.pc: gamma_1 = 1.0; break;
            case GasdynamicUpdate.midpoint: gamma_1 = 0.5; break;
            case GasdynamicUpdate.classic_rk3: gamma_1 = 0.5; break;
            case GasdynamicUpdate.tvd_rk3: gamma_1 = 1.0; break;
            case GasdynamicUpdate.denman_rk3: gamma_1 = 8.0/15.0; break;
            case GasdynamicUpdate.rkl1:
            case GasdynamicUpdate.rkl2: assert(false, "invalid option");
            }
        }
        U1.mass = U0.mass + dt*gamma_1*dUdt0.mass;
        // Side note:
        // It would be convenient (codewise) for the updates of these Vector3 quantities to
        // be done with the Vector3 arithmetic operators but I suspect that the implementation
        // of those oerators is such that a whole lot of Vector3 temporaries would be created.
        U1.momentum.set(U0.momentum.x + dt*gamma_1*dUdt0.momentum.x,
                        U0.momentum.y + dt*gamma_1*dUdt0.momentum.y,
                        U0.momentum.z + dt*gamma_1*dUdt0.momentum.z);
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U1.B.set(U0.B.x + dt*gamma_1*dUdt0.B.x,
                         U0.B.y + dt*gamma_1*dUdt0.B.y,
                         U0.B.z + dt*gamma_1*dUdt0.B.z);
                U1.divB = dUdt0.divB;
                if (myConfig.divergence_cleaning) {
                    U1.psi = U0.psi + dt*gamma_1*dUdt0.psi;
                    U1.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U1.B.clear();
                U1.psi = 0.0;
                U1.divB = 0.0;
            }
        }
        U1.total_energy = U0.total_energy + dt * gamma_1 * dUdt0.total_energy;
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U1.rhoturb[i] = U0.rhoturb[i] + dt * gamma_1 * dUdt0.rhoturb[i];
                U1.rhoturb[i] = fmax(U1.rhoturb[i], U0.mass * myConfig.turb_model.turb_limits(i));
            }
                // ...assuming a minimum value of 1.0 for omega
                // It may occur (near steps in the wall) that a large flux of romega
                // through one of the cell interfaces causes romega within the cell
                // to drop rapidly.
                // The large values of omega come from Menter's near-wall correction that may be
                // applied outside the control of this finite-volume core code.
                // These large values of omega will be convected along the wall and,
                // if they are convected past a corner with a strong expansion,
                // there will be an unreasonably-large flux out of the cell.
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. myConfig.n_species) {
                U1.massf[isp] = U0.massf[isp] + dt*gamma_1*dUdt0.massf[isp];
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U1.energies.length) {
                U1.energies[imode] = U0.energies[imode] + dt*gamma_1*dUdt0.energies[imode];
            }
        }
        return;
    } // end stage_1_update_for_flow_on_fixed_grid()

    @nogc
    void stage_2_update_for_flow_on_fixed_grid(double dt, bool with_local_time_stepping)
    {
        // use the local-time step
        if (with_local_time_stepping) dt = this.dt_local;

        ConservedQuantities dUdt0 = dUdt[0];
        ConservedQuantities dUdt1 = dUdt[1];
        ConservedQuantities U_old = U[0];
        if (myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) U_old = U[1];
        ConservedQuantities U2 = U[2];
        double gamma_1 = 0.5; // Presume predictor-corrector.
        double gamma_2 = 0.5;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.moving_grid_1_stage: assert(false, "invalid for 1-stage update.");
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc: gamma_1 = 0.5, gamma_2 = 0.5; break;
        case GasdynamicUpdate.midpoint: gamma_1 = 0.0; gamma_2 = 1.0; break;
        case GasdynamicUpdate.classic_rk3: gamma_1 = -1.0; gamma_2 = 2.0; break;
        case GasdynamicUpdate.tvd_rk3: gamma_1 = 0.25; gamma_2 = 0.25; break;
        case GasdynamicUpdate.denman_rk3: gamma_1 = -17.0/60.0; gamma_2 = 5.0/12.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        }
        U2.mass = U_old.mass + dt*(gamma_1*dUdt0.mass + gamma_2*dUdt1.mass);
        U2.momentum.set(U_old.momentum.x + dt*(gamma_1*dUdt0.momentum.x + gamma_2*dUdt1.momentum.x),
                        U_old.momentum.y + dt*(gamma_1*dUdt0.momentum.y + gamma_2*dUdt1.momentum.y),
                        U_old.momentum.z + dt*(gamma_1*dUdt0.momentum.z + gamma_2*dUdt1.momentum.z));
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U2.B.set(U_old.B.x + dt*(gamma_1*dUdt0.B.x + gamma_2*dUdt1.B.x),
                         U_old.B.y + dt*(gamma_1*dUdt0.B.y + gamma_2*dUdt1.B.y),
                         U_old.B.z + dt*(gamma_1*dUdt0.B.z + gamma_2*dUdt1.B.z));
                U2.divB = dUdt0.divB;
                if (myConfig.divergence_cleaning) {
                    U2.psi = U_old.psi + dt*(gamma_1*dUdt0.psi + gamma_2*dUdt1.psi);
                    U2.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U2.B.clear();
                U2.psi = 0.0;
                U2.divB = 0.0;
            }
        }
        U2.total_energy = U_old.total_energy + dt*(gamma_1*dUdt0.total_energy + gamma_2*dUdt1.total_energy);
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U2.rhoturb[i] = U_old.rhoturb[i] + dt*(gamma_1*dUdt0.rhoturb[i] + gamma_2*dUdt1.rhoturb[i]);
                U2.rhoturb[i] = fmax(U2.rhoturb[i], U_old.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. myConfig.n_species) {
                U2.massf[isp] = U_old.massf[isp] + dt*(gamma_1*dUdt0.massf[isp] + gamma_2*dUdt1.massf[isp]);
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U2.energies.length) {
                U2.energies[imode] = U_old.energies[imode] + dt*(gamma_1*dUdt0.energies[imode] +
                                                                 gamma_2*dUdt1.energies[imode]);
            }
        }
        return;
    } // end stage_2_update_for_flow_on_fixed_grid()

    @nogc
    void stage_3_update_for_flow_on_fixed_grid(double dt, bool with_local_time_stepping)
    {
        // use the local-time step
        if (with_local_time_stepping) dt = this.dt_local;

        ConservedQuantities dUdt0 = dUdt[0];
        ConservedQuantities dUdt1 = dUdt[1];
        ConservedQuantities dUdt2 = dUdt[2];
        ConservedQuantities U_old = U[0];
        if (myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) U_old = U[2];
        ConservedQuantities U3 = U[3];
        double gamma_1 = 1.0/6.0; // presume TVD_RK3 scheme.
        double gamma_2 = 1.0/6.0;
        double gamma_3 = 4.0/6.0;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc:
        case GasdynamicUpdate.midpoint:
            assert(false, "invalid for 2-stage update.");
        case GasdynamicUpdate.classic_rk3: gamma_1 = 1.0/6.0; gamma_2 = 4.0/6.0; gamma_3 = 1.0/6.0; break;
        case GasdynamicUpdate.tvd_rk3: gamma_1 = 1.0/6.0; gamma_2 = 1.0/6.0; gamma_3 = 4.0/6.0; break;
            // FIX-ME: Check that we have Andrew Denman's scheme ported correctly.
        case GasdynamicUpdate.denman_rk3: gamma_1 = 0.0; gamma_2 = -5.0/12.0; gamma_3 = 3.0/4.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        }
        U3.mass = U_old.mass + dt * (gamma_1*dUdt0.mass + gamma_2*dUdt1.mass + gamma_3*dUdt2.mass);
        U3.momentum.set(U_old.momentum.x + dt*(gamma_1*dUdt0.momentum.x + gamma_2*dUdt1.momentum.x +
                                               gamma_3*dUdt2.momentum.x),
                        U_old.momentum.y + dt*(gamma_1*dUdt0.momentum.y + gamma_2*dUdt1.momentum.y +
                                               gamma_3*dUdt2.momentum.y),
                        U_old.momentum.z + dt*(gamma_1*dUdt0.momentum.z + gamma_2*dUdt1.momentum.z +
                                               gamma_3*dUdt2.momentum.z));
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U3.B.set(U_old.B.x + dt*(gamma_1*dUdt0.B.x + gamma_2*dUdt1.B.x + gamma_3*dUdt2.B.x),
                         U_old.B.y + dt*(gamma_1*dUdt0.B.y + gamma_2*dUdt1.B.y + gamma_3*dUdt2.B.y),
                         U_old.B.z + dt*(gamma_1*dUdt0.B.z + gamma_2*dUdt1.B.z + gamma_3*dUdt2.B.z));
                if (myConfig.divergence_cleaning) {
                    U3.psi = U_old.psi + dt*(gamma_1*dUdt0.psi + gamma_2*dUdt1.psi + gamma_3*dUdt2.psi);
                    U3.psi *= divergence_damping_factor(dt, myConfig.c_h, myConfig.divB_damping_length);
                }
            } else {
                U3.B.clear();
                U3.psi = 0.0;
            }
        }
        U3.total_energy = U_old.total_energy +
            dt*(gamma_1*dUdt0.total_energy + gamma_2*dUdt1.total_energy + gamma_3*dUdt2.total_energy);
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U3.rhoturb[i] = U_old.rhoturb[i] + dt*(gamma_1*dUdt0.rhoturb[i] + gamma_2*dUdt1.rhoturb[i] + gamma_3*dUdt2.rhoturb[i]);
                U3.rhoturb[i] = fmax(U3.rhoturb[i], U_old.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. myConfig.n_species) {
                U3.massf[isp] = U_old.massf[isp] +
                    dt*(gamma_1*dUdt0.massf[isp] + gamma_2*dUdt1.massf[isp] + gamma_3*dUdt2.massf[isp]);
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U3.energies.length) {
                U3.energies[imode] = U_old.energies[imode] +
                    dt*(gamma_1*dUdt0.energies[imode] + gamma_2*dUdt1.energies[imode] +
                        gamma_3*dUdt2.energies[imode]);
            }
        }
        return;
    } // end stage_3_update_for_flow_on_fixed_grid()

    @nogc
    void stage_1_update_for_flow_on_moving_grid(double dt, bool with_local_time_stepping)
    {
        // use the local-time step
        if (with_local_time_stepping) dt = this.dt_local;

        ConservedQuantities dUdt0 = dUdt[0];
        ConservedQuantities U0 = U[0];
        ConservedQuantities U1 = U[1];
        double gamma_1 = 1.0;
        number vr = volume[0] / volume[1];
        U1.mass = vr*(U0.mass + dt*gamma_1*dUdt0.mass);
        U1.momentum.set(vr*(U0.momentum.x + dt*gamma_1*dUdt0.momentum.x),
                        vr*(U0.momentum.y + dt*gamma_1*dUdt0.momentum.y),
                        vr*(U0.momentum.z + dt*gamma_1*dUdt0.momentum.z));
        version(MHD) {
            // Magnetic field
            if (myConfig.MHD) {
                U1.B.set(vr*(U0.B.x + dt*gamma_1*dUdt0.B.x),
                         vr*(U0.B.y + dt*gamma_1*dUdt0.B.y),
                         vr*(U0.B.z + dt*gamma_1*dUdt0.B.z));
            } else {
                U1.B.clear();
            }
        }
        U1.total_energy = vr*(U0.total_energy + dt*gamma_1*dUdt0.total_energy);
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U1.rhoturb[i] = vr*(U0.rhoturb[i] + dt*gamma_1*dUdt0.rhoturb[i]);
                U1.rhoturb[i] = fmax(U1.rhoturb[i], U0.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. myConfig.n_species) {
                U1.massf[isp] = vr*(U0.massf[isp] + dt*gamma_1*dUdt0.massf[isp]);
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U1.energies.length) {
                U1.energies[imode] = vr*(U0.energies[imode] + dt*gamma_1*dUdt0.energies[imode]);
            }
        }
        return;
    } // end stage_1_update_for_flow_on_moving_grid()

    @nogc
    void stage_2_update_for_flow_on_moving_grid(double dt, bool with_local_time_stepping)
    {
        // use the local-time step
        if (with_local_time_stepping) dt = this.dt_local;

        ConservedQuantities dUdt0 = dUdt[0];
        ConservedQuantities dUdt1 = dUdt[1];
        ConservedQuantities U0 = U[0];
        ConservedQuantities U2 = U[2];
        number gamma_2 = 0.5;
        number gamma_1 = 0.5;
        number v_old = volume[0];
        number vol_inv = 1.0 / volume[2];
        gamma_1 *= volume[0]; gamma_2 *= volume[1]; // Roll-in the volumes for convenience below.
        //
        U2.mass = vol_inv * (v_old * U0.mass + dt * (gamma_1 * dUdt0.mass + gamma_2 * dUdt1.mass));
        U2.momentum.set(vol_inv*(v_old*U0.momentum.x + dt*(gamma_1*dUdt0.momentum.x +
                                                           gamma_2*dUdt1.momentum.x)),
                        vol_inv*(v_old*U0.momentum.y + dt*(gamma_1*dUdt0.momentum.y +
                                                           gamma_2*dUdt1.momentum.y)),
                        vol_inv*(v_old*U0.momentum.z + dt*(gamma_1*dUdt0.momentum.z +
                                                           gamma_2*dUdt1.momentum.z)));
        version(MHD) {
            if (myConfig.MHD) {
                // Magnetic field
                U2.B.set(vol_inv*(v_old*U0.B.x + dt*(gamma_1*dUdt0.B.x + gamma_2*dUdt1.B.x)),
                         vol_inv*(v_old*U0.B.y + dt*(gamma_1*dUdt0.B.y + gamma_2*dUdt1.B.y)),
                         vol_inv*(v_old*U0.B.z + dt*(gamma_1*dUdt0.B.z + gamma_2*dUdt1.B.z)));
            } else {
                U2.B.clear();
            }
        }
        U2.total_energy = vol_inv*(v_old*U0.total_energy +
                                   dt*(gamma_1*dUdt0.total_energy + gamma_2*dUdt1.total_energy));
        version(turbulence) {
            foreach(i; 0 .. myConfig.turb_model.nturb){
                U2.rhoturb[i] = vol_inv*(v_old*U0.rhoturb[i] + dt*(gamma_1*dUdt0.rhoturb[i] + gamma_2*dUdt1.rhoturb[i]));
                U2.rhoturb[i] = fmax(U2.rhoturb[i], U0.mass * myConfig.turb_model.turb_limits(i));
            }
        }
        version(multi_species_gas) {
            foreach(isp; 0 .. myConfig.n_species) {
                U2.massf[isp] = vol_inv*(v_old*U0.massf[isp] +
                                         dt*(gamma_1*dUdt0.massf[isp] + gamma_2*dUdt1.massf[isp]));
            }
        }
        version(multi_T_gas) {
            foreach(imode; 0 .. U2.energies.length) {
                U2.energies[imode] = vol_inv*(v_old*U0.energies[imode] +
                                              dt*(gamma_1*dUdt0.energies[imode] +
                                                  gamma_2*dUdt1.energies[imode]));
            }
        }
        return;
    } // end stage_2_update_for_flow_on_moving_grid()

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
            myConfig.thermochemUpdate(fs.gas, dt, dt_chem, dt_therm, params);
            if (myConfig.ignition_zone_active) {
                // Restore actual gas temperature
                fs.gas.T = T_save;
            }
        } catch(ThermochemicalReactorUpdateException err) {
            // It's probably worth one more try but setting dt_chem = -1.0 to give
            // the ODE solver a fresh chance to find a good timestep.
            dt_chem = -1.0;
            try {
                 myConfig.thermochemUpdate(fs.gas, dt, dt_chem, dt_therm, params);
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
        version(multi_species_gas) {
            // Species densities: mass of species isp per unit volume.
            foreach(isp; 0 .. fs.gas.massf.length) {
                U[0].massf[isp] = fs.gas.rho * fs.gas.massf[isp];
            }
        }
        version(multi_T_gas) {
            // Independent energies energy: Joules per unit volume.
            foreach(imode; 0 .. U[0].energies.length) {
                U[0].energies[imode] = fs.gas.rho * fs.gas.u_modes[imode];
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
        fs.mu_t = myConfig.turb_model.turbulent_viscosity(fs, grad, pos[0].y, dwall);
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
        Q.mass = 0.0;
        Q.momentum.clear();
        Q.total_energy = 0.0;
        version(MHD) {
            Q.B.clear();
        }
        version(turbulence) {
            foreach(ref rt; Q.rhoturb) { rt = 0.0; }
        }
        version(multi_species_gas) {
            foreach(ref elem; Q.massf) { elem = 0.0; }
        }
        version(multi_T_gas) {
            foreach(ref elem; Q.energies) { elem = 0.0; }
        }
        Q_rE_rad = 0.0;
    } // end clear_source_vector()

    @nogc
    void add_inviscid_source_vector(int gtl, double omegaz=0.0)
    // Add the components of the source vector, Q, for inviscid flow.
    //
    // Currently, the axisymmetric equations include the
    // pressure contribution to the y-momentum equation
    // here rather than in the boundary fluxes.
    // By default, assume 2D-planar, or 3D-Cartesian flow.
    {
        if (omegaz != 0.0) {
            // Rotating frame.
            number rho = fs.gas.rho;
            number x = pos[gtl].x;
            number y = pos[gtl].y;
            number wx = fs.vel.x;
            number wy = fs.vel.y;
            // Coriolis and centrifugal forces contribute to momenta.
            Q.momentum.refx += rho * (omegaz*omegaz*x + 2.0*omegaz*wy);
            Q.momentum.refy += rho * (omegaz*omegaz*y - 2.0*omegaz*wx);
            // There is no contribution to the energy equation in the rotating frame
            // because it is implicit in the use of rothalpy as the conserved quantity.
        }
        if (myConfig.axisymmetric) {
            // For axisymmetric flow:
            // pressure contribution from the Front and Back (radial) interfaces.
            Q.momentum.refy += fs.gas.p * areaxy[gtl] / volume[gtl];
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
            Q.total_energy += Q_rE_rad;
            version(multi_T_gas) {
                Q.energies.back() += Q_rE_rad; // FIX-ME old C++ code
            }
        }
        return;
    } // end add_inviscid_source_vector()

    @nogc
    void add_viscous_source_vector()
    {
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
            Q.momentum.refy -= tau_00 * areaxy[0] / volume[0];
        } // end if ( myConfig.axisymmetric )

        version(turbulence) {
            if (in_turbulent_zone) {
                myConfig.turb_model.source_terms(fs, grad, pos[0].y, dwall, L_min, L_max, Q.rhoturb);
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
    void add_chemistry_source_vector()
    {
        version(multi_species_gas){
        rmech.eval_source_terms(myConfig.gmodel, fs.gas, chem_conc, chem_rates, chem_source);
        foreach(sp, ref elem; Q.massf) { elem += chem_source[sp]; }
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
        grad.copy_values_from(vtx[0].grad);
        foreach (i; 1 .. vtx.length) grad.accumulate_values_from(vtx[i].grad);
        grad.scale_values_by(to!number(1.0/vtx.length));
    } // end average_vertex_deriv_values()

    @nogc
    void average_interface_deriv_values()
    {
        grad.copy_values_from(iface[0].grad);
        foreach (i; 1 .. iface.length) grad.accumulate_values_from(iface[i].grad);
        grad.scale_values_by(to!number(1.0/iface.length));
    } // end average_interface_deriv_values()

    @nogc
    void lusgs_startup_iteration(number dtInv, double omega, ref number[] dU, number[] R)
    {
        // Compute LHS Jacobian diagonal (D) and evaluate dU0 = D^{-1} * R
        number lambda = 0.0;
        foreach (f; iface) {
            lambda += f.area[0]*f.spectral_radius(omega);
        }

        scalar_diag_inv[] = (dtInv + 0.5*lambda/volume[0])^^(-1.0);
        dU[] = scalar_diag_inv[]*R[];
    } // end lusgs_startup_iteration()

    @nogc
    void lusgs_relaxation_iteration(double omega, bool matrix_based, ref number[] dU, number[] R)
    {
        // Compute a relaxation subiteration dU^{k+1} = D^{-1} * (R - 0.5*LU)

        // Make a stack-local copy of conserved quantities info
        size_t nConserved = myConfig.cqi.nConservedQuantities;
        size_t MASS = myConfig.cqi.mass;
        size_t X_MOM = myConfig.cqi.xMom;
        size_t Y_MOM = myConfig.cqi.yMom;
        size_t Z_MOM = myConfig.cqi.zMom;
        size_t TOT_ENERGY = myConfig.cqi.totEnergy;
        size_t TKE = myConfig.cqi.tke;
        size_t nturb = myConfig.turb_model.nturb;

        LU[] = to!number(0.0);
        // loop through neighbouring cells and approximate off-diagonal terms (L+U)
        foreach (i; 1..cell_cloud.length) {
            FVCell nc = cell_cloud[i]; FVInterface f = iface[i-1];
            number lambda = f.spectral_radius(omega);
            if (matrix_based) {
                nc.roeFluxJacobian(f);
                dot(nc.dFdU, nc.dUk[0..nConserved], nc.dF[0..nConserved]);
            } else { // matrix free flux increment
                nc.evalMatrixFreeFluxIncrement(f);
            }
            LU[MASS] += (nc.dF[MASS]*outsign[i-1] - lambda*nc.dUk[MASS])*f.area[0];
            LU[X_MOM] += (nc.dF[X_MOM]*outsign[i-1] - lambda*nc.dUk[X_MOM])*f.area[0];
            LU[Y_MOM] += (nc.dF[Y_MOM]*outsign[i-1] - lambda*nc.dUk[Y_MOM])*f.area[0];
            if (myConfig.dimensions == 3)
                LU[Z_MOM] += (nc.dF[Z_MOM]*outsign[i-1] - lambda*nc.dUk[Z_MOM])*f.area[0];
            LU[TOT_ENERGY] += (nc.dF[TOT_ENERGY]*outsign[i-1] - lambda*nc.dUk[TOT_ENERGY])*f.area[0];
            foreach(it; 0 .. nturb) {
                LU[TKE+it] += (nc.dF[TKE+it]*outsign[i-1] - lambda*nc.dUk[TKE+it])*f.area[0];
            }
        }
        LU[] *= 0.5/volume[0];

        dU[] = R[] - LU[];
        dU[] *= scalar_diag_inv[];
    }

    @nogc
    void evalMatrixFreeFluxIncrement(FVInterface f)
    {
        // Matrix-Free Flux vector increment
        //
        // As defined on right column, pg 4 of
        // Rieger and Jameson (1988),
        // Solution of Steady Three-Dimensional Compressible Euler and Navier-Stokes Equations by an Implicit LU Scheme
        // AIAA conference paper
        //
        // Uses Roe's split flux scheme for LHS Jacobian as per
        // Luo, Baum, and Lohner (1998)
        // A Fast, Matrix-free Implicit Method for Compressible Flows on Unstructured Grids,
        // Journal of computational physics
        // 
        
        // Make a stack-local copy of conserved quantities info
        size_t nConserved = myConfig.cqi.nConservedQuantities;
        size_t MASS = myConfig.cqi.mass;
        size_t X_MOM = myConfig.cqi.xMom;
        size_t Y_MOM = myConfig.cqi.yMom;
        size_t Z_MOM = myConfig.cqi.zMom;
        size_t TOT_ENERGY = myConfig.cqi.totEnergy;
        size_t TKE = myConfig.cqi.tke;
        
        size_t nturb = myConfig.turb_model.nturb;

        // make sure cells have conserved quantities filled
        encode_conserved(0, 0, 0.0);
        
        // peturb conserved quantities by approximation of dU
        U[1].copy_values_from(U[0]);
        U[1].mass += dUk[MASS];
        U[1].momentum.refx += dUk[X_MOM];
        U[1].momentum.refy += dUk[Y_MOM];
        if (GlobalConfig.dimensions == 3 )
            U[1].momentum.refz += dUk[Z_MOM];
        U[1].total_energy += dUk[TOT_ENERGY];
        foreach(it; 0 .. nturb) {
            U[1].rhoturb[it] += dUk[TKE+it];
        }
        
        // update primitive variables
        decode_conserved(0, 1, 0.0);          
        
        // Peturbed state flux
        number rho = fs.gas.rho;
        number velx = fs.vel.dot(f.n);
        number vely = fs.vel.dot(f.t1);
        number velz = fs.vel.dot(f.t2);
        number p = fs.gas.p;
        auto gmodel = myConfig.gmodel;
        number e = gmodel.internal_energy(fs.gas);
        
        dF[MASS]= rho*velx;
        dF[X_MOM] = p + rho*velx*velx;
        dF[Y_MOM] = rho*velx*vely;
        if (myConfig.dimensions == 3)
            dF[Z_MOM] = rho*velx*velz;
        dF[TOT_ENERGY] = (rho*e + rho*(velx^^2 + vely^^2 + velz^^2)/2.0 + p)*velx;
        foreach(it; 0 .. nturb) {
            dF[TKE+it] = rho*velx*fs.turb[it];
        }
        
        // reset primitive variables to unperturbed state
        decode_conserved(0, 0, 0.0);          
        
        // original state flux
        rho = fs.gas.rho;
        velx = fs.vel.dot(f.n);
        vely = fs.vel.dot(f.t1);
        velz = fs.vel.dot(f.t2);
        p = fs.gas.p;
        e = gmodel.internal_energy(fs.gas);
        
        // flux vector increment
        dF[MASS] -= rho*velx;
        dF[X_MOM] -= p + rho*velx*velx;
        dF[Y_MOM] -= rho*velx*vely;
        if (myConfig.dimensions == 3)
            dF[Z_MOM] -= rho*velx*velz;

        number global_mom_x = dF[X_MOM]*f.n.x + dF[Y_MOM]*f.t1.x; // global-x
        number global_mom_y = dF[X_MOM]*f.n.y + dF[Y_MOM]*f.t1.y; // global-y
        number global_mom_z;
        if (myConfig.dimensions == 3) {
            global_mom_x += dF[Z_MOM]*f.t2.x;
            global_mom_y += dF[Z_MOM]*f.t2.y;
            global_mom_z = dF[X_MOM]*f.n.z + dF[Y_MOM]*f.t1.z + dF[Z_MOM]*f.t2.z; // global-z
        }
        dF[X_MOM] = global_mom_x;
        dF[Y_MOM] = global_mom_y;
        if (myConfig.dimensions == 3)
            dF[Z_MOM] = global_mom_z;
        
        dF[TOT_ENERGY] -= (rho*e + rho*(velx^^2 + vely^^2 + velz^^2)/2.0 + p)*velx;
        foreach(it; 0 .. nturb) {
            dF[TKE+it] -= rho*velx*fs.turb[it];
        }
    } // end evalMatrixFreeFluxIncrement()

    @nogc
    void roeFluxJacobian(FVInterface f)
    {
        // Hand differentiation of Roe's split flux scheme for LHS Jacobian as per
        // Luo, Baum, and Lohner (1998)
        // A Fast, Matrix-free Implicit Method for Compressible Flows on Unstructured Grids,
        // Journal of computational physics
        // 
        
        // Make a stack-local copy of conserved quantities info
        size_t nConserved = myConfig.cqi.nConservedQuantities;
        size_t MASS = myConfig.cqi.mass;
        size_t X_MOM = myConfig.cqi.xMom;
        size_t Y_MOM = myConfig.cqi.yMom;
        size_t Z_MOM = myConfig.cqi.zMom;
        size_t TOT_ENERGY = myConfig.cqi.totEnergy;
        size_t TKE = myConfig.cqi.tke;

        // primitive variables
        auto gmodel = myConfig.gmodel;
        number gam = gmodel.gamma(fs.gas);
        number rho = fs.gas.rho;
        // rotate velocity into interface reference frame
        number u = fs.vel.dot(f.n);
        number v = fs.vel.dot(f.t1);
        number w = fs.vel.dot(f.t2);
        number p = fs.gas.p;
        number e = gmodel.internal_energy(fs.gas);
        
        // conserved variables
        number U1 = rho; 
        number U2 = rho*u;
        number U3 = rho*v;
        number U4 = rho*w;
        number U5 = rho*e + rho*(u^^2 + v^^2 + w^^2)/2.0;
        
        // approximate flux Jacobian based on Roe's approximate split flux scheme
        dFdU[MASS,MASS] = to!number(0.0);
        dFdU[MASS,X_MOM] = to!number(1.0);
        dFdU[MASS,Y_MOM] = to!number(0.0);
        if (myConfig.dimensions == 3) { dFdU[X_MOM,Z_MOM] = to!number(0.0); }
        dFdU[MASS,TOT_ENERGY] = to!number(0.0);
        
        dFdU[X_MOM,MASS] = -(U2*U2)/(U1*U1) + (gam-1.0)*(U2*U2+U3*U3+U4*U4)/(2.0*U1*U1);
        dFdU[X_MOM,X_MOM] = (3.0-gam)*(U2/U1);
        dFdU[X_MOM,Y_MOM] = (1.0-gam)*(U3/U1);
        if (myConfig.dimensions == 3) { dFdU[X_MOM,Z_MOM] = (1.0-gam)*(U4/U1); }
        dFdU[X_MOM,TOT_ENERGY] = (gam-1.0);
        
        dFdU[Y_MOM,MASS] = -(U2*U3)/(U1*U1);
        dFdU[Y_MOM,X_MOM] = U3/U1;
        dFdU[Y_MOM,Y_MOM] = U2/U1;
        if (myConfig.dimensions == 3) { dFdU[Y_MOM,Z_MOM] = to!number(0.0); }
        dFdU[Y_MOM,TOT_ENERGY] = to!number(0.0);
        
        if (myConfig.dimensions == 3) {
            dFdU[Z_MOM,MASS] = -(U2*U4)/(U1*U1);
            dFdU[Z_MOM,X_MOM] = U4/U1;
            dFdU[Z_MOM,Y_MOM] = to!number(0.0);
            dFdU[Z_MOM,Z_MOM] = U2/U1;
            dFdU[Z_MOM,TOT_ENERGY] = to!number(0.0);
        }
        
        dFdU[TOT_ENERGY,MASS] = -gam*(U5*U2)/(U1*U1) + (gam-1.0)*(U2*U2*U2+U2*U3*U3+U2*U4*U4)/(U1*U1*U1);
        dFdU[TOT_ENERGY,X_MOM] = gam*(U5/U1) + (1.0-gam)*(3*U2*U2+U3*U3+U4*U4)/(2*U1*U1);
        dFdU[TOT_ENERGY,Y_MOM] = (1.0-gam)*(U3*U2)/(U1*U1);
        if (myConfig.dimensions == 3) { dFdU[TOT_ENERGY,Z_MOM] = (1.0-gam)*(U4*U2)/(U1*U1); }
        dFdU[TOT_ENERGY,TOT_ENERGY] = gam*(U2/U1); 

        size_t nturb = myConfig.turb_model.nturb;
        foreach(it; 0 .. nturb) {
            dFdU[TKE+it, MASS] = -rho*fs.turb[it]*U2/(U1*U1);
            dFdU[TKE+it, X_MOM] = rho*fs.turb[it]/U1;
            dFdU[TKE+it, TKE+it] = U2/U1;
        }
        
        // rotate matrix back into the global reference frame
        dot(f.Tinv, dFdU, dFdU_rotated);
        dot(dFdU_rotated, f.T, dFdU);

    } // end roeFluxJacobian()

    void gather_residual_stencil_lists(size_t spatial_order_of_jacobian)
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
        bool nearest_neighbours_only = false;
        if ( (spatial_order_of_jacobian == 1 && include_viscous_effects == false) ||
             spatial_order_of_jacobian == 0) { nearest_neighbours_only = true; }

        if (nearest_neighbours_only) {
            // this first order stencil includes the cell and its nearest neighbours

            // gather cells
            size_t[] cell_ids;
            unordered_cell_list ~= cell_cloud[0];
            cell_pos_array[cell_cloud[0].id] = unordered_cell_list.length-1;
            cell_ids ~= cell_cloud[0].id;
            foreach (cell; cell_cloud) {
                bool cell_exists = cell_ids.canFind(cell.id);
                if (!cell_exists && cell.id < 1_000_000_000) {
                    unordered_cell_list ~= cell;
                    cell_pos_array[cell.id] = unordered_cell_list.length-1;
                    cell_ids ~= cell.id;
                }
            } // finished gathering cells

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
        } else {
            // second order (&/or viscous, or first order with viscous effects) stencil includes the cell
            // and its nearest neighbours as well as the nearest neighbours of the nearest neighbours

            // gather cells
            size_t[] cell_ids;
            unordered_cell_list ~= cell_cloud[0];
            cell_pos_array[cell_cloud[0].id] = unordered_cell_list.length-1;
            cell_ids ~= cell_cloud[0].id;
            foreach (icell; 1 .. cell_cloud.length) {
                foreach (cell; cell_cloud[icell].cell_cloud) {
                    bool cell_exists = cell_ids.canFind(cell.id);
                    if (!cell_exists && cell.id < 1_000_000_000) {
                        unordered_cell_list ~= cell;
                        cell_pos_array[cell.id] = unordered_cell_list.length-1;
                        cell_ids ~= cell.id;
                    }
                }
            } // finished gathering cells

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
        }
    } // end gather_residual_stencil_lists()

    void gather_residual_stencil_lists_for_ghost_cells(size_t spatial_order_of_jacobian, FVCell[] neighbour_cell_cloud)
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

        if (nearest_neighbours_only) {
            // this first order stencil includes the ghost cells nearest neighbours
            cell_list ~= neighbour_cell_cloud[0]; // this is the interior cell that shares an interface
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
        } else {
            // second order (&/or viscous, or first order with viscous effects) stencil includes the cell
            // and its nearest neighbours as well as the nearest neighbours of the nearest neighbours

            // gather cells
            foreach (c; neighbour_cell_cloud) {
                if ( c.id != id && c.is_interior_to_domain) {
                    cell_list ~= c;
                }
            } // finished gathering cells

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
        }
    } // end gather_residual_stencil_lists_for_ghost_cells()

} // end class FVCell

//--------------------------------------------------------------------------------
// The following functions define the written fixed-formats for the cell data.
// Other input and output functions should delegate their work to these functions.
//--------------------------------------------------------------------------------

string cell_data_as_string(ref const(Vector3) pos, number volume, ref const(FlowState) fs,
                           number Q_rad_org, number f_rad_org, number Q_rE_rad,
                           bool with_local_time_stepping, double dt_local, double dt_chem, double dt_therm,
                           bool include_quality,
                           bool MHDflag, bool divergence_cleaning,
                           bool radiation, size_t nturb)
{
    // We'll treat this function as the master definition of the data format.
    auto writer = appender!string();
    version(complex_numbers) {
        // For complex_numbers, we presently write out only the real parts.
        // [TODO] Maybe we should write full complex numbers.
        formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
                       pos.x.re, pos.y.re, pos.z.re, volume.re, fs.gas.rho.re,
                       fs.vel.x.re, fs.vel.y.re, fs.vel.z.re);
        version(MHD) {
            if (MHDflag) { formattedWrite(writer, " %.18e %.18e %.18e %.18e", fs.B.x.re, fs.B.y.re, fs.B.z.re, fs.divB.re); }
            if (MHDflag && divergence_cleaning) { formattedWrite(writer, " %.18e", fs.psi.re); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { formattedWrite(writer, " %.18e", fs.gas.quality.re); }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.gas.p.re, fs.gas.a.re, fs.gas.mu.re);
        formattedWrite(writer, " %.18e", fs.gas.k.re);
        version(multi_T_gas) {
            foreach (kvalue; fs.gas.k_modes) { formattedWrite(writer, " %.18e", kvalue.re); }
        }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.mu_t.re, fs.k_t.re, fs.S.re);
        if (radiation) { formattedWrite(writer, " %.18e %.18e %.18e", Q_rad_org.re, f_rad_org.re, Q_rE_rad.re); }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                formattedWrite(writer, " %.18e", fs.turb[it].re);
            }
        }
        version(multi_species_gas) {
            foreach (massfvalue; fs.gas.massf) { formattedWrite(writer, " %.18e", massfvalue.re); } 
            if (fs.gas.massf.length > 1) { formattedWrite(writer, " %.18e", dt_chem); }
        } else {
            formattedWrite(writer, " %.18e", 1.0); // single-species mass fraction
        }
        formattedWrite(writer, " %.18e %.18e", fs.gas.u.re, fs.gas.T.re);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                formattedWrite(writer, " %.18e %.18e", fs.gas.u_modes[imode].re, fs.gas.T_modes[imode].re);
            }
            if (fs.gas.u_modes.length > 0) { formattedWrite(writer, " %.18e", dt_therm); }
        }
        if (with_local_time_stepping) formattedWrite(writer, " %.18e", dt_local);
    } else {
        // version double_numbers
        formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
                       pos.x, pos.y, pos.z, volume, fs.gas.rho,
                       fs.vel.x, fs.vel.y, fs.vel.z);
        version(MHD) {
            if (MHDflag) { formattedWrite(writer, " %.18e %.18e %.18e %.18e", fs.B.x, fs.B.y, fs.B.z, fs.divB); }
            if (MHDflag && divergence_cleaning) { formattedWrite(writer, " %.18e", fs.psi); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { formattedWrite(writer, " %.18e", fs.gas.quality); }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.gas.p, fs.gas.a, fs.gas.mu);
        formattedWrite(writer, " %.18e", fs.gas.k);
        version(multi_T_gas) {
            foreach (kvalue; fs.gas.k_modes) { formattedWrite(writer, " %.18e", kvalue); }
        }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.mu_t, fs.k_t, fs.S);
        if (radiation) { formattedWrite(writer, " %.18e %.18e %.18e", Q_rad_org, f_rad_org, Q_rE_rad); }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                formattedWrite(writer, " %.18e", fs.turb[it]);
            }
        }
        version(multi_species_gas) {
            foreach (massfvalue; fs.gas.massf) { formattedWrite(writer, " %.18e", massfvalue); }
            if (fs.gas.massf.length > 1) { formattedWrite(writer, " %.18e", dt_chem); }
        } else {
            formattedWrite(writer, " %.18e", 1.0); // single-species mass fraction
        }
        formattedWrite(writer, " %.18e %.18e", fs.gas.u, fs.gas.T);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                formattedWrite(writer, " %.18e %.18e", fs.gas.u_modes[imode], fs.gas.T_modes[imode]);
            }
            if (fs.gas.u_modes.length > 0) { formattedWrite(writer, " %.18e", dt_therm); }
        }
        if (with_local_time_stepping) { formattedWrite(writer, " %.18e", dt_local); }
    } // end version double_numbers
    return writer.data;
} // end cell_data_as_string()

void cell_data_to_raw_binary(ref File fout,
                             ref const(Vector3) pos, number volume, ref const(FlowState) fs,
                             number Q_rad_org, number f_rad_org, number Q_rE_rad,
                             bool with_local_time_stepping, double dt_local, double dt_chem, double dt_therm,
                             bool include_quality, bool MHDflag, bool divergence_cleaning,
                             bool radiation, size_t nturb)
{
    // This function should match function cell_data_as_string()
    // which is considered the master definition of the data format.
    // We have tried to keep the same code layout.  There is some history.
    // 2017-09-02:
    // Change to using all double values so that the FlowSolution reader becomes simpler.
    //
    version(complex_numbers) {
        // For complex_numbers, we presently write out only the real parts.
        // [TODO] Maybe we should write full complex numbers.
        // Fixed-length buffers to hold data for sending to binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        //
        dbl4[0] = pos.x.re; dbl4[1] = pos.y.re; dbl4[2] = pos.z.re; dbl4[3] = volume.re;
        fout.rawWrite(dbl4);
        dbl4[0] = fs.gas.rho.re; dbl4[1] = fs.vel.x.re; dbl4[2] = fs.vel.y.re; dbl4[3] = fs.vel.z.re;
        fout.rawWrite(dbl4);
        version(MHD) {
            if (MHDflag) {
                dbl4[0] = fs.B.x.re; dbl4[1] = fs.B.y.re; dbl4[2] = fs.B.z.re; dbl4[3] = fs.divB.re;
                fout.rawWrite(dbl4);
            }
            if (MHDflag && divergence_cleaning) { dbl1[0] = fs.psi.re; fout.rawWrite(dbl1); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { dbl1[0] = fs.gas.quality.re; fout.rawWrite(dbl1); }
        dbl4[0] = fs.gas.p.re; dbl4[1] = fs.gas.a.re; dbl4[2] = fs.gas.mu.re; dbl4[3] = fs.gas.k.re;
        fout.rawWrite(dbl4);
        foreach (kvalue; fs.gas.k_modes) { dbl1[0] = kvalue.re; fout.rawWrite(dbl1); }
        dbl2[0] = fs.mu_t.re; dbl2[1] = fs.k_t.re; fout.rawWrite(dbl2);
        dbl1[0] = fs.S.re; fout.rawWrite(dbl1);
        if (radiation) {
            dbl3[0] = Q_rad_org.re; dbl3[1] = f_rad_org.re; dbl3[2] = Q_rE_rad.re;
            fout.rawWrite(dbl3);
        }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                dbl1[0] = fs.turb[it].re; fout.rawWrite(dbl1);
            }
        }
        version(multi_species_gas) {
            foreach (mf; fs.gas.massf) { dbl1[0] = mf.re; fout.rawWrite(dbl1); }
            if (fs.gas.massf.length > 1) { dbl1[0] = dt_chem; fout.rawWrite(dbl1); }
        } else {
            dbl1[0] = 1.0; fout.rawWrite(dbl1); // single-species mass fraction
        }
        dbl2[0] = fs.gas.u.re; dbl2[1] = fs.gas.T.re; fout.rawWrite(dbl2);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                dbl2[0] = fs.gas.u_modes[imode].re; dbl2[1] = fs.gas.T_modes[imode].re;
                fout.rawWrite(dbl2);
            }
            if (fs.gas.u_modes.length > 0) { dbl1[0] = dt_therm; fout.rawWrite(dbl1); }
        }
        if (with_local_time_stepping) { dbl1[0] = dt_local; fout.rawWrite(dbl1); }
    } else {
        // version double_numbers
        // Fixed-length buffers to hold data for sending to binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        //
        dbl4[0] = pos.x; dbl4[1] = pos.y; dbl4[2] = pos.z; dbl4[3] = volume;
        fout.rawWrite(dbl4);
        dbl4[0] = fs.gas.rho; dbl4[1] = fs.vel.x; dbl4[2] = fs.vel.y; dbl4[3] = fs.vel.z;
        fout.rawWrite(dbl4);
        version(MHD) {
            if (MHDflag) {
                dbl4[0] = fs.B.x; dbl4[1] = fs.B.y; dbl4[2] = fs.B.z; dbl4[3] = fs.divB;
                fout.rawWrite(dbl4);
            }
            if (MHDflag && divergence_cleaning) { dbl1[0] = fs.psi; fout.rawWrite(dbl1); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { dbl1[0] = fs.gas.quality; fout.rawWrite(dbl1); }
        dbl4[0] = fs.gas.p; dbl4[1] = fs.gas.a; dbl4[2] = fs.gas.mu; dbl4[3] = fs.gas.k;
        fout.rawWrite(dbl4);
        version(multi_species_gas) {
            foreach (kvalue; fs.gas.k_modes) { dbl1[0] = kvalue; fout.rawWrite(dbl1); }
        }
        dbl2[0] = fs.mu_t; dbl2[1] = fs.k_t; fout.rawWrite(dbl2);
        dbl1[0] = fs.S; fout.rawWrite(dbl1);
        if (radiation) {
            dbl3[0] = Q_rad_org; dbl3[1] = f_rad_org; dbl3[2] = Q_rE_rad;
            fout.rawWrite(dbl3);
        }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                dbl1[0] = fs.turb[it]; fout.rawWrite(dbl1);
            }
        }
        version(multi_species_gas) {
            fout.rawWrite(fs.gas.massf);
            if (fs.gas.massf.length > 1) { dbl1[0] = dt_chem; fout.rawWrite(dbl1); }
        } else {
            dbl1[0] = 1.0; fout.rawWrite(dbl1); // single-species mass fraction
        }
        dbl2[0] = fs.gas.u; dbl2[1] = fs.gas.T; fout.rawWrite(dbl2);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                dbl2[0] = fs.gas.u_modes[imode]; dbl2[1] = fs.gas.T_modes[imode];
                fout.rawWrite(dbl2);
            }
            if (fs.gas.u_modes.length > 0) { dbl1[0] = dt_therm; fout.rawWrite(dbl1); }
        }
        if (with_local_time_stepping) { dbl1[0] = dt_local; fout.rawWrite(dbl1); }
    } // end version double_numbers
    return;
} // end cell_data_to_raw_binary()

void scan_cell_data_from_fixed_order_string
(string buffer,
 ref Vector3 pos, ref number volume, ref FlowState fs,
 ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
 bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
 bool include_quality, bool MHDflag, bool divergence_cleaning, bool radiation, size_t nturb)
{
    // This function needs to be kept consistent with cell_data_as_string() above.
    auto items = split(buffer);
    version(complex_numbers) {
        // For complex_numbers, we presently set only the real parts.
        // [TODO] Maybe we should read full complex numbers.
        pos.refx = Complex!double(items.front); items.popFront();
        pos.refy = Complex!double(items.front); items.popFront();
        pos.refz = Complex!double(items.front); items.popFront();
        volume = Complex!double(items.front); items.popFront();
        fs.gas.rho = Complex!double(items.front); items.popFront();
        fs.vel.refx = Complex!double(items.front); items.popFront();
        fs.vel.refy = Complex!double(items.front); items.popFront();
        fs.vel.refz = Complex!double(items.front); items.popFront();
        version(MHD) {
            if (MHDflag) {
                fs.B.refx = Complex!double(items.front); items.popFront();
                fs.B.refy = Complex!double(items.front); items.popFront();
                fs.B.refz = Complex!double(items.front); items.popFront();
                fs.divB = Complex!double(items.front); items.popFront();
                if (divergence_cleaning) {
                    fs.psi = Complex!double(items.front); items.popFront();
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fs.gas.quality = Complex!double(items.front); items.popFront();
        } else {
            fs.gas.quality = 1.0;
        }
        fs.gas.p = Complex!double(items.front); items.popFront();
        fs.gas.a = Complex!double(items.front); items.popFront();
        fs.gas.mu = Complex!double(items.front); items.popFront();
        fs.gas.k = Complex!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fs.gas.k_modes[i] = Complex!double(items.front); items.popFront();
            }
        }
        fs.mu_t = Complex!double(items.front); items.popFront();
        fs.k_t = Complex!double(items.front); items.popFront();
        fs.S = Complex!double(items.front); items.popFront();
        if (radiation) {
            Q_rad_org = Complex!double(items.front); items.popFront();
            f_rad_org = Complex!double(items.front); items.popFront();
            Q_rE_rad = Complex!double(items.front); items.popFront();
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(it; 0 .. nturb) {
                fs.turb[it] = Complex!double(items.front); items.popFront();
            }
        }
        version(multi_species_gas) {
            foreach(i; 0 .. fs.gas.massf.length) {
                fs.gas.massf[i] = Complex!double(items.front); items.popFront();
            }
            if (fs.gas.massf.length > 1) {
                dt_chem = to!double(items.front); items.popFront();
            }
        } else {
            items.popFront(); // discard the single-species mass fraction, assumed 1.0
        }
        fs.gas.u = Complex!double(items.front); items.popFront();
        fs.gas.T = Complex!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fs.gas.u_modes[i] = Complex!double(items.front); items.popFront();
                fs.gas.T_modes[i] = Complex!double(items.front); items.popFront();
            }
            if (fs.gas.u_modes.length > 0) {
                dt_therm = to!double(items.front); items.popFront();
            }
        }
        if (with_local_time_stepping) { dt_local = to!double(items.front); items.popFront(); }
    } else {
        // version double_numbers
        pos.refx = to!double(items.front); items.popFront();
        pos.refy = to!double(items.front); items.popFront();
        pos.refz = to!double(items.front); items.popFront();
        volume = to!double(items.front); items.popFront();
        fs.gas.rho = to!double(items.front); items.popFront();
        fs.vel.refx = to!double(items.front); items.popFront();
        fs.vel.refy = to!double(items.front); items.popFront();
        fs.vel.refz = to!double(items.front); items.popFront();
        version(MHD) {
            if (MHDflag) {
                fs.B.refx = to!double(items.front); items.popFront();
                fs.B.refy = to!double(items.front); items.popFront();
                fs.B.refz = to!double(items.front); items.popFront();
                fs.divB = to!double(items.front); items.popFront();
                if (divergence_cleaning) {
                    fs.psi = to!double(items.front); items.popFront();
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fs.gas.quality = to!double(items.front); items.popFront();
        } else {
            fs.gas.quality = 1.0;
        }
        fs.gas.p = to!double(items.front); items.popFront();
        fs.gas.a = to!double(items.front); items.popFront();
        fs.gas.mu = to!double(items.front); items.popFront();
        fs.gas.k = to!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fs.gas.k_modes[i] = to!double(items.front); items.popFront();
            }
        }
        fs.mu_t = to!double(items.front); items.popFront();
        fs.k_t = to!double(items.front); items.popFront();
        fs.S = to!double(items.front); items.popFront();
        if (radiation) {
            Q_rad_org = to!double(items.front); items.popFront();
            f_rad_org = to!double(items.front); items.popFront();
            Q_rE_rad = to!double(items.front); items.popFront();
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(it; 0 .. nturb) {
                fs.turb[it] = to!double(items.front); items.popFront();
            }
        }
        version(multi_species_gas) {
            foreach(i; 0 .. fs.gas.massf.length) {
                fs.gas.massf[i] = to!double(items.front); items.popFront();
            }
            if (fs.gas.massf.length > 1) {
                dt_chem = to!double(items.front); items.popFront();
            }
        } else {
            items.popFront(); // discard single-species mass fraction
        }
        fs.gas.u = to!double(items.front); items.popFront();
        fs.gas.T = to!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fs.gas.u_modes[i] = to!double(items.front); items.popFront();
                fs.gas.T_modes[i] = to!double(items.front); items.popFront();
            }
            if (fs.gas.u_modes.length > 0) {
                dt_therm = to!double(items.front); items.popFront();
            }
        }
        if (with_local_time_stepping) { dt_local = to!double(items.front); items.popFront(); }
    } // end version double_numbers
} // end scan_values_from_fixed_order_string()

void scan_cell_data_from_variable_order_string
(string buffer, string[] varNameList, GasModel gmodel, const ref TurbulenceModel tm,
 ref Vector3 pos, ref number volume, ref FlowState fs,
 ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
 bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
 bool include_quality, bool MHDflag, bool divergence_cleaning, bool radiation)
{
    // This function uses the list of variable names read from the file
    // to work out which data item to assign to each variable.
    version(complex_numbers) {
        Complex!double[] values;
        // Note that we expect only the real part of each item to be in the string.
        foreach (item; buffer.strip().split()) { values ~= to!(Complex!double)(item); }
    } else {
        double[] values;
        foreach (item; buffer.strip().split()) { values ~= to!double(item); }
    }
    pos.refx = values[countUntil(varNameList, flowVarName(FlowVar.pos_x))];
    pos.refy = values[countUntil(varNameList, flowVarName(FlowVar.pos_y))];
    pos.refz = values[countUntil(varNameList, flowVarName(FlowVar.pos_z))];
    volume = values[countUntil(varNameList, flowVarName(FlowVar.volume))];
    fs.gas.rho = values[countUntil(varNameList, flowVarName(FlowVar.rho))];
    fs.vel.refx = values[countUntil(varNameList, flowVarName(FlowVar.vel_x))];
    fs.vel.refy = values[countUntil(varNameList, flowVarName(FlowVar.vel_y))];
    fs.vel.refz = values[countUntil(varNameList, flowVarName(FlowVar.vel_z))];
    version(MHD) {
        if (MHDflag) {
            fs.B.refx = values[countUntil(varNameList, flowVarName(FlowVar.B_x))];
            fs.B.refy = values[countUntil(varNameList, flowVarName(FlowVar.B_y))];
            fs.B.refz = values[countUntil(varNameList, flowVarName(FlowVar.B_z))];
            fs.divB = values[countUntil(varNameList, flowVarName(FlowVar.divB))];
            if (divergence_cleaning) {
                fs.psi = values[countUntil(varNameList, flowVarName(FlowVar.psi))];
            } else {
                fs.psi = 0.0;
            }
        } else {
            fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
        }
    }
    if (include_quality) {
        fs.gas.quality = values[countUntil(varNameList, flowVarName(FlowVar.quality))];
    } else {
        fs.gas.quality = 1.0;
    }
    fs.gas.p = values[countUntil(varNameList, flowVarName(FlowVar.p))];
    fs.gas.a = values[countUntil(varNameList, flowVarName(FlowVar.a))];
    fs.gas.mu = values[countUntil(varNameList, flowVarName(FlowVar.mu))];
    fs.gas.k = values[countUntil(varNameList, flowVarName(FlowVar.k))];
    version(multi_T_gas) {
        foreach(i; 0 .. fs.gas.k_modes.length) {
            fs.gas.k_modes[i] = values[countUntil(varNameList, k_modesName(to!int(i)))];
        }
    }
    fs.mu_t = values[countUntil(varNameList, flowVarName(FlowVar.mu_t))];
    fs.k_t = values[countUntil(varNameList, flowVarName(FlowVar.k_t))];
    fs.S = values[countUntil(varNameList, flowVarName(FlowVar.S))];
    if (radiation) {
        Q_rad_org = values[countUntil(varNameList, flowVarName(FlowVar.Q_rad_org))];
        f_rad_org = values[countUntil(varNameList, flowVarName(FlowVar.f_rad_org))];
        Q_rE_rad = values[countUntil(varNameList, flowVarName(FlowVar.Q_rE_rad))];
    } else {
        Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
    }
    version(turbulence) {
        foreach(i; 0 .. tm.nturb) {
            fs.turb[i] = values[countUntil(varNameList, tm.primitive_variable_name(i))];
        }
    }
    version(multi_species_gas) {
        foreach(i; 0 .. fs.gas.massf.length) {
            fs.gas.massf[i] = values[countUntil(varNameList, massfName(gmodel, to!int(i)))];
        }
        if (fs.gas.massf.length > 1) {
            dt_chem = values[countUntil(varNameList, flowVarName(FlowVar.dt_chem))].re;
        }
    }
    fs.gas.u = values[countUntil(varNameList, flowVarName(FlowVar.u))];
    fs.gas.T = values[countUntil(varNameList, flowVarName(FlowVar.T))];
    version(multi_T_gas) {
        foreach(i; 0 .. fs.gas.u_modes.length) {
            fs.gas.u_modes[i] = values[countUntil(varNameList, u_modesName(to!int(i)))];
            fs.gas.T_modes[i] = values[countUntil(varNameList, T_modesName(to!int(i)))];
        }
        if (fs.gas.u_modes.length > 0) {
            dt_therm = values[countUntil(varNameList, flowVarName(FlowVar.dt_therm))].re;
        }
    }
    if (with_local_time_stepping) { dt_local = values[countUntil(varNameList, flowVarName(FlowVar.dt_local))].re; }
} // end scan_values_from_variable_order_string()

void raw_binary_to_cell_data(ref File fin,
                             ref Vector3 pos, ref number volume, ref FlowState fs,
                             ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
                             bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
                             bool include_quality, bool MHDflag, bool divergence_cleaning,
                             bool radiation, size_t nturb)
{
    // This function needs to be kept consistent with cell_data_to_raw_binary() above.
    //
    version(complex_numbers) {
        // For complex_numbers, we presently set only the real parts.
        // [TODO] Maybe we should read full complex numbers.
        // Fixed-length buffers to hold data while reading binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        fin.rawRead(dbl4);
        pos.set(dbl4[0], dbl4[1], dbl4[2]);
        volume = dbl4[3];
        fin.rawRead(dbl4);
        fs.gas.rho = dbl4[0];
        fs.vel.set(dbl4[1], dbl4[2], dbl4[3]);
        version(MHD) {
            if (MHDflag) {
                fin.rawRead(dbl4);
                fs.B.set(dbl4[0], dbl4[1], dbl4[2]);
                fs.divB = dbl4[3];
                if (divergence_cleaning) {
                    fin.rawRead(dbl1); fs.psi = dbl1[0];
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fin.rawRead(dbl1); fs.gas.quality = dbl1[0];
        } else {
            fs.gas.quality = 1.0;
        }
        fin.rawRead(dbl4);
        fs.gas.p = dbl4[0]; fs.gas.a = dbl4[1]; fs.gas.mu = dbl4[2]; fs.gas.k = dbl4[3];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fin.rawRead(dbl1); fs.gas.k_modes[i] = dbl1[0];
            }
        }
        fin.rawRead(dbl2); fs.mu_t = dbl2[0]; fs.k_t = dbl2[1];
        fin.rawRead(dbl1); fs.S = dbl1[0];
        if (radiation) {
            fin.rawRead(dbl3);
            Q_rad_org = dbl3[0]; f_rad_org = dbl3[1]; Q_rE_rad = dbl3[2];
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        fin.rawRead(dbl2); // tke, omega
        version(turbulence) {
            foreach(i; 0 .. nturb){
                fin.rawRead(dbl1); fs.turb[i] = dbl1[0];
            }
        }
        version(multi_species_gas) {
            foreach (i; 0 .. fs.gas.massf.length) {
                fin.rawRead(dbl1); fs.gas.massf[i] = dbl1[0];
            }
            if (fs.gas.massf.length > 1) { fin.rawRead(dbl1); dt_chem = dbl1[0]; }
        } else {
            fin.rawRead(dbl1); // single-species mass fraction discarded
        }
        fin.rawRead(dbl2);
        fs.gas.u = dbl2[0]; fs.gas.T = dbl2[1];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fin.rawRead(dbl2); fs.gas.u_modes[i] = dbl2[0]; fs.gas.T_modes[i] = dbl2[1];
            }
            if (fs.gas.u_modes.length > 0) { fin.rawRead(dbl1); dt_therm = dbl1[0]; }
        }
        if (with_local_time_stepping) { fin.rawRead(dbl1); dt_local = dbl1[0]; }
    } else {
        // double_numbers
        // Fixed-length buffers to hold data while reading binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        fin.rawRead(dbl4);
        pos.set(dbl4[0], dbl4[1], dbl4[2]);
        volume = dbl4[3];
        fin.rawRead(dbl4);
        fs.gas.rho = dbl4[0];
        fs.vel.set(dbl4[1], dbl4[2], dbl4[3]);
        version(MHD) {
            if (MHDflag) {
                fin.rawRead(dbl4);
                fs.B.set(dbl4[0], dbl4[1], dbl4[2]);
                fs.divB = dbl4[3];
                if (divergence_cleaning) {
                    fin.rawRead(dbl1); fs.psi = dbl1[0];
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fin.rawRead(dbl1); fs.gas.quality = dbl1[0];
        } else {
            fs.gas.quality = 1.0;
        }
        fin.rawRead(dbl4);
        fs.gas.p = dbl4[0]; fs.gas.a = dbl4[1]; fs.gas.mu = dbl4[2]; fs.gas.k = dbl4[3];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fin.rawRead(dbl1); fs.gas.k_modes[i] = dbl1[0];
            }
        }
        fin.rawRead(dbl2); fs.mu_t = dbl2[0]; fs.k_t = dbl2[1];
        fin.rawRead(dbl1); fs.S = dbl1[0];
        if (radiation) {
            fin.rawRead(dbl3);
            Q_rad_org = dbl3[0]; f_rad_org = dbl3[1]; Q_rE_rad = dbl3[2];
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(i; 0 .. nturb){
                fin.rawRead(dbl1); fs.turb[i] = dbl1[0];
            }
        }
        version(multi_species_gas) {
            fin.rawRead(fs.gas.massf);
            if (fs.gas.massf.length > 1) { fin.rawRead(dbl1); dt_chem = dbl1[0]; }
        } else {
            fin.rawRead(dbl1); // single-species mass fraction discarded
        }
        fin.rawRead(dbl2);
        fs.gas.u = dbl2[0]; fs.gas.T = dbl2[1];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fin.rawRead(dbl2); fs.gas.u_modes[i] = dbl2[0]; fs.gas.T_modes[i] = dbl2[1];
            }
            if (fs.gas.u_modes.length > 0) { fin.rawRead(dbl1); dt_therm = dbl1[0]; }
        }
        if (with_local_time_stepping) { fin.rawRead(dbl1); dt_local = dbl1[0]; }
    } // end version double_numbers
} // end raw_binary_to_cell_data()
