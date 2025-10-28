/**
 * fvinterface.d
 * Finite-volume cell-interface class, for use in the CFD codes.
 * Fluxes of conserved quantities are transported (between cells) across cell interfaces.

 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 *          2015-02-13: Keep an eye on the future of the moving_grid option.
 *          2015-05-04: keep references to adjoining cells and defining vertices.
 */

module lmr.fvinterface;

import std.conv;
import std.format;
import std.math;
import std.stdio;

import gas;
import geom;
import nm.bbla;
import nm.number;
import ntypes.complex;

import lmr.conservedquantities;
import lmr.coredata;
import lmr.flowgradients;
import lmr.flowstate;
import lmr.fluidfvcell;
import lmr.fvvertex;
import lmr.globalconfig;
import lmr.lsqinterp;
import lmr.mass_diffusion;

enum IndexDirection {i=0, j, k, none=666}; // Needed for StructuredGrid interpolation.

class FVInterface {
public:
    int id;
    char logical_dir; // logical index direction
    bool is_on_boundary = false;  // by default, assume not on boundary
    size_t bc_id;  // if the face is on a block boundary, which one
    size_t i_bndry; // if the face is on a block boundary, store index into the array of faces attached to bc
    bool use_wall_function_shear_and_heat_flux = false; // for use in viscous_flux_calc()
    bool in_suppress_reconstruction_zone; // if true, we no do reconstruction at this face
    bool in_suppress_viscous_stresses_zone; // if true, we have zero viscous stresses at this face
    //
    // Geometry
    IndexDirection idir;   // For StructuredGrid: in which index-direction is face pointing?
    Vector3 pos;           // position of the (approx) midpoint
    Vector3* gvel;         // grid velocity at interface, m/s
    number Ybar;           // Y-coordinate of the mid-point
    number length;         // Interface length in the x,y-plane
    number[] area;         // Area m**2 for each grid-time-level.
                           // Area per radian in axisymmetric geometry
    Vector3 n;             // Direction cosines for unit normal
    Vector3 t1;            // tangent vector 1 (aka p)
    Vector3 t2;            // tangent vector 2 (aka q)
    FVVertex[] vtx;        // references to vertices for line (2D), quadrilateral and triangular (3D) faces
    //
    // Adjoining cells.
    // These are references to either active cells or ghost cells.
    // The reference may be nil if no cell has been assigned,
    // maybe for a boundary without ghost cells.
    FluidFVCell left_cell;      // interface normal points out of this adjoining cell
    FluidFVCell right_cell;     // interface normal points into this adjoining cell
    // For structured grids, the structure of the grid allows us to carry
    // references for several left- and right- cells.
    FluidFVCell[] left_cells;
    FluidFVCell[] right_cells;
    //
    // Flow
    FVInterfaceData* fvid;
    FlowState* fs;          // Flow properties
    ConservedQuantities F; // Flux conserved quantity per unit area
    number tau_wall_x, tau_wall_y, tau_wall_z; // shear at face (used by wall-function BCs)
    number q;              // heat-flux across face (used by wall-function BCs)
    //
    // Shock-detector-related quantities.
    int[] nbr_id; // list of neighbour ids
    FlowState*[] nbr_fs; // list of neighbouring flow states
    number[] nbr_dist; // distance to neighbour
    //
    // Viscous-flux-related quantities.
    FlowGradients* grad;
    WLSQGradWorkspace* ws_grad;
    Vector3*[] cloud_pos; // Positions of flow points for gradients calculation.
    FlowState*[] cloud_fs; // References to flow states at those points.
    number[] jx; // diffusive mass flux in x
    number[] jy; // diffusive mass flux in y
    number[] jz; // diffusive mass flux in z
    number[] hs; // enthalpies for diffusive flux
    number q_diffusion;
    number q_conduction;
    //
    // Shape sensitivity calculator workspace.
    string global_id;
    version(shape_sensitivity) {
        //string global_id;
	number[][] dFdU;
        // arrays used to temporarily store data intended for the neighbouring block
        // during construction of the external portion of the flow Jacobian.
        size_t[] idList;
        number[] aa;
    }
private:
    LocalConfig myConfig;

public:
    @disable this();
    
    this(LocalConfig myConfig,
         IndexDirection idir,
         FVInterfaceData* fvid,
         int id_init,
         char dir=' ')
    {
        this.myConfig = myConfig;
        this.idir = idir;
        id = id_init;
        logical_dir = dir;
        auto gmodel = myConfig.gmodel;
        uint n_species = myConfig.n_species;
        size_t neq = myConfig.cqi.n;
        size_t ngtl = myConfig.n_grid_time_levels;

        this.fvid = fvid;
        this.fs = &(fvid.flowstates[id]);
        this.grad = &(fvid.gradients[id]);
        this.ws_grad = &(fvid.workspaces[id]);
        this.F = fvid.fluxes[id*neq .. (id+1)*neq];
        this.area = fvid.areas[id*ngtl .. (id+1)*ngtl];
        this.gvel = &(fvid.grid_velocities[id]);
        F.clear();

        version(multi_species_gas) {
            jx.length = n_species;
            jy.length = n_species;
            jz.length = n_species;
            hs.length = n_species;
        }
        version(shape_sensitivity) {
            dFdU.length = 7; // number of conserved variables; FIX-ME for versions
            foreach (ref a; dFdU) a.length = 7;
            foreach (i; 0..dFdU.length) {
                foreach (j; 0..dFdU[i].length) {
                    dFdU[i][j] = 0.0;
                }
            }
        }
        q_diffusion = to!number(0.0);
        q_conduction = to!number(0.0);
    }

    this(FVInterface other, GasModel gm) // not const; see note below
    {
        id = other.id;
        idir = other.idir;
        // We are sort-of promising not to alter the myConfig object,
        // so the rest of the code had better honour that deal...
        myConfig = cast(LocalConfig)other.myConfig;
        is_on_boundary = other.is_on_boundary;
        bc_id = other.bc_id;
        use_wall_function_shear_and_heat_flux = other.use_wall_function_shear_and_heat_flux;
        pos = other.pos;
        gvel = other.gvel;
        Ybar = other.Ybar;
        length = other.length;
        area = other.area.dup;
        n = other.n;
        t1 = other.t1;
        t2 = other.t2;
        tau_wall_x = other.tau_wall_x;
        tau_wall_y = other.tau_wall_y;
        tau_wall_z = other.tau_wall_z;
        q = other.q;
        size_t neq = myConfig.cqi.n;
        this.fvid = other.fvid;
        this.fs = &(fvid.flowstates[id]);
        this.grad = &(fvid.gradients[id]);
        this.ws_grad = &(fvid.workspaces[id]);
        this.F = fvid.fluxes[id*neq .. (id+1)*neq];
        // Because we copy the following pointers and references,
        // we cannot have const (or "in") qualifier on other.
        cloud_pos = other.cloud_pos.dup();
        cloud_fs = other.cloud_fs.dup();
        version(multi_species_gas) {
            jx = other.jx.dup();
            jy = other.jy.dup();
            jz = other.jz.dup();
            hs = other.hs.dup();
        }
        // FIX-ME - KYLE -- there's a typo here. Not sure the intent.
        version(steadystate) {
            dFdU_L.length = 5; // number of conserved variables; FIX-ME for versions
            foreach (ref a; dFdU_L) a.length = 5;
            dFdU_R.length = 5;
            foreach (ref a; dFdU_R) a.length = 5;
        }
        q_diffusion = other.q_diffusion;
        q_conduction = other.q_conduction;
    }

    @nogc
    void setLocalConfig(LocalConfig lcfg)
    {
        myConfig = lcfg;
    }
    @nogc
    void copy_values_from(in FVInterface other, uint type_of_copy)
    {
        switch (type_of_copy) {
        case CopyDataOption.minimal_flow:
        case CopyDataOption.all_flow:
            fs.copy_values_from(other.fs);
            F.copy_values_from(other.F);
            tau_wall_x = other.tau_wall_x;
            tau_wall_y = other.tau_wall_y;
            tau_wall_z = other.tau_wall_z;
            q = other.q;
            q_diffusion = other.q_diffusion;
            break;
        case CopyDataOption.grid:
            pos.set(other.pos);
            gvel.set(other.gvel);
            Ybar = other.Ybar;
            length = other.length;
            area[] = other.area[];
            n.set(other.n); t1.set(other.t1); t2.set(other.t2);
            break;
        case CopyDataOption.all:
        default:
            id = other.id;
            idir = other.idir;
            // We are sort-of promising not to alter the myConfig object,
            // so the rest of the code had better honour that deal...
            myConfig = cast(LocalConfig)other.myConfig;
            pos.set(other.pos);
            gvel.set(other.gvel);
            Ybar = other.Ybar;
            length = other.length;
            area[] = other.area[];
            n.set(other.n); t1.set(other.t1); t2.set(other.t2);
            fs.copy_values_from(other.fs);
            F.copy_values_from(other.F);
            tau_wall_x = other.tau_wall_x;
            tau_wall_y = other.tau_wall_y;
            tau_wall_z = other.tau_wall_z;
            q = other.q;
            grad.copy_values_from(*(other.grad));
            q_diffusion = other.q_diffusion;
            // omit scratch workspace ws_grad
        } // end switch
    }

    @nogc
    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
        area[to_level] = area[from_level];
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "FVInterface(";
        repr ~= "id=" ~ to!string(id);
        repr ~= ", idir=" ~ to!string(idir);
        repr ~= ", universe_blk_id=" ~ to!string(myConfig.universe_blk_id);
        repr ~= ", pos=" ~ to!string(pos);
        repr ~= ", vtx_ids=[";
        foreach (v; vtx) { repr ~= format("%d,", v.id); }
        repr ~= "]";
        repr ~= format(", left_cell_id=%d", left_cell ? left_cell.id : -1);
        repr ~= format(", right_cell_id=%d", right_cell ? right_cell.id : -1);
        repr ~= ", gvel=" ~ to!string(gvel);
        repr ~= ", Ybar=" ~ to!string(Ybar);
        repr ~= ", length=" ~ to!string(length);
        repr ~= ", area=" ~ to!string(area);
        repr ~= ", n=" ~ to!string(n);
        repr ~= ", t1=" ~ to!string(t1);
        repr ~= ", t2=" ~ to!string(t2);
        repr ~= ", fs=" ~ to!string(fs);
        repr ~= ", tau_wall_x=" ~ to!string(tau_wall_x);
        repr ~= ", tau_wall_y=" ~ to!string(tau_wall_y);
        repr ~= ", tau_wall_z=" ~ to!string(tau_wall_z);
        repr ~= ", q=" ~ to!string(q);
        repr ~= ", F=" ~ to!string(F);
        repr ~= ", E=" ~ to!string(fs.electric_field);
        repr ~= ", grad=" ~ to!string(grad);
        repr ~= ", cloud_pos=[";
        // Because the positions are held as pointers to Vector3 objects,
        // we will get more interesting information by printing the objects
        // that they point to.
        foreach (i, vtxptr; cloud_pos) {
            if (i > 0) { repr ~= ", "; }
            repr ~= to!string(*vtxptr);
        }
        repr ~= "]";
        repr ~= ", cloud_fs=" ~ to!string(cloud_fs);
        repr ~= ")";
        return to!string(repr);
    }

    @nogc
    void update_2D_geometric_data(size_t gtl, bool axisymmetric)
    {
        number xA = vtx[0].pos[gtl].x;
        number yA = vtx[0].pos[gtl].y;
        number xB = vtx[1].pos[gtl].x;
        number yB = vtx[1].pos[gtl].y;
        number LAB = sqrt((xB-xA)*(xB-xA) + (yB-yA)*(yB-yA));
        // Direction cosines for the unit normal and two tangential directions.
        if (LAB > 1.0e-12) {
            // Normal is purely in the xy-plane, pointing to the "right"
            // as we sit at A, looking toward B.
            n.set((yB-yA)/LAB, -(xB-xA)/LAB, to!number(0.0));
            t2 = Vector3(0.0, 0.0, 1.0);
            cross(t1, n, t2);
            length = LAB; // Length in the XY-plane.
        } else {
            // A and B coincide.
            n = Vector3(1.0, 0.0, 0.0); // Arbitrary direction
            t2 = Vector3(0.0, 0.0, 1.0);
            t1 = Vector3(0.0, 1.0, 0.0);
            length = 0.0; // Zero length in the xy-plane
        }
        // Mid-point and surface area.
        number Xbar = 0.5*(xA+xB);
        Ybar = 0.5*(yA+yB);
        if (axisymmetric) {
            area[gtl] = length * Ybar; // Face area per radian.
        } else {
            area[gtl] = length; // Assume unit depth in the Z-direction.
        }
        pos.set(Xbar, Ybar, to!number(0.0));
        fvid.normals[id] = n;
        fvid.tangents1[id] = t1;
        fvid.tangents2[id] = t2;
        fvid.positions[id] = pos;
    } // end update_2D_geometric_data()

    @nogc
    void update_3D_geometric_data(size_t gtl)
    {
        switch (vtx.length) {
        case 3:
            triangle_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                                vtx[2].pos[gtl],
                                pos, n, t1, t2, area[gtl]);
            length = sqrt(area[gtl]);
            break;
        case 4:
            quad_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                            vtx[2].pos[gtl], vtx[3].pos[gtl],
                            pos, n, t1, t2, area[gtl]);
            length = sqrt(area[gtl]);
            break;
        default:
            string msg = "FVInterface.update_3D_geometric_data(): Unhandled number of vertices: ";
            debug { msg ~= format("%d", vtx.length); }
            throw new FlowSolverException(msg);
        } // end switch
        fvid.normals[id] = n;
        fvid.tangents1[id] = t1;
        fvid.tangents2[id] = t2;
        fvid.positions[id] = pos;
    } // end update_3D_geometric_data()

    @nogc
    number upwind_weighting(number M)
    {
        // Weighting function defined on page 78 of Ian Johnston's thesis.
        // M is Mach number in direction pointing toward the location
        // at which we are evaluating the weighted sum.
        if (M > 1.0) return to!number(1.0);
        number Mp1 = M + 1.0;
        return (Mp1 > 0.0) ? (0.25*Mp1^^2) : to!number(0.0);
    }

    @nogc
    void average_vertex_deriv_values()
    {
        number[4] uwfs;
        number uwf_sum = 0.0;
        if (myConfig.upwind_vertex_gradients) {
            // Upwind weighting of the gradient values.
            if (!left_cell && !right_cell) {
                throw new Exception("Oops! This face does not have at least one cell attached.");
            }
            FluidFVCell c = (right_cell && right_cell.is_interior_to_domain) ? right_cell : left_cell;
            Vector3 tangent;
            foreach (i; 0 .. vtx.length) {
                tangent.set(pos); tangent -= vtx[i].pos[0]; tangent.normalize();
                number M = c.fs.vel.dot(tangent)/c.fs.gas.a;
                number uwf = upwind_weighting(M);
                uwfs[i] = uwf;
                uwf_sum += uwf;
            }
        } else {
            // Equal weighting of the gradient values.
            foreach (i; 0 .. vtx.length) {
                uwfs[i] = 1.0;
                uwf_sum += 1.0;
            }
        }
        grad.copy_values_from(*(vtx[0].grad));
        grad.scale_values_by(uwfs[0]);
        foreach (i; 1 .. vtx.length) {
            grad.accumulate_values_from(*(vtx[i].grad), uwfs[i]);
        }
        grad.scale_values_by(1.0/uwf_sum);
    } // end average_vertex_deriv_values()

    @nogc
    void average_cell_deriv_values(int gtl)
    {
        if (!left_cell && !right_cell) {
            throw new Exception("Oops! This face does not have at least one cell attached.");
        }

        // helper function for augmenting the averaged face gradients
        pragma(inline, true)
        @safe @nogc nothrow pure
        void augment_gradients(scope const ref number[3] avg_grad, in number u_i, in number u_j,
                               scope const ref number[3] n_tilde, scope const ref number[3] dp_ij_hat,
                               in number dp_ij_mag, ref number[3] aug_grad) {
            // compute component of gradient in the direction of dp_ij
            const number avg_grad_dot_dp_ij_hat = avg_grad[0]*dp_ij_hat[0] + avg_grad[1]*dp_ij_hat[1] + avg_grad[2]*dp_ij_hat[2];
            // compute finite difference gradient in direction of dp_ij
            const number dudp = (u_j - u_i)/dp_ij_mag;
            // augment gradients
            aug_grad[] = avg_grad[] - (avg_grad_dot_dp_ij_hat - dudp)*n_tilde[];
        }

        // determine if there is a valid gradient in the ghost cell
        bool one_interior_cell_exists = false;
        if ( (left_cell && !right_cell) ||
             (!left_cell && right_cell) ) {
            one_interior_cell_exists = true;
        } else { // both a left cell and right cell exists => check if they are interior cells
            if ( (left_cell.is_interior_to_domain && !right_cell.is_interior_to_domain) ||
                 (!left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) ) {
                one_interior_cell_exists = true;
            }
        }

        const bool apply_boundary_correction = myConfig.include_boundary_faces_in_spatial_deriv_correction;
        if (one_interior_cell_exists) {
            // The interface has only one cell interior to the domain and that doesn't have a mapping to a cell in a neighbouring block.
            FluidFVCell c = (right_cell && right_cell.is_interior_to_domain) ? right_cell : left_cell; // i --> the boundary interface is j

            if (apply_boundary_correction) {
                // compute geometric terms
                number[3] dp_ij, dp_ij_hat, n_tilde;
                dp_ij[0] = pos.x - c.pos[gtl].x;
                dp_ij[1] = pos.y - c.pos[gtl].y;
                dp_ij[2] = pos.z - c.pos[gtl].z;
                const number dp_ij_mag = sqrt(dp_ij[0]*dp_ij[0] + dp_ij[1]*dp_ij[1] + dp_ij[2]*dp_ij[2]);
                dp_ij_hat[] = dp_ij[]/dp_ij_mag;
                const number n_dot_dp_ij_hat = n.x*dp_ij_hat[0] + n.y*dp_ij_hat[1] + n.z*dp_ij_hat[2];
                n_tilde[0] = n.x/n_dot_dp_ij_hat;
                n_tilde[1] = n.y/n_dot_dp_ij_hat;
                n_tilde[2] = n.z/n_dot_dp_ij_hat;

                // x-velocity
                augment_gradients(c.grad.vel[0], c.fs.vel.x, fs.vel.x, n_tilde, dp_ij_hat, dp_ij_mag, grad.vel[0]);
                // y-velocity
                augment_gradients(c.grad.vel[1], c.fs.vel.y, fs.vel.y, n_tilde, dp_ij_hat, dp_ij_mag, grad.vel[1]);
                // z-velocity
                augment_gradients(c.grad.vel[2], c.fs.vel.z, fs.vel.z, n_tilde, dp_ij_hat, dp_ij_mag, grad.vel[2]);
                // massf
                version(multi_species_gas) {
                    foreach (isp; 0 .. myConfig.n_species) {
                        augment_gradients(c.grad.massf[isp], c.fs.gas.massf[isp], fs.gas.massf[isp], n_tilde, dp_ij_hat, dp_ij_mag, grad.massf[isp]);
                    }
                }
                // temperature
                augment_gradients(c.grad.T, c.fs.gas.T, fs.gas.T, n_tilde, dp_ij_hat, dp_ij_mag, grad.T);
                // thermal modes
                version(multi_T_gas) {
                    uint n_modes = myConfig.n_modes;
                    foreach (imode; 0 .. n_modes) {
                        augment_gradients(c.grad.T_modes[imode], c.fs.gas.T_modes[imode], fs.gas.T_modes[imode], n_tilde, dp_ij_hat, dp_ij_mag, grad.T_modes[imode]);
                    }
                }
                // turbulence
                version(turbulence) {
                    foreach(i; 0 .. myConfig.turb_model.nturb) {
                        augment_gradients(c.grad.turb[i], c.fs.turb[i], fs.turb[i], n_tilde, dp_ij_hat, dp_ij_mag, grad.turb[i]);
                    }
                }
            } else { // just copy the cell gradients
                // velocities
                grad.vel[0][] = c.grad.vel[0][];
                grad.vel[1][] = c.grad.vel[1][];
                grad.vel[2][] = c.grad.vel[2][];
                // species mass fractions
                version(multi_species_gas) {
                    foreach (isp; 0 .. myConfig.n_species) {
                        grad.massf[isp][] = c.grad.massf[isp][];
                    }
                }
                // temperature
                grad.T[] = c.grad.T[];
                // multi-temperature modes
                version(multi_T_gas) {
                    uint n_modes = myConfig.n_modes;
                    foreach (imode; 0 .. n_modes) {
                        grad.T_modes[imode][] = c.grad.T_modes[imode][];
                    }
                }
                // turbulence variables
                version(turbulence) {
                    foreach (i; 0 .. myConfig.turb_model.nturb) {
                        grad.turb[i][] = c.grad.turb[i][];
                    }
                }
            }
        } else {
            // With two attached cells, we are at a face that is internal to the domain and so we can proceed to compute the average of the gradient values.
            FluidFVCell cL0 = left_cell;  // i
            FluidFVCell cR0 = right_cell; // j

            // compute geometric terms
            number[3] dp_ij, dp_ij_hat, n_tilde;
            dp_ij[0] = cR0.pos[gtl].x - cL0.pos[gtl].x;
            dp_ij[1] = cR0.pos[gtl].y - cL0.pos[gtl].y;
            dp_ij[2] = cR0.pos[gtl].z - cL0.pos[gtl].z;
            const number dp_ij_mag = sqrt(dp_ij[0]*dp_ij[0] + dp_ij[1]*dp_ij[1] + dp_ij[2]*dp_ij[2]);
            dp_ij_hat[] = dp_ij[]/dp_ij_mag;
            const number n_dot_dp_ij_hat = n.x*dp_ij_hat[0] + n.y*dp_ij_hat[1] + n.z*dp_ij_hat[2];
            n_tilde[0] = n.x/n_dot_dp_ij_hat;
            n_tilde[1] = n.y/n_dot_dp_ij_hat;
            n_tilde[2] = n.z/n_dot_dp_ij_hat;

            number[3] avg_grad;
            number HALF = 0.5;
            // x-velocity
            avg_grad[] = HALF*(cL0.grad.vel[0][]+cR0.grad.vel[0][]);
            augment_gradients(avg_grad, cL0.fs.vel.x, cR0.fs.vel.x, n_tilde, dp_ij_hat, dp_ij_mag, grad.vel[0]);
            // y-velocity
            avg_grad[] = HALF*(cL0.grad.vel[1][]+cR0.grad.vel[1][]);
            augment_gradients(avg_grad, cL0.fs.vel.y, cR0.fs.vel.y, n_tilde, dp_ij_hat, dp_ij_mag, grad.vel[1]);
            // z-velocity
            avg_grad[] = HALF*(cL0.grad.vel[2][]+cR0.grad.vel[2][]);
            augment_gradients(avg_grad, cL0.fs.vel.z, cR0.fs.vel.z, n_tilde, dp_ij_hat, dp_ij_mag, grad.vel[2]);
            // massf
            version(multi_species_gas) {
                foreach (isp; 0 .. myConfig.n_species) {
                    avg_grad[] = HALF*(cL0.grad.massf[isp][]+cR0.grad.massf[isp][]);
                    augment_gradients(avg_grad, cL0.fs.gas.massf[isp], cR0.fs.gas.massf[isp], n_tilde, dp_ij_hat, dp_ij_mag, grad.massf[isp]);
                }
            }
            // temperature
            avg_grad[] = HALF*(cL0.grad.T[]+cR0.grad.T[]);
            augment_gradients(avg_grad, cL0.fs.gas.T, cR0.fs.gas.T, n_tilde, dp_ij_hat, dp_ij_mag, grad.T);
            // thermal modes
            version(multi_T_gas) {
                uint n_modes = myConfig.n_modes;
                foreach (imode; 0 .. n_modes) {
                    avg_grad[] = HALF*(cL0.grad.T_modes[imode][]+cR0.grad.T_modes[imode][]);
                    augment_gradients(avg_grad, cL0.fs.gas.T_modes[imode], cR0.fs.gas.T_modes[imode], n_tilde, dp_ij_hat, dp_ij_mag, grad.T_modes[imode]);
                }
            }
            // turbulence
            version(turbulence) {
                foreach(i; 0 .. myConfig.turb_model.nturb) {
                    avg_grad[] = HALF*(cL0.grad.turb[i][]+cR0.grad.turb[i][]);
                    augment_gradients(avg_grad, cL0.fs.turb[i], cR0.fs.turb[i], n_tilde, dp_ij_hat, dp_ij_mag, grad.turb[i]);
                }
            }
        }
    } // end average_cell_deriv_values()

    @nogc
    void average_electric_field()
    {
        if (left_cell && right_cell) {
            if (!left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
                fs.electric_field[0] = right_cell.electric_field[0];
                fs.electric_field[1] = right_cell.electric_field[1];
            }
            else if (left_cell.is_interior_to_domain && !right_cell.is_interior_to_domain) {
                fs.electric_field[0] = left_cell.electric_field[0];
                fs.electric_field[1] = left_cell.electric_field[1];
            }
            else {
                if (isNaN(left_cell.electric_field[0])) {
                    fs.electric_field[0] = right_cell.electric_field[0];
                    fs.electric_field[1] = right_cell.electric_field[1];
                }
                else if (isNaN(right_cell.electric_field[0])) {
                    fs.electric_field[0] = left_cell.electric_field[0];
                    fs.electric_field[1] = left_cell.electric_field[1];
                }
                else {
                    fs.electric_field[0] = 0.5 * (left_cell.electric_field[0] +
                                               right_cell.electric_field[0]);
                    fs.electric_field[1] = 0.5 * (left_cell.electric_field[1] +
                                               right_cell.electric_field[1]);               
                }
            }

        }
        else if (left_cell && !right_cell) {
            fs.electric_field[0] = left_cell.electric_field[0];
            fs.electric_field[1] = left_cell.electric_field[1];
        }
        else if (!left_cell && right_cell) {
            fs.electric_field[0] = right_cell.electric_field[0];
            fs.electric_field[1] = right_cell.electric_field[1];
        }
        return;
        // throw new Exception("Oops! This face does not have at least one cell attached.");
    }

    @nogc
    void average_turbulent_transprops()
    {
        if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
            fs.k_t = 0.5*(left_cell.fs.k_t+right_cell.fs.k_t);
            fs.mu_t = 0.5*(left_cell.fs.mu_t+right_cell.fs.mu_t);
        } else if (left_cell && left_cell.is_interior_to_domain) {
            fs.k_t = left_cell.fs.k_t;
            fs.mu_t = left_cell.fs.mu_t;
        } else if (right_cell && right_cell.is_interior_to_domain) {
            fs.k_t = right_cell.fs.k_t;
            fs.mu_t = right_cell.fs.mu_t;
        } else {
            assert(0, "Oops, don't seem to have a cell available.");
        }
    }

    @nogc
    void viscous_flux_calc()
    // Unified 2D and 3D viscous-flux calculation.
    // Note that the gradient values need to be in place before calling this procedure.
    // Note, also, that the viscous fluxes are added to the flux-vector components.
    {
        if (in_suppress_viscous_stresses_zone) {
            // We wish to ignore the viscous fluxes here.
            return;
        }
        auto gmodel = myConfig.gmodel;
        uint n_species = myConfig.n_species;
        uint n_modes = myConfig.n_modes;
        double viscous_factor = myConfig.viscous_factor;
        number k_eff = fs.gas.k + fs.k_t ;
        number mu_eff = fs.gas.mu + fs.mu_t;
        number lmbda;
        lmbda = -2.0/3.0 * mu_eff;

        number local_pressure = fs.gas.p;
        number shear_stress_limit = myConfig.shear_stress_relative_limit * local_pressure;
        number heat_transfer_limit = (mu_eff > 0.0) ? k_eff/mu_eff*shear_stress_limit : to!number(0.0);

        // Species diffusion: Changed by NNG on 22/01/18.
        // We now apply both laminar and turbulent diffusion additively, to prevent artificially low
        // diffusion in areas with a small turbulent viscosity.
        version(multi_species_gas) {
            if (myConfig.mass_diffusion_model != MassDiffusionModel.none) {
                myConfig.massDiffusion.update_mass_fluxes(*fs, *grad, jx, jy, jz);
                foreach (isp; 0 .. n_species) {
                    jx[isp] *= viscous_factor;
                    jy[isp] *= viscous_factor;
                    jz[isp] *= viscous_factor;
                }
            } else if (myConfig.turb_model.isTurbulent) { // Turbulent but no mass diffusion model
                foreach (isp; 0 .. n_species) {
                    jx[isp] = to!number(0.0);
                    jy[isp] = to!number(0.0);
                    jz[isp] = to!number(0.0);
                }
            }
        }

        if (myConfig.turb_model.isTurbulent) {
            double Sc_t = myConfig.turbulence_schmidt_number;
            number D_t = fs.mu_t / (fs.gas.rho * Sc_t);
            version(multi_species_gas) {
                foreach (isp; 0 .. n_species) {
                    jx[isp] -= fs.gas.rho * D_t * grad.massf[isp][0];
                    jy[isp] -= fs.gas.rho * D_t * grad.massf[isp][1];
                    jz[isp] -= fs.gas.rho * D_t * grad.massf[isp][2];
                }
            }
        }
        number tau_xx = 0.0;
        number tau_yy = 0.0;
        number tau_zz = 0.0;
        number tau_xy = 0.0;
        number tau_xz = 0.0;
        number tau_yz = 0.0;
        if (myConfig.dimensions == 3) {
            number dudx = grad.vel[0][0];
            number dudy = grad.vel[0][1];
            number dudz = grad.vel[0][2];
            number dvdx = grad.vel[1][0];
            number dvdy = grad.vel[1][1];
            number dvdz = grad.vel[1][2];
            number dwdx = grad.vel[2][0];
            number dwdy = grad.vel[2][1];
            number dwdz = grad.vel[2][2];
            // 3-dimensional planar stresses.
            tau_xx = 2.0*mu_eff*dudx + lmbda*(dudx + dvdy + dwdz);
            tau_yy = 2.0*mu_eff*dvdy + lmbda*(dudx + dvdy + dwdz);
            tau_zz = 2.0*mu_eff*dwdz + lmbda*(dudx + dvdy + dwdz);
            tau_xy = mu_eff * (dudy + dvdx);
            tau_xz = mu_eff * (dudz + dwdx);
            tau_yz = mu_eff * (dvdz + dwdy);
        } else {
            // 2D
            number dudx = grad.vel[0][0];
            number dudy = grad.vel[0][1];
            number dvdx = grad.vel[1][0];
            number dvdy = grad.vel[1][1];
            if (myConfig.axisymmetric) {
                // Viscous stresses at the mid-point of the interface.
                // Axisymmetric terms no longer include the radial multiplier
                // as that has been absorbed into the interface area calculation.
                number ybar = Ybar;
                if (ybar > 1.0e-10) { // something very small for a cell height
                    tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                    tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                } else {
                    tau_xx = 0.0;
                    tau_yy = 0.0;
                }
                tau_xy = mu_eff * (dudy + dvdx);
            } else {
                // 2-dimensional-planar stresses.
                tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy);
                tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy);
                tau_xy = mu_eff * (dudy + dvdx);
            }
        }
        // Thermal conductivity (NOTE: q is total energy flux)
        number qx = k_eff * grad.T[0];
        number qy = k_eff * grad.T[1];
        number qz = k_eff * grad.T[2];
        version(multi_T_gas) {
            foreach (imode; 0 .. n_modes) {
                qx += viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][0];
                qy += viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][1];
                qz += viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][2];
            }
        }
        q_conduction = (qx*n.x + qy*n.y + qz*n.z);
        version(multi_species_gas) {
            if (myConfig.turb_model.isTurbulent ||
                myConfig.mass_diffusion_model != MassDiffusionModel.none ) {
                q_diffusion = to!number(0.0);
                gmodel.enthalpies(fs.gas, hs);
                foreach (isp; 0 .. n_species) {
                    qx -= jx[isp] * hs[isp];
                    qy -= jy[isp] * hs[isp];
                    qz -= jz[isp] * hs[isp];
                    q_diffusion -= (jx[isp]*hs[isp]*n.x + jy[isp]*hs[isp]*n.y + jz[isp]*hs[isp]*n.z);
                }
            }
        }
        version(turbulence) {
            if ( myConfig.turb_model.isTurbulent &&
                 !(myConfig.axisymmetric && (Ybar <= 1.0e-10)) ) {
                // Turbulence contribution to the shear stresses.
                number tke = myConfig.turb_model.turbulent_kinetic_energy(*fs);
                tau_xx -= 2.0/3.0 * fs.gas.rho * tke;
                tau_yy -= 2.0/3.0 * fs.gas.rho * tke;
                if (myConfig.dimensions == 3) { tau_zz -= 2.0/3.0 * fs.gas.rho * tke; }

                // Turbulent transport of turbulent kinetic energy
                number[3] qtke = myConfig.turb_model.turbulent_kinetic_energy_transport(*fs, *grad);
                qx += qtke[0];
                qy += qtke[1];
                if (myConfig.dimensions == 3) { qz += qtke[2]; }
            }
        }
        if (myConfig.apply_shear_stress_relative_limit) {
            version(complex_numbers) {
                // Do not try to limit the component values.
                // Something in this limiting plays havoc with the complex derivatives.
            } else {
                // Apply limits to the component values.
                tau_xx = copysign(fmin(fabs(tau_xx),shear_stress_limit), tau_xx);
                tau_yy = copysign(fmin(fabs(tau_yy),shear_stress_limit), tau_yy);
                tau_zz = copysign(fmin(fabs(tau_zz),shear_stress_limit), tau_zz);
                tau_xy = copysign(fmin(fabs(tau_xy),shear_stress_limit), tau_xy);
                tau_xz = copysign(fmin(fabs(tau_xz),shear_stress_limit), tau_xz);
                tau_yz = copysign(fmin(fabs(tau_yz),shear_stress_limit), tau_yz);
                qx = copysign(fmin(fabs(qx),heat_transfer_limit), qx);
                qy = copysign(fmin(fabs(qy),heat_transfer_limit), qy);
                qz = copysign(fmin(fabs(qz),heat_transfer_limit), qz);
            }
        } // end if apply_shear_stress_relative_limit
        //
        // Combine into fluxes: store as the dot product (F.n).
        number nx = n.x;
        number ny = n.y;
        number nz = n.z;
        auto cqi = myConfig.cqi;
        // In some cases, the shear and heat fluxes have been previously
        // computed by the wall functions in the boundary condition call.
        if (use_wall_function_shear_and_heat_flux) {
            // Mass flux -- NO CONTRIBUTION, unless there's diffusion (below)
            // [TODO] As per Jason's recommendation, we need to do something
            // to correct for corner cells.
            // [TODO] Currently implemented for 2D; need to extend to 3D.
            F[cqi.xMom] -= tau_xx*nx + tau_wall_x;
            F[cqi.yMom] -= tau_yy*ny + tau_wall_y;
            if (cqi.threeD) { F[cqi.zMom] -= tau_zz*nz + tau_wall_z; }
            F[cqi.totEnergy] -=
                tau_xx*fs.vel.x*nx + tau_yy*fs.vel.y*ny + tau_zz*fs.vel.z*nz +
                tau_wall_x*fs.vel.x + tau_wall_y*fs.vel.y + tau_wall_z*fs.vel.z + q;
        }
        else { // proceed with locally computed shear and heat flux
            // Mass flux -- NO CONTRIBUTION, unless there's diffusion (below)
            F[cqi.xMom] -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
            F[cqi.yMom] -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
            if (cqi.threeD) { F[cqi.zMom] -= tau_xz*nx + tau_yz*ny + tau_zz*nz; }
            F[cqi.totEnergy] -=
                (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx)*nx +
                (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy)*ny +
                (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz)*nz;
        } // end if
        version(multi_T_gas) {
            foreach (imode; 0 .. n_modes) {
                F[cqi.modes+imode] -= viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][0] * nx;
                F[cqi.modes+imode] -= viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][1] * ny;
                F[cqi.modes+imode] -= viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][2] * nz;
                // Species diffusion contribution to the energy modes (added by NNG, 2022/01/26)
                version(multi_species_gas) {
                    if (myConfig.turb_model.isTurbulent || (myConfig.mass_diffusion_model != MassDiffusionModel.none)) {
                        foreach (isp; 0 .. n_species) {
                            number hMode = gmodel.enthalpyPerSpeciesInMode(fs.gas, cast(int)isp, cast(int)imode);
                            // The sign here needs to be opposite to the thermal conduction, hence +=
                            F[cqi.modes+imode] += viscous_factor * hMode *(jx[isp]*n.x
                                                                             + jy[isp]*n.y
                                                                             + jz[isp]*n.z);
                        }
                    }
                }
            }
        }
        version(turbulence) {
            if ( myConfig.turb_model.isTurbulent &&
                 !(myConfig.axisymmetric && (Ybar <= 1.0e-10)) ) {
                //
                // Turbulence transport of the turbulence properties themselves.
                foreach(i; 0 .. myConfig.turb_model.nturb){
                    number tau_tx = 0.0;
                    number tau_ty = 0.0;
                    number tau_tz = 0.0;
                    //
                    number mu_effective = myConfig.turb_model.viscous_transport_coeff(*fs, i);
                    // Apply a limit on mu_effective in the same manner as that applied to mu_t.
                    mu_effective = fmin(mu_effective, myConfig.max_mu_t_factor * fs.gas.mu);
                    tau_tx = mu_effective * grad.turb[i][0];
                    tau_ty = mu_effective * grad.turb[i][1];
                    if (myConfig.dimensions == 3) { tau_tz = mu_effective * grad.turb[i][2]; }
                    //
                    F[cqi.rhoturb+i] -= tau_tx * nx + tau_ty * ny + tau_tz * nz;
                }
            }
        }
        version(multi_species_gas) {
            if (myConfig.turb_model.isTurbulent ||
                myConfig.mass_diffusion_model != MassDiffusionModel.none) {
                if (cqi.n_species > 1) {
                    foreach (isp; 0 .. cqi.n_species) {
                        F[cqi.species+isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
                    }
                }
            }
        }
    } // end viscous_flux_calc()

} // end of class FV_Interface

@nogc
void navier_stokes_viscous_fluxes(in FlowState fs, in FlowGradients grad, LocalConfig myConfig, size_t n_modes, size_t nturb,
                                  bool isTurbulent, bool axisymmetric, bool is3d, Vector3 n, double Ybar,
                                  ConservedQuantities F)
// Unified 2D and 3D viscous-flux calculation.
// Note that the gradient values need to be in place before calling this procedure.
// Note, also, that the viscous fluxes are added to the flux-vector components.
{
    number k_eff = fs.gas.k + fs.k_t ;
    number mu_eff = fs.gas.mu + fs.mu_t;
    number lmbda;
    lmbda = -2.0/3.0 * mu_eff;

    number tau_xx = 0.0;
    number tau_yy = 0.0;
    number tau_zz = 0.0;
    number tau_xy = 0.0;
    number tau_xz = 0.0;
    number tau_yz = 0.0;
    if (is3d) {
        number dudx = grad.vel[0][0];
        number dudy = grad.vel[0][1];
        number dudz = grad.vel[0][2];
        number dvdx = grad.vel[1][0];
        number dvdy = grad.vel[1][1];
        number dvdz = grad.vel[1][2];
        number dwdx = grad.vel[2][0];
        number dwdy = grad.vel[2][1];
        number dwdz = grad.vel[2][2];
        // 3-dimensional planar stresses.
        tau_xx = 2.0*mu_eff*dudx + lmbda*(dudx + dvdy + dwdz);
        tau_yy = 2.0*mu_eff*dvdy + lmbda*(dudx + dvdy + dwdz);
        tau_zz = 2.0*mu_eff*dwdz + lmbda*(dudx + dvdy + dwdz);
        tau_xy = mu_eff * (dudy + dvdx);
        tau_xz = mu_eff * (dudz + dwdx);
        tau_yz = mu_eff * (dvdz + dwdy);
    } else {
        // 2D
        number dudx = grad.vel[0][0];
        number dudy = grad.vel[0][1];
        number dvdx = grad.vel[1][0];
        number dvdy = grad.vel[1][1];
        if (axisymmetric) {
            // Viscous stresses at the mid-point of the interface.
            // Axisymmetric terms no longer include the radial multiplier
            // as that has been absorbed into the interface area calculation.
            number ybar = Ybar;
            if (ybar > 1.0e-10) { // something very small for a cell height
                tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy + fs.vel.y / ybar);
            } else {
                tau_xx = 0.0;
                tau_yy = 0.0;
            }
            tau_xy = mu_eff * (dudy + dvdx);
        } else {
            // 2-dimensional-planar stresses.
            tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy);
            tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy);
            tau_xy = mu_eff * (dudy + dvdx);
        }
    }
    // Thermal conductivity (NOTE: q is total energy flux)
    number qx = k_eff * grad.T[0];
    number qy = k_eff * grad.T[1];
    number qz = k_eff * grad.T[2];
    version(multi_T_gas) {
        foreach (imode; 0 .. n_modes) {
            qx += fs.gas.k_modes[imode] * grad.T_modes[imode][0];
            qy += fs.gas.k_modes[imode] * grad.T_modes[imode][1];
            qz += fs.gas.k_modes[imode] * grad.T_modes[imode][2];
        }
    }
    version(turbulence) {
        if ( isTurbulent ) {
            // Turbulence contribution to the shear stresses.
            number tke = myConfig.turb_model.turbulent_kinetic_energy(fs);
            tau_xx -= 2.0/3.0 * fs.gas.rho * tke;
            tau_yy -= 2.0/3.0 * fs.gas.rho * tke;
            if (is3d) { tau_zz -= 2.0/3.0 * fs.gas.rho * tke; }

            // Turbulent transport of turbulent kinetic energy
            number[3] qtke = myConfig.turb_model.turbulent_kinetic_energy_transport(fs, grad);
            qx += qtke[0];
            qy += qtke[1];
            if (is3d) { qz += qtke[2]; }
        }
    }
    //
    // Combine into fluxes: store as the dot product (F.n).
    auto cqi = myConfig.cqi;
    number nx = n.x;
    number ny = n.y;
    number nz = n.z;
    F[cqi.xMom] -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
    F[cqi.yMom] -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
    if (cqi.threeD) { F[cqi.zMom] -= tau_xz*nx + tau_yz*ny + tau_zz*nz; }
    F[cqi.totEnergy] -=
        (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx)*nx +
        (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy)*ny +
        (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz)*nz;

    version(multi_T_gas) {
        foreach (imode; 0 .. n_modes) {
            F[cqi.modes+imode] -= fs.gas.k_modes[imode] * grad.T_modes[imode][0] * nx;
            F[cqi.modes+imode] -= fs.gas.k_modes[imode] * grad.T_modes[imode][1] * ny;
            F[cqi.modes+imode] -= fs.gas.k_modes[imode] * grad.T_modes[imode][2] * nz;
        }
    }
    version(turbulence) {
        if ( isTurbulent ) {
            // Turbulence transport of the turbulence properties themselves.
            foreach(i; 0 .. nturb){
                number tau_tx = 0.0;
                number tau_ty = 0.0;
                number tau_tz = 0.0;
                //
                number mu_effective = myConfig.turb_model.viscous_transport_coeff(fs, i);
                tau_tx = mu_effective * grad.turb[i][0];
                tau_ty = mu_effective * grad.turb[i][1];
                if (is3d) { tau_tz = mu_effective * grad.turb[i][2]; }
                //
                F[cqi.rhoturb+i] -= tau_tx * nx + tau_ty * ny + tau_tz * nz;
            }
        }
    }
} // end viscous_flux_calc()

@nogc
void diffusion_viscous_fluxes(in FlowState fs, in FlowGradients grad, LocalConfig myConfig, size_t n_species, size_t n_modes,
                              double Sc_t, bool laminarDiffusion, bool isTurbulent, Vector3 n, number[] jx, number[] jy, number[] jz, number[] hs,
                              ConservedQuantities F) {
version(multi_species_gas) {
    auto gmodel = myConfig.gmodel;

    // Species diffusion: Changed by NNG on 22/01/18.
    // We now apply both laminar and turbulent diffusion additively, to prevent artificially low
    // diffusion in areas with a small turbulent viscosity.
    if (laminarDiffusion) {
        myConfig.massDiffusion.update_mass_fluxes(fs, grad, jx, jy, jz);
    } else if (isTurbulent) { // Turbulent but no mass diffusion model
        foreach (isp; 0 .. n_species) {
            jx[isp] = 0.0;
            jy[isp] = 0.0;
            jz[isp] = 0.0;
        }
    }

    version(turbulence) {
        if (isTurbulent) {
            number D_t = fs.mu_t / (fs.gas.rho * Sc_t);
            foreach (isp; 0 .. n_species) {
                jx[isp] -= fs.gas.rho * D_t * grad.massf[isp][0];
                jy[isp] -= fs.gas.rho * D_t * grad.massf[isp][1];
                jz[isp] -= fs.gas.rho * D_t * grad.massf[isp][2];
            }
        }
    }
    number qx = 0.0;
    number qy = 0.0;
    number qz = 0.0;
    gmodel.enthalpies(fs.gas, hs);
    for(int isp=0; isp<n_species; isp++) {
        qx -= jx[isp] * hs[isp];
        qy -= jy[isp] * hs[isp];
        qz -= jz[isp] * hs[isp];
    }

    // Combine into fluxes: store as the dot product (F.n).
    number nx = n.x;
    number ny = n.y;
    number nz = n.z;
    auto cqi = myConfig.cqi;
    F[cqi.totEnergy] -= qx*nx + qy*ny + qz*nz;

    version(multi_T_gas) {
        for (int imode=0; imode<n_modes; imode++) {
            // Species diffusion contribution to the energy modes (added by NNG, 2022/01/26)
            for (int isp=0; isp<n_species; isp++) {
                number hMode = gmodel.enthalpyPerSpeciesInMode(fs.gas, isp, imode);
                // The sign here needs to be opposite to the thermal conduction, hence +=
                F[cqi.modes+imode] += hMode*(jx[isp]*n.x + jy[isp]*n.y + jz[isp]*n.z);
            }
        }
    }
    foreach (isp; 0 .. cqi.n_species) {
        F[cqi.species+isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
    }
} // end version(multi_species_gas)
} // end viscous_flux_calc()
