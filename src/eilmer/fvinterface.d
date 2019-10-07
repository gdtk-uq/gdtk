/**
 * fvinterface.d
 * Finite-volume cell-interface class, for use in the CFD codes.
 * Fluxes of conserved quantities are transported (between cells) across cell interfaces.

 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 *          2015-02-13: Keep an eye on the future of the moving_grid option.
 *          2015-05-04: keep references to adjoining cells and defining vertices.
 */

module fvinterface;

import std.stdio;

import std.conv;
import std.format;
import std.math;
import nm.complex;
import nm.number;

import geom;
import gas;
import fvcore;
import fvvertex;
import fvcell;
import flowstate;
import flowgradients;
import conservedquantities;
import globalconfig;
import lsqinterp;
import mass_diffusion;

class FVInterface {
public:
    int id;
    bool is_on_boundary = false;  // by default, assume not on boundary
    size_t bc_id;  // if the face is on a block boundary, which one
    size_t i_bndry; // if the face is on a block boundary, store index into the array of faces attached to bc
    bool use_wall_function_shear_and_heat_flux = false; // for use in viscous_flux_calc()
    bool in_suppress_reconstruction_zone; // if true, we no do reconstruction at this face
    //
    // Geometry
    Vector3 pos;           // position of the (approx) midpoint
    Vector3 gvel;          // grid velocity at interface, m/s
    number Ybar;           // Y-coordinate of the mid-point
    number length;         // Interface length in the x,y-plane
    number[] area;         // Area m**2 for each grid-time-level.
                           // Area per radian in axisymmetric geometry
    Vector3 n;             // Direction cosines for unit normal
    Vector3 t1;            // tangent vector 1 (aka p)
    Vector3 t2;            // tangent vector 2 (aka q)
    FVVertex[] vtx;        // references to vertices for line (2D) and quadrilateral (3D) faces
    //
    // Adjoining cells.
    // These are references to either active cells or ghost cells.
    // The reference may be nil if no cell has been assigned,
    // maybe for a boundary without ghost cells.
    FVCell left_cell;      // interface normal points out of this adjoining cell
    FVCell right_cell;     // interface normal points into this adjoining cell
    //
    // Flow
    FlowState fs;          // Flow properties
    ConservedQuantities F; // Flux conserved quantity per unit area
    number tau_wall_x, tau_wall_y, tau_wall_z; // shear at face (used by wall-function BCs)
    number q;              // heat-flux across face (used by wall-function BCs)
    //
    // Viscous-flux-related quantities.
    FlowGradients grad;
    WLSQGradWorkspace ws_grad;
    Vector3*[] cloud_pos; // Positions of flow points for gradients calculation.
    FlowState[] cloud_fs; // References to flow states at those points.
    number[] jx; // diffusive mass flux in x
    number[] jy; // diffusive mass flux in y
    number[] jz; // diffusive mass flux in z
    number q_diffusion; 
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
    this(LocalConfig myConfig,
         bool allocate_spatial_deriv_lsq_workspace,
         int id_init=-1)
    {
        this.myConfig = myConfig;
        id = id_init;
        area.length = myConfig.n_grid_time_levels;
        gvel = Vector3(0.0,0.0,0.0); // default to fixed grid
        auto gmodel = myConfig.gmodel;
        uint n_species = myConfig.n_species;
        uint n_modes = myConfig.n_modes;
        double T = 300.0;
        double[] T_modes; foreach(i; 0 .. n_modes) { T_modes ~= 300.0; }
        fs = new FlowState(gmodel, 100.0e3, T, T_modes, Vector3(0.0,0.0,0.0));
        F = new ConservedQuantities(n_species, n_modes);
        F.clear();
        grad = new FlowGradients(myConfig);
        if (allocate_spatial_deriv_lsq_workspace) {
            ws_grad = new WLSQGradWorkspace();
        }
        version(multi_species_gas) {
            jx.length = n_species;
            jy.length = n_species;
            jz.length = n_species;
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
    }

    this(FVInterface other, GasModel gm)
    {
        id = other.id;
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
        fs = new FlowState(other.fs, gm);
        F = new ConservedQuantities(other.F);
        tau_wall_x = other.tau_wall_x;
        tau_wall_y = other.tau_wall_y;
        tau_wall_z = other.tau_wall_z;
        q = other.q;
        grad = new FlowGradients(other.grad);
        if (other.ws_grad) ws_grad = new WLSQGradWorkspace(other.ws_grad);
        // Because we copy the following pointers and references,
        // we cannot have const (or "in") qualifier on other.
        cloud_pos = other.cloud_pos.dup();
        cloud_fs = other.cloud_fs.dup();
        version(multi_species_gas) {
            jx = other.jx.dup();
            jy = other.jy.dup();
            jz = other.jz.dup();
        }
        version(steadystate) {
            dFdU_L.length = 5; // number of conserved variables; FIX-ME for versions
            foreach (ref a; dFdU_L) a.length = 5;
            dFdU_R.length = 5;
            foreach (ref a; dFdU_R) a.length = 5;
        }
        q_diffusion = other.q_diffusion;
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
            grad.copy_values_from(other.grad);
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
    } // end update_2D_geometric_data()

    @nogc
    void update_3D_geometric_data(size_t gtl)
    {
        switch (vtx.length) {
        case 3:
            triangle_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                                vtx[2].pos[gtl],
                                pos, n, t1, t2, area[gtl]);
            break;
        case 4:
            quad_properties(vtx[0].pos[gtl], vtx[1].pos[gtl],
                            vtx[2].pos[gtl], vtx[3].pos[gtl],
                            pos, n, t1, t2, area[gtl]);
            break;
        default:
            string msg = "FVInterface.update_3D_geometric_data(): Unhandled number of vertices: ";
            debug { msg ~= format("%d", vtx.length); }
            throw new FlowSolverException(msg);
        } // end switch     
    } // end update_3D_geometric_data()

    @nogc
    void average_vertex_deriv_values()
    {
        grad.copy_values_from(vtx[0].grad);
        foreach (i; 1 .. vtx.length) grad.accumulate_values_from(vtx[i].grad);
        grad.scale_values_by(to!number(1.0/vtx.length));
    } // end average_vertex_deriv_values()

    @nogc
    void average_cell_deriv_values(int gtl)
    {
        // if the interface is along a boundary that doesn't have a mapping to
        // a cell in a neighbouring block, i.e. the interface is along a domain
        // boundary, then we just copy the gradient from the cell-center to the interface.
        if (left_cell.is_interior_to_domain == false ||
            right_cell.is_interior_to_domain == false) {
            FVCell c;
            if (left_cell.is_interior_to_domain == false) c = right_cell;
            else c = left_cell;
            // vel-x
            grad.vel[0][0] = c.grad.vel[0][0];
            grad.vel[0][1] = c.grad.vel[0][1];
            grad.vel[0][2] = c.grad.vel[0][2];
            
            // vel-y
            grad.vel[1][0] = c.grad.vel[1][0];
            grad.vel[1][1] = c.grad.vel[1][1];
            grad.vel[1][2] = c.grad.vel[1][2];
            
            // vel-z
            grad.vel[2][0] = c.grad.vel[2][0];
            grad.vel[2][1] = c.grad.vel[2][1];
            grad.vel[2][2] = c.grad.vel[2][2];
            
            // massf
            version(multi_species_gas) {
                uint nsp = (myConfig.sticky_electrons) ? myConfig.n_heavy : myConfig.n_species;
                foreach (isp; 0 .. nsp) {
                    grad.massf[isp][0] = c.grad.massf[isp][0];
                    grad.massf[isp][1] = c.grad.massf[isp][1];
                    grad.massf[isp][2] = c.grad.massf[isp][2];
                }
            }
            
            // T
            grad.T[0] = c.grad.T[0];
            grad.T[1] = c.grad.T[1];
            grad.T[2] = c.grad.T[2];
            
            version(komega) {
                // tke
                grad.tke[0] = c.grad.tke[0];
                grad.tke[1] = c.grad.tke[1];
                grad.tke[2] = c.grad.tke[2];
                
                // omega
                grad.omega[0] = c.grad.omega[0];
                grad.omega[1] = c.grad.omega[1];
                grad.omega[2] = c.grad.omega[2];
            }
        } else {
            number qL; number qR;
            FVCell cL0 = left_cell; // i
            FVCell cR0 = right_cell; // j
            // interface normal
            number nx = n.x;
            number ny = n.y;
            number nz = n.z;
            // vector from left-cell-centre to face midpoint
            number rLx = pos.x - cL0.pos[gtl].x;
            number rLy = pos.y - cL0.pos[gtl].y;
            number rLz = pos.z - cL0.pos[gtl].z;
            number rRx = pos.x - cR0.pos[gtl].x;
            number rRy = pos.y - cR0.pos[gtl].y;
            number rRz = pos.z - cR0.pos[gtl].z;
            // vector from left-cell-centre to right-cell-centre
            number ex = cR0.pos[gtl].x - cL0.pos[gtl].x;
            number ey = cR0.pos[gtl].y - cL0.pos[gtl].y;
            number ez = cR0.pos[gtl].z - cL0.pos[gtl].z;                
            // ehat
            number emag = sqrt(ex*ex + ey*ey + ez*ez);
            number ehatx = ex/emag;
            number ehaty = ey/emag;
            number ehatz = ez/emag;                
            // ndotehat
            number ndotehat = nx*ehatx + ny*ehaty + nz*ehatz;
            number avgdotehat;
            number jump;
            
            // vel-x
            avgdotehat = 0.5*(cL0.grad.vel[0][0]+cR0.grad.vel[0][0])*ehatx +
                0.5*(cL0.grad.vel[0][1]+cR0.grad.vel[0][1])*ehaty +
                0.5*(cL0.grad.vel[0][2]+cR0.grad.vel[0][2])*ehatz;
            jump = avgdotehat - (cR0.fs.vel.x - cL0.fs.vel.x)/emag;
            grad.vel[0][0] = 0.5*(cL0.grad.vel[0][0]+cR0.grad.vel[0][0]) - jump*(nx/ndotehat);
            grad.vel[0][1] = 0.5*(cL0.grad.vel[0][1]+cR0.grad.vel[0][1]) - jump*(ny/ndotehat);
            grad.vel[0][2] = 0.5*(cL0.grad.vel[0][2]+cR0.grad.vel[0][2]) - jump*(nz/ndotehat);
            
            // vel-y
            avgdotehat = 0.5*(cL0.grad.vel[1][0]+cR0.grad.vel[1][0])*ehatx +
                0.5*(cL0.grad.vel[1][1]+cR0.grad.vel[1][1])*ehaty +
                0.5*(cL0.grad.vel[1][2]+cR0.grad.vel[1][2])*ehatz;
            jump = avgdotehat - (cR0.fs.vel.y - cL0.fs.vel.y)/emag;
            grad.vel[1][0] = 0.5*(cL0.grad.vel[1][0]+cR0.grad.vel[1][0]) - jump*(nx/ndotehat);
            grad.vel[1][1] = 0.5*(cL0.grad.vel[1][1]+cR0.grad.vel[1][1]) - jump*(ny/ndotehat);
            grad.vel[1][2] = 0.5*(cL0.grad.vel[1][2]+cR0.grad.vel[1][2]) - jump*(nz/ndotehat);
            
            // vel-z
            avgdotehat = 0.5*(cL0.grad.vel[2][0]+cR0.grad.vel[2][0])*ehatx +
                0.5*(cL0.grad.vel[2][1]+cR0.grad.vel[2][1])*ehaty +
                0.5*(cL0.grad.vel[2][2]+cR0.grad.vel[2][2])*ehatz;
            jump = avgdotehat - (cR0.fs.vel.z - cL0.fs.vel.z)/emag;
            grad.vel[2][0] = 0.5*(cL0.grad.vel[2][0]+cR0.grad.vel[2][0]) - jump*(nx/ndotehat);
            grad.vel[2][1] = 0.5*(cL0.grad.vel[2][1]+cR0.grad.vel[2][1]) - jump*(ny/ndotehat);
            grad.vel[2][2] = 0.5*(cL0.grad.vel[2][2]+cR0.grad.vel[2][2]) - jump*(nz/ndotehat);
            
            // massf
            version(multi_species_gas) {
                uint nsp = (myConfig.sticky_electrons) ? myConfig.n_heavy : myConfig.n_species;
                foreach (isp; 0 .. nsp) {
                    avgdotehat = 0.5*(cL0.grad.massf[isp][0]+cR0.grad.massf[isp][0])*ehatx +
                        0.5*(cL0.grad.massf[isp][1]+cR0.grad.massf[isp][1])*ehaty +
                        0.5*(cL0.grad.massf[isp][2]+cR0.grad.massf[isp][2])*ehatz;
                    jump = avgdotehat - (cR0.fs.gas.massf[isp] - cL0.fs.gas.massf[isp])/emag;
                    grad.massf[isp][0] = 0.5*(cL0.grad.massf[isp][0]+cR0.grad.massf[isp][0]) - jump*(nx/ndotehat);
                    grad.massf[isp][1] = 0.5*(cL0.grad.massf[isp][1]+cR0.grad.massf[isp][1]) - jump*(ny/ndotehat);
                    grad.massf[isp][2] = 0.5*(cL0.grad.massf[isp][2]+cR0.grad.massf[isp][2]) - jump*(nz/ndotehat);
                }
            }
            
            // T
            avgdotehat = 0.5*(cL0.grad.T[0]+cR0.grad.T[0])*ehatx +
                0.5*(cL0.grad.T[1]+cR0.grad.T[1])*ehaty +
                0.5*(cL0.grad.T[2]+cR0.grad.T[2])*ehatz;
            jump = avgdotehat - (cR0.fs.gas.T - cL0.fs.gas.T)/emag;
            grad.T[0] = 0.5*(cL0.grad.T[0]+cR0.grad.T[0]) - jump*(nx/ndotehat);
            grad.T[1] = 0.5*(cL0.grad.T[1]+cR0.grad.T[1]) - jump*(ny/ndotehat);
            grad.T[2] = 0.5*(cL0.grad.T[2]+cR0.grad.T[2]) - jump*(nz/ndotehat);
            
            version(komega) {
                // tke
                avgdotehat = 0.5*(cL0.grad.tke[0]+cR0.grad.tke[0])*ehatx +
                    0.5*(cL0.grad.tke[1]+cR0.grad.tke[1])*ehaty +
                    0.5*(cL0.grad.tke[2]+cR0.grad.tke[2])*ehatz;
                jump = avgdotehat - (cR0.fs.tke - cL0.fs.tke)/emag;
                grad.tke[0] = 0.5*(cL0.grad.tke[0]+cR0.grad.tke[0]) - jump*(nx/ndotehat);
                grad.tke[1] = 0.5*(cL0.grad.tke[1]+cR0.grad.tke[1]) - jump*(ny/ndotehat);
                grad.tke[2] = 0.5*(cL0.grad.tke[2]+cR0.grad.tke[2]) - jump*(nz/ndotehat);
                
                // omega
                avgdotehat = 0.5*(cL0.grad.omega[0]+cR0.grad.omega[0])*ehatx +
                    0.5*(cL0.grad.omega[1]+cR0.grad.omega[1])*ehaty +
                    0.5*(cL0.grad.omega[2]+cR0.grad.omega[2])*ehatz;
                jump = avgdotehat - (cR0.fs.omega - cL0.fs.omega)/emag;
                grad.omega[0] = 0.5*(cL0.grad.omega[0]+cR0.grad.omega[0]) - jump*(nx/ndotehat);
                grad.omega[1] = 0.5*(cL0.grad.omega[1]+cR0.grad.omega[1]) - jump*(ny/ndotehat);
                grad.omega[2] = 0.5*(cL0.grad.omega[2]+cR0.grad.omega[2]) - jump*(nz/ndotehat);
            }
        }
    } // end average_cell_spatial_derivs()
            
    @nogc
    void viscous_flux_calc()
    // Unified 2D and 3D viscous-flux calculation.
    // Note that the gradient values need to be in place before calling this procedure.
    {
        auto gmodel = myConfig.gmodel;
        uint n_species = (myConfig.sticky_electrons) ? myConfig.n_heavy : myConfig.n_species;
        uint n_modes = myConfig.n_modes;
        double viscous_factor = myConfig.viscous_factor;
        number k_laminar = fs.gas.k;
        number mu_laminar = fs.gas.mu;
        if (myConfig.use_viscosity_from_cells) {
            // Emulate Eilmer3 behaviour by using the viscous transport coefficients
            // from the cells either side of the interface.
            if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
                k_laminar = 0.5*(left_cell.fs.gas.k+right_cell.fs.gas.k);
                mu_laminar = 0.5*(left_cell.fs.gas.mu+right_cell.fs.gas.mu);
            } else if (left_cell && left_cell.is_interior_to_domain) {
                k_laminar = left_cell.fs.gas.k;
                mu_laminar = left_cell.fs.gas.mu;
            } else if (right_cell && right_cell.is_interior_to_domain) {
                k_laminar = right_cell.fs.gas.k;
                mu_laminar = right_cell.fs.gas.mu;
            } else {
                assert(0, "Oops, don't seem to have a cell available.");
            }
        }
        number k_eff;
        number mu_eff;
        number lmbda;
        // we would like to use the most up to date turbulent properties, so take averages of the neighbouring cell values
        if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
            k_eff = viscous_factor * (k_laminar + 0.5*(left_cell.fs.k_t+right_cell.fs.k_t));
            mu_eff = viscous_factor * (mu_laminar + 0.5*(left_cell.fs.mu_t+right_cell.fs.mu_t));
        } else if (left_cell && left_cell.is_interior_to_domain) {
            k_eff = viscous_factor * (k_laminar + left_cell.fs.k_t);
            mu_eff = viscous_factor * (mu_laminar + left_cell.fs.mu_t);
        } else if (right_cell && right_cell.is_interior_to_domain) {
            k_eff = viscous_factor * (k_laminar + right_cell.fs.k_t);
            mu_eff = viscous_factor * (mu_laminar + right_cell.fs.mu_t);
        } else {
            assert(0, "Oops, don't seem to have a cell available.");
        }
        lmbda = -2.0/3.0 * mu_eff;
        //
        number local_pressure;
        if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
            local_pressure = 0.5*(left_cell.fs.gas.p+right_cell.fs.gas.p);
        } else if (left_cell && left_cell.is_interior_to_domain) {
            local_pressure = left_cell.fs.gas.p;
        } else if (right_cell && right_cell.is_interior_to_domain) {
            local_pressure = right_cell.fs.gas.p;
        } else {
            assert(0, "Oops, don't seem to have a cell available.");
        }
        number shear_stress_limit = myConfig.shear_stress_relative_limit * local_pressure;
        number heat_transfer_limit = (mu_eff > 0.0) ? k_eff/mu_eff*shear_stress_limit : to!number(0.0);

        if (myConfig.spatial_deriv_from_many_points) {
            // Viscous fluxes are constructed with gradients from many points.
            // The alternative, below, is just to use two points.
            //
            // We separate diffusion based on laminar or turbulent
            // and treat the differently.
            if (myConfig.turbulence_model != TurbulenceModel.none) {
                double Sc_t = myConfig.turbulence_schmidt_number;
                number D_t; // = fs.mu_t / (fs.gas.rho * Sc_t)
                // we would like to use the most up to date turbulent properties,
                // so take averages of the neighbouring cell values
                if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
                    D_t = 0.5*(left_cell.fs.mu_t+right_cell.fs.mu_t) / (fs.gas.rho * Sc_t);
                } else if (left_cell && left_cell.is_interior_to_domain) {
                    D_t = left_cell.fs.mu_t / (fs.gas.rho * Sc_t);
                } else if (right_cell && right_cell.is_interior_to_domain) {
                    D_t = right_cell.fs.mu_t / (fs.gas.rho * Sc_t);
                } else {
                    assert(0, "Oops, don't seem to have a cell available.");
                }
                version(multi_species_gas) {
                    foreach (isp; 0 .. n_species) {
                        jx[isp] = -fs.gas.rho * D_t * grad.massf[isp][0];
                        jy[isp] = -fs.gas.rho * D_t * grad.massf[isp][1];
                        jz[isp] = -fs.gas.rho * D_t * grad.massf[isp][2];
                    }
                }
            }
            else { // apply molecular diffusion instead of turbulence in laminar flows
                version(multi_species_gas) {
                    if (myConfig.mass_diffusion_model != MassDiffusionModel.none) {
                        myConfig.massDiffusion.update_mass_fluxes(fs, grad, jx, jy, jz);
                        foreach (isp; 0 .. n_species) {
                            jx[isp] *= viscous_factor;
                            jy[isp] *= viscous_factor;
                            jz[isp] *= viscous_factor;
                        }
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
            version(multi_species_gas) {
                if (myConfig.turbulence_model != TurbulenceModel.none ||
                    myConfig.mass_diffusion_model != MassDiffusionModel.none ) {
                    q_diffusion = to!number(0.0);
                    foreach (isp; 0 .. n_species) {
                        number h = gmodel.enthalpy(fs.gas, cast(int)isp);
                        qx -= jx[isp] * h;
                        qy -= jy[isp] * h;
                        qz -= jz[isp] * h;
                        q_diffusion -= (jx[isp]*h*n.x + jy[isp]*h*n.y + jz[isp]*h*n.z);
                        version(multi_T_gas) {
                            foreach (imode; 0 .. n_modes) {
                                number hMode = gmodel.enthalpyPerSpeciesInMode(fs.gas, cast(int)isp, cast(int)imode);
                                q_diffusion -= (jx[isp]*hMode*n.x + jy[isp]*hMode*n.y + jz[isp]*hMode*n.z);
                            }
                        } // end multi_T_gas
                    }
                } // multi_species_gas
            }
            version(komega) {
                number tau_kx = 0.0;
                number tau_ky = 0.0;
                number tau_kz = 0.0;
                number tau_wx = 0.0;
                number tau_wy = 0.0;
                number tau_wz = 0.0;
                if ( myConfig.turbulence_model == TurbulenceModel.k_omega &&
                     !(myConfig.axisymmetric && (Ybar <= 1.0e-10)) ) {
                    // Turbulence contribution to the shear stresses.
                    tau_xx -= 2.0/3.0 * fs.gas.rho * fs.tke;
                    tau_yy -= 2.0/3.0 * fs.gas.rho * fs.tke;
                    if (myConfig.dimensions == 3) { tau_zz -= 2.0/3.0 * fs.gas.rho * fs.tke; }
                    // Turbulence contribution to heat transfer.
                    number sigma_star = 0.6;
                    number mu_effective = fs.gas.mu + sigma_star * fs.gas.rho * fs.tke / fs.omega;
                    // Apply a limit on mu_effective in the same manner as that applied to mu_t.
                    mu_effective = fmin(mu_effective, myConfig.max_mu_t_factor * fs.gas.mu);
                    qx += mu_effective * grad.tke[0];
                    qy += mu_effective * grad.tke[1];
                    if (myConfig.dimensions == 3) { qz += mu_effective * grad.tke[2]; }
                    // Turbulence transport of the turbulence properties themselves.
                    tau_kx = mu_effective * grad.tke[0]; 
                    tau_ky = mu_effective * grad.tke[1];
                    if (myConfig.dimensions == 3) { tau_kz = mu_effective * grad.tke[2]; }
                    number sigma = 0.5;
                    mu_effective = fs.gas.mu + sigma * fs.gas.rho * fs.tke / fs.omega;
                    // Apply a limit on mu_effective in the same manner as that applied to mu_t.
                    mu_effective = fmin(mu_effective, myConfig.max_mu_t_factor * fs.gas.mu);
                    tau_wx = mu_effective * grad.omega[0]; 
                    tau_wy = mu_effective * grad.omega[1]; 
                    if (myConfig.dimensions == 3) { tau_wz = mu_effective * grad.omega[2]; } 
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
            // In some cases, the shear and heat fluxes have been previously
            // computed by the wall functions in the boundary condition call.
            if (use_wall_function_shear_and_heat_flux) {
                // Mass flux -- NO CONTRIBUTION, unless there's diffusion (below)
                // [TODO] As per Jason's recommendation, we need to do something 
                // to correct for corner cells. 
                // [TODO] Currently implemented for 2D; need to extend to 3D.
                F.momentum.refx -= tau_xx*nx + tau_wall_x;
                F.momentum.refy -= tau_yy*ny + tau_wall_y;
                F.momentum.refz -= tau_zz*nz + tau_wall_z;
                F.total_energy -=
                    tau_xx*fs.vel.x*nx + tau_yy*fs.vel.y*ny + tau_zz*fs.vel.z*nz +
                    tau_wall_x*fs.vel.x + tau_wall_y*fs.vel.y + tau_wall_z*fs.vel.z + q;
            }
            else { // proceed with locally computed shear and heat flux
                // Mass flux -- NO CONTRIBUTION, unless there's diffusion (below)
                F.momentum.refx -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
                F.momentum.refy -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
                F.momentum.refz -= tau_xz*nx + tau_yz*ny + tau_zz*nz;
                F.total_energy -=
                    (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx)*nx +
                    (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy)*ny +
                    (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz)*nz;
            } // end if
            version(multi_T_gas) {
                foreach (imode; 0 .. n_modes) {
                    F.energies[imode] -= viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][0] * nx;
                    F.energies[imode] -= viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][1] * ny;
                    F.energies[imode] -= viscous_factor * fs.gas.k_modes[imode] * grad.T_modes[imode][2] * nz;
                }
            }
            version(komega) {
                if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    F.tke -= tau_kx * nx + tau_ky * ny + tau_kz * nz;
                    F.omega -= tau_wx * nx + tau_wy * ny + tau_wz * nz;
                }
            }
            version(multi_species_gas) {
                if (myConfig.turbulence_model != TurbulenceModel.none ||
                    myConfig.mass_diffusion_model != MassDiffusionModel.none) {
                    foreach (isp; 0 .. n_species) {
                        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
                    }
                }
            }
        } // end of viscous-flux calculation with gradients from many points
        else {
            // Compute viscous fluxes with gradients come just from two points,
            // as suggested by Paul Petrie-Repar long ago.
            //
            // First, select the 2 points.
            number x0, x1, y0, y1, z0, z1;
            number velx0, velx1, vely0, vely1, velz0, velz1, T0, T1;
            if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
                x0 = left_cell.pos[0].x; x1 = right_cell.pos[0].x;
                y0 = left_cell.pos[0].y; y1 = right_cell.pos[0].y;
                z0 = left_cell.pos[0].z; z1 = right_cell.pos[0].z;
                velx0 = left_cell.fs.vel.x; velx1 = right_cell.fs.vel.x;
                vely0 = left_cell.fs.vel.y; vely1 = right_cell.fs.vel.y;
                velz0 = left_cell.fs.vel.z; velz1 = right_cell.fs.vel.z;
                T0 = left_cell.fs.gas.T; T1 = right_cell.fs.gas.T;
            } else if (left_cell && left_cell.is_interior_to_domain) {
                x0 = left_cell.pos[0].x; x1 = pos.x;
                y0 = left_cell.pos[0].y; y1 = pos.y;
                z0 = left_cell.pos[0].z; z1 = pos.z;
                velx0 = left_cell.fs.vel.x; velx1 = fs.vel.x;
                vely0 = left_cell.fs.vel.y; vely1 = fs.vel.y;
                velz0 = left_cell.fs.vel.z; velz1 = fs.vel.z;
                T0 = left_cell.fs.gas.T; T1 = fs.gas.T;
            } else if (right_cell && right_cell.is_interior_to_domain) {
                x0 = pos.x; x1 = right_cell.pos[0].x;
                y0 = pos.y; y1 = right_cell.pos[0].y;
                z0 = pos.z; z1 = right_cell.pos[0].z;
                velx0 = fs.vel.x; velx1 = right_cell.fs.vel.x;
                vely0 = fs.vel.y; vely1 = right_cell.fs.vel.y;
                velz0 = fs.vel.z; velz1 = right_cell.fs.vel.z;
                T0 = fs.gas.T; T1 = right_cell.fs.gas.T;
            } else {
                assert(0, "Oops, don't seem to have a cell available.");
            }
            // Distance in direction of face normal is what we will use
            // in the finite-difference approximation.
            number nx = n.x;
            number ny = n.y;
            number nz = n.z;
            // Distance between points, in face-normal direction.
            number deln = (x1-x0)*nx + (y1-y0)*ny;
            // Velocity in face-normal direction.
            number veln0 = velx0*nx + vely0*ny;
            number veln1 = velx1*nx + vely1*ny;
            if (myConfig.dimensions == 3) {
                deln += (z1-z0)*nz;
                veln0 += velz0*nz;
                veln1 += velz1*nz;
            } 
            number veln_face;
            if (left_cell && right_cell && left_cell.is_interior_to_domain && right_cell.is_interior_to_domain) {
                veln_face = 0.5*(veln0+veln1);
            } else if (left_cell && left_cell.is_interior_to_domain) {
                veln_face = veln1;
            } else if (right_cell && right_cell.is_interior_to_domain) {
                veln_face = veln0;
            }
            number ke0 = 0.5*(velx0^^2 + vely0^^2 + velz0^^2);
            number ke1 = 0.5*(velx1^^2 + vely1^^2 + velz1^^2);
            // Derivatives normal to face.
            number dvelxdn = (velx1 - velx0)/deln; // gradient of momentum per mass
            number dvelydn = (vely1 - vely0)/deln;
            number dvelzdn = (myConfig.dimensions == 3) ? (velz1 - velz0)/deln : to!number(0.0);
            number dkedn = (ke1 - ke0)/deln; // gradient of kinetic energy per mass
            number dTdn = (T1 - T0)/deln; // gradient of temperature
            number dvelndn = (veln1 - veln0)/deln; // for velocity dilatation
            //
            // Finally, compute the actual fluxes.
            //
            if (use_wall_function_shear_and_heat_flux) {
                assert(0, "Oops, not implemented.");
            }
            else { // proceed with locally computed shear and heat flux
                number tau_x = mu_eff * dvelxdn;
                number tau_y = mu_eff * dvelydn;
                number tau_z = mu_eff * dvelzdn;
                number qn = k_eff*dTdn + mu_eff*dkedn + lmbda*dvelndn*veln_face;
                tau_x = copysign(fmin(fabs(tau_x),shear_stress_limit), tau_x);
                tau_y = copysign(fmin(fabs(tau_y),shear_stress_limit), tau_y);
                tau_z = copysign(fmin(fabs(tau_z),shear_stress_limit), tau_z);
                qn = copysign(fmin(fabs(qn),heat_transfer_limit), qn);
                // Mass flux -- NO CONTRIBUTION, unless there's diffusion (below)
                F.momentum.refx -= tau_x;
                F.momentum.refy -= tau_y;
                F.momentum.refz -= tau_z;
                F.total_energy -= qn;
            } // end if (use_wall_function_shear_and_heat_flux)
            version(multi_T_gas) {
                // [TODO] Rowan, Modal energy flux?
            }
            version(komega) {
                if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    number tau_kn = 0.0; // [TODO] PJ, 2018-05-05, talk to Wilson
                    number tau_wn = 0.0;
                    F.tke -= tau_kn;
                    F.omega -= tau_wn;
                }
            }
            version(multi_species_gas) {
                if (myConfig.turbulence_model != TurbulenceModel.none ||
                    myConfig.mass_diffusion_model != MassDiffusionModel.none) {
                    // Mass diffusion done separately.
                    assert(0, "Oops, not implemented.");
                }
            }
        } // end of viscous-flux calculation with gradients from two points
    } // end viscous_flux_calc()

} // end of class FV_Interface
