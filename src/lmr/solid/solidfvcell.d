/**
 * solidfvcell.d
 *
 * A solid finite-volume cell, to be held by SolidBlock objects.
 *
 * Author: Rowan G. Kyle D. and Peter J.
 * Version: 2015-22-04
 */

module solidfvcell;

import std.conv;
import std.string;
import std.array;
import std.format;
import std.math;
import std.algorithm;
import ntypes.complex;
import nm.number;
import geom;
import solidfvinterface;
import solidfvvertex;
import lmr.solid.solidstate;
import lmr.solid.solidthermalmodel;
import std.stdio;
import globalconfig;
import lmr.fvcell : FVCell;

class SolidFVCell : FVCell {
public:
    size_t id;
    // Cell properties
    number volume;
    number areaxy;
    Vector3 pos;
    // Cell state
    SolidState ss;
    number T;
    number[] e;
    number[] dedt;
    number de_prev;
    // Cell source term
    number Q;
    // Connections
    SolidFVInterface[] iface;
    SolidFVVertex[] vtx;
    // Cell-centered gradients
    Vector3[] cloud_pos;
    number*[] cloud_T;
    number[12] wx, wy, wz;
    number dTdx;
    number dTdy;
    number dTdz;
    bool is_ghost = true;
    // list of cells and faces in the residual stencil
    SolidFVCell[] cell_list;
    SolidFVInterface[] face_list;
    // storage for local Jacobian entry
    double dRde;

private:
    LocalConfig myConfig;

public:

    this(LocalConfig myConfig)
    {
        this.myConfig = myConfig;
        e.length = myConfig.n_flow_time_levels;
        dedt.length = myConfig.n_flow_time_levels;
    }

    /**
     * A stripped down initialisation for use in file writing.
     *
     * Authors: RJG
     * Date: 2024-03-03
     */
    this(in Vector3 pos, number volume, double T, SolidThermalModel stm, int id_init=-1)
    // stripped down initialisation
    {
        id = id_init;
        this.pos = pos;
        this.volume = volume;
        this.T = T;
        this.ss = SolidState();
        ss.T = T;
        stm.updateEnergy(ss);
        this.e.length = 1;
        this.e[0] = ss.e;
    }


    @nogc
    void copy_values_from(SolidFVCell other) {
        volume = other.volume;
        areaxy = other.areaxy;
        pos = other.pos;
        ss = other.ss;
        T = other.T;
        foreach (i; 0..e.length) { e[i] = other.e[i] ; }
        foreach (i; 0..dedt.length) { dedt[i] = other.dedt[i] ; }
        de_prev = other.de_prev;
        Q = other.Q;
        dTdx = other.dTdx;
        dTdy = other.dTdy;
        dTdz = other.dTdz;
        is_ghost = other.is_ghost;
    }

    void scanValuesFromString(string buffer)
    {
        auto items = split(buffer);
        pos.x = to!double(items.front); items.popFront();
        pos.y = to!double(items.front); items.popFront();
        pos.z = to!double(items.front); items.popFront();
        volume = to!double(items.front); items.popFront();
        e[0] = to!double(items.front); items.popFront();
        T = to!double(items.front); items.popFront();
        ss.rho = to!double(items.front); items.popFront();
        ss.Cp = to!double(items.front); items.popFront();
        ss.k = to!double(items.front); items.popFront();
    }

    void timeDerivatives(int ftl, int dimensions)
    {
        SolidFVInterface IFn = iface[Face.north];
        SolidFVInterface IFe = iface[Face.east];
        SolidFVInterface IFs = iface[Face.south];
        SolidFVInterface IFw = iface[Face.west];
        SolidFVInterface IFt, IFb;
        if (dimensions == 3) {
            IFt = iface[Face.top];
            IFb = iface[Face.bottom];
        }
        // Cell volume (inverted).
        number volInv = 1.0 / volume;
        number integral;
        // Sum up fluxes (of form q.n)
        integral = -IFe.flux * IFe.area - IFn.flux * IFn.area
            + IFw.flux * IFw.area + IFs.flux * IFs.area;
        if (dimensions == 3) {
            integral += (-IFt.flux * IFt.area + IFb.flux * IFb.area);
        }
        dedt[ftl] = volInv * integral + Q;

        if (isNaN(dedt[ftl])) {
            string msg = "dedt is Not A Number in solid cell ";
            msg ~= format("(id: %d, pos: [%.8f, %.8f, %.8f])", id, pos.x.re, pos.y.re, pos.z.re);
            throw new FlowSolverException(msg);
        }
    }

    void stage1RKL1Update(double dt, int j, int s)
    {
        // RKL1 stage 1 coefficients
        double muj_tilde;
        muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);

        // stage 1 update
        e[1] = e[0] + muj_tilde*dt*dedt[0];

        return;
    } // end stage1RKL1Update()

    void stage2RKL1Update(double dt, int j, int s)
    {
        // stage j coefficients
        double muj; double vuj; double muj_tilde;
        muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);
        muj = (2.0*j-1)/j;
        vuj = (1.0-j)/j;

        // stage j update
        e[2] = muj*e[1] + vuj*e[0] + muj_tilde*dt*dedt[1];

        return;
    } // end stage2RKL1Update()

    void stage1RKL2Update(double dt, int j, int s)
    {
        // stage 1 coefficients
        double muj_tilde;
        muj_tilde = 4.0/(3.0*(s*s+s-2.0));

        // stage 1 update
        e[1] = e[0] + muj_tilde*dt*dedt[0];

        // make a copy of the initial conserved quantities for the stage j updates
        e[3] = e[0];

        return;
    } // end stage1RKL2Update()

    void stage2RKL2Update(double dt, int j, int s)
    {
        //  stage j coefficients
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

        // stage j update
        e[2] = muj*e[1] + vuj*e[0] + (1.0-muj-vuj)*e[3] + muj_tilde*dt*dedt[1] + gam_tilde*dt*dedt[0];

        return;
    } // end stage2RKL2Update()

    void eulerUpdate(double dt)
    {
        e[1] = e[0] + dt*dedt[0];

    }

    void stage1Update(double dt)
    {
        double gamma1 = 1.0; // Assume Euler
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc: gamma1 = 1.0; break;
        case GasdynamicUpdate.midpoint: gamma1 = 0.5; break;
        case GasdynamicUpdate.classic_rk3: gamma1 = 0.5; break;
        case GasdynamicUpdate.tvd_rk3: gamma1 = 1.0; break;
        case GasdynamicUpdate.denman_rk3: gamma1 = 8.0/15.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        case GasdynamicUpdate.classic_rk4: gamma1 = 0.5; break;
        }
        e[1] = e[0] + dt*gamma1*dedt[0];

    }

    void stage2Update(double dt)
    {
        // Assuming predictor-corrector
        double gamma1 = 0.5;
        double gamma2 = 0.5;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.moving_grid_1_stage: assert(false, "invalid for 1-stage update.");
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc: gamma1 = 0.5, gamma2 = 0.5; break;
        case GasdynamicUpdate.midpoint: gamma1 = 0.0; gamma2 = 1.0; break;
        case GasdynamicUpdate.classic_rk3: gamma1 = -1.0; gamma2 = 2.0; break;
        case GasdynamicUpdate.tvd_rk3: gamma1 = 0.25; gamma2 = 0.25; break;
        case GasdynamicUpdate.denman_rk3: gamma1 = -17.0/60.0; gamma2 = 5.0/12.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        case GasdynamicUpdate.classic_rk4: gamma1 = 0.0; gamma2 = 0.5; break;
        }
        e[2] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1]);

    }

    void stage3Update(double dt)
    {
        // Assuming TVD_RK3 scheme as done in flow update
        double gamma1 = 1.0/6.0; // presume TVD_RK3 scheme.
        double gamma2 = 1.0/6.0;
        double gamma3 = 4.0/6.0;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc:
        case GasdynamicUpdate.midpoint:
            assert(false, "invalid for 2-stage update.");
        case GasdynamicUpdate.classic_rk3: gamma1 = 1.0/6.0; gamma2 = 4.0/6.0; gamma3 = 1.0/6.0; break;
        case GasdynamicUpdate.tvd_rk3: gamma1 = 1.0/6.0; gamma2 = 1.0/6.0; gamma3 = 4.0/6.0; break;
            // FIX-ME: Check that we have Andrew Denman's scheme ported correctly.
        case GasdynamicUpdate.denman_rk3: gamma1 = 0.0; gamma2 = -5.0/12.0; gamma3 = 3.0/4.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        case GasdynamicUpdate.classic_rk4: gamma1 = 0.0; gamma2 = 0.0; gamma3 = 1.0; break;
        }
        e[3] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1] + gamma3*dedt[2]);

    }

    void stage4Update(double dt)
    {
        double gamma1 = 1.0/6.0; // presume classic_rk4 scheme.
        double gamma2 = 2.0/6.0;
        double gamma3 = 2.0/6.0;
        double gamma4 = 1.0/6.0;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc:
        case GasdynamicUpdate.midpoint:
        case GasdynamicUpdate.classic_rk3:
        case GasdynamicUpdate.tvd_rk3:
        case GasdynamicUpdate.denman_rk3:
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2:
        case GasdynamicUpdate.classic_rk4: gamma1 = 1.0/6.0; gamma2 = 2.0/6.0; gamma3 = 2.0/6.0; gamma4 = 1.0/6.0; break;
        }
        e[3] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1] + gamma3*dedt[2] + gamma4*dedt[3]);

    }

    void gather_residual_stencil_lists(int spatial_order_of_jacobian)
    {
        /*
          This function gathers the interfaces and cells that make up the residual stencil.
          These stencils are used in constructing the residual/flux Jacobian.
          The stencils can be thought of in terms of what neighbouring cells will have perturbed
          residuals in the event this cells conserved quantities are perturbed.

          Note that we need the cells to be in numerical order according to their local id for
          entry into the flow Jacobian later.

          We expect the spatial_order_of_jacobian will either be:
          0 : include just this cell (i.e. we will later form only the diagonal components of the Jacobian)
          1 : include nearest neighbours (approximate Jacobian)
          2 : include full stencil
         */

        SolidFVCell[] unordered_cell_list;  // TODO: think about possibly pre-sizing this array
        size_t[size_t] cell_pos_array;      // this is used to retrieve a cell from the unordered list

        // always add this cell
        size_t[] cell_ids;
        unordered_cell_list ~= this;
        cell_pos_array[this.id] = unordered_cell_list.length-1;
        cell_ids ~= this.id;
        bool nearest_neighbours = false;
        if ( spatial_order_of_jacobian >= 1) { nearest_neighbours = true; }

        if (nearest_neighbours) {
            foreach(f; iface) {
                foreach (cell; [f.cellLeft, f.cellRight]) {
                    bool cell_exists = cell_ids.canFind(cell.id);
                    if (!cell_exists && cell.id < 1_000_000_000) {
                        unordered_cell_list ~= cell;
                        cell_pos_array[cell.id] = unordered_cell_list.length-1;
                        cell_ids ~= cell.id;
                    }
                }
            }
        } // end if (nearest_neighbours)

        bool extended_neighbours = false;
        if ( nearest_neighbours && spatial_order_of_jacobian > 1) { extended_neighbours = true; }

        if (extended_neighbours) {
            size_t np = unordered_cell_list.length;
            foreach (idx; 1 .. np) {
                foreach(f; unordered_cell_list[idx].iface) {
                    foreach (cell; [f.cellLeft, f.cellRight]) {
                        bool cell_exists = cell_ids.canFind(cell.id);
                        if (!cell_exists && cell.id < 1_000_000_000) {
                            unordered_cell_list ~= cell;
                            cell_pos_array[cell.id] = unordered_cell_list.length-1;
                            cell_ids ~= cell.id;
                        }
                    }
                }
            }
        } // end if (extended_neighbours)

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

    version(complex_numbers) {
        @nogc void clear_imaginary_components()
        // When performing complex-step differentiation the variables used in the flux evaluation accumulate imaginary components,
        // this routine clears out these components so that we start with a clean slate, so to speak.
        {
            foreach (i; 0 .. e.length) { e[i].im = 0.0; }
            foreach (i; 0 .. dedt.length) { dedt[i].im = 0.0; }
            T.im = 0.0;
            dTdx.im = 0.0;
            dTdy.im = 0.0;
            dTdz.im = 0.0;
            foreach (face; iface) {
                face.T.im = 0.0;
                face.e.im = 0.0;
                face.flux.im = 0.0;
            }
        } // end clear_imaginary_components()
    } // end version(complex)
}

string[] varListForSolidCell()
{
    string[] list;
    list ~= ["pos.x", "pos.y", "pos.z", "volume"];
    list ~= ["e", "T"];
    list ~= ["rho", "Cp", "k"];
    list ~= ["k11", "k12", "k13"];
    list ~= ["k21", "k22", "k23"];
    list ~= ["k31", "k32", "k33"];
    return list;
}
