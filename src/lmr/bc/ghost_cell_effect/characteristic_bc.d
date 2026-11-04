// characteristic_bc.d
// 
// JSW & RWJ 2026-10-01
// 
// Formulation of characteristic boundary condition at the free stream
//
// From paper Motheau, E., Almgren, A. and Bell, J.B., 2017.
// Navier–Stokes characteristic boundary conditions using ghost cells. 
// AIAA Journal, 55(10), pp.3399-3408.
// Only working for structured mesh
//
// Only working for 2 layers ghost cells
// With Same dx and dy near this boundary
//
// In this stage, only test for east

module lmr.bc.ghost_cell_effect.characteristic_bc;

import std.conv;
import std.json;
import std.math;
import std.stdio;
import std.string;
import nm.number;
import ntypes.complex;


import gas;
import geom;

import lmr.bc;
import lmr.flowstate;
import lmr.fluidblock;
import lmr.fluidfvcell;
import lmr.fvinterface;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.sfluidblock;

// Stencil near boundary (eg. east boundary)
// 
//            +-----+
//            | j+1 |
// +----+-----+-----+------+-----+
// | i-2| i-1 |  i  | gc0  | gc1 |
// +----+-----+-----+------+-----+
//            | j-1 |
//            +-----+
//

class GhostCellCharacteristic : GhostCellEffect {
public:

   // double p_target, l_x;

     this(int id, int boundary)
     {
         super(id, boundary, "characteristic_bc");
         // this.p_target = ptarget;
         // this.l_x = lx;
     }

    //From this moment, we giving the p_target and l_x.
   // double p_target = 1.0e4;
   // double l_x = 0.5;
    
    
    // @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
    throw new Error("This Characteristic BC only working for structure grid.");
    }
    
    // @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
    throw new Error("This Characteristic BC only working for structure grid.");
    }

    // @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        
        throw new Error("This Characteristic BC only working for structure grid.");
     }

    // @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t[3] ijk;
        // size_t f_tmp_idx;
        FVInterface f_tmp;
        FluidFVCell c_i, c_im1, c_im2, c_jp1, c_jm1, gc0, gc1;
        charac_deriv cd;
        auto gmodel = blk.myConfig.gmodel;
        number gamma;
        
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        assert(blk.n_ghost_cell_layers == 2, "Oops, the ghost cell layers should be 2 for this characteristic BC");
        BoundaryCondition bc = blk.bc[which_boundary];
        size_t nfaces = bc.faces.length;
        foreach (f_idx, f; bc.faces) {
            final switch (which_boundary) { 
            case Face.east:
                c_i = f.left_cells[0];
                ijk = blk.to_ijk_indices_for_cell(c_i.id);
                f_tmp = blk.get_ifi(ijk[0], ijk[1], ijk[2]);
                c_im1 = f_tmp.left_cells[0];
                ijk = blk.to_ijk_indices_for_cell(c_im1.id);
                f_tmp = blk.get_ifi(ijk[0], ijk[1], ijk[2]);
                c_im2 = f_tmp.left_cells[0];
                gamma = gmodel.gamma(c_i.fs.gas);
                if (f_idx == 0 && f_idx == (nfaces-1)) {
                    // do something special at bottom corner
                    gc0.fs.copy_values_from(c_i.fs);
                    gc1.fs.copy_values_from(c_i.fs);

                }
                else {
                    c_jp1 = bc.faces[f_idx+1].left_cells[0];
                    c_jm1 = bc.faces[f_idx-1].left_cells[0];
                    double dx = c_i.pos[gtl].x.re - c_im1.pos[gtl].x.re;
                    double dy = c_i.pos[gtl].y.re - c_jm1.pos[gtl].x.re;
                    gradient(c_i.fs, c_im1.fs, c_im2.fs, c_jp1.fs, c_jm1.fs, gamma, dx, dy, cd);
                    characteristic(c_i.fs, c_im1.fs, gc0.fs, gc1.fs, dx, cd);
                }
                break;
            case Face.north, Face.south, Face.west, Face.top, Face.bottom:
                throw new Error("Characteristic BC not implemented for boundary.");
                break;
            }
                
           }
     }

private:
    // double p_target, l_x;
    // this(int id, int boundary, double ptarget, double lx)
    // {
    //     super(id, boundary, "characteristic_bc");
    //     this.p_target = ptarget;
    //     this.l_x = lx;
    // }

    // override string toString() const
    // {
    //     return "Characteristic BC, p_target= " ~to!string(p_target)~  ",l_x=" ~to!string(l_x);
    // }

    struct charac_deriv
    {
        number Dpdx;
        number Drhodx;
        number Dudx;
        number Dvdx;
    }
    @nogc
    void gradient(FlowState* fsi0j0, FlowState* fsim1j0, FlowState* fsim2j0,
                  FlowState* fsi0j1, FlowState* fsi0jm1, number gamma, double dx, double dy,
                  out charac_deriv cd)
    {
        immutable double sigma = 0.25;
        // From Boundary Conditions for Direct Simulations of
        // Compressible Viscous Flows, Poinsot and Lele, most case the beta can be 1.
        // But here we first follow the paper from Motheau et el.
        immutable double beta = 0;
        
        // define p_target and l_x before or get the data from input file;
        
        // The only one unknown wave going to domain, named L1 from Eq. 36
        // This L1 only for subsonic, but here we are going to implement this method for supersonic,
        // even the local region is subsonic
        // auto T1 = fsi0j0.vel.y * (((fsi0j1.gas.p - fsi0jm1.gas.p) / (2*dy)) 
        //   - fsi0j0.gas.rho * fsi0j0.gas.a * ((fsi0j1.vel.x - fsi0jm1.vel.x) / (2*dy)))
        //   + gamma * fsi0j0.gas.p * ((fsi0j1.vel.y - fsi0jm1.vel.y) / (2*dy));
        // auto K = sigma * fsi0j0.gas.a * (1 - (sqrt(fsi0j0.vel.x^^2 + fsi0j0.vel.y^^2) / fsi0j0.gas.a)^^2) / l_x;
        // auto L1 = K * (fsi0j0.gas.p - p_target) - (1 - beta) * T1;
        
        // L2 ~ L4
        auto drhodx = (fsim2j0.gas.rho - 4 * fsim1j0.gas.rho + 3 * fsi0j0.gas.rho) / (2 * dx);
        auto dpdx = (fsim2j0.gas.p - 4 * fsim1j0.gas.p + 3 * fsi0j0.gas.p) / (2 * dx);
        auto dudx = (fsim2j0.vel.x - 4 * fsim1j0.vel.x + 3 * fsi0j0.vel.x) / (2 * dx);
        auto dvdx = (fsim2j0.vel.y - 4 * fsim1j0.vel.y + 3 * fsi0j0.vel.y) / (2 * dx);

        auto L1 = (fsi0j0.vel.x - fsi0j0.gas.a) * (dpdx - fsi0j0.gas.rho * fsi0j0.gas.a * dudx);
        auto L2 = fsi0j0.vel.x * (fsi0j0.gas.a^^2 * drhodx - dpdx);
        auto L3 = fsi0j0.vel.x * dvdx;
        auto L4 = (fsi0j0.vel.x + fsi0j0.gas.a) * (dpdx + fsi0j0.gas.rho * fsi0j0.gas.a * dudx);

        // gradient for p, rho, u
        cd.Dpdx = 0.5 * ((L1/(fsi0j0.vel.x-fsi0j0.gas.a)) + (L4/(fsi0j0.vel.x+fsi0j0.gas.a)));
        cd.Drhodx = (1/fsi0j0.gas.a^^2) * ((L2/fsi0j0.vel.x) +
         0.5*(L1/(fsi0j0.vel.x-fsi0j0.gas.a) + L4/(fsi0j0.vel.x+fsi0j0.gas.a)));
        cd.Dudx = ((L4/(fsi0j0.vel.x+fsi0j0.gas.a)) - (L1/(fsi0j0.vel.x-fsi0j0.gas.a))) 
         * (1 / (2*fsi0j0.gas.rho*fsi0j0.gas.a));  
        cd.Dvdx = L3 / fsi0j0.vel.y;

        return;
    }
    @nogc
    void characteristic(FlowState* fsi0j0, FlowState* fsim1j0,
                   FlowState* gc0, FlowState* gc1, double dx, in charac_deriv cd)
    {
        auto gmodel = blk.myConfig.gmodel;
        gc0.gas.rho = fsim1j0.gas.rho + 2 * dx * cd.Drhodx;
        gc1.gas.rho = -2*fsim1j0.gas.rho - 3*fsi0j0.gas.rho + 6*gc0.gas.rho - 6*dx*cd.Drhodx; 
        gc0.vel.x = fsim1j0.vel.x + 2 * dx * cd.Dudx;
        gc1.vel.x = -2*fsim1j0.vel.x - 3*fsi0j0.vel.x + 6*gc0.vel.x - 6*dx*cd.Dudx; 
        gc0.vel.y = fsim1j0.vel.y + 2 * dx * cd.Dvdx;
        gc1.vel.y = -2*fsim1j0.vel.y - 3*fsi0j0.vel.y + 6*gc0.vel.y - 6*dx*cd.Dvdx; 
        gc0.gas.p = fsim1j0.gas.p + 2 * dx * cd.Dpdx;
        gc1.gas.p = -2*fsim1j0.gas.p - 3*fsi0j0.gas.p + 6*gc0.gas.p - 6*dx*cd.Dpdx;
        gmodel.update_thermo_from_rhou(gc0.gas);
        gmodel.update_thermo_from_rhou(gc1.gas);
   }
}
