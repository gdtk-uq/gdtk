/**
 * Eilmer4 boundary interface loads writing functions.
 *
 * Author: Kyle Damm
 * First code: 2017-03-17
 * Edits by Will Landsberg and Tim Cullen
 * 
 * 2018-01-20 Some code clean up by PJ;
 * [TODO] more is needed for MPI context.
 */

module loads;

import std.stdio;
import std.string;
import std.file;
import std.array;
import std.algorithm;
import std.format;
import std.parallelism;
import std.math;

import fvcore;
import globalconfig;
import globaldata;
import fvcell;
import fileutil;
import solidfvcell;
import solidblock;
import flowstate;
import flowgradients;
import geom;
import fluidblock;
import sfluidblock: SFluidBlock;
import ufluidblock: UFluidBlock;
import fvinterface;
import bc;

string loadsDir = "loads";

void init_loads_times_file()
{
    string fname = loadsDir ~ "/" ~ GlobalConfig.base_file_name ~ "-loads.times";
    std.file.write(fname, "# 1:loads_index 2:sim_time\n");
}

void update_loads_times_file(double sim_time, int current_loads_tindx)
{
    string fname = loadsDir ~ "/" ~ GlobalConfig.base_file_name ~ "-loads.times";
    std.file.append(fname, format("%04d %.18e \n", current_loads_tindx, sim_time));
}

void write_boundary_loads_to_file(double sim_time, int current_loads_tindx) {
    foreach (blk; localFluidBlocks) {
        if (blk.active) {
            final switch (blk.grid_type) {
            case Grid_t.unstructured_grid: 
                auto ublk = cast(UFluidBlock) blk;
                assert(ublk !is null, "Oops, this should be a UFluidBlock object.");
                apply_unstructured_grid(ublk, sim_time, current_loads_tindx);
                break;
            case Grid_t.structured_grid:
                auto sblk = cast(SFluidBlock) blk;
                assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                apply_structured_grid(sblk, sim_time, current_loads_tindx);
            } // end final switch
        } // if (blk.active)
    } // end foreach
} // end write_boundary_loads_to_file()

string generate_boundary_load_file(int current_loads_tindx, double sim_time, string group)
{
    // generate data file -- naming format tindx_groupname.dat
    string fname = format("%s/t%04d-%s.dat", loadsDir, current_loads_tindx, group);
    // check if file exists already, if not create file
    if (!exists(fname)) {
        auto f = File(fname, "w");
        f.writeln("# t = ", sim_time);
        f.writeln("# 1:pos.x 2:pos.y 3:pos.z 4:area 5:q 6:tau 7:l_tau 8:m_tau 9:n_tau 10:sigma "~
                  "11:n.x 12:n.y 13:n.z 14:T 15:Re_wall 16:y+ 17:cellWidthNormalToSurface");
        f.close();
    }
    return fname;
} // end generate_boundary_load_file()

void apply_unstructured_grid(UFluidBlock blk, double sim_time, int current_loads_tindx)
{
    foreach (bndary; blk.bc) {
        if (canFind(bndary.group, GlobalConfig.boundary_group_for_loads)) {
            string fname = generate_boundary_load_file(current_loads_tindx, sim_time, bndary.group);
            foreach (i, iface; bndary.faces) {
                // cell width normal to surface
                double w = (bndary.outsigns[i] == 1) ? iface.left_cell.L_min : iface.right_cell.L_min;
                compute_and_store_loads(iface, w, sim_time, fname);
            }
        }
    }
} // end apply_unstructured_grid()

void apply_structured_grid(SFluidBlock blk, double sim_time, int current_loads_tindx) {
    foreach (bndary; blk.bc) {
        if (canFind(bndary.group, GlobalConfig.boundary_group_for_loads)) {
            string fname = generate_boundary_load_file(current_loads_tindx, sim_time, bndary.group);
            size_t i, j, k;
            FVCell cell;
            FVInterface IFace;
            double cellWidthNormalToSurface;
            final switch (bndary.which_boundary) {
            case Face.north:
                j = blk.jmax;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        cell = blk.get_cell(i,j,k);
                        IFace = cell.iface[Face.north];
                        cellWidthNormalToSurface = cell.jLength;
                        compute_and_store_loads(IFace, cellWidthNormalToSurface, sim_time, fname);
                    } // end i loop
                } // end for k
                break;
            case Face.east:
                i = blk.imax;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        cell = blk.get_cell(i,j,k);
                        IFace = cell.iface[Face.east];
                        cellWidthNormalToSurface = cell.iLength;
                        compute_and_store_loads(IFace, cellWidthNormalToSurface, sim_time, fname);
                    } // end j loop
                } // end for k
                break;
            case Face.south:
                j = blk.jmin;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        cell = blk.get_cell(i,j,k);
                        IFace = cell.iface[Face.south];
                        cellWidthNormalToSurface = cell.jLength;                        
                        compute_and_store_loads(IFace, cellWidthNormalToSurface, sim_time, fname);
                    } // end i loop
                } // end for k
                break;
            case Face.west:
                i = blk.imin;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        cell = blk.get_cell(i,j,k);
                        IFace = cell.iface[Face.west];
                        cellWidthNormalToSurface = cell.iLength;
                        compute_and_store_loads(IFace, cellWidthNormalToSurface, sim_time, fname);
                    } // end j loop
                } // end for k
                break;
            case Face.top:
                k = blk.kmax;
                for (i = blk.imin; i <= blk.imax; ++i) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        cell = blk.get_cell(i,j,k);
                        IFace = cell.iface[Face.top];
                        cellWidthNormalToSurface = cell.kLength;
                        compute_and_store_loads(IFace, cellWidthNormalToSurface, sim_time, fname);
                    } // end j loop
                } // end for i
                break;
            case Face.bottom:
                k = blk.kmin;
                for (i = blk.imin; i <= blk.imax; ++i) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        cell = blk.get_cell(i,j,k);
                        IFace = cell.iface[Face.bottom];
                        cellWidthNormalToSurface = cell.kLength;
                        compute_and_store_loads(IFace, cellWidthNormalToSurface, sim_time, fname);
                    } // end j loop
                } // end for i
                break;
            } // end switch which_boundary
        } // end if
    } // end foreach
} // end apply_structured_grid()

void compute_and_store_loads(FVInterface iface, double cellWidthNormalToSurface, double sim_time, string fname)
{
    FlowState fs = iface.fs;
    FlowGradients grad = iface.grad;
    // iface orientation
    double nx = iface.n.x; double ny = iface.n.y; double nz = iface.n.z;
    double t1x = iface.t1.x; double t1y = iface.t1.y; double t1z = iface.t1.z;
    double t2x = iface.t2.x; double t2y = iface.t2.y; double t2z = iface.t2.z;
    // iface properties
    double mu_wall = fs.gas.mu;
    double k_wall = fs.gas.k;
    double P = fs.gas.p;
    double T_wall = fs.gas.T;
    double a_wall = fs.gas.a;
    double rho_wall = fs.gas.rho;
    double sigma_wall, tau_wall, l_tau, m_tau, n_tau, q, Re_wall, nu_wall, u_star, y_plus;
    if (GlobalConfig.viscous) {
        double dTdx = grad.T[0]; double dTdy = grad.T[1]; double dTdz = grad.T[2]; 
        double dudx = grad.vel[0][0]; double dudy = grad.vel[0][1]; double dudz = grad.vel[0][2];
        double dvdx = grad.vel[1][0]; double dvdy = grad.vel[1][1]; double dvdz = grad.vel[1][2];
        double dwdx = grad.vel[2][0]; double dwdy = grad.vel[2][1]; double dwdz = grad.vel[2][2];
        // compute heat load
        double dTdn = dTdx*nx + dTdy*ny + dTdz*nz; // dot product
        q = k_wall * dTdn; // heat load (positive sign means heat flows to the wall)
        // compute stress tensor at interface in global reference frame
        double lmbda = -2.0/3.0 * mu_wall;
        double tau_xx = 2.0*mu_wall*dudx + lmbda*(dudx + dvdy + dwdz);
        double tau_yy = 2.0*mu_wall*dvdy + lmbda*(dudx + dvdy + dwdz);
        double tau_zz = 2.0*mu_wall*dwdz + lmbda*(dudx + dvdy + dwdz);
        double tau_xy = mu_wall * (dudy + dvdx);
        double tau_xz = mu_wall * (dudz + dwdx);
        double tau_yz = mu_wall * (dvdz + dwdy);
        double sigma_x = P + tau_xx;
        double sigma_y = P + tau_yy;
        double sigma_z = P + tau_zz;
        // compute direction cosines for interface normal
        // [TODO] Kyle, shouldn't the nx, ny and nz values already be direction cosines? PJ 2017-09-09
        double l = nx / sqrt(nx*nx + ny*ny + nz*nz);
        double m = ny / sqrt(nx*nx + ny*ny + nz*nz);
        double n = nz / sqrt(nx*nx + ny*ny + nz*nz);
        // transform stress tensor -- we only need stress on a single surface
        // we can avoid performing the entire transformation by following
        // the procedure in Roark's Formulas for Stress and Strain (Young and Budynas) pg. 21.
        sigma_wall = sigma_x*l*l+sigma_y*m*m+sigma_z*n*n+2.0*tau_xy*l*m+2.0*tau_yz*m*n+2.0*tau_xz*n*l;
        tau_wall = sqrt((sigma_x*l+tau_xy*m+tau_xz*n)*(sigma_x*l+tau_xy*m+tau_xz*n)
                        +(tau_xy*l+sigma_y*m+tau_yz*n)*(tau_xy*l+sigma_y*m+tau_yz*n)
                        +(tau_xz*l+tau_yz*m+sigma_z*n)*(tau_xz*l+tau_yz*m+sigma_z*n)
                        -sigma_wall*sigma_wall);
        // tau_wall directional cosines
        l_tau = 1.0/tau_wall * ((sigma_x - sigma_wall)*l+tau_xy*m+tau_xz*n);
        m_tau = 1.0/tau_wall * (tau_xy*l+(sigma_y - sigma_wall)*m+tau_yz*n);
        n_tau = 1.0/tau_wall * (tau_xz*l+tau_yz*m+(sigma_z-sigma_wall)*n);
        // compute y+
        nu_wall = mu_wall / rho_wall;
        u_star = sqrt(tau_wall / rho_wall);
        y_plus = u_star*(cellWidthNormalToSurface/2.0) / nu_wall;
        // compute cell Reynolds number
        Re_wall = rho_wall * a_wall * cellWidthNormalToSurface / mu_wall;
    } else {
        // For an inviscid simulation, we have only pressure.
        sigma_wall = P;
        tau_wall = 0.0;
        l_tau = 0.0; m_tau = 0.0; n_tau = 0.0;
        q = 0.0;
        Re_wall = 0.0;
        y_plus = 0.0;
    }
    // store in file
    auto writer = format("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                         iface.pos.x, iface.pos.y, iface.pos.z, iface.area[0], q, tau_wall, l_tau, m_tau, n_tau, sigma_wall, nx, ny, nz, T_wall, Re_wall, y_plus, cellWidthNormalToSurface);
    std.file.append(fname, writer);    
} // end compute_and_store_loads()
