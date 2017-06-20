/**
 * Eilmer4 boundary interface loads writing functions.
 *
 * Author: Kyle Damm
 * First code: 2017-03-17
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
import grid;
import block;
import fvinterface;
import bc;

string loadsDir = "loads";

void init_loads_dir()
{
    ensure_directory_is_present(loadsDir);
}

void write_boundary_loads_to_file(double sim_time, int current_tindx) {
    foreach (blk; gasBlocks) { //foreach (blk; parallel(gasBlocks,1)) {
	if (blk.active) {
	    final switch (blk.grid_type) {
	    case Grid_t.unstructured_grid: 
		apply_unstructured_grid(blk, sim_time, current_tindx);
		break;
	    case Grid_t.structured_grid:
		apply_structured_grid(blk, sim_time, current_tindx);
	    } // end final switch
	} // if (blk.active)
    } // end foreach
} // end write_boundary_loads_to_file

string generate_boundary_load_file(int current_tindx, double sim_time, string group) {
    // generate data file -- naming format tindx_groupname.dat
    string fname = format("%s/t%d-%s.dat", 
			  loadsDir, current_tindx, group);
    // check if file exists already, if not create file
    if (!exists(fname)) {
	auto f = File(fname, "w");
	f.writeln("# t = ", sim_time);
	f.write("# 1:pos.x 2:pos.y 3:pos.z 4:area 5:q 6:tau 7:l_tau 8:m_tau 9:n_tau 10:sigma 11:n.x 12:n.y 13:n.z \n");
	f.close();
    }
    return fname;
}

void apply_unstructured_grid(Block blk, double sim_time, int current_tindx) {
    foreach (bndary; blk.bc) {
	if (canFind(bndary.group, GlobalConfig.boundary_group_for_loads)) {
	    string fname = generate_boundary_load_file(current_tindx, sim_time, bndary.group);
		foreach (iface; bndary.faces) {
		    compute_and_store_loads(iface, sim_time, current_tindx, fname);
		} // end foreach face
	} // end foreach (iface; bndary.faces)
    } // end if (bndary.group != "")
} // foreach (bndary; blk.bc)

void apply_structured_grid(Block blk, double sim_time, int current_tindx) {
    foreach (bndary; blk.bc) {
	if (canFind(bndary.group, GlobalConfig.boundary_group_for_loads)) {
	    string fname = generate_boundary_load_file(current_tindx, sim_time, bndary.group);
	    size_t i, j, k;
	    FVCell cell;
	    FVInterface IFace;
	    final switch (bndary.which_boundary) {
	    case Face.north:
		j = blk.jmax;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (i = blk.imin; i <= blk.imax; ++i) {
			cell = blk.get_cell(i,j,k);
			IFace = cell.iface[Face.north];
			compute_and_store_loads(IFace, sim_time, current_tindx, fname);
		    } // end i loop
		} // end for k
		break;
	    case Face.east:
		i = blk.imax;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			cell = blk.get_cell(i,j,k);
			IFace = cell.iface[Face.east];
			compute_and_store_loads(IFace, sim_time, current_tindx, fname);
		    } // end j loop
		} // end for k
		break;
	    case Face.south:
		j = blk.jmin;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (i = blk.imin; i <= blk.imax; ++i) {
			cell = blk.get_cell(i,j,k);
			IFace = cell.iface[Face.south];
			compute_and_store_loads(IFace, sim_time, current_tindx, fname);
		    } // end i loop
		} // end for k
		break;
	    case Face.west:
		i = blk.imin;
		for (k = blk.kmin; k <= blk.kmax; ++k) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			cell = blk.get_cell(i,j,k);
			IFace = cell.iface[Face.west];
			compute_and_store_loads(IFace, sim_time, current_tindx, fname);
		    } // end j loop
		} // end for k
		break;
	    case Face.top:
		k = blk.kmax;
		for (i = blk.imin; i <= blk.imax; ++i) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			cell = blk.get_cell(i,j,k);
			IFace = cell.iface[Face.top];
			compute_and_store_loads(IFace, sim_time, current_tindx, fname);
		    } // end j loop
		} // end for i
		break;
	    case Face.bottom:
		k = blk.kmin;
		for (i = blk.imin; i <= blk.imax; ++i) {
		    for (j = blk.jmin; j <= blk.jmax; ++j) {
			cell = blk.get_cell(i,j,k);
			IFace = cell.iface[Face.bottom];
			compute_and_store_loads(IFace, sim_time, current_tindx, fname);
		    } // end j loop
		} // end for i
		break;
	    } // end switch which_boundary
	}
    }
}

void compute_and_store_loads(FVInterface iface, double sim_time, int current_tindx, string fname)
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
    double dTtrdx = grad.Ttr[0]; double dTtrdy = grad.Ttr[1]; double dTtrdz = grad.Ttr[2]; 
    double dudx = grad.vel[0][0]; double dudy = grad.vel[0][1]; double dudz = grad.vel[0][2];
    double dvdx = grad.vel[1][0]; double dvdy = grad.vel[1][1]; double dvdz = grad.vel[1][2];
    double dwdx = grad.vel[2][0]; double dwdy = grad.vel[2][1]; double dwdz = grad.vel[2][2];
    // compute heat load
    double dTdn = dTtrdx*nx + dTtrdy*ny + dTtrdz*nz; // dot product
    double q = k_wall * dTdn; // heat load (positive sign means heat flows to the wall)
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
    double l = nx / sqrt(nx*nx + ny*ny + nz*nz);
    double m = ny / sqrt(nx*nx + ny*ny + nz*nz);
    double n = nz / sqrt(nx*nx + ny*ny + nz*nz);
    // transform stress tensor -- we only need stress on a single surface
    // we can avoid performing the entire transformation by following
    // the procedure in Roark's Formulas for Stress and Strain (Young and Budynas) pg. 21.
    double sigma_wall = sigma_x*l*l+sigma_y*m*m+sigma_z*n*n+2.0*tau_xy*l*m+2.0*tau_yz*m*n+2.0*tau_xz*n*l;
    double tau_wall = sqrt((sigma_x*l+tau_xy*m+tau_xz*n)*(sigma_x*l+tau_xy*m+tau_xz*n)
		    +(tau_xy*l+sigma_y*m+tau_yz*n)*(tau_xy*l+sigma_y*m+tau_yz*n)
		    +(tau_xz*l+tau_yz*m+sigma_z*n)*(tau_xz*l+tau_yz*m+sigma_z*n)-sigma_wall*sigma_wall);
    // tau_wall directional cosines
    double l_tau = 1.0/tau_wall * ((sigma_x - sigma_wall)*l+tau_xy*m+tau_xz*n);
    double m_tau = 1.0/tau_wall * (tau_xy*l+(sigma_y - sigma_wall)*m+tau_yz*n);
    double n_tau = 1.0/tau_wall * (tau_xz*l+tau_yz*m+(sigma_z-sigma_wall)*n);
    // store in file
    auto writer = format("%f %f %f %f %f %f %f %f %f %f %f %f %f \n", iface.pos.x, iface.pos.y, iface.pos.z, iface.area[0], q, tau_wall, l_tau, m_tau, n_tau, sigma_wall, nx, ny, nz);
    append(fname, writer);    
} // end compute_and_store_loads()
