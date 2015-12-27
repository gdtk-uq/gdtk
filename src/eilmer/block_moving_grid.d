// block_moving_grid.d
// Module (collection) of functions for implementing a moving grid in Eilmer4.
// Kyle D. original implmentation (udf moving grid, shock fitting) Nov 2015.

module block_moving_grid;

import std.conv;
import std.stdio;
import std.math;
import std.json;
import util.lua;
import util.lua_service;
import std.string;

import kinetics;
import globalconfig;
import globaldata;
import fvcore;
import fvvertex;
import fvinterface;
import fvcell;
import flowgradients;
import bc;
import block;
import sblock;
import geom;

int set_gcl_interface_properties(Block blk, size_t gtl, double dt) {
    size_t i, j, k;
    FVVertex vtx1, vtx2;
    FVInterface IFace;
    Vector3 pos1, pos2, temp;
    Vector3 averaged_ivel;
    k = blk.kmin;
    // loop over i-interfaces and compute interface velocity wif'.
    for (j = blk.jmin; j <= blk.jmax; ++j) {
	for (i = blk.imin; i <= blk.imax+1; ++i) {
	    vtx1 = blk.get_vtx(i,j,k);
	    vtx2 = blk.get_vtx(i,j+1,k);
	    IFace = blk.get_ifi(i,j,k);   	
	    pos1 = vtx1.pos[gtl] - vtx2.pos[0];
	    pos2 = vtx2.pos[gtl] - vtx1.pos[0];
	    averaged_ivel = (vtx1.vel[0] + vtx2.vel[0]) / 2.0;
	    // Use effective edge velocity
	    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
	    // Full Potential and Euler solutions for transonic unsteady flow
	    // Aeronautical Journal November 1994 Eqn 25
	    temp = 0.5 * cross( pos1, pos2 ) / ( dt * IFace.area[0] );
	    // temp is the interface velocity (W_if) from the GCL
	    // interface area determined at gtl 0 since GCL formulation
	    // recommends using initial interfacial area in calculation.
	    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    IFace.gvel.refx = temp.z;
	    IFace.gvel.refy = averaged_ivel.y;
	    IFace.gvel.refz = averaged_ivel.z;
	    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);	    
	    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);	    			 
	}
    }
    // loop over j-interfaces and compute interface velocity wif'.
    for (j = blk.jmin; j <= blk.jmax+1; ++j) {
	for (i = blk.imin; i <= blk.imax; ++i) {
	    vtx1 = blk.get_vtx(i,j,k);
	    vtx2 = blk.get_vtx(i+1,j,k);
	    IFace = blk.get_ifj(i,j,k);
	    pos1 = vtx2.pos[gtl] - vtx1.pos[0];
	    pos2 = vtx1.pos[gtl] - vtx2.pos[0];
	    averaged_ivel = (vtx1.vel[0] + vtx2.vel[0]) / 2.0;
	    temp = 0.5 * cross( pos1, pos2 ) / ( dt * IFace.area[0] );
	    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);	    
	    IFace.gvel.refx = temp.z;
	    IFace.gvel.refy = averaged_ivel.y;
	    IFace.gvel.refz = averaged_ivel.z;	    
	    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);	    
	    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);	  		       	
	}
    }
    return 0;
}

void predict_vertex_positions(Block blk, size_t dimensions, double dt) {
    size_t krangemax = ( dimensions == 2 ) ? blk.kmax : blk.kmax+1;
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
	for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
	    for ( size_t i = blk.imin; i <= blk.imax+1; ++i ) {
		FVVertex vtx = blk.get_vtx(i,j,k);
		vtx.pos[1] = vtx.pos[0] + dt *  vtx.vel[0];
	    }
	}
    }
    return;
}

