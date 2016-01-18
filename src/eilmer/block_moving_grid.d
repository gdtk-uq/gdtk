// block_moving_grid.d
// Module for implementing a moving grid in Eilmer4.
// Kyle D. original implmentation (moving grid, shock fitting) Nov 2015.

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
import boundary_flux_effect;

immutable double SHOCK_DETECT_THRESHOLD = 0.2;

int set_gcl_interface_properties(Block blk, size_t gtl, double dt) {
    size_t i, j, k;
    FVVertex vtx1, vtx2;
    FVInterface IFace;
    Vector3 pos1, pos2, temp;
    Vector3 averaged_ivel, vol;
    k = blk.kmin;
    // loop over i-interfaces and compute interface velocity wif'.
    for (j = blk.jmin; j <= blk.jmax; ++j) {
	for (i = blk.imin; i <= blk.imax+1; ++i) {//  i <= blk.imax+1
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
	    if ( blk.myConfig.axisymmetric == false ) vol = 0.5 * cross( pos1, pos2 );
	    else vol = 0.5 * cross( pos1, pos2 ) * ( ( vtx1.pos[gtl].y + vtx1.pos[0].y + vtx2.pos[gtl].y + vtx2.pos[0].y ) /4.0 );
	    temp = vol / ( dt * IFace.area[0] );
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
    for (j = blk.jmin; j <= blk.jmax+1; ++j) { //j <= blk.jmax+
	for (i = blk.imin; i <= blk.imax; ++i) {
	    vtx1 = blk.get_vtx(i,j,k);
	    vtx2 = blk.get_vtx(i+1,j,k);
	    IFace = blk.get_ifj(i,j,k);
	    pos1 = vtx2.pos[gtl] - vtx1.pos[0];
	    pos2 = vtx1.pos[gtl] - vtx2.pos[0];
	    averaged_ivel = (vtx1.vel[0] + vtx2.vel[0]) / 2.0;
	    if ( blk.myConfig.axisymmetric == false ) vol = 0.5 * cross( pos1, pos2 );
	    else vol = 0.5 * cross( pos1, pos2 ) * ( ( vtx1.pos[gtl].y + vtx1.pos[0].y + vtx2.pos[gtl].y + vtx2.pos[0].y ) /4.0 );
	    temp = vol / ( dt * IFace.area[0] );
	    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
	    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);	    
	    IFace.gvel.refx = temp.z;
	    IFace.gvel.refy = averaged_ivel.y;
	    IFace.gvel.refz = averaged_ivel.z;
	    if ( blk.myConfig.axisymmetric && j == blk.jmin ) {
		// For axis symmetric cases the cells along the axis of symmetry have 0 interface area,
		// this is a problem for determining Wif, so we have to catch the NaN from dividing by 0.
		// We choose to set the y and z directions to 0, but take an averaged value for the
		// x-direction so as to not force the grid to be stationary, defeating the moving grid's purpose.
		IFace.gvel.refx = averaged_ivel.x; IFace.gvel.refy = 0.0; IFace.gvel.refz = 0.0;
	    }
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

void shock_fitting_vertex_velocities(Block blk, size_t dimensions, int step, double sim_time) {
    /++ for a given block, loop through cell vertices and update the vertex
      + velocities. The boundary vertex velocities are set via the methodology laid out
      + in Ian Johnston's thesis available on the cfcfd website. The internal velocities
      + are then assigned based off the boundary velocity and a user chosen weighting.
      + NB: the current implementation is hard coded for a moving WEST boundary
      + NB: Only for Euler stepping currently (i.e. gtl = 0 for velocities/positions)
      ++/
    size_t krangemax = ( dimensions == 2 ) ? blk.kmax : blk.kmax+1;
    // make sure all the vertices are given a velocity to begin with
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
	for ( size_t j = blk.jmin; j <= blk.jmax+2; ++j ) {    
	    for ( size_t i = blk.imin; i <= blk.imax+2; ++i ) {
		FVVertex vtx = blk.get_vtx(i,j,k);
		vtx.vel[0] =  Vector3(0.0, 0.0, 0.0);
	    }
	}
    }

    // let's give the shock some time to form before searching for it
    if (sim_time < GlobalConfig.shock_fitting_delay || blk.bc[Face.west].type != "\"inflow_shock_fitting\"") return;
    
    // inflow is a reference to a duplicated flow state, it should be treated as read-only (do not alter)
    auto constFluxEffect = cast(BFE_ConstFlux) blk.bc[Face.west].postConvFluxAction[0];
    auto inflow = constFluxEffect.fstate;
    
    // #####-------------------------- SHOCK SEARCH --------------------------##### //
    FVVertex vtx, vtx_left, vtx_right;
    FVInterface iface_neighbour;
    FVCell cell_toprght, cell_botrght, cell;
    Vector3 temp_vel;
    Vector3 unit_d;    // this will be a unit vector which points along a radial line toward the EAST boundary
    double rho, shock_detect, temp1, temp2, ws1, ws2, rho_lft, rho_rght,  p_lft, p_rght, M;
    Vector3[4] interface_ws;
    double[4] w;
    Vector3 u_lft, u_rght, ns, tav;
    // First update all the WEST boundary vertex velocities (these are the masters)
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
	for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
	    size_t i = blk.imin;
	    vtx = blk.get_vtx(i,j,k);
	    vtx_left = blk.get_vtx(i-1,j,k);
	    vtx_right = blk.get_vtx(i+1,j,k);
	    cell_toprght = blk.get_cell(i, j, k);            
	    cell_botrght = blk.get_cell(i, j-1, k);
	    rho = 0.5*(cell_toprght.fs.gas.rho+cell_botrght.fs.gas.rho);
	    if (j == blk.jmin && blk.bc[Face.south].type != "\"exchange_over_full_face\"") rho = cell_toprght.fs.gas.rho;   // if on an outer boundary just take
	    if (j == blk.jmax+1 && blk.bc[Face.north].type != "\"exchange_over_full_face\"") rho = cell_botrght.fs.gas.rho; // the first internal cell
	    shock_detect = abs(inflow.gas.rho - rho)/fmax(inflow.gas.rho, rho);
	    if (shock_detect < SHOCK_DETECT_THRESHOLD) { // no shock across boundary set vertex velocity to wave speed (u+a)
		temp_vel.refx = (inflow.vel.x + inflow.gas.a);
		temp_vel.refy = -1.0*(inflow.vel.y + inflow.gas.a);
		temp_vel.refz = 0.0;
	    }
	    else { // shock_detect > 0.2  shock detected across boundary
		// loop over cells which neighbour current vertex and calculate wave speed at interfaces
		// left of the boundary is the constant flux condition
		rho_lft = inflow.gas.rho;       // density in left cell
		u_lft = inflow.vel;             // velocity vector in left cell
		p_lft = inflow.gas.p;           // pressure in left cell
		// loop over neighbouring interfaces (above vtx first, then below vtx)
		for (int face = 0; face < 2; face++) {
		    cell =  blk.get_cell(i, j-face, k);
		    iface_neighbour = blk.get_ifi(i,j-face,k);
		    if (j == blk.jmin && blk.bc[Face.south].type =="\"exchange_over_full_face\"" && face == 1)  cell = gasBlocks[0].get_cell(gasBlocks[0].imin, gasBlocks[0].jmax+1-1, k);
		    if (j == blk.jmax+1 && blk.bc[Face.north].type=="\"exchange_over_full_face\"" && face == 0) cell = gasBlocks[1].get_cell(gasBlocks[1].imin, gasBlocks[1].jmin, k);
		    if (j == blk.jmin && blk.bc[Face.south].type =="\"exchange_over_full_face\"" && face == 1)  iface_neighbour = gasBlocks[0].get_ifi(gasBlocks[0].imin, gasBlocks[0].jmax+1-1, k);
		    if (j == blk.jmax+1 && blk.bc[Face.north].type=="\"exchange_over_full_face\"" && face == 0) iface_neighbour = gasBlocks[1].get_ifi(gasBlocks[1].imin, gasBlocks[1].jmin+1, k);
		    rho_rght = cell.fs.gas.rho;         // density in top right cell
		    u_rght = cell.fs.vel;               // velocity vector in right cell
		    p_rght = cell.fs.gas.p;             // pressure in right cell 
		    ns =  cell.iface[Face.west].n;      // normal to the shock front (taken to be WEST face normal of right cell)
		    ws1 = (rho_lft*dot(u_lft, ns) - rho_rght*dot(u_rght, ns))/(rho_lft - rho_rght);
		    temp1 = ((p_rght - p_lft)/(abs(p_rght - p_lft)))/rho_lft; // just need the sign of p_rght - p_lft
		    temp2 =  sqrt(abs((p_rght - p_lft)/(1/rho_lft - 1/rho_rght)));
		    ws2 = dot(u_lft, ns) - temp1 * temp2;
		    interface_ws[face] = (0.5*ws1 + (1-0.5)*ws2)*ns;
		    // tav is a unit vector which points from the neighbouring interface to the current vertex
		    tav = (vtx.pos[0]-iface_neighbour.pos)/sqrt(dot(vtx.pos[0] - iface_neighbour.pos, vtx.pos[0] - iface_neighbour.pos));
		    M = dot(0.5*(u_rght+u_lft), tav)/cell.fs.gas.a; // equation explained in Ian Johnston's thesis on page 77, note...
		    // we are currently just using the right cell (i.e. first non-ghost cell) as the "post-shock" value, for higher accuracy
		    // we will need to update this with the right hand side reconstructed value.
		    //w[face] = ( M + abs(M) ) / 2;  // alternate weighting 
		    if (M <= 1.0) w[face] = ((M+1)*(M+1)+(M+1)*abs(M+1))/8.0;
		    else w[face] = M;
		    if (j == blk.jmin && blk.bc[Face.south].type != "\"exchange_over_full_face\"") w[1] = 0.0, interface_ws[1] = Vector3(0.0, 0.0, 0.0); // south boundary vertex has no bottom neighbour
		    if (j == blk.jmax+1 && blk.bc[Face.north].type != "\"exchange_over_full_face\"") w[0] = 0.0, interface_ws[0] = Vector3(0.0, 0.0, 0.0); // north boundary vertex has no top neighbour
		}
		// now that we have the surrounding interface velocities, let's combine them to approximate the central vertex velocity
		if (w[0] == 0.0 && w[1] == 0.0) w[0] = 1.0, w[1] = 1.0; // prevents a division by zero. Reverts back to unweighted average
		temp_vel = 0.8*(w[0] * interface_ws[0] + w[1] * interface_ws[1]) / (w[0] + w[1] ); // this is the vertex velocity, dampened to 80% for stabilit
	    }
	    unit_d = correct_direction(unit_d, vtx.pos[0], vtx_left.pos[0], vtx_right.pos[0], i, blk.imin, blk.imax);
	    temp_vel = dot(temp_vel, unit_d)*unit_d;
	    vtx.vel[0] = temp_vel;
	}
    }
    // Next update the internal vertex velocities (these are slaves dervied from the WEST boundary master velocities) 
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
	for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
	    for ( size_t i = blk.imin+1; i <= blk.imax+1; ++i ) {
		vtx = blk.get_vtx(i,j,k);
		vtx_left = blk.get_vtx(i-1,j,k);
		vtx_right = blk.get_vtx(i+1,j,k);
		temp_vel = weighting_function(blk.get_vtx(blk.imin, j, k).vel[0],blk.imax+1,i);
		unit_d = correct_direction(unit_d, vtx.pos[0], vtx_left.pos[0], vtx_right.pos[0], i, blk.imin, blk.imax);
		temp_vel = dot(temp_vel, unit_d)*unit_d;
		vtx.vel[0] = temp_vel;
	    }
	}
    }
    return;
}

Vector3 correct_direction(Vector3 unit_d, Vector3 pos, Vector3 left_pos, Vector3 right_pos, size_t i, size_t imin, size_t imax) {
    // as Ian Johnston recommends in his thesis, we force the vertices to move along "rails" which are the radial lines
    // originating from the EAST boundary spanning to the west
    Vector3 lft_temp;
    Vector3 rght_temp;
    lft_temp = (pos-left_pos)/sqrt(dot(left_pos - pos, left_pos - pos));
    rght_temp = (right_pos-pos)/sqrt(dot(right_pos-pos, right_pos-pos));
    unit_d = 0.5 * (lft_temp + rght_temp); // take an average of the left and right unit vectors
    if (i == imin) unit_d = rght_temp;     // west boundary vertex has no left neighbour
    if (i == imax+1) unit_d = lft_temp;    // east boundary vertex has no right neighbour
    return unit_d;
}

Vector3 weighting_function(Vector3 vel_max, size_t imax, size_t i) {
    // TODO: user chooses weighting function. Currently just linear.
    Vector3 vel = -(vel_max/(imax-2)) * (i-2) + vel_max; // y = mx + c
    return vel;
}
