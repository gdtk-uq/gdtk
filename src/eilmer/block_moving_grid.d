// block_moving_grid.d
// Module for implementing a moving grid in Eilmer4.
// Kyle D. original implmentation (moving grid, shock fitting) Nov 2015.

module block_moving_grid;

import std.conv;
import std.stdio;
import std.math;
import nm.complex;
import nm.number;
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
import fluidblock;
import sfluidblock;
import geom;

@nogc
int set_gcl_interface_properties(SFluidBlock blk, size_t gtl, double dt) {
    size_t i, j, k;
    FVInterface IFace;
    Vector3 pos1, pos2, temp;
    Vector3 averaged_ivel, vol;
    if (blk.myConfig.dimensions == 2) {
    FVVertex vtx1, vtx2;
    k = blk.kmin;
    // loop over i-interfaces and compute interface velocity wif'.
    for (j = blk.jmin; j <= blk.jmax; ++j) {
        for (i = blk.imin; i <= blk.imax+1; ++i) {//  i <= blk.imax+1
            vtx1 = blk.get_vtx(i,j,k);
            vtx2 = blk.get_vtx(i,j+1,k);
            IFace = blk.get_ifi(i,j,k);         
            pos1 = vtx1.pos[gtl];
            pos1 -= vtx2.pos[0];
            pos2 = vtx2.pos[gtl];
            pos2 -= vtx1.pos[0];
            averaged_ivel = vtx1.vel[0];
            averaged_ivel += vtx2.vel[0];
            averaged_ivel.scale(0.5);
            // Use effective edge velocity
            // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
            // Full Potential and Euler solutions for transonic unsteady flow
            // Aeronautical Journal November 1994 Eqn 25
            cross(vol, pos1, pos2);
            if (blk.myConfig.axisymmetric == false) {
                // vol = 0.5*cross(pos1, pos2);
                vol.scale(0.5);
            } else {
                // vol=0.5*cross(pos1, pos2)*((vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y)/4.0);
                vol.scale(0.125*(vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y));
            }
            temp = vol; temp /= dt*IFace.area[0];
            // temp is the interface velocity (W_if) from the GCL
            // interface area determined at gtl 0 since GCL formulation
            // recommends using initial interfacial area in calculation.
            IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            IFace.gvel.set(temp.z, averaged_ivel.y, averaged_ivel.z);
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
            pos1 = vtx2.pos[gtl]; pos1 -= vtx1.pos[0];
            pos2 = vtx1.pos[gtl]; pos2 -= vtx2.pos[0];
            averaged_ivel = vtx1.vel[0]; averaged_ivel += vtx2.vel[0]; averaged_ivel.scale(0.5);
            cross(vol, pos1, pos2);
            if (blk.myConfig.axisymmetric == false) {
                // vol=0.5*cross( pos1, pos2 );
                vol.scale(0.5);
            } else {
                // vol=0.5*cross(pos1, pos2)*((vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y)/4.0);
                vol.scale(0.125*(vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y));
            }
            IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);  
            if (blk.myConfig.axisymmetric && j == blk.jmin && IFace.area[0] == 0.0) {
                // For axi-symmetric cases the cells along the axis of symmetry have 0 interface area,
                // this is a problem for determining Wif, so we have to catch the NaN from dividing by 0.
                // We choose to set the y and z directions to 0, but take an averaged value for the
                // x-direction so as to not force the grid to be stationary, defeating the moving grid's purpose.
                IFace.gvel.set(averaged_ivel.x, to!number(0.0), to!number(0.0));
            } else {
                temp = vol; temp /= dt*IFace.area[0];
                IFace.gvel.set(temp.z, averaged_ivel.y, averaged_ivel.z);
            }
            averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
            IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } 
    }

    // Do 3-D cases, where faces move and create a new hexahedron
    if (blk.myConfig.dimensions == 3) {
        FVVertex vtx0, vtx1, vtx2, vtx3;
        Vector3 p0, p1, p2, p3, p4, p5, p6, p7;
        Vector3 centroid_hex, sub_centroid;
        number volume, sub_volume, temp2;
        // loop over i-interfaces and compute interface velocity wif'.
        for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (j = blk.jmin; j <= blk.jmax; ++j) {
                for (i = blk.imin; i <= blk.imax+1; ++i) {//  i <= blk.imax+1
                    // Calculate volume generated by sweeping face 0123 from pos[0] to pos[gtl]
                    vtx0 = blk.get_vtx(i,j  ,k  );
                    vtx1 = blk.get_vtx(i,j+1,k  );
                    vtx2 = blk.get_vtx(i,j+1,k+1);
                    vtx3 = blk.get_vtx(i,j  ,k+1  );
                    p0 = vtx0.pos[0]; p1 = vtx1.pos[0];
                    p2 = vtx2.pos[0]; p3 = vtx3.pos[0];
                    p4 = vtx0.pos[gtl]; p5 = vtx1.pos[gtl];
                    p6 = vtx2.pos[gtl]; p7 = vtx3.pos[gtl];
                    // use 6x pyramid approach as used to calculate internal volume of hex cells
                    centroid_hex.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    volume = 0.0; 
                    pyramid_properties(p6, p7, p3, p2, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p5, p6, p2, p1, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p4, p5, p1, p0, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p4, p0, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p6, p5, p4, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p0, p1, p2, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    //
                    IFace = blk.get_ifi(i,j,k);         
                    averaged_ivel = vtx0.vel[0];
                    averaged_ivel += vtx1.vel[0];
                    averaged_ivel += vtx2.vel[0];
                    averaged_ivel += vtx3.vel[0];
                    averaged_ivel.scale(0.25);
                    // Use effective face velocity, analoguous to edge velocity concept
                    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
                    // Full Potential and Euler solutions for transonic unsteady flow
                    // Aeronautical Journal November 1994 Eqn 25
                    temp2 = volume; temp /= dt*IFace.area[0];
                    // temp2 is the interface velocity (W_if) from the GCL
                    // interface area determined at gtl 0 since GCL formulation
                    // recommends using initial interfacial area in calculation.
                    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    IFace.gvel.set(temp2, averaged_ivel.y, averaged_ivel.z);
                    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
                    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);  
                }
            }
        }
        // loop over j-interfaces and compute interface velocity wif'.
        for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (j = blk.jmin; j <= blk.jmax+1; ++j) { //j <= blk.jmax+
                for (i = blk.imin; i <= blk.imax; ++i) {
                    // Calculate volume generated by sweeping face 0123 from pos[0] to pos[gtl]
                    vtx0 = blk.get_vtx(i  ,j,k  );
                    vtx1 = blk.get_vtx(i  ,j,k+1);
                    vtx2 = blk.get_vtx(i+1,j,k+1);
                    vtx3 = blk.get_vtx(i+1,j,k  );
                    p0 = vtx0.pos[0]; p1 = vtx1.pos[0];
                    p2 = vtx2.pos[0]; p3 = vtx3.pos[0];
                    p4 = vtx0.pos[gtl]; p5 = vtx1.pos[gtl];
                    p6 = vtx2.pos[gtl]; p7 = vtx3.pos[gtl];
                    // use 6x pyramid approach as used to calculate internal volume of hex cells
                    centroid_hex.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    volume = 0.0; 
                    pyramid_properties(p6, p7, p3, p2, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p5, p6, p2, p1, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p4, p5, p1, p0, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p4, p0, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p6, p5, p4, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p0, p1, p2, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    // 
                    IFace = blk.get_ifj(i,j,k);         
                    averaged_ivel = vtx0.vel[0];
                    averaged_ivel += vtx1.vel[0];
                    averaged_ivel += vtx2.vel[0];
                    averaged_ivel += vtx3.vel[0];
                    averaged_ivel.scale(0.25);
                    // Use effective face velocity, analoguous to edge velocity concept
                    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
                    // Full Potential and Euler solutions for transonic unsteady flow
                    // Aeronautical Journal November 1994 Eqn 25
                    temp2 = volume; temp /= dt*IFace.area[0];
                    // temp2 is the interface velocity (W_if) from the GCL
                    // interface area determined at gtl 0 since GCL formulation
                    // recommends using initial interfacial area in calculation.
                    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    IFace.gvel.set(temp2, averaged_ivel.y, averaged_ivel.z);
                    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
                    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);  
                }
            } 
        }
        // loop over k-interfaces and compute interface velocity wif'.
        for (k = blk.kmin; k <= blk.kmax+1; ++k) {
            for (j = blk.jmin; j <= blk.jmax; ++j) { //j <= blk.jmax+
                for (i = blk.imin; i <= blk.imax; ++i) {
                    // Calculate volume generated by sweeping face 0123 from pos[0] to pos[gtl]
                    vtx0 = blk.get_vtx(i  ,j  ,k);
                    vtx1 = blk.get_vtx(i+1,j  ,k);
                    vtx2 = blk.get_vtx(i+1,j+1,k);
                    vtx3 = blk.get_vtx(i  ,j+1,k);
                    p0 = vtx0.pos[0]; p1 = vtx1.pos[0];
                    p2 = vtx2.pos[0]; p3 = vtx3.pos[0];
                    p4 = vtx0.pos[gtl]; p5 = vtx1.pos[gtl];
                    p6 = vtx2.pos[gtl]; p7 = vtx3.pos[gtl];
                    // use 6x pyramid approach as used to calculate internal volume of hex cells
                    centroid_hex.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    volume = 0.0; 
                    pyramid_properties(p6, p7, p3, p2, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p5, p6, p2, p1, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p4, p5, p1, p0, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p4, p0, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p6, p5, p4, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p0, p1, p2, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    //
                    IFace = blk.get_ifk(i,j,k);         
                    averaged_ivel = vtx0.vel[0];
                    averaged_ivel += vtx1.vel[0];
                    averaged_ivel += vtx2.vel[0];
                    averaged_ivel += vtx3.vel[0];
                    averaged_ivel.scale(0.25);
                    // Use effective face velocity, analoguous to edge velocity concept
                    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
                    // Full Potential and Euler solutions for transonic unsteady flow
                    // Aeronautical Journal November 1994 Eqn 25
                    temp2 = volume; temp2 /= dt*IFace.area[0];
                    // temp2 is the interface velocity (W_if) from the GCL
                    // interface area determined at gtl 0 since GCL formulation
                    // recommends using initial interfacial area in calculation.
                    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    IFace.gvel.set(temp2, averaged_ivel.y, averaged_ivel.z);
                    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
                    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);  
                }
            } 
        }
    } 
    return 0;
}

void predict_vertex_positions(SFluidBlock blk, size_t dimensions, double dt, int gtl) {
    size_t krangemax = ( dimensions == 2 ) ? blk.kmax : blk.kmax+1;
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
            for ( size_t i = blk.imin; i <= blk.imax+1; ++i ) {
                FVVertex vtx = blk.get_vtx(i,j,k);
                if (gtl == 0) // if this is the predictor/sole step then update grid
                    vtx.pos[1] = vtx.pos[0] + dt *  vtx.vel[0];
                else   // else if this is the corrector step keep then grid fixed
                    vtx.pos[2] = vtx.pos[1];
            }
        }
    }
    return;
}

void shock_fitting_vertex_velocities(SFluidBlock blk, int step, double sim_time) {
    /++ for a given block, loop through cell vertices and update the vertex
      + velocities. The boundary vertex velocities are set via the methodology laid out
      + in Ian Johnston's thesis available on the cfcfd website. The internal velocities
      + are then assigned based on the boundary velocity and a user chosen weighting.
      + NB: Implementation is hard coded for a moving WEST boundary
      ++/
    FVVertex vtx, vtx_left, vtx_right;
    FVInterface iface_neighbour;
    FVCell cell_toprght, cell_botrght, cell, top_cell_R1, top_cell_R2, bot_cell_R1, bot_cell_R2, cell_R1, cell_R2;
    Vector3 temp_vel, unit_d, u_lft, u_rght, ns, tav; 
    number rho, shock_detect, temp1, temp2, ws1, ws2, rho_lft, rho_rght, p_lft, p_rght, M, rho_recon_top, rho_recon_bottom;
    Vector3[4] interface_ws;
    number[4] w;
    int interpolation_order = blk.myConfig.shock_fitting_interpolation_order;
    immutable double SHOCK_DETECT_THRESHOLD =  0.2;
    immutable double VTX_VEL_SCALING_FACTOR = blk.myConfig.shock_fitting_scale_factor;
    size_t krangemax = ( blk.myConfig.dimensions == 2 ) ? blk.kmax : blk.kmax+1;

    // make sure all the vertices are given a velocity to begin with
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+2; ++j ) {    
            for ( size_t i = blk.imin; i <= blk.imax+2; ++i ) {
                vtx = blk.get_vtx(i,j,k);
                vtx.vel[0].clear();
            }
        }
    }
   
    // let's give the shock some time to form before searching for it
    if (sim_time < GlobalConfig.shock_fitting_delay || blk.bc[Face.west].type != "inflow_shock_fitting") return;
    // inflow is a reference to a duplicated flow state, it should be treated as read-only (do not alter)
    auto constFluxEffect = cast(BFE_ConstFlux) blk.bc[Face.west].postConvFluxAction[0];
    auto inflow = constFluxEffect.fstate;
    number U_plus_a = geom.abs(inflow.vel) + inflow.gas.a; // inflow gas wave speed
    // #####-------------------------- SHOCK SEARCH --------------------------##### //
    // First update all the WEST boundary vertex velocities
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
            size_t i = blk.imin;
            vtx = blk.get_vtx(i,j,k);
            vtx_left = blk.get_vtx(i-1,j,k);
            vtx_right = blk.get_vtx(i+1,j,k);
            /++ the next several lines are necessarily bulky, what we are doing here is estimating
             ++ a density value (rho) to compare with the inflow density for the shock detector.
             ++ This estimation can take on many different forms depending on whether the current
             ++ vertex is 1. on the flow domain boundary, 2. on a block edge, 3. neither (typical vertex). 
             ++ The following code determines the vertex position and then applies the necessary method.
             ++/
            // if vtx is on block edge then grab cell data from neighbour block
            if (j == blk.jmin && blk.bc[Face.south].type =="exchange_over_full_face") {
                auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.south].preReconAction[0];
                int neighbourBlock =  ffeBC.neighbourBlock.id;
                auto nbblk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
                cell_botrght = nbblk.get_cell(nbblk.imin, nbblk.jmax, k);
            }
            else // else grab data from neighbour cell in current block
                cell_botrght = blk.get_cell(i, j-1, k);
            // if vtx is on block edge then grab cell data from neighbour block
            if (j == blk.jmax+1 && blk.bc[Face.north].type=="exchange_over_full_face") {
                auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.north].preReconAction[0];
                int neighbourBlock =  ffeBC.neighbourBlock.id;
                auto nbblk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
                cell_toprght = nbblk.get_cell(nbblk.imin, nbblk.jmin, k);
            }
            else // else grab data from neighbour cell in current block
                cell_toprght = blk.get_cell(i, j, k);            
            if (interpolation_order == 2) { 
                if (j == blk.jmin && blk.bc[Face.south].type =="exchange_over_full_face") {
                    // if reconsturction is true and vtx is on block edge grab rhs cell data from neighbour block
                    auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.south].preReconAction[0];
                    int neighbourBlock =  ffeBC.neighbourBlock.id;
                    auto nbblk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                    assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
                    bot_cell_R1 = nbblk.get_cell(nbblk.imin+1, nbblk.jmax, k);
                    bot_cell_R2 = nbblk.get_cell(nbblk.imin+2, nbblk.jmax, k);
                }
                else { // else if reconstruction is true and vtx is not on block edge
                    bot_cell_R1 = blk.get_cell(i+1, j-1, k);
                    bot_cell_R2 = blk.get_cell(i+2, j-1, k);
                }
                if (j == blk.jmax+1 && blk.bc[Face.north].type=="exchange_over_full_face") {
                    // if reconsturction is true and vtx is on block edge grab rhs cell data from neighbour block
                    auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.north].preReconAction[0];
                    int neighbourBlock =  ffeBC.neighbourBlock.id;
                    auto nbblk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                    assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
                    top_cell_R1 = nbblk.get_cell(nbblk.imin+1, nbblk.jmin, k);
                    top_cell_R2 = nbblk.get_cell(nbblk.imin+2, nbblk.jmin, k);
                }
                else { // else if reconstruction is true and vtx is not on block edge
                    top_cell_R1 = blk.get_cell(i+1, j, k);
                    top_cell_R2 = blk.get_cell(i+2, j, k);
                }
                // perform reconstruction
                rho_recon_top = scalar_reconstruction(cell_toprght.fs.gas.rho,  top_cell_R1.fs.gas.rho,
                                                      top_cell_R2.fs.gas.rho, cell_toprght.iLength,
                                                      top_cell_R1.iLength,  top_cell_R2.iLength, inflow.gas.rho);
                rho_recon_bottom = scalar_reconstruction(cell_botrght.fs.gas.rho,  bot_cell_R1.fs.gas.rho,
                                                         bot_cell_R2.fs.gas.rho, cell_botrght.iLength,
                                                         bot_cell_R1.iLength,  bot_cell_R2.iLength, inflow.gas.rho);
            }
            
            // using stored data estimate a value for rho
            // NB: if on flow domain boundary then just use the first internal cell to estimate rho
            if (interpolation_order == 2) { // use reconstruction
                if (j == blk.jmin && blk.bc[Face.south].type != "exchange_over_full_face")
                    rho = rho_recon_top;
                else if (j == blk.jmax+1 && blk.bc[Face.north].type != "exchange_over_full_face")
                    rho = rho_recon_bottom;
                else
                    rho = 0.5*(rho_recon_top + rho_recon_bottom);
            }
            else { 
                if (j == blk.jmin && blk.bc[Face.south].type != "exchange_over_full_face")
                    rho = cell_toprght.fs.gas.rho;
                else if (j == blk.jmax+1 && blk.bc[Face.north].type != "exchange_over_full_face")
                    rho = cell_botrght.fs.gas.rho;
                else
                    rho = 0.5*(cell_toprght.fs.gas.rho+cell_botrght.fs.gas.rho);
            }
            shock_detect = abs(inflow.gas.rho - rho)/fmax(inflow.gas.rho, rho);
            if (shock_detect < SHOCK_DETECT_THRESHOLD) { // no shock across boundary set vertex velocity to wave speed (u+a)
                temp_vel.set(inflow.vel.x + inflow.gas.a, -1.0*(inflow.vel.y+inflow.gas.a), to!number(0.0));
            }
            else { // shock detected across boundary
                // loop over cells which neighbour current vertex and calculate wave speed at interfaces
                // left of the boundary is the constant flux condition
                rho_lft = inflow.gas.rho;       // density
                u_lft = inflow.vel;             // velocity vector
                p_lft = inflow.gas.p;           // pressure
                // loop over neighbouring interfaces (above vtx first, then below vtx)
                foreach (int jOffSet; [0, 1]) {
                    cell =  blk.get_cell(i, j-jOffSet, k);
                    iface_neighbour = blk.get_ifi(i,j-jOffSet,k);
                    // For the case where two blocks connect (either on the south or north faces) we need to reach across to the neighbour block and
                    // retrieve the neighbouring vertex and interface to ensure the vertices on the connecting interface move together and
                    // do not move independently
                    if (j == blk.jmin && blk.bc[Face.south].type =="exchange_over_full_face" && jOffSet == 1) {
                        auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.south].preReconAction[0];
                        int neighbourBlock =  ffeBC.neighbourBlock.id;
                        auto nbblk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                        assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
                        cell = nbblk.get_cell(nbblk.imin, nbblk.jmax, k);
                        iface_neighbour = nbblk.get_ifi(nbblk.imin, nbblk.jmax, k);
                    }
                    if (j == blk.jmax+1 && blk.bc[Face.north].type=="exchange_over_full_face" && jOffSet == 0) {
                        auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.north].preReconAction[0];
                        int neighbourBlock =  ffeBC.neighbourBlock.id;
                        auto nbblk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                        assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
                        cell = nbblk.get_cell(nbblk.imin, nbblk.jmin, k);
                        iface_neighbour = nbblk.get_ifi(nbblk.imin, nbblk.jmin, k);
                    }
                    if (interpolation_order == 2) {
                        cell_R1 = blk.get_cell(i+1, j-jOffSet, k);
                        cell_R2 =  blk.get_cell(i+2, j-jOffSet, k);
                        // As suggested in onedinterp.d we are transforming all cells velocities to a local reference frame relative to the
                        // interface where the reconstruction is taking place
                        cell.fs.vel.transform_to_local_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                        cell_R1.fs.vel.transform_to_local_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                        cell_R2.fs.vel.transform_to_local_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                        rho_rght = scalar_reconstruction(cell.fs.gas.rho,  cell_R1.fs.gas.rho,
                                                         cell_R2.fs.gas.rho, cell.iLength,
                                                         cell_R1.iLength,  cell_R2.iLength, inflow.gas.rho);
                        u_rght.refx = scalar_reconstruction(cell.fs.vel.x,  cell_R1.fs.vel.x,
                                                            cell_R2.fs.vel.x, cell.iLength,
                                                            cell_R1.iLength, cell_R2.iLength, inflow.vel.x*cell.iface[Face.west].n.x);
                        u_rght.refy = scalar_reconstruction(cell.fs.vel.y,  cell_R1.fs.vel.y,
                                                            cell_R2.fs.vel.y, cell.iLength,
                                                            cell_R1.iLength,  cell_R2.iLength, inflow.vel.y*cell.iface[Face.west].n.y);
                        u_rght.refz = scalar_reconstruction(cell.fs.vel.z,  cell_R1.fs.vel.z,
                                                            cell_R2.fs.vel.z, cell.iLength,
                                                            cell_R1.iLength, cell_R2.iLength, inflow.vel.z*cell.iface[Face.west].n.z);
                        p_rght = scalar_reconstruction(cell.fs.gas.p, cell_R1.fs.gas.p,
                                                       cell_R2.fs.gas.p, cell.iLength,
                                                       cell_R1.iLength, cell_R2.iLength, inflow.gas.p);
                        // here we are transforming the cell velocities back to the global reference frame
                        u_rght.transform_to_global_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                        cell.fs.vel.transform_to_global_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                        cell_R1.fs.vel.transform_to_global_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                        cell_R2.fs.vel.transform_to_global_frame(cell.iface[Face.west].n, cell.iface[Face.west].t1, cell.iface[Face.west].t2);
                    }
                    else {
                        rho_rght = cell.fs.gas.rho;         // density in top right cell
                        u_rght = cell.fs.vel;               // velocity vector in right cell
                        p_rght = cell.fs.gas.p;             // pressure in right cell 
                    }
                    ns =  cell.iface[Face.west].n;      // normal to the shock front (taken to be WEST face normal of right cell)
                    ws1 = (rho_lft*dot(u_lft, ns) - rho_rght*dot(u_rght, ns))/(rho_lft - rho_rght);
                    temp1 = sgn(p_rght.re - p_lft.re)/rho_lft; // just need the sign of p_rght - p_lft
                    temp2 =  sqrt(abs((p_rght - p_lft)/(1/rho_lft - 1/rho_rght)));
                    ws2 = dot(u_lft, ns) - temp1 * temp2;
                    interface_ws[jOffSet] = (0.5*ws1 + (1-0.5)*ws2)*ns;
                    // tav is a unit vector which points from the neighbouring interface to the current vertex
                    Vector3 del1 = vtx.pos[0] - iface_neighbour.pos;
                    tav = (del1)/geom.abs(del1);
                    M = dot(u_rght, tav)/cell.fs.gas.a; // equation explained in Ian Johnston's thesis on page 77, note...
                    // we are currently just using the right cell (i.e. first non-ghost cell) as the "post-shock" value, for higher accuracy
                    // we will need to update this with the right hand side reconstructed value.
                    //w[jOffSet] = ( M + abs(M) ) / 2;  // alternate weighting 
                    if (M <= 1.0) w[jOffSet] = ((M+1)*(M+1)+(M+1)*abs(M+1))/8.0;
                    else w[jOffSet] = M;
                }
                if (j == blk.jmin && blk.bc[Face.south].type != "exchange_over_full_face") {
                    // south boundary vertex has no bottom neighbour
                    w[1] = 0.0;
                    interface_ws[1].clear();
                }
                if (j == blk.jmax+1 && blk.bc[Face.north].type != "exchange_over_full_face") {
                    // north boundary vertex has no top neighbour
                    w[0] = 0.0;
                    interface_ws[0].clear();
                }
                // now that we have the surrounding interface velocities, let's combine them to approximate the central vertex velocity
                if (abs(w[0]) < 1.0e-10 && abs(w[1]) < 1.0e-10) {
                    w[0] = 1.0;
                    w[1] = 1.0;
                    // prevents a division by zero. Reverts back to unweighted average
                }
                // this is the vertex velocity, 80% for stability               
                temp_vel =  (w[0] * interface_ws[0] + w[1] * interface_ws[1]) / (w[0] + w[1] );
                if (geom.abs(temp_vel) > U_plus_a) {
                    // safety catch: if an extreme velocity has been assigned let's just limit vel to  u+a
                    temp_vel.set(inflow.vel.x + inflow.gas.a, -1.0*(inflow.vel.y+inflow.gas.a), to!number(0.0));
                }
                
            }
            unit_d = correct_direction(unit_d, vtx.pos[0], vtx_left.pos[0], vtx_right.pos[0], i, blk.imin, blk.imax);
            temp_vel = dot(temp_vel, unit_d)*unit_d;
            vtx.vel[0] = VTX_VEL_SCALING_FACTOR*temp_vel;
        }
    }
    // Next update the internal vertex velocities (these are slaves dervied from the WEST boundary master velocities) 
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
            Vector3 vel_max = blk.get_vtx(blk.imin, j, k).vel[0];
            // Set up a linear weighting on fraction of distance from east boundary back to west boundary.
            number[300] distance;
            assert(300 > blk.imax+1, "my distance array is not big enough");
            distance[blk.imax+1] = 0.0;
            for (size_t i = blk.imax; i >= blk.imin; --i) {
                vtx =  blk.get_vtx(i,j,k);
                vtx_right = blk.get_vtx(i+1,j,k);
                Vector3 delta = vtx.pos[0] - vtx_right.pos[0];
                distance[i] = geom.abs(delta) + distance[i+1];
            }
            number west_distance = distance[blk.imin];
            for (size_t i = blk.imax; i >= blk.imin; --i) {
                distance[i] /= west_distance;
            }
            for ( size_t i = blk.imin+1; i <= blk.imax; ++i ) {
                vtx = blk.get_vtx(i,j,k);
                vtx_left = blk.get_vtx(i-1,j,k);
                vtx_right = blk.get_vtx(i+1,j,k);
                temp_vel = distance[i]*vel_max; // we used to call weighting function here
                unit_d = correct_direction(unit_d, vtx.pos[0], vtx_left.pos[0], vtx_right.pos[0], i, blk.imin, blk.imax);
                temp_vel = dot(temp_vel, unit_d)*unit_d;
                vtx.vel[0] = temp_vel;
            }
            vtx = blk.get_vtx(blk.imax+1,j,k); // east boundary vertex is fixed
            vtx.vel[0].clear();
        }
    }
    return;
}

Vector3 correct_direction(Vector3 unit_d, Vector3 pos, Vector3 left_pos, Vector3 right_pos, size_t i, size_t imin, size_t imax) {
    // as Ian Johnston recommends in his thesis, we force the vertices to move along "rails" which are the radial lines
    // originating from the EAST boundary spanning to the west
    Vector3 lft_temp;
    Vector3 rght_temp;
    Vector3 delL = pos-left_pos;
    Vector3 delR = right_pos-pos;
    lft_temp = delL/geom.abs(delL);
    rght_temp = delR/geom.abs(delR);
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

number scalar_reconstruction(number x1, number x2, number x3, number h1, number h2, number h3, number g1) {
    bool johnston_reconstruction = false;
    double eps = 1.0e-12;
    number reconstructed_value;
    if (johnston_reconstruction) {
        // This is a special one sided reconstruction presented in Ian Johnston's thesis. 
        number delta_theta = 2*(x2-x1)/(h1+h2);
        number delta_cross = 2*(x3-x2)/(h2+h3);
        number kappa = 0.0; // blending parameter
        number s = (2*(delta_cross)*(delta_theta)+eps)/((delta_cross)*(delta_cross)+(delta_theta)*(delta_theta)+eps);
        number delta_1 = s/2 * ((1-s*kappa)*delta_cross + (1+s*kappa)*delta_theta);
        reconstructed_value =  x1 - (delta_1 * 0.5*h1);
    }
    else {
        // linear one-sided reconstruction function 
        number delta = 2.0*(x2-x1)/(h1+h2);
        number r = (x1-g1+eps)/(x2-x1+eps);
        //number phi = (r*r + r)/(r*r+1); // van albada 1
        //number phi = (2*r)/(r*r+1);        // van albada 2
        number phi = (r + abs(r))/(1+abs(r));   // van leer
        reconstructed_value =  x1 - phi*(delta * 0.5*h1);
    }
    return reconstructed_value;
}
