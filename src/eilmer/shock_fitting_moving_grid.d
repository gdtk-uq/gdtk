// shock_fitting_moving_grid.d
// Module for implementing shock-fitting with a moving grid in Eilmer4.
// Kyle D. original implmentation (moving grid, shock fitting) Nov 2015.

module shock_fitting_moving_grid;

import std.conv;
import std.stdio;
import std.math;
import nm.complex;
import nm.number;
import std.json;
import util.lua;
import util.lua_service;
import std.string;
import std.algorithm;
import std.datetime;
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

import json_helper;
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
import grid_motion;
import bc;


// Begin MPI_Wait_a_while- copied from PJs program used in the ghost cell data transfer. Prevents the program from locking up if something goes wrong in a process
version(mpi_parallel) {
    int MPI_Wait_a_while(MPI_Request* request, MPI_Status *status) {
        int ierr = 0;
        number timeout = 10000;
        SysTime startTime = Clock.currTime();
        int flag = 0;
        while (!flag) {
            ierr = MPI_Test(request, &flag, status);
            number elapsedtime = (Clock.currTime() - startTime).total!"msecs"();
            if (elapsedtime > timeout) {
                throw new Exception("MPI has timed out");
            }
        }
        return ierr;
    }
}

// The mpi version of the program that will assign a normalised radial distance to each vertex. Could be more efficient, but only runs once so not really an issue.
// Start assign_radial_dist_mpi
version(mpi_parallel) {
    void assign_radial_dist_mpi(SFluidBlock blk) {
        size_t krange = (blk.myConfig.dimensions == 3) ? blk.kmax+1 : blk.kmax;
        MPI_Status MPI_incoming_dist_status, MPI_incoming_partner_status, MPI_incoming_chained_status;
        MPI_Request MPI_incoming_dist_request, MPI_incoming_partner_request, MPI_incoming_chained_request;
        FVVertex vtx, vtx_prev;

        number[] local_dist = new number[(blk.jmax - blk.jmin + 2) * (krange - blk.kmin + 1)], total_dist = new number[local_dist.length];

        // Populate the list with zeros as we will be incrementing the elements- can't leave as nans
        for (size_t i = 0; i < local_dist.length; i++) {
            local_dist[i] = 0.0;
        }
        
        bool requested_data = false, sent_data = false;
        int no_partners, east_neighbour_rank, west_neighbour_rank, east_neighbour, west_neighbour;

        // Assign absolute radial distance within the block
        for (size_t k = blk.kmin; k <= krange; k++) {
            for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                for (size_t i = blk.imax; i >= blk.imin; i--) {
                    vtx = blk.get_vtx!()(i, j, k);
                    vtx_prev = blk.get_vtx!()(i+1, j, k);
                    Vector3 delta = vtx.pos[0] - vtx_prev.pos[0];
                    local_dist[(blk.jmax - blk.jmin + 2) * (k - blk.kmin) + (j - blk.jmin)] += geom.abs(delta);
                    vtx.radial_pos_norm = local_dist[(blk.jmax - blk.jmin + 2) * (k - blk.kmin) + (j - blk.jmin)];
                }
            }
        }
        // Each block needs to go do its own bit of calculation first so we'll post a non-blocking receive for 2 things- the total radial distance up to this block, and the number of blocks
        // run through up to this point.
        if (blk.bc[Face.east].type == "exchange_over_full_face") {
            GhostCellFullFaceCopy ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.east].preReconAction[0];
            east_neighbour = ffeBC.neighbourBlock.id;
            if (find(GlobalConfig.localBlockIds, east_neighbour).empty) {
                
                //Block is non-local in the mpi context, request some data

                requested_data = true;

                int mpi_dist_recv_tag = make_mpi_tag(east_neighbour, 4, 0), mpi_partner_recv_tag = make_mpi_tag(east_neighbour, 4, 1);
                east_neighbour_rank = GlobalConfig.mpi_rank_for_block[east_neighbour];
                int ne = to!int(total_dist.length);
                
                MPI_Irecv(total_dist.ptr, ne, MPI_DOUBLE, east_neighbour_rank, mpi_dist_recv_tag, MPI_COMM_WORLD, &MPI_incoming_dist_request);

                MPI_Irecv(&no_partners, 1, MPI_INT, east_neighbour_rank, mpi_partner_recv_tag, MPI_COMM_WORLD, &MPI_incoming_partner_request);
            }
        }
        else {

            // Block is local in the mpi context 
            for (size_t i = 0; i < local_dist.length; i++) {
                total_dist[i] = 0.0;
            }
        }

        // Wait for the receive request for the distance array to be completed
        if (requested_data) {
            MPI_Wait(&MPI_incoming_dist_request, &MPI_incoming_dist_status);
        }
        // This is where we get the proper radial position of all the vertices by adding the radial distance of the previous blocks to the local vertices 
        for (size_t k = blk.kmin; k <= krange; k++) {
            for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                for (size_t i = blk.imin; i <= blk.imax+1; i++) {
                    vtx = blk.get_vtx!()(i, j, k);
                    vtx.radial_pos_norm += total_dist[(blk.jmax - blk.jmin + 2) * (k - blk.kmin) + (j - blk.jmin)];
                }
            }
        }

        // Create a new array that will be sent away to the next block which contains the radial distance up to the west edge of this block
        number[] total_dist_send = new number[(blk.jmax - blk.jmin + 2) * (krange - blk.kmin + 1)];
        for (size_t m = 0; m < total_dist_send.length; m++) {
            total_dist_send[m] = total_dist[m] + local_dist[m];
        }

        // Wait for the receive request for the size of the chained_to array to be completed
        if (requested_data) {MPI_Wait_a_while(&MPI_incoming_partner_request, &MPI_incoming_partner_status);}

        // Now we know how big the chained_to array is, call for it. We can't request it until we know how many elements it should have
        int[] chained_to = new int[no_partners];
        if (requested_data) {
            int mpi_chained_recv_tag = make_mpi_tag(east_neighbour, 4, 2);
            MPI_Recv(chained_to.ptr, no_partners, MPI_INT, east_neighbour_rank, mpi_chained_recv_tag, MPI_COMM_WORLD, &MPI_incoming_chained_status);
        }

        // Make its changes to the chained_to array
        chained_to.length += 1;
        chained_to[$-1] = blk.id;

        // Increment the number of inflow partners
        no_partners++;

        // Ensure that the number of partners is the same length as the chained_to array
        assert(to!int(chained_to.length) == no_partners);


        // We want to make the send request as soon as possible to other blocks can continue their work
        if (blk.bc[Face.west].type == "exchange_over_full_face") {
            GhostCellFullFaceCopy ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.west].preReconAction[0];
            west_neighbour = ffeBC.neighbourBlock.id;
            if (find(GlobalConfig.localBlockIds, west_neighbour).empty) {
                sent_data = true;
                int ne = to!int(total_dist_send.length);
                int mpi_dist_send_tag = make_mpi_tag(blk.id, 4, 0), mpi_partner_send_tag = make_mpi_tag(blk.id, 4, 1), mpi_chained_send_tag = make_mpi_tag(blk.id, 4, 2);
                west_neighbour_rank = GlobalConfig.mpi_rank_for_block[west_neighbour];

                // Send away all the data; the total distance to the edge of this block, the number of blocks partnered to this one, 
                MPI_Send(total_dist_send.ptr, ne, MPI_DOUBLE, west_neighbour_rank, mpi_dist_send_tag, MPI_COMM_WORLD);

                MPI_Send(&no_partners, 1, MPI_INT, west_neighbour_rank, mpi_partner_send_tag, MPI_COMM_WORLD);
                
                MPI_Send(chained_to.ptr, no_partners, MPI_INT, west_neighbour_rank, mpi_chained_send_tag, MPI_COMM_WORLD);
                
            }
        }

        /++ At this point, every vertex in the simulation should have a non-normalised value for their radial position, and the outer-most block should have the total radial distance along
            each line and the list of all its slaves. Now lets transfer that information back along the radial direction to assign a normalised radial position and give every block a copy
            of those it is chained to.
        ++/

        if (sent_data) {
            // If I sent data on the way out, I must receive data on the way back
            int ne = to!int(total_dist_send.length);
            int mpi_dist_recv_tag = make_mpi_tag(west_neighbour, 3, 3), mpi_partner_recv_tag = make_mpi_tag(west_neighbour, 3, 4);
            MPI_Irecv(total_dist.ptr, ne, MPI_DOUBLE, west_neighbour_rank, mpi_dist_recv_tag, MPI_COMM_WORLD, &MPI_incoming_dist_request);
            MPI_Irecv(&no_partners, 1, MPI_INT, west_neighbour_rank, mpi_partner_recv_tag, MPI_COMM_WORLD, &MPI_incoming_partner_request);
            MPI_Wait(&MPI_incoming_dist_request, &MPI_incoming_dist_status);
        }
        else {
            total_dist = total_dist_send;
        }


        for (size_t k = blk.kmin; k <= krange; k++) {
            for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                for (size_t i = blk.imin; i <= blk.imax+1; i++) {
                    vtx = blk.get_vtx!()(i, j, k);
                    vtx.radial_pos_norm /= total_dist[(blk.jmax - blk.jmin + 2) * (k - blk.kmin) + (j - blk.jmin)];
                    if (j == blk.jmin) {
                    }
                }
            }
        }
        
        if (sent_data) {
            MPI_Wait(&MPI_incoming_partner_request, &MPI_incoming_partner_status);
            int mpi_chained_recv_tag = make_mpi_tag(west_neighbour, 3, 5);
            chained_to.length = no_partners;
            MPI_Recv(chained_to.ptr, no_partners, MPI_INT, west_neighbour_rank, mpi_chained_recv_tag, MPI_COMM_WORLD, &MPI_incoming_chained_status);
        }
         
        if (requested_data) {
            // If I requested data on the way out, I must send data on the way back
            int ne = to!int(total_dist_send.length);
            int mpi_dist_send_tag = make_mpi_tag(blk.id, 3, 3), mpi_partner_send_tag = make_mpi_tag(blk.id, 3, 4), mpi_chained_send_tag = make_mpi_tag(blk.id, 3, 5);
            MPI_Send(total_dist.ptr, ne, MPI_DOUBLE, east_neighbour_rank, mpi_dist_send_tag, MPI_COMM_WORLD);
            MPI_Send(&no_partners, 1, MPI_INT, east_neighbour_rank, mpi_partner_send_tag, MPI_COMM_WORLD);
            MPI_Send(chained_to.ptr, no_partners, MPI_INT, east_neighbour_rank, mpi_chained_send_tag, MPI_COMM_WORLD);
        }
        blk.inflow_partners = chained_to;
    }
}

// End assign_radial_dist_mpi 

//-------------------------------------------------------------------------------------------------------------------------------

// Assigning radial distances in the shared memory version
// Start assign_radial_dist 
void assign_radial_dist(SFluidBlock blk) {
    /++ This will assign all the vertices a normalised radial distance for use
      + with the shock fitting internal velocities.
    ++/
    //Write the radial distance data for the current block
    size_t krange = ( blk.myConfig.dimensions == 3 ) ? blk.kmax+1 : blk.kmax;
    FVVertex vtx, vtx_next;
    Vector3 delta;
    SFluidBlock master_block = blk;
    number radial_dist, block_crossover_value;
    bool last_block;
    for ( size_t k = blk.kmin; k <= krange; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
            radial_dist = 0;
            last_block = false;
            while (last_block == false) {
                for ( size_t i = blk.imin; i <= blk.imax; ++i ) {
                    vtx = blk.get_vtx!()(i, j, k);
                    vtx_next = blk.get_vtx!()(i+1, j, k);
                    delta = vtx_next.pos[0] - vtx.pos[0];
                    radial_dist += geom.abs(delta);
                }
                // We want to know which blocks are slaved to the inflow shock fitted block later
                if (k == blk.kmin && j == blk.jmin) {   //Make sure this only happens once per block
                    master_block.inflow_partners.length += 1;
                    master_block.inflow_partners[$-1] = blk.id;
                }
                
                if (blk.bc[Face.east].type == "exchange_over_full_face") {
                    // There is another block to go through
                    auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.east].preReconAction[0];
                    int neighbourBlock = ffeBC.neighbourBlock.id;
                    auto next_blk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                    blk = next_blk;
                }
                else {
                    // No more blocks to go through
                    last_block = true;
                }
                // Somewhat messy implementation- is there a cleaner way?
            }

            /++ Now we want to go and assign the normalised radial distances; we will run from the
              + east-most boundary to the western shock fitting boundary
            ++/
            block_crossover_value = 0; //It is important we initialise this to 0 for the very first vertex
            last_block = false;
            while (last_block == false) {
                for ( size_t i = blk.imax+1; i >= blk.imin; --i ) {
                    /++ Some special consideration needs to happen at the boundaries here, for both
                      + concerns of information transfer between blocks and consistency between block
                      + boundaries. By setting the east-most vertex of the current block to the same
                      + value as the west-most boundary of the last block, we ensure the boundaries
                      + should never move apart due to having different velocities.
                    ++/
                    if (i == blk.imax+1) {
                        blk.get_vtx!()(i, j, k).radial_pos_norm = block_crossover_value;
                    }
                    else {
                        vtx = blk.get_vtx!()(i, j, k);
                        vtx_next = blk.get_vtx!()(i+1, j, k);
                        delta = vtx.pos[0] - vtx_next.pos[0];
                        vtx.radial_pos_norm = geom.abs(delta) / radial_dist + vtx_next.radial_pos_norm;
                    }
                    if (i == blk.imin) block_crossover_value = blk.get_vtx!()(i, j, k).radial_pos_norm;
                }
                if (blk.bc[Face.west].type == "exchange_over_full_face") {
                    auto ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.west].preReconAction[0];
                    int neighbourBlock =  ffeBC.neighbourBlock.id;
                    auto next_blk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];
                    blk = next_blk;
                }
                else {
                    last_block = true;
                }
            }
        }
    }
}

// End assign_radial_dist
//-----------------------------------------------------------------------------------------------------------------------------

// The meat of the shock fitting program
Vector3[] shock_fitting_vertex_velocities(SFluidBlock blk) {
    /++ for a given block, loop through cell vertices and update the vertex
      + velocities. The boundary vertex velocities are set via the methodology laid out
      + in Ian Johnston's thesis available on the cfcfd website. The internal velocities
      + are then assigned based on the boundary velocity and a user chosen weighting.
      + NB: Implementation is hard coded for a moving WEST boundary
      ++/

    // Lets be verbose about where each of these are used
    FVVertex vtx;        // The vertex we are calculating velocity at
    Vector3[] vertex_master_velocities;
    int interpolation_order = blk.myConfig.shock_fitting_interpolation_order;
    immutable number SHOCK_DETECT_THRESHOLD =  0.2;
    immutable number vtx_vel_stability_factor = blk.myConfig.shock_fitting_scale_factor;
    size_t krangemax = ( blk.myConfig.dimensions == 3 ) ? blk.kmax+1 : blk.kmax;
    // Clear the velocities of all the vertices within the block- this shouldn't be necessary. Take this out once its working and check still works.
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {    
            for ( size_t i = blk.imin; i <= blk.imax+1; ++i ) {
                vtx = blk.get_vtx!()(i,j,k);
                vtx.vel[0].clear();
            }
        }
    }
    // Define the size of the array required- the number of vertices on the western face. The extra +1 is due to indexing starting at 0 but length starting at 1
    vertex_master_velocities.length = ( blk.myConfig.dimensions == 3) ? ((blk.jmax - blk.jmin + 1) * (blk.kmax - blk.kmin + 1) + 1) : (blk.jmax - blk.jmin + 2);
	/++ Layout the naming configuration for the cells and vertices used
			+-------------------------------------------------------------------------------+
			|					|					|					|					|
			|					|					|					|					|
			|					|					|					|					|
			|		inflow		|	    cell    	|	   cell_R1  	|	    cell_R2 	|
			|		(ghost)		|					|		    		|		        	|
			|					|					|					|					|
			|					|					|					|					|
			|					|					|					|					|
	  vtx(i-1,j,k)---------vtx(i,j,k)---------vtx(i+1,j,k)--------------+-------------------+
			|					|					|					|					|
			|					|					|					|					|
			|					|					|					|					|
			|		inflow		|	    cell    	|	   cell_R1	    |	   cell_R2  	|
			|		(ghost)		|					|		        	|		        	|
			|					|					|					|					|
			|					|					|					|					|
			|					|					|					|					|
			+-------------------------------------------------------------------------------+
    	
	  + This shows the reconstruction process in relation to the vertex being shock fitted vtx(i, j, k). 					
	++/


    // let's give the shock some time to form before searching for it. Move this to simcore- the function shouldn't be called
    // if we aren't shock fitting.

    populate_ghost_cell_interface_geometry(blk);

    // inflow is a reference to a duplicated flow state, it should be treated as read-only (do not alter)
    auto constFluxEffect = cast(BFE_ConstFlux) blk.bc[Face.west].postConvFluxAction[0];
    auto inflow = constFluxEffect.fstate;

    // #####-------------------------- SHOCK SEARCH --------------------------##### //

    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
            size_t i = blk.imin;
            vtx = blk.get_vtx!()(i,j,k);

            number rho_east = 0.0, p_east;   
            Vector3 u_east, numerator;    
            // Create the arrays used to store the cell data around the vertex.
            FVCell[4] cell, cell_R1, cell_R2;
            Vector3[4] interface_ws;
            number[4] w;
            number denom;
            Vector3 exact_vel;
            size_t numcells = 0;

            // loop over neighbouring interfaces (above vtx first, then below vtx)
            foreach (int jOffSet; [0, 1]) {
                foreach (kOffSet; (blk.myConfig.dimensions == 3) ? [0, 1] : [0]) {
                    if ((j-jOffSet) < blk.jmin && blk.bc[Face.south].type != "exchange_over_full_face") {} // Do nothing
                    else if ((j-jOffSet) > blk.jmax && blk.bc[Face.north].type != "exchange_over_full_face") {} // Do nothing
                    else {
                        // Normal cell  
                        cell[numcells] = blk.get_cell!()(i, j-jOffSet, k-kOffSet);
                        cell_R1[numcells] = blk.get_cell!()(i+1, j-jOffSet, k-kOffSet);
                        cell_R2[numcells++] = blk.get_cell!()(i+2, j-jOffSet, k-kOffSet);
                    }
                }
            }

            for(i = 0; i < numcells; ++i) {
                rho_east += scalar_reconstruction(cell[i].fs.gas.rho,  cell_R1[i].fs.gas.rho,
                                                cell_R2[i].fs.gas.rho, cell[i].iLength,
                                                cell_R1[i].iLength,  cell_R2[i].iLength, inflow.gas.rho);
            }

            number shock_detect = abs(inflow.gas.rho - (rho_east/numcells))/fmax(inflow.gas.rho, (rho_east/numcells));
            if (shock_detect < SHOCK_DETECT_THRESHOLD) { 
                // We don't think there's a shock here, so set vertex velocity to the maximum, U+a.
                exact_vel.set(inflow.vel.x + inflow.gas.a, -1.0*(inflow.vel.y+inflow.gas.a), (blk.myConfig.dimensions == 3) ? inflow.vel.z+inflow.gas.a: to!number(0.0));
            }
            else {
                // We do have a shock here
                for (i = 0; i < numcells; ++i) {
                    // Reconstruct the density and pressure values at the interface
                    rho_east = scalar_reconstruction(cell[i].fs.gas.rho,  cell_R1[i].fs.gas.rho,
                                                cell_R2[i].fs.gas.rho, cell[i].iLength,
                                                cell_R1[i].iLength,  cell_R2[i].iLength, inflow.gas.rho);
                    p_east = scalar_reconstruction(cell[i].fs.gas.p, cell_R1[i].fs.gas.p,
                                                cell_R2[i].fs.gas.p, cell[i].iLength,
                                                cell_R1[i].iLength, cell_R2[i].iLength, inflow.gas.p);

                    u_east.refx = scalar_reconstruction(cell[i].fs.vel.x,  cell_R1[i].fs.vel.x,
                                                cell_R2[i].fs.vel.x, cell[i].iLength,
                                                cell_R1[i].iLength, cell_R2[i].iLength, inflow.vel.x);
                    u_east.refy = scalar_reconstruction(cell[i].fs.vel.y,  cell_R1[i].fs.vel.y,
                                                cell_R2[i].fs.vel.y, cell[i].iLength,
                                                cell_R1[i].iLength,  cell_R2[i].iLength, inflow.vel.y);
                    u_east.refz = scalar_reconstruction(cell[i].fs.vel.z,  cell_R1[i].fs.vel.z,
                                                cell_R2[i].fs.vel.z, cell[i].iLength,
                                                cell_R1[i].iLength, cell_R2[i].iLength, inflow.vel.z);

                    //Find the shock speeds ws1 and ws2 as outlined on page 73 of Ian Johnston's thesis
                    Vector3 ns, t_dim, t_unit;
                    number sign_ws2, sqrt_ws2, M, mixing_param = 0.5, ws1, ws2;

                    ns =  cell[i].iface[Face.west].n;
                    ws1 = (inflow.gas.rho * dot(inflow.vel, ns) - rho_east * dot(u_east, ns))/(inflow.gas.rho - rho_east);
                    sign_ws2 = sgn(p_east.re - inflow.gas.rho.re)/inflow.gas.rho;
                    sqrt_ws2 =  sqrt(abs((p_east - inflow.gas.p)/(1/inflow.gas.rho - 1/rho_east)));
                    ws2 = dot(inflow.vel, ns) - sign_ws2 * sqrt_ws2;
                    interface_ws[i] = (mixing_param * ws1 + (1-mixing_param) * ws2) * ns;
                    // Create a unit vector which points from the centre of the interface to the vertex we are interested in
                    t_dim = vtx.pos[0] - cell[i].iface[Face.west].pos;
                    t_unit = (t_dim) / geom.abs(t_dim);

                    // Find the post shock Mach number in the direction of this pointing vector for use in the weighting function
                    M = dot(u_east, t_unit)/cell[i].fs.gas.a;
                    
                    // Weighting functions defined on page 78
                    if (M <= 1.0) {
                        w[i] = ((M+1)*(M+1)+(M+1)*abs(M+1))/8.0;
                    }
                    else {
                        w[i] = M;
                    }
                }
                // All interface velocities and weightings are now found- we will do a quick check of the weights to prevent divide by 0
                // Reset wsum to 0
                number wsum = 0;
                for (i = 0; i < numcells; ++i) {
                    wsum += w[i];
                }
                if (wsum < 1e-7) {
                    for (i = 0; i < numcells; ++i) {
                        w[i] = 1.0;
                    }
                }
                denom = 0.0;
                numerator.clear();

                for (i = 0; i < numcells; ++i) {
                    numerator += w[i] * interface_ws[i];
                    denom += w[i];
                }
                // The theoretically correct vertex velocity- we still need to do a bit more for numerical stability               
                exact_vel = numerator / denom;
                // Catch if the velocity is too large- set to unshocked velocity
                if (geom.abs(exact_vel) > (geom.abs(inflow.vel) + inflow.gas.a)) {
                    exact_vel.set(inflow.vel.x + inflow.gas.a, -1.0*(inflow.vel.y+inflow.gas.a), (blk.myConfig.dimensions == 3) ? inflow.vel.z+inflow.gas.a : to!number(0.0));
                }
            }
            /++ Fill in the correct slot in the master array of vertex velocities. Do k-kmin to find how many rows have passed, multiply by jmax-jmin+1
              + (+1 as we want vertices instead of cells) for the number of cells in each row, then add the current j-jmin for the cell in the current row.
            ++/
            vertex_master_velocities[(blk.jmax - blk.jmin + 1) * (k - blk.kmin) + (j - blk.jmin)] = vtx_vel_stability_factor * exact_vel;
        
        }
    }
    return vertex_master_velocities;
}

// end shock_fitting_vertex_velocities
//--------------------------------------------------------------------------------------------------------------------------------------------
// Begin various helper functions called within shock_fitting_vertex_velocities
//--------------------------------------------------------------------------------------------------------------------------------------------

// Assigns velocities to all the inner vertices from the master array of velocities

void assign_slave_velocities(SFluidBlock blk, Vector3[] velocity_array) {
    size_t krangemax = ( blk.myConfig.dimensions == 3 ) ? blk.kmax+1 : blk.kmax;
    FVVertex vtx, vtx_left, vtx_right;
    Vector3 unit_d;
    get_ghost_vertex_positions(blk);
    for (size_t k = blk.kmin; k <= krangemax; ++k) {
        for (size_t j = blk.jmin; j <= blk.jmax+1; ++j) {
            for (size_t i = blk.imin; i <= blk.imax+1; ++i) {
                           
                vtx = blk.get_vtx(i, j, k);
                vtx_left = blk.get_vtx(i-1, j, k);
                vtx_right = blk.get_vtx(i+1, j, k);
                
                // Scale it by its normalised radial position
                vtx.vel[0] = vtx.radial_pos_norm * velocity_array[(blk.jmax - blk.jmin + 1) * (k - blk.kmin) + (j - blk.jmin)];

                // Find its radial direction and point it in that direction
                unit_d = correct_direction([vtx_left.pos[0], vtx.pos[0], vtx_right.pos[0]]);
                vtx.vel[0] = unit_d * dot(vtx.vel[0], unit_d);
            }
        }
    }
}

//-----------------------------------------------------------------------------------------------------

// Sets the direction for vertex movement
Vector3 correct_direction(Vector3[] pos) {
    // as Ian Johnston recommends in his thesis, we force the vertices to move along "rails" which are the radial lines
    // originating from the EAST boundary spanning to the west
    Vector3 unit_d;
    Vector3 left_pointer = pos[1] - pos[0];
    Vector3 right_pointer = pos[2] - pos[1];
    unit_d = 0.5 * (left_pointer / geom.abs(left_pointer) + right_pointer / geom.abs(right_pointer));
    return unit_d;
}
// ----------------------------------------------------------------------------------------------------
// Begin populate_ghost_cell_interface_geometry- used to get some extra ghost cell information required

void populate_ghost_cell_interface_geometry(SFluidBlock blk) {
    size_t krangemax = ( blk.myConfig.dimensions == 3 ) ? blk.kmax+1 : blk.kmax;

    if (blk.bc[Face.south].type == "exchange_over_full_face") {
        GhostCellFullFaceCopy ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.south].preReconAction[0];
        int neighbour = ffeBC.neighbourBlock.id;
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, neighbour).empty) {

                MPI_Request MPI_incoming_request;
                MPI_Status MPI_incoming_status;
                int mpi_recv_tag, mpi_send_tag, ne, neighbour_rank;
                number[] incoming_cells, outgoing_cells;
                // Neighbour block is non-local in the mpi context
                neighbour_rank = GlobalConfig.mpi_rank_for_block[neighbour];

                mpi_recv_tag = make_mpi_tag(blk.id, 1, 0);
                ne = to!int((krangemax - blk.kmin) + 1) * 6;
                incoming_cells.length = ne;
                outgoing_cells.length = ne;
                // Make the receive request to prepare for receiving info

                MPI_Irecv(incoming_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_recv_tag, MPI_COMM_WORLD, &MPI_incoming_request);

                // Construct the outgoing array
                size_t i = 0;
                for (size_t k = blk.kmin; k <= krangemax; k++) {
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmin, k).iface[Face.west].n.x;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmin, k).iface[Face.west].n.y;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmin, k).iface[Face.west].n.z;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmin, k).iface[Face.west].pos.x;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmin, k).iface[Face.west].pos.y;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmin, k).iface[Face.west].pos.z;
                }

                mpi_send_tag = make_mpi_tag(neighbour, 2, 0);
                MPI_Send(outgoing_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_send_tag, MPI_COMM_WORLD);
                if (blk.id == 1) {
                }
                MPI_Wait(&MPI_incoming_request, &MPI_incoming_status);
                i = 0;
                for (size_t k = blk.kmin; k <= krangemax; k++) {
                    blk.get_cell!()(blk.imin, blk.jmin-1, k).iface[Face.west].n.set(incoming_cells[i], incoming_cells[i+1], incoming_cells[i+2]);
                    blk.get_cell!()(blk.imin, blk.jmin-1, k).iface[Face.west].pos.set(incoming_cells[i+3], incoming_cells[i+4], incoming_cells[i+5]);
                    i += 6;
                }
            }

            else { // Neighbour block is local in the mpi context
                SFluidBlock neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
                for (size_t k = blk.kmin; k <= krangemax; k++) {
                blk.get_cell!()(blk.imin, blk.jmin-1, k) = neighbour_blk.get_cell!()(neighbour_blk.imin, neighbour_blk.jmax, k);
                }
            }
        }
        else {
            SFluidBlock neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
            for (size_t k = blk.kmin; k <= krangemax; k++) {
                blk.get_cell!()(blk.imin, blk.jmin-1, k) = neighbour_blk.get_cell!()(neighbour_blk.imin, neighbour_blk.jmax, k);
            }
        }
    }

    if (blk.bc[Face.north].type == "exchange_over_full_face") {
        GhostCellFullFaceCopy ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.north].preReconAction[0];
        int neighbour = ffeBC.neighbourBlock.id;
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, neighbour).empty) {
                MPI_Request MPI_incoming_request;
                MPI_Status MPI_incoming_status;
                int mpi_recv_tag, mpi_send_tag, ne, neighbour_rank;
                number[] incoming_cells, outgoing_cells;
                // Neighbour block is non-local in the mpi context
                neighbour_rank = GlobalConfig.mpi_rank_for_block[neighbour];

                mpi_recv_tag = make_mpi_tag(blk.id, 2, 0);
                ne = to!int((krangemax - blk.kmin) + 1) * 6;
                incoming_cells.length = ne;
                outgoing_cells.length = ne;
                // Make the receive request to prepare for receiving info
                MPI_Irecv(incoming_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_recv_tag, MPI_COMM_WORLD, &MPI_incoming_request);

                // Construct the outgoing array
                size_t i = 0;
                for (size_t k = blk.kmin; k <= krangemax; k++) {
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmax, k).iface[Face.west].n.x;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmax, k).iface[Face.west].n.y;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmax, k).iface[Face.west].n.z;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmax, k).iface[Face.west].pos.x;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmax, k).iface[Face.west].pos.y;
                    outgoing_cells[i++] = blk.get_cell!()(blk.imin, blk.jmax, k).iface[Face.west].pos.z;
                }

                mpi_send_tag = make_mpi_tag(neighbour, 1, 0);
                MPI_Send(outgoing_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_send_tag, MPI_COMM_WORLD);

                MPI_Wait(&MPI_incoming_request, &MPI_incoming_status);
                i = 0;
                for (size_t k = blk.kmin; k <= krangemax; k++) {
                    blk.get_cell!()(blk.imin, blk.jmax+1, k).iface[Face.west].n.set(incoming_cells[i], incoming_cells[i+1], incoming_cells[i+2]);
                    blk.get_cell!()(blk.imin, blk.jmax+1, k).iface[Face.west].pos.set(incoming_cells[i+3], incoming_cells[i+4], incoming_cells[i+5]);
                    i += 6;
                }
            }

            // Neighbour block is local in the mpi context
            else {
                SFluidBlock neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
                for (size_t k = blk.kmin; k <= krangemax; k++) {
                    blk.get_cell!()(blk.imin, blk.jmax+1, k) = neighbour_blk.get_cell!()(neighbour_blk.imin, neighbour_blk.jmin, k);
                }
            }
        }
        else {
            SFluidBlock neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
            for (size_t k = blk.kmin; k <= krangemax; k++) {
                blk.get_cell!()(blk.imin, blk.jmax+1, k) = neighbour_blk.get_cell!()(neighbour_blk.imin, neighbour_blk.jmin, k);
            }
        }
    }
}
// End populate_ghost_cell_geometry
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// Begin get_ghost_vertex_positions- used to get the ghost vertex locations for determine vertex direction

void get_ghost_vertex_positions(SFluidBlock blk) {
    
    size_t krange = (GlobalConfig.dimensions == 3) ? blk.kmax+1 : blk.kmax;
    // For each boundary in the inflow direction that is an exchange_over_full_face we need to do an exchange of information
    if (blk.bc[Face.west].type == "exchange_over_full_face") {
        GhostCellFullFaceCopy ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.west].preReconAction[0];
        int neighbour = ffeBC.neighbourBlock.id;

        // Check if the block is local in the mpi sense
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, neighbour).empty) {
                // We must exchange information using MPI
                MPI_Request MPI_incoming_request;
                MPI_Status MPI_incoming_status;
                int mpi_recv_tag, mpi_send_tag, ne, neighbour_rank;
                number[] incoming_cells, outgoing_cells;

                ne = to!int(3 * ((blk.jmax - blk.jmin + 2) * (krange - blk.kmin + 1)));
                neighbour_rank = GlobalConfig.mpi_rank_for_block[neighbour];

                mpi_recv_tag = make_mpi_tag(blk.id, 3, 0);
                
                incoming_cells.length = ne;
                outgoing_cells.length = ne;

                // Make a non-blocking receive to prepare to receive information
                MPI_Irecv(incoming_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_recv_tag, MPI_COMM_WORLD, &MPI_incoming_request);
            
                size_t i = 0;
                // Unpack the vertex positions
                for (size_t k = blk.kmin; k <= krange; k++) {
                    for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                        Vector3 vtx_pos = blk.get_vtx!()(blk.imin+1, j, k).pos[0];
                        outgoing_cells[i++] = vtx_pos.x;
                        outgoing_cells[i++] = vtx_pos.y;
                        outgoing_cells[i++] = vtx_pos.z;
                    }
                }

                // Send away the vertex position data
                mpi_send_tag = make_mpi_tag(neighbour, 4, 0);
                MPI_Send(outgoing_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_send_tag, MPI_COMM_WORLD);

                // Wait until the receive request is complete
                MPI_Wait(&MPI_incoming_request, &MPI_incoming_status);
                
                i = 0;
                for (size_t k = blk.kmin; k <= krange; k++) {
                    for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                        blk.get_vtx!()(blk.imin-1, j, k).pos[0].set(incoming_cells[i++], incoming_cells[i++], incoming_cells[i++]);
                    }
                }
            }

            else {  // Block is local in the mpi context- we can just grab the neighbours info
                auto neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
                for (size_t k = blk.kmin; k <= krange; k++) {
                    for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                        blk.get_vtx!()(blk.imin-1, j, k).pos[0] = neighbour_blk.get_vtx!()(neighbour_blk.imax, j, k).pos[0];
                    }
                }
            }
        }

        else {
            auto neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
            for (size_t k = krange; k <= krange; k++) {
                for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                    blk.get_vtx!()(blk.imin-1, j, k).pos[0] = neighbour_blk.get_vtx!()(neighbour_blk.imax, j, k).pos[0];
                }
            }
        }
    }
    else { // There's no real vertex here- we'll just extrapolate it out based on the last 2 vertex locations.
        for (size_t k = krange; k <= krange; k++) {
            for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                blk.get_vtx!()(blk.imin-1, j, k).pos[0] = blk.get_vtx!()(blk.imin, j, k).pos[0] + (blk.get_vtx!()(blk.imin, j, k).pos[0] - blk.get_vtx!()(blk.imin+1, j, k).pos[0]);
            }
        }
    }

    // Now we do the same process on the east boundary
    if (blk.bc[Face.east].type == "exchange_over_full_face") {
        GhostCellFullFaceCopy ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.east].preReconAction[0];
        int neighbour = ffeBC.neighbourBlock.id;

        // Check if the block is local in the mpi sense
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, neighbour).empty) {
                // We must exchange information using MPI
                MPI_Request MPI_incoming_request;
                MPI_Status MPI_incoming_status;
                int mpi_recv_tag, mpi_send_tag, ne, neighbour_rank;
                number[] incoming_cells, outgoing_cells;

                ne = to!int(3 * ((blk.jmax - blk.jmin + 2) * (krange - blk.kmin + 1)));
                neighbour_rank = GlobalConfig.mpi_rank_for_block[neighbour];

                mpi_recv_tag = make_mpi_tag(blk.id, 4, 0);
                
                incoming_cells.length = ne;
                outgoing_cells.length = ne;

                // Make a non-blocking receive to prepare to receive information
                MPI_Irecv(incoming_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_recv_tag, MPI_COMM_WORLD, &MPI_incoming_request);
            
                size_t i = 0;
                // Unpack the vertex positions
                for (size_t k = blk.kmin; k <= krange; k++) {
                    for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                        Vector3 vtx_pos = blk.get_vtx!()(blk.imax, j, k).pos[0];
                        outgoing_cells[i++] = vtx_pos.x;
                        outgoing_cells[i++] = vtx_pos.y;
                        outgoing_cells[i++] = vtx_pos.z;
                    }
                }

                // Send away the vertex position data
                mpi_send_tag = make_mpi_tag(neighbour, 3, 0);
                MPI_Send(outgoing_cells.ptr, ne, MPI_DOUBLE, neighbour_rank, mpi_send_tag, MPI_COMM_WORLD);

                // Wait until the receive request is complete
                MPI_Wait(&MPI_incoming_request, &MPI_incoming_status);
                
                i = 0;
                for (size_t k = blk.kmin; k <= krange; k++) {
                    for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                        blk.get_vtx!()(blk.imax+2, j, k).pos[0].set(incoming_cells[i++], incoming_cells[i++], incoming_cells[i++]);
                    }
                }
            }

            else {  // Block is local in the mpi context- we can just grab the neighbours info
                auto neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
                for (size_t k = blk.kmin; k <= krange; k++) {
                    for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                        blk.get_vtx!()(blk.imax+2, j, k).pos[0] = neighbour_blk.get_vtx!()(neighbour_blk.imin+1, j, k).pos[0];
                    }
                }
            }
        }

        else {
            auto neighbour_blk = cast(SFluidBlock) globalFluidBlocks[neighbour];
            for (size_t k = blk.kmin; k <= krange; k++) {
                for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                    blk.get_vtx!()(blk.imax+2, j, k).pos[0] = neighbour_blk.get_vtx!()(neighbour_blk.imin+1, j, k).pos[0];
                }
            }
        }
    }
    else { // There's no real vertex here- we'll just extrapolate it out based on the last 2 vertex locations
        for (size_t k = blk.kmin; k <= krange; k++) {
            for (size_t j = blk.jmin; j <= blk.jmax+1; j++) {
                blk.get_vtx!()(blk.imax+2, j, k).pos[0] = blk.get_vtx!()(blk.imax+1, j, k).pos[0] + (blk.get_vtx!()(blk.imax+1, j, k).pos[0] - blk.get_vtx!()(blk.imax, j, k).pos[0]);
            }
        }
    }
}

// End get_ghost_vertex_positions
//------------------------------------------------------------------------------------------------------------------

// Makes a (hopefully) unique integer for identifying MPI messages
@nogc
int make_mpi_tag(int blk_id, int bndry_id, int seq) {
    // unique integer tag for communication processes
    assert(seq < 10, "mpi-oops, too many differents messages");
    assert(bndry_id < 100, "mpi-oops, too many boundaries");
    assert(blk_id < 100000, "mpi_oops, too many blocks");
    return blk_id * 1000 + bndry_id * 10 + seq;
}

//-------------------------------------------------------------------------------------------------------------------

// Next two are actually referred to by simcore- need to unpack the velocities from Vector3s to an array of numbers
number[] unpack_vertex_velocities(Vector3[] vertex_velocities_as_vector) {
    size_t ne = vertex_velocities_as_vector.length * 3, i = 0;
    number[] unpacked_vertex_velocities;
    unpacked_vertex_velocities.length = ne;

    foreach(j, vector; vertex_velocities_as_vector) {
        unpacked_vertex_velocities[i++] = vector.x;
        unpacked_vertex_velocities[i++] = vector.y;
        unpacked_vertex_velocities[i++] = vector.z;
    }
    return unpacked_vertex_velocities;
}

// Pack the array of numbers describing the velocities back into Vector3s
Vector3[] pack_vertex_velocities(number[] vertex_velocities_as_elements) {
    size_t ne = vertex_velocities_as_elements.length / 3, i = 0;
    Vector3[] packed_vertex_velocities;
    packed_vertex_velocities.length = ne;
    // For some reason a foreach did not work as expected here
    for (size_t j = 0; j < ne; j++) {
        packed_vertex_velocities[j].set(vertex_velocities_as_elements[i++], vertex_velocities_as_elements[i++], vertex_velocities_as_elements[i++]);
    }
    return packed_vertex_velocities;
}

//--------------------------------------------------------------------------------------------------------------------

// Reconstruct the right quantities for shock detection and movement
number scalar_reconstruction(number x1, number x2, number x3, number h1, number h2, number h3, number g1) {
    bool johnston_reconstruction = true;
    number eps = 1.0e-12;
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

//-------------------------------------------------------------------------------------------------------------------
