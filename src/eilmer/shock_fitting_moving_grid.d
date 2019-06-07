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

/++ Still to do- parallelise the process of running through the inner blocks and assigning the internal vertex velocities; this includes the task of assigning each shock fitted block a list of 'partners'
which are tied to it. This assigning should be done here in assign_radial_dist, but thinking needs to done regarding what to attach this list to; perhaps just the shock fitted block as an attribute?

Other to do is find a place to implement this assign_radial_dist function, essentially, where do the blocks become SFluidBlocks in the initialisation, should I ask for some other type of input which fits better for initialisation? ++/

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

void shock_fitting_vertex_velocities(SFluidBlock blk, int step, double sim_time) {
    /++ for a given block, loop through cell vertices and update the vertex
      + velocities. The boundary vertex velocities are set via the methodology laid out
      + in Ian Johnston's thesis available on the cfcfd website. The internal velocities
      + are then assigned based on the boundary velocity and a user chosen weighting.
      + NB: Implementation is hard coded for a moving WEST boundary
      ++/

    // Lets be verbose about where each of these are used
    FVVertex vtx;        // The vertex we are calculating velocity at
    FVInterface[] iface_neighbour;      // Used in calculating interface speed weights
    FVCell[] cell, cell_R1, cell_R2;
    SFluidBlock new_blk;
    Vector3 exact_vel, directional_vel, unit_d, u_west, u_east, ns, t_dim, t_unit, numerator; 
    number rho, shock_detect, sign_ws2, sqrt_ws2, ws1, ws2, rho_west, rho_east, p_west, p_east, M, rho_northbot, rho_southbot, rho_northtop, rho_southtop;
    Vector3[] interface_ws;
    number[] w;
    int array_size;
    number wsum, denom;
    Vector3[] vertex_master_velocities;
	bool isSouth, isNorth, isTop, isBot;
    int interpolation_order = blk.myConfig.shock_fitting_interpolation_order;
    immutable double SHOCK_DETECT_THRESHOLD =  0.2;
    immutable double vtx_vel_stability_factor = blk.myConfig.shock_fitting_scale_factor;
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
    if (sim_time < GlobalConfig.shock_fitting_delay || blk.bc[Face.west].type != "inflow_shock_fitting") return;
    // inflow is a reference to a duplicated flow state, it should be treated as read-only (do not alter)
    auto constFluxEffect = cast(BFE_ConstFlux) blk.bc[Face.west].postConvFluxAction[0];
    auto inflow = constFluxEffect.fstate;
    // #####-------------------------- SHOCK SEARCH --------------------------##### //
    for ( size_t k = blk.kmin; k <= krangemax; ++k ) {
        for ( size_t j = blk.jmin; j <= blk.jmax+1; ++j ) {
            size_t i = blk.imin;
            vtx = blk.get_vtx!()(i,j,k);
            // West of the vertex is the constant flux condition
            rho_west = inflow.gas.rho;       // density
            u_west = inflow.vel;             // velocity vector
            p_west = inflow.gas.p;           // pressure
            rho_east = 0.0;       
            // Clear the arrays used to store reconstruction cells.
            cell = [];
            iface_neighbour = [];
            cell_R1 = [];
            cell_R2 = [];
            interface_ws = [];
            w = [];

            // loop over neighbouring interfaces (above vtx first, then below vtx)
            foreach (int jOffSet; [0, 1]) {
                foreach (kOffSet; (blk.myConfig.dimensions == 3) ? [0, 1] : [0]) {
                    if ((j-jOffSet) >= blk.jmin && (j-jOffSet) <= blk.jmax && (k-kOffSet) >= blk.kmin && (k-kOffSet) <= krangemax) {
                        // This is the regular case where the cell is on the normal block i.e. not a ghost cell.
                        // Add another index to the cell arrays
                        cell.length += 1;
                        iface_neighbour.length += 1;
                        cell_R1.length += 1;
                        cell_R2.length += 1;
                        w.length += 1;
                        interface_ws.length += 1;

                        // Call the current cell
                        cell[$-1] = blk.get_cell!()(i, j-jOffSet, k-kOffSet);
                        iface_neighbour[$-1] = blk.get_ifi!()(i,j-jOffSet, k-kOffSet);
                        cell_R1[$-1] = blk.get_cell!()(i+1, j-jOffSet, k-kOffSet);
                        cell_R2[$-1] = blk.get_cell!()(i+2, j-jOffSet, k-kOffSet);
                        /++ Unfortunately we need to pull data from other blocks here rather than just relying on ghost cells,
                          + As ghost cells do not contain any information about their orientation which is required for shock
                          + speed estimation in the form of the normal vector ns. Can't think of any way to do this elegantly;
                          + will just manually go through every situation. Corners in 3D... Yikes.
                        ++/
                    }
                        /++ The order in which we make these border checks is relevant in 3D. We need to check for blocks touching
                          + at the corner first, as half of its logic conditions for existence will be true for any edge cell, so
                          + with the usage of else if statements a corner cell will just run as a regular edge cell case. E.g. If
                          + we have a cell off the bottom-south corner, the function will see "our j/k index is less than jmin/kmin
                          + for the block, so run the south/bottom edge case then continue to interface speed calculation" if we
                          + order our checks incorrectly.
                        ++/
                    else if ((j-jOffSet) < blk.jmin && (k-kOffSet) < blk.kmin && blk.bc[Face.south].type == "exchange_over_full_face"
                                && blk.bc[Face.bottom].type == "exchange_over_full_face") {
                        // We're in the bottom-south corner block, but we can't increment the size of our cell array yet; the south
                        // and bottom blocks may border eachother, rather than there being a 4th block there.
                        SFluidBlock south_blk = call_neighbour_block("south", blk);
                        SFluidBlock bottom_blk = call_neighbour_block("bottom", blk);

                        if (south_blk.bc[Face.bottom].type == "exchange_over_full_face" && bottom_blk.bc[Face.south].type == "exchange_over_full_face") {
                            SFluidBlock south_bottom_blk = call_neighbour_block("south", bottom_blk);
                            SFluidBlock bottom_south_blk = call_neighbour_block("south", south_blk);
                            if (south_bottom_blk.id == bottom_south_blk.id) {
                                // There is a corner block here
                                cell.length += 1;
                                iface_neighbour.length += 1;
                                cell_R1.length += 1;
                                cell_R2.length += 1;
                                w.length += 1;
                                interface_ws.length += 1;
                                    
                                cell[$-1] = south_bottom_blk.get_cell!()(i, south_bottom_blk.jmax, south_bottom_blk.kmax);
                                iface_neighbour[$-1] = south_bottom_blk.get_ifi!()(i, south_bottom_blk.jmax, south_bottom_blk.kmax);
                                cell_R1[$-1] = south_bottom_blk.get_cell!()(i+1, south_bottom_blk.jmax, south_bottom_blk.kmax);
                                cell_R2[$-1] = south_bottom_blk.get_cell!()(i+2, south_bottom_blk.jmax, south_bottom_blk.kmax);
                            }
                        }
                    }    
                    else if ((j-jOffSet) > blk.jmax && (k-kOffSet) < blk.kmin && blk.bc[Face.north].type == "exchange_over_full_face"
                                && blk.bc[Face.bottom].type == "exchange_over_full_face") {
                        // We're in the bottom-north corner block, but we can't increment the size of our cell array yet; the north
                        // and bottom blocks may border eachother, rather than there being a 4th block there.
                        SFluidBlock north_blk = call_neighbour_block("north", blk);
                        SFluidBlock bottom_blk = call_neighbour_block("bottom", blk);

                        if (north_blk.bc[Face.bottom].type == "exchange_over_full_face" && bottom_blk.bc[Face.north].type == "exchange_over_full_face") {
                            SFluidBlock north_bottom_blk = call_neighbour_block("north", bottom_blk);
                            SFluidBlock bottom_north_blk = call_neighbour_block("north", north_blk);
                            if (north_bottom_blk.id == bottom_north_blk.id) {
                                // There is a corner block here
                                cell.length += 1;
                                iface_neighbour.length += 1;
                                cell_R1.length += 1;
                                cell_R2.length += 1;
                                w.length += 1;
                                interface_ws.length += 1;
                                    
                                cell[$-1] = north_bottom_blk.get_cell!()(i, north_bottom_blk.jmin, north_bottom_blk.kmax);
                                iface_neighbour[$-1] = north_bottom_blk.get_ifi!()(i, north_bottom_blk.jmax, north_bottom_blk.kmax);
                                cell_R1[$-1] = north_bottom_blk.get_cell!()(i+1, north_bottom_blk.jmin, north_bottom_blk.kmax);
                                cell_R2[$-1] = north_bottom_blk.get_cell!()(i+2, north_bottom_blk.jmin, north_bottom_blk.kmax);
                            }
                        }
                    }    
                    else if ((j-jOffSet) < blk.jmin && (k-kOffSet) > krangemax && blk.bc[Face.south].type == "exchange_over_full_face"
                                && blk.bc[Face.top].type == "exchange_over_full_face") {
                        // We're in the top-south corner block, but we can't increment the size of our cell array yet; the south
                        // and top blocks may border eachother, rather than there being a 4th block there.
                        SFluidBlock south_blk = call_neighbour_block("south", blk);
                        SFluidBlock top_blk = call_neighbour_block("top", blk);

                        if (south_blk.bc[Face.top].type == "exchange_over_full_face" && top_blk.bc[Face.south].type == "exchange_over_full_face") {
                            SFluidBlock south_top_blk = call_neighbour_block("south", top_blk);
                            SFluidBlock top_south_blk = call_neighbour_block("south", south_blk);
                            if (south_top_blk.id == top_south_blk.id) {
                                // There is a corner block here
                                cell.length += 1;
                                iface_neighbour.length += 1;
                                cell_R1.length += 1;
                                cell_R2.length += 1;
                                w.length += 1;
                                interface_ws.length += 1;
                                    
                                cell[$-1] = south_top_blk.get_cell!()(i, south_top_blk.jmax, south_top_blk.kmin);
                                iface_neighbour[$-1] = south_top_blk.get_ifi!()(i, south_top_blk.jmax, south_top_blk.kmin);
                                cell_R1[$-1] = south_top_blk.get_cell!()(i+1, south_top_blk.jmax, south_top_blk.kmin);
                                cell_R2[$-1] = south_top_blk.get_cell!()(i+2, south_top_blk.jmax, south_top_blk.kmin);
                            }
                        }
                    }    
                    else if ((j-jOffSet) > blk.jmax && (k-kOffSet) > krangemax && blk.bc[Face.north].type == "exchange_over_full_face"
                                && blk.bc[Face.top].type == "exchange_over_full_face") {
                        // We're in the top-north corner block, but we can't increment the size of our cell array yet; the north
                        // and top blocks may border eachother, rather than there being a 4th block there.
                        SFluidBlock north_blk = call_neighbour_block("north", blk);
                        SFluidBlock top_blk = call_neighbour_block("top", blk);

                        if (north_blk.bc[Face.top].type == "exchange_over_full_face" && top_blk.bc[Face.north].type == "exchange_over_full_face") {
                            SFluidBlock north_top_blk = call_neighbour_block("north", top_blk);
                            SFluidBlock top_north_blk = call_neighbour_block("north", north_blk);
                            if (north_top_blk.id == top_north_blk.id) {
                                // There is a corner block here
                                cell.length += 1;
                                iface_neighbour.length += 1;
                                cell_R1.length += 1;
                                cell_R2.length += 1;
                                w.length += 1;
                                interface_ws.length += 1;
                                    
                                cell[$-1] = north_top_blk.get_cell!()(i, north_top_blk.jmin, north_top_blk.kmin);
                                iface_neighbour[$-1] = north_top_blk.get_ifi!()(i, north_top_blk.jmin, north_top_blk.kmin);
                                cell_R1[$-1] = north_top_blk.get_cell!()(i+1, north_top_blk.jmin, north_top_blk.kmin);
                                cell_R2[$-1] = north_top_blk.get_cell!()(i+2, north_top_blk.jmin, north_top_blk.kmin);
                            }
                        }
                    }
                    /++ This should be all the general corner cases done. There are fringe cases where this will fail; the only one
                      + I can think of currently is the case where there is a 'hanging' block e.g. you have a block to the south of
                      + your bottom block, but no south block (see the pretty picture below). While this may happen relatively often
                writeln("shock");
                      + around geometries, I can't see it happening at an inflow boundary. If it does and it's important for your
                      + problem, talk to Lachlan Whyborn on level 5 Mansergh Shaw.
                    -------------------------------------------------
                    |                       |                       |
                    |                       |                       |
                    |                       |                       |
                    |                       |                       |
                    |     North Block       |     Current Block     |
                    |                       |                       |
                    |                       |                       |
                    |                       |                       |
                    |                       |                       |
                    -------------------------------------------------
                    |                       |
                    |                       |   Imagine the flow going into the page in this diagram.
                    |                       |
                    |                       |
                    |  North-Bottom Block   |
                    |                       |
                    |                       |
                    |                       |
                    |                       |
                    -------------------------    
                    ++/

                    // Now we consider the much simpler edge-sharing blocks
                    else if ((j-jOffSet) < blk.jmin && blk.bc[Face.south].type == "exchange_over_full_face") {
                        // We're in the block to the south of our current block
                        // Add another index to the cell arrays
                        cell.length += 1;
                        iface_neighbour.length += 1;
                        cell_R1.length += 1;
                        cell_R2.length += 1;
                        w.length += 1;
                        interface_ws.length += 1;

                        new_blk = call_neighbour_block("south", blk);
                        cell[$-1] = new_blk.get_cell!()(i, new_blk.jmax, k-kOffSet);
                        iface_neighbour[$-1] = new_blk.get_ifi!()(i, new_blk.jmax, k-kOffSet);
                        cell_R1[$-1] = new_blk.get_cell!()(i+1, new_blk.jmax, k-kOffSet);
                        cell_R2[$-1] = new_blk.get_cell!()(i+2, new_blk.jmax, k-kOffSet);
                    }
                    else if ((j-jOffSet) > blk.jmax && blk.bc[Face.north].type == "exchange_over_full_face" ) {
                        // We're in the block to the north of our current block
                        // Add another index to the cell arrays
                        cell.length += 1;
                        iface_neighbour.length += 1;
                        cell_R1.length += 1;
                        cell_R2.length += 1;
                        w.length += 1;
                        interface_ws.length += 1;

                        new_blk = call_neighbour_block("north", blk);
                        cell[$-1] = new_blk.get_cell!()(i, new_blk.jmin, k-kOffSet);
                        iface_neighbour[$-1] = new_blk.get_ifi!()(i, new_blk.jmin, k-kOffSet);
                        cell_R1[$-1] = new_blk.get_cell!()(i+1, new_blk.jmin, k-kOffSet);
                        cell_R2[$-1] = new_blk.get_cell!()(i+2, new_blk.jmin, k-kOffSet);
                    }
                    else if ((k-kOffSet) < blk.kmin && blk.bc[Face.bottom].type == "exchange_over_full_face" ) {
                        // We're in the block to the north of our current block
                        // Add another index to the cell arrays
                        cell.length += 1;
                        iface_neighbour.length += 1;
                        cell_R1.length += 1;
                        cell_R2.length += 1;
                        w.length += 1;
                        interface_ws.length += 1;

                        new_blk = call_neighbour_block("bottom", blk);
                        cell[$-1] = new_blk.get_cell!()(i, j-jOffSet, new_blk.kmax);
                        iface_neighbour[$-1] = new_blk.get_ifi!()(i, j-jOffSet, new_blk.kmax);
                        cell_R1[$-1] = new_blk.get_cell!()(i+1, j-jOffSet, new_blk.kmax);
                        cell_R2[$-1] = new_blk.get_cell!()(i+2, j-jOffSet, new_blk.kmax);
                    }
                    
                    else if ((k-kOffSet) > krangemax && blk.bc[Face.top].type == "exchange_over_full_face" ) {
                        // We're in the block to the north of our current block
                        // Add another index to the cell arrays
                        cell.length += 1;
                        iface_neighbour.length += 1;
                        cell_R1.length += 1;
                        cell_R2.length += 1;
                        w.length += 1;
                        interface_ws.length += 1;

                        new_blk = call_neighbour_block("top", blk);
                        cell[$-1] = new_blk.get_cell!()(i, j-jOffSet, new_blk.kmin);
                        iface_neighbour[$-1] = new_blk.get_ifi!()(i, j-jOffSet, new_blk.kmin);
                        cell_R1[$-1] = new_blk.get_cell!()(i+1, j-jOffSet, new_blk.kmin);
                        cell_R2[$-1] = new_blk.get_cell!()(i+2, j-jOffSet, new_blk.kmin);
                    }
                    else {}
                        // This is just here for good practice- if it hits this, it means there is no cell for the 
                        // given j-jOffSet, k-kOffset and nothing should be done.
                }
            }
            for(i = 0; i < cell.length; ++i) {
                rho_east += scalar_reconstruction(cell[i].fs.gas.rho,  cell_R1[i].fs.gas.rho,
                                                cell_R2[i].fs.gas.rho, cell[i].iLength,
                                                cell_R1[i].iLength,  cell_R2[i].iLength, inflow.gas.rho);
            }
            shock_detect = abs(inflow.gas.rho - (rho_east/cell.length))/fmax(inflow.gas.rho, (rho_east/cell.length));
            if (shock_detect < SHOCK_DETECT_THRESHOLD) { 
                // We don't think there's a shock here, so set vertex velocity to the maximum, U+a.
                exact_vel.set(inflow.vel.x + inflow.gas.a, -1.0*(inflow.vel.y+inflow.gas.a), (blk.myConfig.dimensions == 3) ? inflow.vel.z+inflow.gas.a: to!number(0.0));
            }
            else {
                // We do have a shock here
                for (i = 0; i < cell.length; ++i) {
                    // Reconstruct the density and pressure values at the interface
                    cell[i].fs.vel.transform_to_local_frame(cell[i].iface[Face.west].n, cell[i].iface[Face.west].t1, cell[i].iface[Face.west].t2);
                    cell_R1[i].fs.vel.transform_to_local_frame(cell_R1[i].iface[Face.west].n, cell_R1[i].iface[Face.west].t1, cell_R1[i].iface[Face.west].t2);
                    cell_R2[i].fs.vel.transform_to_local_frame(cell_R2[i].iface[Face.west].n, cell_R2[i].iface[Face.west].t1, cell_R2[i].iface[Face.west].t2);
                    
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
                    
                    u_east.transform_to_global_frame(cell[i].iface[Face.west].n, cell[i].iface[Face.west].t1, cell[i].iface[Face.west].t2);
                    cell[i].fs.vel.transform_to_global_frame(cell[i].iface[Face.west].n, cell[i].iface[Face.west].t1, cell[i].iface[Face.west].t2);
                    cell_R1[i].fs.vel.transform_to_global_frame(cell_R1[i].iface[Face.west].n, cell_R1[i].iface[Face.west].t1, cell_R1[i].iface[Face.west].t2);
                    cell_R2[i].fs.vel.transform_to_global_frame(cell_R2[i].iface[Face.west].n, cell_R2[i].iface[Face.west].t1, cell_R2[i].iface[Face.west].t2);

                    //Find the shock speeds ws1 and ws2 as outlined on page 73 of Ian Johnston's thesis
                    ns =  cell[i].iface[Face.west].n;
                    ws1 = (rho_west * dot(u_west, ns) - rho_east * dot(u_east, ns))/(rho_west - rho_east);
                    sign_ws2 = sgn(p_east.re - p_west.re)/rho_west;
                    sqrt_ws2 =  sqrt(abs((p_east - p_west)/(1/rho_west - 1/rho_east)));
                    ws2 = dot(u_west, ns) - sign_ws2 * sqrt_ws2;
                    double mixing_param = 0.5;
                    interface_ws[i] = (mixing_param * ws1 + (1-mixing_param) * ws2) * ns;
                    // Create a unit vector which points from the centre of the interface to the vertex we are interested in
                    t_dim = vtx.pos[0] - iface_neighbour[i].pos;
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
                wsum = 0;
                for (i = 0; i < cell.length; ++i) {
                    wsum += w[i];
                }
                if (wsum < 1e-7) {
                    for (i = 0; i < cell.length; ++i) {
                        w[i] = 1.0;
                    }
                }
                denom = 0.0;
                numerator.clear();

                for (i = 0; i < cell.length; ++i) {
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
    // Now we assign all the slave (interior) vertex velocities
    foreach (indx; blk.inflow_partners) {
        SFluidBlock slave_block = cast(SFluidBlock) globalFluidBlocks[indx];
        assign_slave_velocities(slave_block, vertex_master_velocities);
    }
}

void assign_slave_velocities(SFluidBlock blk, Vector3[] velocity_array) {
    size_t krangemax = ( blk.myConfig.dimensions == 3 ) ? blk.kmax+1 : blk.kmax;
    FVVertex vtx, vtx_left, vtx_right;
    Vector3 unit_d;
    for (size_t k = blk.kmin; k <= krangemax; ++k) {
        for (size_t j = blk.jmin; j <= blk.jmax+1; ++j) {
            for (size_t i = blk.imin; i <= blk.imax+1; ++i) {
                vtx = blk.get_vtx(i, j, k);
                vtx_left = blk.get_vtx(i-1, j, k);
                vtx_right = blk.get_vtx(i+1, j, k);
                vtx.vel[0] = vtx.radial_pos_norm * velocity_array[(blk.jmax - blk.jmin + 1) * (k - blk.kmin) + (j - blk.jmin)];
                if (i == blk.imin) {
                    unit_d = correct_direction([vtx.pos[0], vtx_right.pos[0]]);
                }
                else if (i == blk.imax+1) {
                    unit_d = correct_direction([vtx_left.pos[0], vtx.pos[0]]);
                }
                else {
                    unit_d = correct_direction([vtx_left.pos[0], vtx.pos[0], vtx_right.pos[0]]);
                }
            vtx.vel[0] = unit_d * dot(vtx.vel[0], unit_d);
            }
        }
    }
}
    
Vector3 correct_direction(Vector3[] pos) {
    // as Ian Johnston recommends in his thesis, we force the vertices to move along "rails" which are the radial lines
    // originating from the EAST boundary spanning to the west
    Vector3 unit_d;
    if (pos.length == 3) {
        Vector3 left_pointer = pos[1] - pos[0];
        Vector3 right_pointer = pos[2] - pos[1];
        // We've been given 3 vertex positions- assume its in order of inflow position
        unit_d = 0.5 * (left_pointer / geom.abs(left_pointer) + right_pointer / geom.abs(right_pointer));
    }
    else {
        // We've only been given 2 vertex positions- must be on a boundary. Use single pointing vector
        Vector3 pointer = pos[1] - pos[0];
        unit_d = pointer / geom.abs(pointer);
    }
    return unit_d;
}

SFluidBlock call_neighbour_block(string edge, SFluidBlock blk) {
    // Can't find a clean way to do this
    // Give ffeBC a dummy value
    GhostCellFullFaceCopy ffeBC;
    if (edge == "south") ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.south].preReconAction[0];
    else if (edge == "north") ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.north].preReconAction[0];
    else if (edge == "bottom") ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.bottom].preReconAction[0];
    else if (edge == "top") ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.top].preReconAction[0];
    else if (edge == "east") ffeBC = cast(GhostCellFullFaceCopy) blk.bc[Face.east].preReconAction[0];
        
    int neighbourBlock = ffeBC.neighbourBlock.id;
    auto new_blk = cast(SFluidBlock) globalFluidBlocks[neighbourBlock];

    return new_blk;
}

Vector3 weighting_function(Vector3 vel_max, size_t imax, size_t i) {
    // TODO: user chooses weighting function. Currently just linear.
    Vector3 vel = -(vel_max/(imax-2)) * (i-2) + vel_max; // y = mx + c
    return vel;
}

number scalar_reconstruction(number x1, number x2, number x3, number h1, number h2, number h3, number g1) {
    bool johnston_reconstruction = true;
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
