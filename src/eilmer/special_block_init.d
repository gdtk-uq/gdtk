/** special_block_init.d
 * 
 * This module houses functions that do some special initialisation
 * of blocks at runtime.
 *
 * Authors: Rowan G. and Peter J.
 *
 */

module special_block_init;

import std.stdio;
import std.algorithm.searching;

import globalconfig;
import block;
import fvcell;
import fvinterface;

void diffuseWallConditionsIntoBlock(Block blk, int nsteps)
{
    size_t[] cellsAlongWalls;
    bool[size_t] cellsInDiffusionZone;
    size_t[] cellsAddedLastStep;
    size_t[] cellsAddedThisStep;

    // Determine which walls if any are no-slip walls
    bool[size_t] noSlipWalls;
    foreach (bcId, bc; blk.bc) {
	writeln("bcId: ", bcId, " is_wall= ", bc.is_wall);
	if (bc.is_wall) {
	    noSlipWalls[bcId] = true;
	}
    }
    
    // Determine which cells are against a no-slip wall.
    foreach (bcId; noSlipWalls.byKey()) {
	foreach (i, face; blk.bc[bcId].faces) {
	    if (blk.bc[bcId].outsigns[i] == 1) {
		cellsAlongWalls ~= face.left_cell.id;
	    }
	    else {
		cellsAlongWalls ~= face.right_cell.id;
	    }
	    //writeln("Just added wall cell: ", cellsAlongWalls[$-1]);
	}
    }

    // Cells along a wall get the properties at the wall.
    foreach (bcId; noSlipWalls.byKey()) {
	blk.bc[bcId].applyPreSpatialDerivAction(0.0, 0, 0);
    }
    foreach (cellId; cellsAlongWalls) {
	auto cell = blk.get_cell_from_array(cellId);
	foreach (face; cell.iface) {
	    if (face.is_on_boundary && (face.bc_id in noSlipWalls)) {
		cell.fs.gas.Ttr = face.fs.gas.Ttr;
		cell.fs.vel.set(face.fs.vel);
		cell.fs.tke = face.fs.tke;
		cell.fs.omega = face.fs.omega;
		cell.fs.mu_t = face.fs.mu_t;
		cell.fs.k_t = face.fs.k_t;
		blk.myConfig.gmodel.update_thermo_from_pT(cell.fs.gas);
		cell.encode_conserved(0, 0, 0.0);
		cell.decode_conserved(0, 0, 0.0);
	    } 
	}
    }
    
    // Gather the cells adjacent to the wall cells
    // as the first layer in the diffusion zone.
    foreach (cellId; cellsAlongWalls) {
	auto cell = blk.get_cell_from_array(cellId);
	foreach (face; cell.iface) {
	    if (face.is_on_boundary) {
		// The only cells this interface touches
		// are ghost cells or wall cells.
		// We do not want these in the diffusion zone.
		continue;
	    }
	    if (!cellsAlongWalls.canFind(face.left_cell.id)) {
		//writeln("Adding cell to diffusion zone: ", face.left_cell.id);
		cellsInDiffusionZone[face.left_cell.id] = true;
		if (!cellsAddedLastStep.canFind(face.left_cell.id)) {
		    cellsAddedLastStep ~= face.left_cell.id;
		}
	    }
	    if (!cellsAlongWalls.canFind(face.right_cell.id)) {
		//writeln("Adding cell to diffusion zone: ", face.right_cell.id);
		cellsInDiffusionZone[face.right_cell.id] = true;
		if (!cellsAddedLastStep.canFind(face.right_cell.id)) {
		    cellsAddedLastStep ~= face.right_cell.id;
		}
	    }
	}
    }

    // Now perform averaging of flow properties for cells
    // in the diffusion zone.
    foreach (step; 0 .. nsteps) {
	// Do the averaging
	foreach (cellId; cellsInDiffusionZone.byKey()) {
	    auto cell = blk.get_cell_from_array(cellId);
	    int nNbrCells = 0;
	    double T_avg = 0.0;
	    double velx_avg = 0.0;
	    double vely_avg = 0.0;
	    double velz_avg = 0.0;
	    double tke_avg = 0.0;
	    double omega_avg = 0.0;
	    double mu_t_avg = 0.0;
	    double k_t_avg = 0.0;
	    foreach (face; cell.iface) {
		if (face.is_on_boundary) continue;
		if (face.left_cell.id == cell.id) {
		    // Then right cell must be a neighbour
		    ++nNbrCells;
		    T_avg += face.right_cell.fs.gas.Ttr;
		    velx_avg += face.right_cell.fs.vel.x;
		    vely_avg += face.right_cell.fs.vel.y;
		    velz_avg += face.right_cell.fs.vel.z;
		    tke_avg += face.right_cell.fs.tke;
		    omega_avg += face.right_cell.fs.omega;
		    mu_t_avg += face.right_cell.fs.mu_t;
		    k_t_avg += face.right_cell.fs.k_t;
		}
		else {
		    // The left cell must be a neighbour;
		    ++nNbrCells;
		    T_avg += face.left_cell.fs.gas.Ttr;
		    velx_avg += face.left_cell.fs.vel.x;
		    vely_avg += face.left_cell.fs.vel.y;
		    velz_avg += face.left_cell.fs.vel.z;
		    tke_avg += face.left_cell.fs.tke;
		    omega_avg += face.left_cell.fs.omega;
		    mu_t_avg += face.left_cell.fs.mu_t;
		    k_t_avg += face.left_cell.fs.k_t;
		}
	    }
	    // Place the averaged value in cell.
	    cell.fs.gas.Ttr = T_avg / nNbrCells;
	    cell.fs.vel.refx = velx_avg / nNbrCells;
	    cell.fs.vel.refy = vely_avg / nNbrCells;
	    cell.fs.vel.refz = velz_avg / nNbrCells;
	    cell.fs.tke = tke_avg / nNbrCells;
	    cell.fs.omega = omega_avg / nNbrCells;
	    cell.fs.mu_t = mu_t_avg / nNbrCells;
	    cell.fs.k_t = k_t_avg / nNbrCells;
	    blk.myConfig.gmodel.update_thermo_from_pT(cell.fs.gas);
	    cell.encode_conserved(0, 0, 0.0);
	    cell.decode_conserved(0, 0, 0.0);
	}

	// Add a new layer of cells for diffusion step
	if (step != (nsteps-1)) {
	    cellsAddedThisStep.length = 0;
	    foreach (cellId; cellsAddedLastStep) {
		auto cell = blk.get_cell_from_array(cellId);
		foreach (face; cell.iface) {
		    if (face.is_on_boundary) {
			// The only cells this interface touches
			// are ghost cells or wall cells.
			// We do not want these in the diffusion zone.
			continue;
		    }
		    if (face.left_cell.id == cell.id) {
			// Interested in right_cell.
			if ( !(face.right_cell.id in cellsInDiffusionZone) &&
			     !cellsAlongWalls.canFind(face.right_cell.id) ) {
			    //writeln("Adding cell to diffusion zone: ", face.right_cell.id);
			    cellsInDiffusionZone[face.right_cell.id] = true;
			    cellsAddedThisStep ~= face.right_cell.id;
			}
		    }
		    else {
			// Interested in left_cell.
			if ( !(face.left_cell.id in cellsInDiffusionZone) &&
			     !cellsAlongWalls.canFind(face.left_cell.id) ) {
			    //writeln("Adding cell to diffusion zone: ", face.left_cell.id);
			    cellsInDiffusionZone[face.left_cell.id] = true;
			    cellsAddedThisStep ~= face.left_cell.id;
			}
		    }
		}
	    }
	    cellsAddedLastStep = cellsAddedThisStep.dup;
	}
    }
}
