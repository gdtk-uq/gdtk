/** shape_sensitivity_core.d
 * 
 * Eilmer4 shape sensitivity calculator core functions.
 *
 * Author: Kyle D.
**/

module shape_sensitivity_core;

import core.memory;
import std.stdio;
import std.math;
import std.process;
import std.algorithm;
import std.string;
import std.file;
import std.parallelism;
import std.conv;

import util.lua;
import util.lua_service;

import nm.bbla;
import nm.smla;
import nm.luabbla;
import nm.complex;
import nm.number;

import special_block_init;
import steadystate_core;
import fluidblock;
import sfluidblock;
import ufluidblock;
import fvcell;
import fvinterface;
import fvvertex;
import globaldata;
import globalconfig;
import bc;
import onedinterp;
import lsqinterp;
import grid_deform;
import fvcore;
import fileutil;
import geom;
import geom.luawrap;
import lua_helper;
import simcore;
import fluxcalc;
import user_defined_source_terms;
import flowgradients;

shared enum ghost_cell_start_id = 1_000_000_000;
shared immutable double ESSENTIALLY_ZERO = 1.0e-50;
// Module-local, global memory arrays and matrices for GMRES
number[] g0;
number[] g1;
number[] h;
number[] hR;
Matrix!number H0;
Matrix!number H1;
Matrix!number Gamma;
Matrix!number Q0;
Matrix!number Q1;

private lua_State* L; // module-local Lua interpreter

/**********************/
/* Frechet Derivative */
/**********************/
void evalPrimitiveJacobianVecProd(FluidBlock blk, size_t nPrimitive, number[] v, ref number[] p, number EPS) {

    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;
    
    with_k_omega = (blk.myConfig.turbulence_model == TurbulenceModel.k_omega);
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    int cellCount = 0;
    foreach (cell; blk.cells) {
        cell.fs.gas.rho += EPS*v[cellCount+blk.MASS];
        cell.fs.vel.refx += EPS*v[cellCount+blk.X_MOM];
        cell.fs.vel.refy += EPS*v[cellCount+blk.Y_MOM];
        cell.fs.gas.p += EPS*v[cellCount+blk.TOT_ENERGY];
        if ( blk.myConfig.dimensions == 3 )
            cell.fs.vel.refz += EPS*v[cellCount+blk.Z_MOM];
        if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
            cell.fs.tke += EPS*v[cellCount+blk.TKE];
            cell.fs.omega += EPS*v[cellCount+blk.OMEGA];
        }
                
        blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
        blk.myConfig.gmodel.update_trans_coeffs(cell.fs.gas);
        blk.myConfig.gmodel.update_sound_speed(cell.fs.gas);
        //foreach(isp; 0 .. blk.myConfig.gmodel.n_species) cell.fs.gas.massf[isp] = (cell.U[0].massf[isp] * (1.0/cell.fs.gas.rho));
        cellCount += nPrimitive;
    }    
    steadystate_core.evalRHS(0.0, 0);
    cellCount = 0;
    foreach (cell; blk.cells) {
        p[cellCount+blk.MASS] = -cell.dUdt[0].mass.im/EPS.im; 
        p[cellCount+blk.X_MOM] = -cell.dUdt[0].momentum.x.im/EPS.im;
        p[cellCount+blk.Y_MOM] = -cell.dUdt[0].momentum.y.im/EPS.im;
        if ( blk.myConfig.dimensions == 3 )
            p[cellCount+blk.Z_MOM] = -cell.dUdt[0].momentum.z.im/EPS.im;
        p[cellCount+blk.TOT_ENERGY] = -cell.dUdt[0].total_energy.im/EPS.im;
        if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
            p[cellCount+blk.TKE] = -cell.dUdt[0].tke.im/EPS.im;
            p[cellCount+blk.OMEGA] = -cell.dUdt[0].omega.im/EPS.im;
        }
        cellCount += nPrimitive;
    }
}

void evalConservativeJacobianVecProd(FluidBlock blk, size_t nConserved, number[] v, ref number[] p, number EPS) {

    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;

    // We perform a Frechet derivative to evaluate J*D^(-1)v
    with_k_omega = (blk.myConfig.turbulence_model == TurbulenceModel.k_omega);
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    int cellCount = 0;
    foreach (cell; blk.cells) {
        cell.U[1].copy_values_from(cell.U[0]);
        cell.U[1].mass += EPS*v[cellCount+blk.MASS];
        cell.U[1].momentum.refx += EPS*v[cellCount+blk.X_MOM];
        cell.U[1].momentum.refy += EPS*v[cellCount+blk.Y_MOM];
        if ( blk.myConfig.dimensions == 3 )
            cell.U[1].momentum.refz += EPS*v[cellCount+blk.Z_MOM];
        cell.U[1].total_energy += EPS*v[cellCount+blk.TOT_ENERGY];
        if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
            cell.U[1].tke += EPS*v[cellCount+blk.TKE];
            cell.U[1].omega += EPS*v[cellCount+blk.OMEGA];
        }
        cell.decode_conserved(0, 1, 0.0);
        cellCount += nConserved;
    }
    steadystate_core.evalRHS(0.0, 1);
    cellCount = 0;
    foreach (cell; blk.cells) {
        p[cellCount+blk.MASS] = -cell.dUdt[1].mass.im/EPS.im;
        p[cellCount+blk.X_MOM] = -cell.dUdt[1].momentum.x.im/EPS.im;
        p[cellCount+blk.Y_MOM] = -cell.dUdt[1].momentum.y.im/EPS.im;
        if ( blk.myConfig.dimensions == 3 )
            p[cellCount+blk.Z_MOM] = -cell.dUdt[1].momentum.z.im/EPS.im;
        p[cellCount+blk.TOT_ENERGY] = -cell.dUdt[1].total_energy.im/EPS.im;
        if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
            p[cellCount+blk.TKE] = -cell.dUdt[1].tke.im/EPS.im;
            p[cellCount+blk.OMEGA] = -cell.dUdt[1].omega.im/EPS.im;
        }
        cellCount += nConserved;
    }
}


/***************************/
/* FLOW JACOBIAN FUNCTIONS */
/***************************/
void initialisation(ref FluidBlock blk, size_t nPrimitive) {
    blk.cellSave = new FVCell(blk.myConfig);
    foreach(i; 0..blk.MAX_PERTURBED_INTERFACES) blk.ifaceP[i] = new FVInterface(blk.myConfig, false);
}

void ghost_cell_connectivity_for_gradients(ref FluidBlock blk) {
    
    // if a mapped cell has two ghost cell positions (i.e. the same mapped cell attachs to two interior cells) then remove one of the ghostcells.
    foreach (bc; blk.bc) {
	if (bc.type == "exchange_using_mapped_cells") {
	    foreach (i, face0; bc.faces) {
		FVCell ghostcell0;
		if (bc.outsigns[i] == 1) {
		    ghostcell0 = face0.right_cell;
                } else {
		    ghostcell0 = face0.left_cell;
                }
		foreach (j, face1; bc.faces) {
		    FVCell ghostcell1;
		    if (bc.outsigns[j] == 1) {
			ghostcell1 = face1.right_cell;
			if (ghostcell1.global_id == ghostcell0.global_id && ghostcell1.id != ghostcell0.id) {
			    face1.right_cell = ghostcell0; 
			    face1.left_cell.cell_cloud = [];
			    face1.left_cell.cell_cloud ~= face1.left_cell;
			    foreach(face; face1.left_cell.iface) {
				if(face.left_cell.id != face1.left_cell.id) face1.left_cell.cell_cloud ~= face.left_cell; 
				if(face.right_cell.id != face1.left_cell.id) face1.left_cell.cell_cloud ~= face.right_cell; 
			    }
                            ghostcell1.global_id = ghostcell1.id;
			    ghostcell1.is_interior_to_domain = false;
			}
		    } else {
			ghostcell1 = face1.left_cell;
			if (ghostcell1.global_id == ghostcell0.global_id && ghostcell1.id != ghostcell0.id)  {
			    face1.left_cell = ghostcell0;
			    face1.right_cell.cell_cloud = [];
			    face1.right_cell.cell_cloud ~= face1.right_cell;
			    foreach(face; face1.right_cell.iface) {
				if(face.left_cell.id != face1.right_cell.id) face1.right_cell.cell_cloud ~= face.left_cell; 
				if(face.right_cell.id != face1.right_cell.id) face1.right_cell.cell_cloud ~= face.right_cell; 
			    }
                            ghostcell1.global_id = ghostcell1.id;
			    ghostcell1.is_interior_to_domain = false;
                        }
                    }
		}
	    }
	}
    }


    // log the ids of current ghostcells and local cells to make sure we don't add them again
    size_t[] ghost_cell_global_id_list = [];
    foreach (bc; blk.bc) {
        if (bc.type == "exchange_using_mapped_cells") {
            foreach (gc; bc.ghostcells) {
		if(!ghost_cell_global_id_list.canFind(gc.global_id)) {
		    ghost_cell_global_id_list ~= gc.global_id;
		}
	    }
	}
    }

    size_t[] interior_cell_global_id_list;
    foreach(cell; blk.cells) {
        interior_cell_global_id_list ~= cell.global_id;
    }
    
    // collect the extra ring of ghost cells (we will collect the cells in each mapped cells cell cloud (this should be sufficient)
    foreach (bcond; blk.bc) {
	foreach (gce; bcond.preReconAction) {
            auto mygce = cast(GhostCellMappedCellCopy)gce;
            if (mygce && !blk.myConfig.in_mpi_context) {
		foreach (i, bface; bcond.faces) {
		    auto mapped_cell = mygce.get_mapped_cell(i); FVCell ghostcell;
		    if (bcond.outsigns[i] == 1) {
			ghostcell = bface.right_cell;
		    } else {
			ghostcell = bface.left_cell;
		    }
		    ghostcell.outsign = mapped_cell.outsign;
		    ghostcell.in_turbulent_zone = mapped_cell.in_turbulent_zone;
		    foreach(cell; mapped_cell.cell_cloud) {
			if(!ghost_cell_global_id_list.canFind(cell.global_id) && !interior_cell_global_id_list.canFind(cell.global_id)) {
			    // make a new cell and copy the neighbour blocks cell information
			    FVCell new_cell = new FVCell(blk.myConfig);
			    new_cell.copy_values_from(cell, CopyDataOption.all);
			    new_cell.global_id = cell.global_id;
			    new_cell.is_interior_to_domain = cell.is_interior_to_domain;
			    new_cell.in_turbulent_zone = cell.in_turbulent_zone;
			    auto nsp = blk.myConfig.gmodel.n_species;
			    auto nmodes = blk.myConfig.gmodel.n_modes;
			    new_cell.gradients = new LSQInterpGradients(nsp, nmodes);
			    new_cell.gradients.copy_values_from(cell.gradients);
			    new_cell.dqdQ = cell.dqdQ;
			    // add the new cell to the boundary ghostcell list
			    bcond.ghostcells ~= new_cell;
			    ghost_cell_global_id_list ~= new_cell.global_id;
			}
		    }
		}
		//writeln("ghost_cells: ", blk.id, ", ", ghost_cell_global_id_list);
	    }
	}
    }
    // collect the interfaces of the nearest ghostcells
    // make a note of the interface along the boundary - we don't want to copy them from the neighbour block
    string[] ghost_interface_global_id_list;
    foreach (bcond; blk.bc) {
        if (bcond.type == "exchange_using_mapped_cells") {
	    foreach (bface; bcond.faces) {
		ghost_interface_global_id_list ~= bface.global_id;
	    }
	}
    }
    
    foreach (bcond; blk.bc) {
	bool found_mapped_cell_bc = false;
        foreach (gce; bcond.preReconAction) {
            auto mygce = cast(GhostCellMappedCellCopy)gce;
            if (mygce && !blk.myConfig.in_mpi_context) {
                found_mapped_cell_bc = true;
		foreach (i, bface; bcond.faces) {
                    auto mapped_cell = mygce.get_mapped_cell(i);
                    FVCell ghost_cell; FVCell interior_cell;
                    if (bcond.outsigns[i] == 1) {
                        ghost_cell = bface.right_cell;
                        interior_cell = bface.left_cell;
                    } else {
                        ghost_cell = bface.left_cell;
                        interior_cell = bface.right_cell;
                    }
                    foreach(iface; mapped_cell.iface) {
                        if(!ghost_interface_global_id_list.canFind(iface.global_id)) {
			    FVInterface new_face = new FVInterface(blk.myConfig, false);
                            new_face.copy_values_from(iface, CopyDataOption.all);
                            new_face.global_id = iface.global_id;
                            bcond.neighbour_block_faces ~= new_face;
                            ghost_interface_global_id_list ~= new_face.global_id;
                        }
                    }
                }
		//writeln("ghost_ifaces: ", blk.id, ", ", ghost_interface_global_id_list);
	    }
	}
    }
                
    // fill left_cell && right_cell for the added interfaces
    foreach (bcond; blk.bc) {
	bool found_mapped_cell_bc = false;
        foreach (gce; bcond.preReconAction) {
            auto mygce = cast(GhostCellMappedCellCopy)gce;
            if (mygce && !blk.myConfig.in_mpi_context) {
                found_mapped_cell_bc = true;
		string[] processed_interface_global_id;
                foreach (i, bface; bcond.faces) {
                    auto mapped_cell = mygce.get_mapped_cell(i);
                    FVCell ghost_cell; FVCell interior_cell;
                    if (bcond.outsigns[i] == 1) {
                        ghost_cell = bface.right_cell;
                        interior_cell = bface.left_cell;
                    } else {
                        ghost_cell = bface.left_cell;
                        interior_cell = bface.right_cell;
                    }
                    foreach(iface_ext; mapped_cell.iface) {
                        if(iface_ext.global_id != bface.global_id) { // we already have the bface left_cell & right_cell
                            size_t lid = iface_ext.left_cell.global_id;
                            size_t rid = iface_ext.right_cell.global_id;
                            foreach(cell; bcond.ghostcells) {
                                if(cell.global_id == lid) {
                                    foreach(iface_int; bcond.neighbour_block_faces) {
                                        if(iface_ext.global_id == iface_int.global_id) {
                                            iface_int.left_cell = cell;
                                        }
                                    }
                                }
                                if(cell.global_id == rid) {
                                    foreach(iface_int; bcond.neighbour_block_faces) {
                                        if(iface_ext.global_id == iface_int.global_id) {
                                            iface_int.right_cell = cell;
                                        }                                                
                                    }
                                }
                            }
                        }
                    }
                }
	    }
	}
    }
     
    // attach the interfaces to first ring of ghost cells
    foreach (bcond; blk.bc) {
	bool found_mapped_cell_bc = false;
        foreach (gce; bcond.preReconAction) {
            auto mygce = cast(GhostCellMappedCellCopy)gce;
            if (mygce && !blk.myConfig.in_mpi_context) {
                found_mapped_cell_bc = true;
		foreach (i, bface; bcond.faces) {
                    auto mapped_cell = mygce.get_mapped_cell(i);
                    FVCell ghost_cell; FVCell interior_cell;
                    if (bcond.outsigns[i] == 1) {
                        ghost_cell = bface.right_cell;
                        interior_cell = bface.left_cell;
                    } else {
                        ghost_cell = bface.left_cell;
                        interior_cell = bface.right_cell;
                    }
		    if(ghost_cell.iface.length < 1) { // otherwise duplicates will collect each interface twice
			foreach(j, iface_ext; mapped_cell.iface) {
			    if(bface.global_id == iface_ext.global_id) {
				ghost_cell.iface ~= bface;
				ghost_cell.outsign[j] = mapped_cell.outsign[j];
			    } else {
				foreach(iface_int; bcond.neighbour_block_faces) {
				    if(iface_ext.global_id == iface_int.global_id) {
					ghost_cell.iface ~= iface_int;
					ghost_cell.outsign[j] = mapped_cell.outsign[j];
				    }
				}
				// we might have a boundary face which isn't the current bface...
				foreach (iface_int; bcond.faces) {
				    if(iface_ext.global_id == iface_int.global_id) {
					ghost_cell.iface ~= iface_int;
					ghost_cell.outsign[j] = mapped_cell.outsign[j];
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

    // attach interfaces to outer ghost cells
    foreach (bcond; blk.bc) {
        if (bcond.type == "exchange_using_mapped_cells") {
	    foreach (cell; bcond.ghostcells) {
		if(cell.iface.length < 1) {
		    foreach(face; bcond.neighbour_block_faces) {
			if(face.left_cell.global_id == cell.global_id) cell.iface ~= face;
			if(face.right_cell.global_id == cell.global_id) cell.iface ~= face;
		    }
		    //writeln("ghost_cell ifaces: ", blk.id, ", ", cell.global_id, ", ", cell.iface.length);
		}
	    }
	}
    }
    
    // collect the cell cloud for the inviscid gradients
    foreach (bcond; blk.bc) {
	if (bcond.type == "exchange_using_mapped_cells") {
	    foreach (i, bface; bcond.faces) {
		FVCell ghost_cell; FVCell interior_cell;
		if (bcond.outsigns[i] == 1) {
		    ghost_cell = bface.right_cell;
		    interior_cell = bface.left_cell;
		} else {
		    ghost_cell = bface.left_cell;
		    interior_cell = bface.right_cell;
		}
		if(ghost_cell.cell_cloud.length < 1) { // otherwise duplicates will collect each cell twice
		    auto nsp = blk.myConfig.gmodel.n_species;
		    auto nmodes = blk.myConfig.gmodel.n_modes;
		    ghost_cell.gradients = new LSQInterpGradients(nsp, nmodes);
		    ghost_cell.ws = new LSQInterpWorkspace();
		    ghost_cell.cell_cloud ~= ghost_cell;
		    foreach(face; ghost_cell.iface) {
			if(face.left_cell.global_id != ghost_cell.global_id) ghost_cell.cell_cloud ~= face.left_cell;
			if(face.right_cell.global_id != ghost_cell.global_id) ghost_cell.cell_cloud ~= face.right_cell;
		    }
		    ghost_cell.ws.assemble_and_invert_normal_matrix(ghost_cell.cell_cloud, blk.myConfig.dimensions, 0);
		}
	    }
	}	
    }

    // copy limiter values incase we have frozen the limiter
    foreach (bcond; blk.bc) {
	if (bcond.ghost_cell_data_available == false) { continue; }
	// Proceed to do some work only if we have ghost cells in which to insert gradients.
	bool found_mapped_cell_bc = false;
	foreach (gce; bcond.preReconAction) {
	    auto mygce = cast(GhostCellMappedCellCopy)gce;
	    if (mygce && !blk.myConfig.in_mpi_context) {
		found_mapped_cell_bc = true;
		// There is a mapped-cell backing the ghost cell, so we can copy its gradients.
		foreach (i, f; bcond.faces) {
		    // Only FVCell objects in an unstructured-grid are expected to have
		    // precomputed gradients.  There will be an initialized reference
		    // in the FVCell object of a structured-grid block, so we need to
		    // test and avoid copying from such a reference.
		    auto mapped_cell_grad = mygce.get_mapped_cell(i).gradients;
		    if (bcond.outsigns[i] == 1) {
			if (mapped_cell_grad) {
			    f.right_cell.gradients.copy_values_from(mapped_cell_grad);
			} else {
			    // Fall back to looking over the face for suitable gradient data.
			    f.right_cell.gradients.copy_values_from(f.left_cell.gradients);
			}
		    } else {
			if (mapped_cell_grad) {
			    f.left_cell.gradients.copy_values_from(mapped_cell_grad);
			} else {
			    f.left_cell.gradients.copy_values_from(f.right_cell.gradients);
			}
		    }
		} // end foreach f
	    } // end if (mygce)
	} // end foreach gce
    }
        
    // collect the cloud for the viscous gradients
    foreach (bcond; blk.bc) {
	if (bcond.type == "exchange_using_mapped_cells") {
	    // add spatial gradient workspace to ghost cells
	    foreach(cell; bcond.ghostcells) {
		if (cell.iface.length >= 3) {
		    cell.grad = new FlowGradients(blk.myConfig);
		    cell.ws_grad = new WLSQGradWorkspace();
		    // add stencil
		    cell.cloud_pos ~= &(cell.pos[0]);
		    cell.cloud_fs ~= cell.fs;
		    foreach(iface; cell.iface) {
			cell.cloud_pos ~= &(iface.pos);
			cell.cloud_fs ~= iface.fs;
		    }
		    cell.grad.set_up_workspace_leastsq(cell.cloud_pos, cell.pos[0], false, cell.ws_grad);
		    cell.grad.gradients_leastsq(cell.cloud_fs, cell.cloud_pos, cell.ws_grad);
		}
	    }
	}
    }
}

string computeGhostCellDerivatives(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "blk.cellSave.copy_values_from(bcells[0], CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "bcells[0].fs."~varName~" += EPS;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(bcells[0].fs.gas);";
    }
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "qP[blk.MASS] = pcell.fs.gas.rho;";
    codeStr ~= "qP[blk.X_MOM] = pcell.fs.vel.x;";
    codeStr ~= "qP[blk.Y_MOM] = pcell.fs.vel.y;";
    codeStr ~= "if (blk.myConfig.dimensions == 3) {";
    codeStr ~= "qP[blk.Z_MOM] = pcell.fs.vel.z;";
    codeStr ~= "}";
    codeStr ~= "qP[blk.TOT_ENERGY] = pcell.fs.gas.p;";
    codeStr ~= "if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {";
    codeStr ~= "qP[blk.TKE] = pcell.fs.tke;";
    codeStr ~= "qP[blk.OMEGA] = pcell.fs.omega;";
    codeStr ~= "}";
    codeStr ~= "bcells[0].copy_values_from(blk.cellSave, CopyDataOption.all);";
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "dqdQ[blk.MASS][" ~ posInArray ~ "] = qP[blk.MASS].im/(EPS.im);";
    codeStr ~= "dqdQ[blk.X_MOM][" ~ posInArray ~ "] = qP[blk.X_MOM].im/(EPS.im);";
    codeStr ~= "dqdQ[blk.Y_MOM][" ~ posInArray ~ "] = qP[blk.Y_MOM].im/(EPS.im);";
    codeStr ~= "if (blk.myConfig.dimensions == 3) {";
    codeStr ~= "dqdQ[blk.Z_MOM][" ~ posInArray ~ "] = qP[blk.Z_MOM].im/(EPS.im);";
    codeStr ~= "}";
    codeStr ~= "dqdQ[blk.TOT_ENERGY][" ~ posInArray ~ "] = qP[blk.TOT_ENERGY].im/(EPS.im);";         
    codeStr ~= "if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {";
    codeStr ~= "dqdQ[blk.TKE][" ~ posInArray ~ "] = qP[blk.TKE].im/(EPS.im);";
    codeStr ~= "dqdQ[blk.OMEGA][" ~ posInArray ~ "] = qP[blk.OMEGA].im/(EPS.im);";
    codeStr ~= "}";
    return codeStr;
}

void fill_boundary_conditions(FluidBlock blk, size_t np, size_t orderOfJacobian, number EPS, bool transformToConserved, bool preconditionMatrix) {
    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;

    // initialise some re-used data objects here
    //number[][] dRdq; number[][] dqdQ; number[][] Aext; number[] qP;
    //qP.length = np; dRdq.length = np; dqdQ.length = np; Aext.length = np; 
    //foreach (ref a; dRdq) a.length = np;
    //foreach (ref a; dqdQ) a.length = np;
    //foreach (ref a; Aext) a.length = np;
    
    foreach ( bndary; blk.bc ) {
        if ( bndary.type != "exchange_using_mapped_cells") {
            foreach ( bi, bface; bndary.faces) {                
		// initialise some re-used data objects here
		number[][] dRdq; number[][] dqdQ; number[][] Aext; number[] qP;
		qP.length = np; dRdq.length = np; dqdQ.length = np; Aext.length = np; 
		foreach (ref a; dRdq) a.length = np;
		foreach (ref a; dqdQ) a.length = np;
		foreach (ref a; Aext) a.length = np;
		
		// collect interior boundary cells (bcells) and exterior ghost cell (pcell)
                FVCell[] bcells; FVCell pcell;
                if (bndary.outsigns[bi] == 1) {
                    bcells ~= bface.left_cell;
                    pcell = bface.right_cell;
                } else {
                    bcells ~= bface.right_cell;
                    pcell = bface.left_cell;
                }
                
                /* form dqdQ - ghost cell derivatives */

                // 0th perturbation: rho
                mixin(computeGhostCellDerivatives("gas.rho", "blk.MASS", true));
                // 1st perturbation: u
                mixin(computeGhostCellDerivatives("vel.refx", "blk.X_MOM", false));
                // 2nd perturbation: v
                mixin(computeGhostCellDerivatives("vel.refy", "blk.Y_MOM", false));
                if ( blk.myConfig.dimensions == 3 )
                    mixin(computeGhostCellDerivatives("vel.refz", "blk.Z_MOM", false));
                // 3rd perturbation: P
                mixin(computeGhostCellDerivatives("gas.p", "blk.TOT_ENERGY", true));
                if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    // 4th perturbation: k
                    mixin(computeGhostCellDerivatives("tke", "blk.TKE", false));
                    // 5th perturbation: omega
                    mixin(computeGhostCellDerivatives("omega", "blk.OMEGA", false));
                }
		
		pcell.dqdQ = dqdQ;
	    }
	}
    }
}

void apply_boundary_conditions(ref SMatrix!number A, FluidBlock blk, size_t np, size_t orderOfJacobian, number EPS, bool transformToConserved, bool preconditionMatrix) {

    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;
    
    // initialise some re-used data objects here
    number[][] dRdq; number[][] dqdQ; number[][] Aext; number[] qP;
    qP.length = np; dRdq.length = np; dqdQ.length = np; Aext.length = np; 
    foreach (ref a; dRdq) a.length = np;
    foreach (ref a; dqdQ) a.length = np;
    foreach (ref a; Aext) a.length = np;

    foreach ( bndary; blk.bc ) {
        if ( bndary.type != "exchange_using_mapped_cells") {
            foreach ( bi, bface; bndary.faces) {                
                // collect interior boundary cells (bcells) and exterior ghost cell (pcell)
                FVCell[] bcells; FVCell pcell;
                if (bndary.outsigns[bi] == 1) {
                    bcells ~= bface.left_cell;
                    pcell = bface.right_cell;
                } else {
                    bcells ~= bface.right_cell;
                    pcell = bface.left_cell;
                }
                
                /* form dqdQ - ghost cell derivatives */

                // 0th perturbation: rho
                mixin(computeGhostCellDerivatives("gas.rho", "blk.MASS", true));
                // 1st perturbation: u
                mixin(computeGhostCellDerivatives("vel.refx", "blk.X_MOM", false));
                // 2nd perturbation: v
                mixin(computeGhostCellDerivatives("vel.refy", "blk.Y_MOM", false));
                if ( blk.myConfig.dimensions == 3 )
                    mixin(computeGhostCellDerivatives("vel.refz", "blk.Z_MOM", false));
                // 3rd perturbation: P
                mixin(computeGhostCellDerivatives("gas.p", "blk.TOT_ENERGY", true));
                if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    // 4th perturbation: k
                    mixin(computeGhostCellDerivatives("tke", "blk.TKE", false));
                    // 5th perturbation: omega
                    mixin(computeGhostCellDerivatives("omega", "blk.OMEGA", false));
                }

		///////
		//pcell.dqdQ = dqdQ;
		/////
                
                /* form dRdq */                
                // TODO: Currently only works for nearest-neighbour reconstruction stencil. Think about the MLP limiter.
                
                if (orderOfJacobian > 1 || ( blk.myConfig.viscous && preconditionMatrix == false ) ) {
                    
                    size_t[] idList;
                    foreach ( face; bcells[0].iface) {
                        FVCell lftCell = face.left_cell;
                        FVCell rghtCell = face.right_cell;
                        if (lftCell.id != bcells[0].id && idList.canFind(lftCell.id) == false && lftCell.is_interior_to_domain) { //lftCell.is_interior) { //lftCell.id < ghost_cell_start_id) 
                            bcells ~= lftCell;
                            idList ~= lftCell.id;
                        }
                        if (rghtCell.id != bcells[0].id && idList.canFind(rghtCell.id) == false && rghtCell.is_interior_to_domain) { //rghtCell.is_interior) { //rghtCell.id < ghost_cell_start_id) {
                            bcells ~= rghtCell; 
                            idList ~= rghtCell.id;
                        }
                    }
		    /*
                    //
                    foreach (cell; bcells) {
                        foreach ( face; cell.iface) {
                            FVCell lftCell = face.left_cell;
                            FVCell rghtCell = face.right_cell;
                            if (lftCell.id != bcells[0].id && idList.canFind(lftCell.id) == false && lftCell.id < ghost_cell_start_id) {
                                bcells ~= lftCell;
                                idList ~= lftCell.id;
                            }
                            if (rghtCell.id != bcells[0].id && idList.canFind(rghtCell.id) == false && rghtCell.id < ghost_cell_start_id) {
                                bcells ~= rghtCell; 
                                idList ~= rghtCell.id;
                            }
                        }
                    }

		    
                    // Turbulent flow needs a larger stencil for 2nd order....
                    foreach (cell; bcells) {
                        foreach ( face; cell.iface) {
                            FVCell lftCell = face.left_cell;
                            FVCell rghtCell = face.right_cell;
                            if (lftCell.id != bcells[0].id && idList.canFind(lftCell.id) == false && lftCell.id < ghost_cell_start_id) {
                                bcells ~= lftCell;
                                idList ~= lftCell.id;
                            }
                            if (rghtCell.id != bcells[0].id && idList.canFind(rghtCell.id) == false && rghtCell.id < ghost_cell_start_id) {
                                bcells ~= rghtCell; 
                                idList ~= rghtCell.id;
                            }
                        }
                    }
                    //

                    //
                    */
                }
                
                pcell.jacobian_cell_stencil ~= bcells;
                
                size_t[] idList;
                foreach ( bcell; bcells) {
                    foreach ( face; bcell.iface) {
                        if ( idList.canFind(face.id) == false ) {
                            pcell.jacobian_face_stencil ~= face;
                            idList ~= face.id;
                        }
                    }
                }
                
                // 0th perturbation: rho
                mixin(computeFluxDerivativesAroundCell("gas.rho", "blk.MASS", true));
                // 1st perturbation: u
                mixin(computeFluxDerivativesAroundCell("vel.refx", "blk.X_MOM", false));
                // 2nd perturbation: v
                mixin(computeFluxDerivativesAroundCell("vel.refy", "blk.Y_MOM", false));
                if ( blk.myConfig.dimensions == 3 )
                    mixin(computeFluxDerivativesAroundCell("vel.refz", "blk.Z_MOM", false));
                // 3rd perturbation: P
                mixin(computeFluxDerivativesAroundCell("gas.p", "blk.TOT_ENERGY", true));
                if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    // 4th perturbation: k
                    mixin(computeFluxDerivativesAroundCell("tke", "blk.TKE", false));
                    // 5th perturbation: omega
                    mixin(computeFluxDerivativesAroundCell("omega", "blk.OMEGA", false));
                } 
  
                foreach(bcell; pcell.jacobian_cell_stencil) {
                    if (bcell.id < ghost_cell_start_id) {
			number integral;
			number volInv = 1.0 / bcell.volume[0];
			for ( size_t ip = 0; ip < np; ++ip ) {
			    for ( size_t jp = 0; jp < np; ++jp ) {
				integral = 0.0;
				foreach(fi, iface; bcell.iface) {
				    if(bcells[0].global_id == 4 && bcell.global_id == 8) {
                                        //writef(" %d    %d    %.6e    %.6e   %d    %d    %.4e    %.12e    %.12e    %.12e \n", bcells[0].global_id, bcell.global_id, iface.pos.x.re, iface.pos.y.re, ip, jp, bcell.outsign[fi].re, iface.area[0].re, iface.dFdU[ip][jp].re, volInv.re);
                                    }
				    
				    integral -= bcell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
				}
				number entry = volInv * integral + bcell.dQdU[ip][jp];                    
				dRdq[ip][jp] = -entry;
			    }
			}

			if(bcells[0].global_id == 4 && bcell.global_id == 8) {
                            //writeln(bcells[0].global_id, " ------ ");
                            //writeln(dqdQ);
                            //writeln(dRdq);
                        }
			
			// perform matrix-matrix multiplication
			for (size_t i = 0; i < np; i++) {
			    for (size_t j = 0; j < np; j++) {
				Aext[i][j] = 0;
				for (size_t k = 0; k < np; k++) {
				    Aext[i][j] += dRdq[i][k]*dqdQ[k][j];
				}
			    }
			}
                    
                      
			// transform the sub-matrix from primitive form to conservative form
			if (transformToConserved) {
			    auto gmodel = blk.myConfig.gmodel;
			    // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
			    number gamma = gmodel.gamma(bcells[0].fs.gas);
			    // form inverse transformation matrix
			    blk.Minv[blk.MASS,blk.MASS] = to!number(1.0);
			    blk.Minv[blk.MASS,blk.X_MOM] = to!number(0.0);
			    blk.Minv[blk.MASS,blk.Y_MOM] = to!number(0.0);
			    blk.Minv[blk.MASS,blk.TOT_ENERGY] = to!number(0.0);
			    // second row
			    blk.Minv[blk.X_MOM,blk.MASS] = -bcells[0].fs.vel.x/bcells[0].fs.gas.rho;
			    blk.Minv[blk.X_MOM,blk.X_MOM] = 1.0/bcells[0].fs.gas.rho;
			    blk.Minv[blk.X_MOM,blk.Y_MOM] = to!number(0.0);
			    blk.Minv[blk.X_MOM,blk.TOT_ENERGY] = to!number(0.0);
			    // third row
			    blk.Minv[blk.Y_MOM,blk.MASS] = -bcells[0].fs.vel.y/bcells[0].fs.gas.rho;
			    blk.Minv[blk.Y_MOM,blk.X_MOM] = to!number(0.0);
			    blk.Minv[blk.Y_MOM,blk.Y_MOM] = 1.0/bcells[0].fs.gas.rho;
			    blk.Minv[blk.Y_MOM,blk.TOT_ENERGY] = to!number(0.0);
			    // fourth row
			    blk.Minv[blk.TOT_ENERGY,blk.MASS] = 0.5*(gamma-1.0)*(bcells[0].fs.vel.x*bcells[0].fs.vel.x+bcells[0].fs.vel.y*bcells[0].fs.vel.y+bcells[0].fs.vel.z*bcells[0].fs.vel.z);
			    blk.Minv[blk.TOT_ENERGY,blk.X_MOM] = -bcells[0].fs.vel.x*(gamma-1);
			    blk.Minv[blk.TOT_ENERGY,blk.Y_MOM] = -bcells[0].fs.vel.y*(gamma-1);
			    blk.Minv[blk.TOT_ENERGY,blk.TOT_ENERGY] = gamma-1.0;
			    
			    if (blk.myConfig.dimensions == 3) {
				blk.Minv[blk.MASS,blk.Z_MOM] = to!number(0.0);
				blk.Minv[blk.X_MOM,blk.Z_MOM] = to!number(0.0);
				blk.Minv[blk.Y_MOM,blk.Z_MOM] = to!number(0.0);
				blk.Minv[blk.TOT_ENERGY,blk.Z_MOM] = -bcells[0].fs.vel.z*(gamma-1);
				
				blk.Minv[blk.Z_MOM,blk.MASS] = -bcells[0].fs.vel.z/bcells[0].fs.gas.rho;
				blk.Minv[blk.Z_MOM,blk.X_MOM] = to!number(0.0);
				blk.Minv[blk.Z_MOM,blk.Y_MOM] = to!number(0.0);
				blk.Minv[blk.Z_MOM,blk.Z_MOM] = 1.0/bcells[0].fs.gas.rho;
				blk.Minv[blk.Z_MOM,blk.TOT_ENERGY] = to!number(0.0);
			    }
			    
			    if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
				blk.Minv[blk.MASS,blk.TKE] = to!number(0.0);
				blk.Minv[blk.MASS,blk.OMEGA] = to!number(0.0);
				// second row
				blk.Minv[blk.X_MOM,blk.TKE] = to!number(0.0);
				blk.Minv[blk.X_MOM,blk.OMEGA] = to!number(0.0);
				// third row
				blk.Minv[blk.Y_MOM,blk.TKE] = to!number(0.0);
				blk.Minv[blk.Y_MOM,blk.OMEGA] = to!number(0.0);
				// fourth row
				blk.Minv[blk.TOT_ENERGY,blk.TKE] = -(gamma-1.0); 
				blk.Minv[blk.TOT_ENERGY,blk.OMEGA] = to!number(0.0);
				// fifth row
				blk.Minv[blk.TKE,blk.MASS] = -bcells[0].fs.tke/bcells[0].fs.gas.rho;
				blk.Minv[blk.TKE,blk.X_MOM] = to!number(0.0);
				blk.Minv[blk.TKE,blk.Y_MOM] = to!number(0.0);
				blk.Minv[blk.TKE,blk.TOT_ENERGY] = to!number(0.0);
				blk.Minv[blk.TKE,blk.TKE] = 1.0/bcells[0].fs.gas.rho;
				blk.Minv[blk.TKE,blk.OMEGA] = to!number(0.0);
				// sixth row
				blk.Minv[blk.OMEGA,blk.MASS] = -bcells[0].fs.omega/bcells[0].fs.gas.rho;
				blk.Minv[blk.OMEGA,blk.X_MOM] = to!number(0.0);
				blk.Minv[blk.OMEGA,blk.Y_MOM] = to!number(0.0);
				blk.Minv[blk.OMEGA,blk.TOT_ENERGY] = to!number(0.0);
				blk.Minv[blk.OMEGA,blk.TKE] = to!number(0.0);
				blk.Minv[blk.OMEGA,blk.OMEGA] = 1.0/bcells[0].fs.gas.rho;
				
				if (blk.myConfig.dimensions == 3) {
				    blk.Minv[blk.Z_MOM,blk.TKE] = to!number(0.0);
				    blk.Minv[blk.Z_MOM,blk.OMEGA] = to!number(0.0);
				    blk.Minv[blk.TKE,blk.Z_MOM] = to!number(0.0);
				    blk.Minv[blk.OMEGA,blk.Z_MOM] = to!number(0.0);
				}
			    }
			    
			    //writeln(blk.Minv);
			    number[][] tmp;
			    tmp.length = np;
			    foreach (ref a; tmp) a.length = np;
			    for (size_t i = 0; i < np; i++) {
				for (size_t j = 0; j < np; j++) {
				    tmp[i][j] = to!number(0.0);
				    for (size_t k = 0; k < np; k++) {
					tmp[i][j] += Aext[i][k]*blk.Minv[k,j];
				    }
				}
			    }
			    
			    foreach (i; 0..np) {
				foreach (j; 0..np) {
				    Aext[i][j] = tmp[i][j];
				    blk.Minv[i,j] = to!number(0.0);
				}
			    }
			}
						
			
			// add correction to boundary entry in Jacobian
			size_t I, J;
			for ( size_t ip = 0; ip < np; ++ip ) {
			    I = bcell.id*np + ip; // column index
			    for ( size_t jp = 0; jp < np; ++jp ) {
				J = bcells[0].id*np + jp; // row index
				if(bcells[0].global_id == 4) {
				    //writef("check: %d    %d    %.12e    %.12e    %.12e \n", ip, jp, A[J,I], Aext[ip][jp], A[J,I] + Aext[ip][jp]);
				}			
				
				A[J,I] = A[J,I] + Aext[ip][jp];
			    }
			}
		    }
                }

                // clear the interface flux Jacobian entries
                foreach (iface; pcell.jacobian_face_stencil) {
                    foreach (i; 0..iface.dFdU.length) {
                        foreach (j; 0..iface.dFdU[i].length) {
                            iface.dFdU[i][j] = 0.0;
                        }
                    }
                }

                // clear the cell source term Jacobian entries
                foreach (cell; pcell.jacobian_cell_stencil) {
                    foreach (i; 0..cell.dQdU.length) {
                        foreach (j; 0..cell.dQdU[i].length) {
                            cell.dQdU[i][j] = 0.0;
                        }
                    }
                }
                
                pcell.jacobian_cell_stencil = [];
                pcell.jacobian_face_stencil = [];
            }
        }
    }
}

void form_external_flow_jacobian_block_phase0(FluidBlock blk, size_t np, int orderOfJacobian, number EPS) {
    
    // construct the cell and face stencils for the primary ghost cells
    size_t[] processed_cell_id_list;
    foreach (bc; blk.bc) {
	if (bc.type == "exchange_using_mapped_cells") {
	    foreach (i, bface; bc.faces) {
		FVCell interior_cell; FVCell ghost_cell;
		if (bc.outsigns[i] == 1) {
                    interior_cell = bface.left_cell;
                    ghost_cell = bface.right_cell;
                } else {
                    interior_cell = bface.right_cell;
                    ghost_cell = bface.left_cell;
                }
		if (!processed_cell_id_list.canFind(ghost_cell.global_id)) {		    
		    processed_cell_id_list ~= ghost_cell.global_id;
		    // cell stencil
		    size_t[] cell_id_list;
		    

		    foreach (icell; ghost_cell.cell_cloud) {
			foreach (jcell; icell.cell_cloud) {
			    if (!cell_id_list.canFind(jcell.global_id) && jcell.is_interior_to_domain && jcell.cell_cloud.length > 1) {
				ghost_cell.jacobian_cell_stencil ~= jcell;
				cell_id_list ~= jcell.global_id;
			    }
			}
		    }

		    // face stencil
		    string[] face_id_list;
		    foreach (cell; ghost_cell.jacobian_cell_stencil) {
			if (cell.is_interior_to_domain) {
			    foreach(face; cell.iface) {
				if (!face_id_list.canFind(face.global_id)) {
				    ghost_cell.jacobian_face_stencil ~= face;
				    face_id_list ~= face.global_id;
				} 
			    }
			}
		    }
		}
	    }
	}
    }
    
    // construct the cell and face stencils for the secondary ghost cells
    foreach (bc; blk.bc) {
	if (bc.type == "exchange_using_mapped_cells") {
	    foreach (gcell; bc.ghostcells) {
		if (!processed_cell_id_list.canFind(gcell.global_id)) {		    
		    processed_cell_id_list ~= gcell.global_id;		    
		    // cell stencil
		    size_t[] cell_id_list;
		    foreach(face; gcell.iface) {
			if(face.left_cell.global_id != gcell.global_id && !cell_id_list.canFind(face.left_cell.global_id) && face.left_cell.is_interior_to_domain && face.left_cell.cell_cloud.length > 1) {
			    gcell.jacobian_cell_stencil ~= face.left_cell;
			    cell_id_list ~= face.left_cell.global_id;
			}
			if(face.right_cell.global_id != gcell.global_id && !cell_id_list.canFind(face.right_cell.global_id) && face.right_cell.is_interior_to_domain && face.right_cell.cell_cloud.length > 1) {
			    gcell.jacobian_cell_stencil ~= face.right_cell;
			    processed_cell_id_list ~= face.right_cell.global_id;
			}
		    }
		    foreach(cell; gcell.jacobian_cell_stencil) {
			foreach(cloud; cell.cell_cloud) {
			    if(!cell_id_list.canFind(cloud.global_id) && cloud.is_interior_to_domain && cloud.cell_cloud.length > 1) {
				gcell.jacobian_cell_stencil ~= cloud;
				cell_id_list ~= cloud.global_id;
			    }
			}
		    }
		    // face stencil
		    string[] face_id_list;
		    foreach(cell; gcell.jacobian_cell_stencil) {
			foreach(face; cell.iface) {
			    if (!face_id_list.canFind(face.global_id)) {
				gcell.jacobian_face_stencil ~= face;
				face_id_list ~= face.global_id;
			    }
			}
		    }
		}
	    }
	}
    }
    
    // compute Jacobian entries
    foreach (bc; blk.bc) {
	if (bc.type == "exchange_using_mapped_cells") {
	    foreach (cell; bc.ghostcells) {
		if (cell.global_id < ghost_cell_start_id) {
		    size_t dim = blk.myConfig.dimensions;
                    construct_flow_jacobian_for_boundary_cells(cell, cell, blk, dim, np, orderOfJacobian, EPS);
		    // clear stencil
		    cell.jacobian_face_stencil = [];
		    cell.jacobian_cell_stencil = [];
		} // end foreach bc.faces
	    }
	} // end if mapped_cells
    } // end foreach blk.bc
    
    foreach (bc; blk.bc) {
	if (bc.type == "exchange_using_mapped_cells") {
	    foreach (cell; bc.ghostcells) {
		//writeln(blk.id, ", ", cell.global_id, ", ", cell.pcell_global_coord_list, ", ", cell.ecell_global_coord_list );
	    }
	}
    }

    // place all values into the mapped cells
    size_t[] processed_global_id_list = [];
    foreach (bc; blk.bc) {
	if (bc.type == "exchange_using_mapped_cells") {
	    foreach (i, bface; bc.faces) {
		FVCell interior_cell; FVCell ghost_cell;
		if (bc.outsigns[i] == 1) {
		    interior_cell = bface.left_cell;
		    ghost_cell = bface.right_cell;
		} else {
		    interior_cell = bface.right_cell;
		    ghost_cell = bface.left_cell;
		}
		foreach(cell; ghost_cell.cell_cloud) {
		    //if (!processed_global_id_list.canFind(cell.global_id) && cell.pcell_global_coord_list.length > 0) {
		    if (cell.pcell_global_coord_list.length > 0) {
			interior_cell.pcell_global_coord_list ~= cell.pcell_global_coord_list;
			interior_cell.ecell_global_coord_list ~= cell.ecell_global_coord_list;
			interior_cell.entry_list ~= cell.entry_list;
			processed_global_id_list ~= cell.global_id;
		    }
		}
	    }
	}
    }
} // end form_external_flow_jacobian_block_phase0()

void form_external_flow_jacobian_block_phase1(size_t np, int orderOfJacobian, number EPS) {
    
    // collect global data for transmission to all blocks
    size_t[] global_pcell_global_coord_list;
    size_t[][] global_ecell_global_coord_list;
    number[][] global_entry_list;
    
    foreach(blk; localFluidBlocks) {
	foreach (bc; blk.bc) {
	    if (bc.type == "exchange_using_mapped_cells") {
		foreach (i, bface; bc.faces) {
		    FVCell interior_cell; FVCell ghost_cell;
		    if (bc.outsigns[i] == 1) {
			interior_cell = bface.left_cell;
			ghost_cell = bface.right_cell;
		    } else {
			interior_cell = bface.right_cell;
			ghost_cell = bface.left_cell;
		    }
		    global_pcell_global_coord_list ~= interior_cell.pcell_global_coord_list;
		    global_ecell_global_coord_list ~= interior_cell.ecell_global_coord_list;
		    global_entry_list ~= interior_cell.entry_list;
		}
	    }
	}
    }
    
    // transmit to all blocks
    foreach(blk; localFluidBlocks) {
	blk.local_pcell_global_coord_list = global_pcell_global_coord_list;
	blk.local_ecell_global_coord_list = global_ecell_global_coord_list;
	blk.local_entry_list = global_entry_list;
    }

} // end form_external_flow_jacobian_block_phase1()

void form_external_flow_jacobian_block_phase2(ref SMatrix!number A, FluidBlock blk, size_t np, int orderOfJacobian, number EPS) {
    // restructure Jacobian data into a sparse matrix
    size_t[] pcell_id; size_t[size_t] pcell_pos_array; size_t[][] ecell_id; number[][] entry_list;

    size_t[] local_cell_ids;
    foreach(cell; blk.cells) {
	local_cell_ids ~= cell.global_id;
    }
    
    foreach(j, id; blk.local_pcell_global_coord_list) {
	if(local_cell_ids.canFind(id)) {
	    if (pcell_id.canFind(id) == false) {
		pcell_id ~= to!size_t(id);
		pcell_pos_array[to!size_t(id)] = to!int(pcell_pos_array.length);
		ecell_id ~= blk.local_ecell_global_coord_list[j];
		entry_list ~= blk.local_entry_list[j];
	    } else {
		size_t pos = pcell_pos_array[id];
		size_t count = 0;
		foreach(k; 0..blk.local_ecell_global_coord_list[j].length) {
		    if (!ecell_id[pos].canFind(blk.local_ecell_global_coord_list[j][k])
			&& !local_cell_ids.canFind(blk.local_ecell_global_coord_list[j][k])) {
			ecell_id[pos] ~= blk.local_ecell_global_coord_list[j][k];
			foreach(l; 0..np*np) {
			    entry_list[pos] ~= blk.local_entry_list[j][l+count];
			}
			count += np*np;
		    } 
		}
	    }
	}   
    }
    
    //writeln("blk: ", blk.id);
    //writeln(pcell_id);
    //writeln(ecell_id);
    //writeln(entry_list);
    //writeln("blk ", blk.id, ": phase 2 end -----------------------------------------------");

    // order the pcells in ascending order
    pcell_id.sort();    
    
    // store in sparse matrix
    //idx = 0;
    A.ia ~= 0;
    foreach (cell; blk.cells) {
	if (pcell_id.canFind(cell.global_id)) {
	    size_t[] ecells; size_t[size_t] pos_ecells; number[][] entries;
	    ecells = ecell_id[pcell_pos_array[cell.global_id]];
	    foreach(i, c; ecells) {
		pos_ecells[c] = to!int(i);
	    }
	    int count = 0;
	    size_t width = np*np;
	    foreach(ecell; ecells) {
		number[] tmp;
		foreach( i; 0..np*np) {
		    tmp ~= entry_list[pcell_pos_array[cell.global_id]][count+i];
		}
		entries ~= tmp;
		count += width;
	    }
	    //writeln("check: ", cell.global_id, ", ", ecells);
	    foreach(i; 0..entries.length) {
		foreach(j; 0..entries[i].length) {
		    //writef("%.12e ", entries[i][j].re);
		}
		//writef("\n");
		//writef("\n");
	    }
	    
	    ecells.sort();
	    for ( size_t ip = 0; ip < np; ++ip ) {
		size_t aaLen = A.aa.length;
		count = 0;
		foreach(ecell; ecells) {
		    for ( size_t jp = 0; jp < np; ++jp ) {
			//size_t idx = ip*np+jp;
			size_t J = ecell*np + jp;
			A.aa ~= entries[pos_ecells[ecell]][ip+jp*np];
			A.ja ~= J;
		    }
		}
		
		
		if (aaLen == A.aa.length) {
		    // no entries have been added for this row, for CRS we must store a dummy 0 entry
		    A.aa ~= to!number(0.0);
		    A.ja ~= 0;
		}
		A.ia ~= A.aa.length;
	    }
	} else { // internal cell (no entries are required this row, for CRS we must store a dummy 0 entry)
            for ( size_t ip = 0; ip < np; ++ip ) {
                A.aa ~= to!number(0.0);
                A.ja ~= 0;
                A.ia ~= A.aa.length;
            }
        } // end else
    }
} // end form_external_flow_jacobian_block_phase1()


void form_external_flow_jacobian_block_phase3(ref SMatrix!number A, FluidBlock blk, size_t np, int orderOfJacobian, number EPS) {
    
    // account for boundary conditions ----------------------------------------------
    // 1: loop around boundaries computing effects for cells
    number[][] dRdq; number[][] dqdQ; number[][] Aext; number[] qP;
    qP.length = np; dRdq.length = np; dqdQ.length = np; Aext.length = np; 
    foreach (ref a; dRdq) a.length = np;
    foreach (ref a; dqdQ) a.length = np;
    foreach (ref a; Aext) a.length = np;
    
    foreach ( bndary; blk.bc ) {
        if ( bndary.type == "exchange_using_mapped_cells") {
            foreach (bi, bface; bndary.faces) {                
                // collect interior boundary cells (bcells) and exterior ghost cell (pcell)
                FVCell interior; FVCell ghost;
                if (bndary.outsigns[bi] == 1) {
                    interior = bface.left_cell;
                    ghost = bface.right_cell;
                } else {
                    interior = bface.right_cell;
                    ghost = bface.left_cell;
                }
                
                bool found = false;
                foreach (cell; interior.cell_cloud) {
                    if (cell.is_interior_to_domain == false) found = true;
                }
                if (found) {
                    FVCell bcell = ghost;
                    FVCell pcell;
		    foreach (cell; interior.cell_cloud) {
			if (cell.is_interior_to_domain == false) pcell = cell;
		    }
                    
                    // form dqdQ - ghost cell derivatives
                    FVCell[] bcells;
                    bcells ~= bcell;
		    dqdQ = pcell.dqdQ;			    
                    
                    // form dRdq                
                    pcell.jacobian_cell_stencil = [];
                    pcell.jacobian_face_stencil = [];
                    pcell.jacobian_cell_stencil ~= bcells;
                    pcell.jacobian_cell_stencil ~= interior;
                    
                    string[] idList;
                    foreach ( cell; pcell.jacobian_cell_stencil) {
                        foreach ( face; cell.iface) {
                            if ( idList.canFind(face.global_id) == false ) {
                                pcell.jacobian_face_stencil ~= face;
                                idList ~= face.global_id;
                            }
                        }
                    }
                    
                    // 0th perturbation: rho
                    mixin(computeFluxDerivativesAroundCell("gas.rho", "blk.MASS", true));
                    // 1st perturbation: u
                    mixin(computeFluxDerivativesAroundCell("vel.refx", "blk.X_MOM", false));
                    // 2nd perturbation: v
                    mixin(computeFluxDerivativesAroundCell("vel.refy", "blk.Y_MOM", false));
                    if ( blk.myConfig.dimensions == 3 )
                        mixin(computeFluxDerivativesAroundCell("vel.refz", "blk.Z_MOM", false));
                    // 3rd perturbation: P
                    mixin(computeFluxDerivativesAroundCell("gas.p", "blk.TOT_ENERGY", true));
                    if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                        // 4th perturbation: k
                        mixin(computeFluxDerivativesAroundCell("tke", "blk.TKE", false));
                        // 5th perturbation: omega
                        mixin(computeFluxDerivativesAroundCell("omega", "blk.OMEGA", false));
                    }
                    
		    

                    foreach(cell; pcell.jacobian_cell_stencil) {
			number integral;
                        number volInv = 1.0 / cell.volume[0];
                        for ( size_t ip = 0; ip < np; ++ip ) {
                            for ( size_t jp = 0; jp < np; ++jp ) {
                                integral = 0.0;
                                foreach(fi, iface; cell.iface) {
				    integral -= cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                                }
                                number entry = volInv * integral + cell.dQdU[ip][jp];                    
                                dRdq[ip][jp] = -entry*-1; // multiply by negative one due to outsign difference
                            }
                        }
                        
                        // perform matrix-matrix multiplication
                        for (size_t i = 0; i < np; i++) {
                            for (size_t j = 0; j < np; j++) {
                                Aext[i][j] = 0;
                                for (size_t k = 0; k < np; k++) {
                                    Aext[i][j] += dRdq[i][k]*dqdQ[k][j];
                                }
                            }
                        }

                        // 3: update entry
                        if(cell.global_id == ghost.global_id) {
			    for ( size_t ip = 0; ip < np; ++ip ) {
				for ( size_t jp = 0; jp < np; ++jp ) {
				    
				    A[interior.id*np+ip,ghost.global_id*np+jp] = A[interior.id*np+ip,ghost.global_id*np+jp] + Aext[jp][ip];
				}
			    }
			}
		    }
		    // clear the interface flux Jacobian entries
		    foreach (iface; pcell.jacobian_face_stencil) {
			foreach (i; 0..iface.dFdU.length) {
			    foreach (j; 0..iface.dFdU[i].length) {
				iface.dFdU[i][j] = 0.0;
			    }
			}
		    }
                    
		    // clear the cell source term Jacobian entries
		    foreach (c; pcell.jacobian_cell_stencil) {
			foreach (i; 0..c.dQdU.length) {
			    foreach (j; 0..c.dQdU[i].length) {
				c.dQdU[i][j] = 0.0;
			    }
			}
		    }
		    pcell.jacobian_cell_stencil = [];
                    pcell.jacobian_face_stencil = [];
                }
            }
        }
    }
}

void construct_flow_jacobian_for_boundary_cells(FVCell pcell, FVCell icell, FluidBlock blk, size_t ndim, size_t np, size_t orderOfJacobian, number EPS) 
{
    // perturb flowstate & compute interface flux sensitivities
    mixin(computeFluxDerivativesAroundCell("gas.rho", "blk.MASS", true));
    mixin(computeFluxDerivativesAroundCell("vel.refx", "blk.X_MOM", false));
    mixin(computeFluxDerivativesAroundCell("vel.refy", "blk.Y_MOM", false));
    if ( blk.myConfig.dimensions == 3 )
        mixin(computeFluxDerivativesAroundCell("vel.refz", "blk.Z_MOM", false));
    mixin(computeFluxDerivativesAroundCell("gas.p", "blk.TOT_ENERGY", true));
    if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
        mixin(computeFluxDerivativesAroundCell("tke", "blk.TKE", false));
        mixin(computeFluxDerivativesAroundCell("omega", "blk.OMEGA", false));
    }

    // compute cell residual sensitivities (flow Jacobian entries) by integrating interface flux sensitivities
    icell.pcell_global_coord_list ~= pcell.global_id;
    size_t[] ecell_coord_list;
    number[] ecell_entry_list;
    foreach(cell; pcell.jacobian_cell_stencil) {
	if (cell.is_interior_to_domain && cell.id < ghost_cell_start_id) { 
	    ecell_coord_list ~= to!int(cell.global_id);
	    number integral;
	    number volInv = 1.0 / cell.volume[0];
	    for ( size_t ip = 0; ip < np; ++ip ) {
		for ( size_t jp = 0; jp < np; ++jp ) {
		    integral = 0.0;
		    
		    if(pcell.global_id == 14 && cell.global_id == 2) {
			//writeln(pcell.global_id, ", ", cell.global_id, ", ", pcell.pos[0].x.re, ", ", pcell.pos[0].y.re, ", ", cell.pos[0].x.re, ", ", cell.pos[0].y.re);
			foreach(fi, iface; cell.iface) {
			    //writef(" %d    %s    %.6e    %.6e   %d    %d    %.4e    %.12e    %.12e    %.12e \n", iface.id, iface.global_id, iface.pos.x.re, iface.pos.y.re, ip, jp, cell.outsign[fi].re, iface.area[0].re, iface.dFdU[ip][jp].re, volInv.re);
			}
		    }



		    foreach(fi, iface; cell.iface) {
			integral -= cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
		    }
		    number JacEntry = volInv * integral + cell.dQdU[ip][jp]; // add source term contribution
		    //if(pcell.global_id == 325 && cell.global_id == 408) writeln(JacEntry, ", ", cell.dQdU[ip][jp]);
		    ecell_entry_list ~= -JacEntry;
		}
                }
	}
    }
    icell.ecell_global_coord_list ~= ecell_coord_list;
    icell.entry_list ~= ecell_entry_list;

    // clear the interface flux Jacobian entries
    foreach (iface; pcell.jacobian_face_stencil) {
        foreach (i; 0..iface.dFdU.length) {
            foreach (j; 0..iface.dFdU[i].length) {
                iface.dFdU[i][j] = 0.0;
            }
        }
    }

    // clear the cell source term Jacobian entries
    foreach (cell; pcell.jacobian_cell_stencil) {
        foreach (i; 0..cell.dQdU.length) {
            foreach (j; 0..cell.dQdU[i].length) {
                cell.dQdU[i][j] = 0.0;
            }
        }
    }
    pcell.jacobian_cell_stencil = [];
    pcell.jacobian_face_stencil = [];
}


void approximate_residual_stencil(FVCell pcell, size_t orderOfJacobian) {
    // used for the steady-state preconditioner (always first order)
    FVCell[] refs_ordered; FVCell[] refs_unordered;
    size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
    size_t[] cell_ids;
    // clear the stencil arrays
    pcell.jacobian_cell_stencil = [];
    pcell.jacobian_face_stencil = [];
    size_t[] face_ids;
    // collect faces
    foreach(face; pcell.iface) {
        pcell.jacobian_face_stencil ~= face;
        face_ids ~= face.id;
    }

    // for each effected face, add the neighbouring cells
    foreach(face; pcell.jacobian_face_stencil) {
        // collect (non-ghost) neighbour cells
        if (cell_ids.canFind(face.left_cell.id) == false && face.left_cell.id < 1_000_000_000) { //.ghost_cell_start_id) {
            refs_unordered ~= face.left_cell;
            pos_array[face.left_cell.id] = refs_unordered.length-1;
            cell_ids ~= face.left_cell.id;
        }
        if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.id < 1_000_000_000) { //ghost_cell_start_id) {
            refs_unordered ~= face.right_cell;
            pos_array[face.right_cell.id] = refs_unordered.length-1;
            cell_ids ~= face.right_cell.id;
        }
        else continue;
    }
    // we collect the outer most faces as well to capture the viscous effects on the 1st order stencil
    foreach (cell; refs_unordered) {
        foreach (face; cell.iface) {
            if( face_ids.canFind(face.id) == false) {
                pcell.jacobian_face_stencil ~= face;
                face_ids ~= face.id;
            }
        }
    }
    // sort ids, and store sorted cell references
    cell_ids.sort();
    foreach(id; cell_ids) refs_ordered ~= refs_unordered[pos_array[id]];
    pcell.jacobian_cell_stencil ~= refs_ordered;
}
 
void residual_stencil(FVCell pcell, size_t orderOfJacobian) {

    if (orderOfJacobian == 0) {
        pcell.jacobian_cell_stencil ~= pcell; 
        foreach(face; pcell.iface) {
            pcell.jacobian_face_stencil ~= face;
        }
    }
    else if (orderOfJacobian == 1 && GlobalConfig.viscous == false) {
        FVCell[] refs_ordered; FVCell[] refs_unordered;
        size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
        size_t[] cell_ids;
        
        // clear the stencil arrays
        pcell.jacobian_cell_stencil = [];
        pcell.jacobian_face_stencil = [];

        size_t[] face_ids;
        // collect faces
        foreach(face; pcell.iface) {
            pcell.jacobian_face_stencil ~= face;
            face_ids ~= face.id;
        }
        
        // for each effected face, add the neighbouring cells
        foreach(face; pcell.jacobian_face_stencil) {
            // collect (non-ghost) neighbour cells
            if (cell_ids.canFind(face.left_cell.id) == false && face.left_cell.is_interior_to_domain) { // face.left_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.left_cell;
                pos_array[face.left_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.left_cell.id;
            }
            if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.is_interior_to_domain) { //face.right_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.right_cell;
                pos_array[face.right_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.right_cell.id;
            }
            else continue;
        }

        // finally collect the faces of the outer most cells for simulations that use source terms (i.e. axisymmetric)
        foreach(cell; refs_unordered) {
            foreach(face; cell.iface) {
                if (face_ids.canFind(face.id) == false) {
                    pcell.jacobian_face_stencil ~= face;
                    face_ids ~= face.id;
                }
             }
        }
        
        // sort ids, and store sorted cell references
        cell_ids.sort();
        foreach(id; cell_ids) refs_ordered ~= refs_unordered[pos_array[id]];
        pcell.jacobian_cell_stencil ~= refs_ordered;
    }
    else { // 2nd order || viscous || 2nd order + viscous
        FVCell[] refs_ordered; FVCell[] refs_unordered;
        size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
        size_t[] cell_ids; string[] face_ids;
        
        foreach(cell; pcell.cell_cloud) {
            // collect faces
            if(cell.id < ghost_cell_start_id) {
		foreach(face; cell.iface) {
		    if (face_ids.canFind(face.global_id) == false) {
			pcell.jacobian_face_stencil ~= face;
			face_ids ~= face.global_id;
		    }
		}
	    }
        }
        
        // for each effected face, add the neighbouring cells
        foreach(face; pcell.jacobian_face_stencil) {
            // collect (non-ghost) neighbour cells
            if (cell_ids.canFind(face.left_cell.id) == false && face.right_cell.is_interior_to_domain && face.left_cell.cell_cloud.length > 1) { //face.right_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.left_cell;
                pos_array[face.left_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.left_cell.id;
            }
            if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.is_interior_to_domain && face.right_cell.cell_cloud.length > 1) { //face.right_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.right_cell;
                pos_array[face.right_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.right_cell.id;
            }
            else continue;
        }

        // finally collect the faces of the outer most cells for simulations that use source terms (i.e. axisymmetric)
        foreach(cell; refs_unordered) {
	    if(cell.id < ghost_cell_start_id) {
		foreach(face; cell.iface) {
		    if (face_ids.canFind(face.global_id) == false) {
			pcell.jacobian_face_stencil ~= face;
			face_ids ~= face.global_id;
		    }
		}
	    }
	}
	
	/*
        // turbulent modelling requires a larger stencil - add more cells and faces....
        // for each effected face, add the neighbouring cells
        foreach(face; pcell.jacobian_face_stencil) {
            // collect (non-ghost) neighbour cells
            if (cell_ids.canFind(face.left_cell.id) == false && face.right_cell.is_interior_to_domain) { //face.right_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.left_cell;
                pos_array[face.left_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.left_cell.id;
            }
            if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.is_interior_to_domain) { //face.right_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.right_cell;
                pos_array[face.right_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.right_cell.id;
            }
            else continue;
        }

        // finally collect the faces of the outer most cells for simulations that use source terms (i.e. axisymmetric)
        foreach(cell; refs_unordered) {
            foreach(face; cell.iface) {
                 if (face_ids.canFind(face.id) == false) {
                    pcell.jacobian_face_stencil ~= face;
                    face_ids ~= face.id;
                }
            }
        }
        
        // for each effected face, add the neighbouring cells
        foreach(face; pcell.jacobian_face_stencil) {
            // collect (non-ghost) neighbour cells
            if (cell_ids.canFind(face.left_cell.id) == false && face.left_cell.is_interior_to_domain) { //face.left_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.left_cell;
                pos_array[face.left_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.left_cell.id;
            }
            if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.is_interior_to_domain) { //face.right_cell.id < ghost_cell_start_id) {
                refs_unordered ~= face.right_cell;
                pos_array[face.right_cell.id] = refs_unordered.length-1;
                cell_ids ~= face.right_cell.id;
            }
            else continue;
        }

        // finally collect the faces of the outer most cells for simulations that use source terms (i.e. axisymmetric)
        foreach(cell; refs_unordered) {
            foreach(face; cell.iface) {
                 if (face_ids.canFind(face.id) == false) {
                    pcell.jacobian_face_stencil ~= face;
                    face_ids ~= face.id;
                }
            }
        }
        */
        
        // finally sort ids, and store sorted cell references
        cell_ids.sort();
        foreach(id; cell_ids) {
            refs_ordered ~= refs_unordered[pos_array[id]];
        }
        pcell.jacobian_cell_stencil ~= refs_ordered;            
        //writeln(pcell.global_id, " ---- ");
        //foreach(cell; pcell.jacobian_cell_stencil) writeln(cell.global_id);
    }

}

void local_flow_jacobian_transpose(ref SMatrix!number A, ref FluidBlock blk, size_t np, size_t orderOfJacobian, number EPS, bool preconditionMatrix = false, bool transformToConserved = false) {
    //writeln("JACOBIAN START: ", blk.id);
    // set the interpolation order to that of the Jacobian
    if (orderOfJacobian < 2) blk.myConfig.interpolation_order = 1;
    else blk.myConfig.interpolation_order = 2;    
    // initialise re-used objects here to prevent memory bloat
    //writeln("PRIOR STENCIL");
    if (preconditionMatrix)
        foreach (cell; blk.cells) approximate_residual_stencil(cell, orderOfJacobian);
    else
        foreach (cell; blk.cells) residual_stencil(cell, orderOfJacobian);
    //writeln("STENCIL SET");
    number[][] aa; size_t[][] ja; size_t ia = 0;
    foreach(cell; blk.cells) {
	//writeln(cell.global_id, " ----");
	//foreach(c; cell.jacobian_cell_stencil) writef("%d    (%.12e    %.12e)    ", c.global_id, c.pos[0].x.re, c.pos[0].y.re);
	foreach(c; cell.jacobian_cell_stencil) //writef("%d    ", c.global_id);
	//writef("\n");
	foreach(f; cell.jacobian_face_stencil) //writef("%s    ", f.global_id);
	//writef("\n");

        aa.length = np; ja.length = np;
        compute_flow_jacobian_rows_for_cell(aa, ja, cell, blk, np, orderOfJacobian, EPS, transformToConserved);
        foreach (i; 0 .. np ) {
            A.aa ~= aa[i];
            A.ja ~= ja[i];
            A.ia ~= ia;
            ia += aa[i].length;
        }
        aa = [][];
        ja = [][];
    }
    A.ia ~= A.aa.length;
    foreach (cell; blk.cells) {
        cell.jacobian_face_stencil = [];
        cell.jacobian_cell_stencil = [];
    }
    apply_boundary_conditions(A, blk, np, orderOfJacobian, EPS, transformToConserved, preconditionMatrix);
    // reset the interpolation order
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;

    foreach (cell; blk.cells) {
        cell.jacobian_face_stencil = [];
        cell.jacobian_cell_stencil = [];
    }
}

void compute_flow_jacobian_rows_for_cell(number[][] aa, size_t[][] ja, FVCell pcell, FluidBlock blk, size_t np, size_t orderOfJacobian, number EPS, bool transformToConserved) {
    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;

    // 0th perturbation: rho
    mixin(computeFluxDerivativesAroundCell("gas.rho", "blk.MASS", true));
    // 1st perturbation: u
    mixin(computeFluxDerivativesAroundCell("vel.refx", "blk.X_MOM", false));
    // 2nd perturbation: v
    mixin(computeFluxDerivativesAroundCell("vel.refy", "blk.Y_MOM", false));
    if ( blk.myConfig.dimensions == 3 )
        mixin(computeFluxDerivativesAroundCell("vel.refz", "blk.Z_MOM", false));
    // 3rd perturbation: P
    mixin(computeFluxDerivativesAroundCell("gas.p", "blk.TOT_ENERGY", true));
    if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
        // 4th perturbation: tke
        mixin(computeFluxDerivativesAroundCell("tke", "blk.TKE", false));
        // 5th perturbation: omega
        mixin(computeFluxDerivativesAroundCell("omega", "blk.OMEGA", false));
    }
    // transform the face flux Jacobians from primitive to conservative form
    if (transformToConserved) {
        auto gmodel = blk.myConfig.gmodel;
        foreach (f; pcell.jacobian_face_stencil) {          
            // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
            number gamma = gmodel.gamma(pcell.fs.gas);
            // form inverse transformation matrix
            blk.Minv[blk.MASS,blk.MASS] = to!number(1.0);
            blk.Minv[blk.MASS,blk.X_MOM] = to!number(0.0);
            blk.Minv[blk.MASS,blk.Y_MOM] = to!number(0.0);
            blk.Minv[blk.MASS,blk.TOT_ENERGY] = to!number(0.0);
            // second row
            blk.Minv[blk.X_MOM,blk.MASS] = -pcell.fs.vel.x/pcell.fs.gas.rho;
            blk.Minv[blk.X_MOM,blk.X_MOM] = 1.0/pcell.fs.gas.rho;
            blk.Minv[blk.X_MOM,blk.Y_MOM] = to!number(0.0);
            blk.Minv[blk.X_MOM,blk.TOT_ENERGY] = to!number(0.0);
            // third row
            blk.Minv[blk.Y_MOM,blk.MASS] = -pcell.fs.vel.y/pcell.fs.gas.rho;
            blk.Minv[blk.Y_MOM,blk.X_MOM] = to!number(0.0);
            blk.Minv[blk.Y_MOM,blk.Y_MOM] = 1.0/pcell.fs.gas.rho;
            blk.Minv[blk.Y_MOM,blk.TOT_ENERGY] = to!number(0.0);
            // fourth row
            blk.Minv[blk.TOT_ENERGY,blk.MASS] = 0.5*(gamma-1.0)*(pcell.fs.vel.x*pcell.fs.vel.x+pcell.fs.vel.y*pcell.fs.vel.y+pcell.fs.vel.z*pcell.fs.vel.z);
            blk.Minv[blk.TOT_ENERGY,blk.X_MOM] = -pcell.fs.vel.x*(gamma-1);
            blk.Minv[blk.TOT_ENERGY,blk.Y_MOM] = -pcell.fs.vel.y*(gamma-1);
            blk.Minv[blk.TOT_ENERGY,blk.TOT_ENERGY] = gamma-1.0;

            if (blk.myConfig.dimensions == 3) {
                blk.Minv[blk.MASS,blk.Z_MOM] = to!number(0.0);
                blk.Minv[blk.X_MOM,blk.Z_MOM] = to!number(0.0);
                blk.Minv[blk.Y_MOM,blk.Z_MOM] = to!number(0.0);
                blk.Minv[blk.TOT_ENERGY,blk.Z_MOM] = -pcell.fs.vel.z*(gamma-1);
                
                blk.Minv[blk.Z_MOM,blk.MASS] = -pcell.fs.vel.z/pcell.fs.gas.rho;
                blk.Minv[blk.Z_MOM,blk.X_MOM] = to!number(0.0);
                blk.Minv[blk.Z_MOM,blk.Y_MOM] = to!number(0.0);
                blk.Minv[blk.Z_MOM,blk.Z_MOM] = 1.0/pcell.fs.gas.rho;
                blk.Minv[blk.Z_MOM,blk.TOT_ENERGY] = to!number(0.0);
            }
            
            if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                blk.Minv[blk.MASS,blk.TKE] = to!number(0.0);
                blk.Minv[blk.MASS,blk.OMEGA] = to!number(0.0);
                // second row
                blk.Minv[blk.X_MOM,blk.TKE] = to!number(0.0);
                blk.Minv[blk.X_MOM,blk.OMEGA] = to!number(0.0);
                // third row
                blk.Minv[blk.Y_MOM,blk.TKE] = to!number(0.0);
                blk.Minv[blk.Y_MOM,blk.OMEGA] = to!number(0.0);
                // fourth row
                blk.Minv[blk.TOT_ENERGY,blk.TKE] = -(gamma-1.0);
                blk.Minv[blk.TOT_ENERGY,blk.OMEGA] = to!number(0.0);
                // fifth row
                blk.Minv[blk.TKE,blk.MASS] = -pcell.fs.tke/pcell.fs.gas.rho;
                blk.Minv[blk.TKE,blk.X_MOM] = to!number(0.0);
                blk.Minv[blk.TKE,blk.Y_MOM] = to!number(0.0);
                blk.Minv[blk.TKE,blk.TOT_ENERGY] = to!number(0.0);
                blk.Minv[blk.TKE,blk.TKE] = 1.0/pcell.fs.gas.rho;
                blk.Minv[blk.TKE,blk.OMEGA] = to!number(0.0);
                // sixth row
                blk.Minv[blk.OMEGA,blk.MASS] = -pcell.fs.omega/pcell.fs.gas.rho;
                blk.Minv[blk.OMEGA,blk.X_MOM] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.Y_MOM] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.TOT_ENERGY] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.TKE] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.OMEGA] = 1.0/pcell.fs.gas.rho;

                if (blk.myConfig.dimensions == 3) {
                    blk.Minv[blk.Z_MOM,blk.TKE] = to!number(0.0);
                    blk.Minv[blk.Z_MOM,blk.OMEGA] = to!number(0.0);
                    blk.Minv[blk.TKE,blk.Z_MOM] = to!number(0.0);
                    blk.Minv[blk.OMEGA,blk.Z_MOM] = to!number(0.0);
                }
            }
            //writeln(blk.Minv);
            // writeln(pcell.id, ", ", blk.Minv);
            //if(pcell.id == 0) writeln(blk.Minv);
            //writeln(blk.Minv);
            number[][] tmp;
            tmp.length = np;
            foreach (ref a; tmp) a.length = np;
            for (size_t i = 0; i < np; i++) {
                for (size_t j = 0; j < np; j++) {
                    tmp[i][j] = to!number(0.0);
                    for (size_t k = 0; k < np; k++) {
                        tmp[i][j] += f.dFdU[i][k]*blk.Minv[k,j];
                    }
                }
            }

            foreach (i; 0..np) {
                foreach (j; 0..np) {
                    f.dFdU[i][j] = tmp[i][j];
                    blk.Minv[i,j] = to!number(0.0);
                }
            }
        }
        foreach (cell; pcell.jacobian_cell_stencil) {          
            // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
            number gamma = gmodel.gamma(pcell.fs.gas);
            // form inverse transformation matrix
            blk.Minv[blk.MASS,blk.MASS] = to!number(1.0);
            blk.Minv[blk.MASS,blk.X_MOM] = to!number(0.0);
            blk.Minv[blk.MASS,blk.Y_MOM] = to!number(0.0);
            blk.Minv[blk.MASS,blk.TOT_ENERGY] = to!number(0.0);
            // second row
            blk.Minv[blk.X_MOM,blk.MASS] = -pcell.fs.vel.x/pcell.fs.gas.rho;
            blk.Minv[blk.X_MOM,blk.X_MOM] = 1.0/pcell.fs.gas.rho;
            blk.Minv[blk.X_MOM,blk.Y_MOM] = to!number(0.0);
            blk.Minv[blk.X_MOM,blk.TOT_ENERGY] = to!number(0.0);
            // third row
            blk.Minv[blk.Y_MOM,blk.MASS] = -pcell.fs.vel.y/pcell.fs.gas.rho;
            blk.Minv[blk.Y_MOM,blk.X_MOM] = to!number(0.0);
            blk.Minv[blk.Y_MOM,blk.Y_MOM] = 1.0/pcell.fs.gas.rho;
            blk.Minv[blk.Y_MOM,blk.TOT_ENERGY] = to!number(0.0);
            // fourth row
            blk.Minv[blk.TOT_ENERGY,blk.MASS] = 0.5*(gamma-1.0)*(pcell.fs.vel.x*pcell.fs.vel.x+pcell.fs.vel.y*pcell.fs.vel.y+pcell.fs.vel.z*pcell.fs.vel.z);
            blk.Minv[blk.TOT_ENERGY,blk.X_MOM] = -pcell.fs.vel.x*(gamma-1);
            blk.Minv[blk.TOT_ENERGY,blk.Y_MOM] = -pcell.fs.vel.y*(gamma-1);
            blk.Minv[blk.TOT_ENERGY,blk.TOT_ENERGY] = gamma-1.0;

            if (blk.myConfig.dimensions == 3) {
                blk.Minv[blk.MASS,blk.Z_MOM] = to!number(0.0);
                blk.Minv[blk.X_MOM,blk.Z_MOM] = to!number(0.0);
                blk.Minv[blk.Y_MOM,blk.Z_MOM] = to!number(0.0);
                blk.Minv[blk.TOT_ENERGY,blk.Z_MOM] = -pcell.fs.vel.z*(gamma-1);
                
                blk.Minv[blk.Z_MOM,blk.MASS] = -pcell.fs.vel.z/pcell.fs.gas.rho;
                blk.Minv[blk.Z_MOM,blk.X_MOM] = to!number(0.0);
                blk.Minv[blk.Z_MOM,blk.Y_MOM] = to!number(0.0);
                blk.Minv[blk.Z_MOM,blk.Z_MOM] = 1.0/pcell.fs.gas.rho;
                blk.Minv[blk.Z_MOM,blk.TOT_ENERGY] = to!number(0.0);
            }
            
            if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                blk.Minv[blk.MASS,blk.TKE] = to!number(0.0);
                blk.Minv[blk.MASS,blk.OMEGA] = to!number(0.0);
                // second row
                blk.Minv[blk.X_MOM,blk.TKE] = to!number(0.0);
                blk.Minv[blk.X_MOM,blk.OMEGA] = to!number(0.0);
                // third row
                blk.Minv[blk.Y_MOM,blk.TKE] = to!number(0.0);
                blk.Minv[blk.Y_MOM,blk.OMEGA] = to!number(0.0);
                // fourth row
                blk.Minv[blk.TOT_ENERGY,blk.TKE] = -(gamma-1.0);
                blk.Minv[blk.TOT_ENERGY,blk.OMEGA] = to!number(0.0);
                // fifth row
                blk.Minv[blk.TKE,blk.MASS] = -pcell.fs.tke/pcell.fs.gas.rho;
                blk.Minv[blk.TKE,blk.X_MOM] = to!number(0.0);
                blk.Minv[blk.TKE,blk.Y_MOM] = to!number(0.0);
                blk.Minv[blk.TKE,blk.TOT_ENERGY] = to!number(0.0);
                blk.Minv[blk.TKE,blk.TKE] = 1.0/pcell.fs.gas.rho;
                blk.Minv[blk.TKE,blk.OMEGA] = to!number(0.0);
                // sixth row
                blk.Minv[blk.OMEGA,blk.MASS] = -pcell.fs.omega/pcell.fs.gas.rho;
                blk.Minv[blk.OMEGA,blk.X_MOM] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.Y_MOM] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.TOT_ENERGY] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.TKE] = to!number(0.0);
                blk.Minv[blk.OMEGA,blk.OMEGA] = 1.0/pcell.fs.gas.rho;

                if (blk.myConfig.dimensions == 3) {
                    blk.Minv[blk.Z_MOM,blk.TKE] = to!number(0.0);
                    blk.Minv[blk.Z_MOM,blk.OMEGA] = to!number(0.0);
                    blk.Minv[blk.TKE,blk.Z_MOM] = to!number(0.0);
                    blk.Minv[blk.OMEGA,blk.Z_MOM] = to!number(0.0);
                }
            }
            //writeln(blk.Minv);
            number[][] tmp;
            tmp.length = np;
            foreach (ref a; tmp) a.length = np;
            for (size_t i = 0; i < np; i++) {
                for (size_t j = 0; j < np; j++) {
                    tmp[i][j] = to!number(0.0);
                    for (size_t k = 0; k < np; k++) {
                        tmp[i][j] += cell.dQdU[i][k]*blk.Minv[k,j];
                    }
                }
            }
            foreach (i; 0..np) {
                foreach (j; 0..np) {
                    cell.dQdU[i][j] = tmp[i][j];
                    blk.Minv[i,j] = to!number(0.0);
                }
            }
        }
    }

    // compute Jacobian rows for perturbed cell
    foreach(cell; pcell.jacobian_cell_stencil) {
	if (cell.id < ghost_cell_start_id) {
	    size_t I, J; // indices in Jacobian matrix
	    number integral;
	    number volInv = 1.0 / cell.volume[0];
	    for ( size_t ip = 0; ip < np; ++ip ) {
		I = cell.id*np + ip; // row index
		for ( size_t jp = 0; jp < np; ++jp ) {
		    integral = 0.0;
		    J = jp; // column index

		    if(pcell.global_id == 6 && cell.global_id == 10) {
			//writeln(pcell.global_id, ", ", cell.global_id, ", ", pcell.pos[0].x.re, ", ", pcell.pos[0].y.re, ", ", cell.pos[0].x.re, ", ", cell.pos[0].y.re);
			foreach(fi, iface; cell.iface) {
			    //writef(" %.6e    %.6e   %d    %d    %.4e    %.12e    %.12e    %.12e \n", iface.pos.x.re, iface.pos.y.re, ip, jp, cell.outsign[fi].re, iface.area[0].re, iface.dFdU[ip][jp].re, volInv.re);
			}
		    }


		    foreach(fi, iface; cell.iface) {
			//if (cell.id == 0 && pcell.id == 1 && iface.id == 1) writeln("dFdU_f1: ", iface.dFdU);
			//if (cell.id == 90 && pcell.id == 91 && iface.id == 100) writeln("dFdU_f100: ", iface.dFdU);
			integral -= cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
		    }
		    number JacEntry = volInv * integral + cell.dQdU[ip][jp];
		    aa[J] ~= -JacEntry;
		    ja[J] ~= I;
		}
	    }
	}
    }
    // clear the interface flux Jacobian entries
    foreach (iface; pcell.jacobian_face_stencil) {
        foreach (i; 0..iface.dFdU.length) {
            foreach (j; 0..iface.dFdU[i].length) {
                iface.dFdU[i][j] = 0.0;
            }
        }
    }

    // clear the cell source term Jacobian entries
    foreach (cell; pcell.jacobian_cell_stencil) {
        foreach (i; 0..cell.dQdU.length) {
            foreach (j; 0..cell.dQdU[i].length) {
                cell.dQdU[i][j] = 0.0;
            }
        }
    }
}

string computeFluxDerivativesAroundCell(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "pcell.encode_conserved(0, 0, 0.0);";
    codeStr ~= "blk.cellSave.copy_values_from(pcell, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "pcell.fs."~varName~" += EPS;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(pcell.fs.gas);";
        codeStr ~= "blk.myConfig.gmodel.update_trans_coeffs(pcell.fs.gas);";
        codeStr ~= "blk.myConfig.gmodel.update_sound_speed(pcell.fs.gas);";
        //codeStr ~= "foreach(isp; 0 .. blk.myConfig.gmodel.n_species) pcell.fs.gas.massf[isp] = (pcell.U[0].massf[isp] * (1.0/pcell.fs.gas.rho));";
    }
    codeStr ~= "compute_flux(pcell, blk, orderOfJacobian, pcell.jacobian_cell_stencil, pcell.jacobian_face_stencil, blk.ifaceP);"; 
    //codeStr ~= "pcell.copy_values_from(cellSave, CopyDataOption.all);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "foreach (i, iface; pcell.jacobian_face_stencil) {";
    codeStr ~= "iface.dFdU[blk.MASS][" ~ posInArray ~ "] = blk.ifaceP[i].F.mass.im/EPS.im;";         
    codeStr ~= "iface.dFdU[blk.X_MOM][" ~ posInArray ~ "] = blk.ifaceP[i].F.momentum.x.im/EPS.im;";
    codeStr ~= "iface.dFdU[blk.Y_MOM][" ~ posInArray ~ "] = blk.ifaceP[i].F.momentum.y.im/EPS.im;";
    codeStr ~= "if (blk.myConfig.dimensions == 3) {";
    codeStr ~= "iface.dFdU[blk.Z_MOM][" ~ posInArray ~ "] = blk.ifaceP[i].F.momentum.z.im/EPS.im;";
    codeStr ~= "}";
    codeStr ~= "iface.dFdU[blk.TOT_ENERGY][" ~ posInArray ~ "] = blk.ifaceP[i].F.total_energy.im/EPS.im;";
    codeStr ~= "if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {";
    codeStr ~= "iface.dFdU[blk.TKE][" ~ posInArray ~ "] = blk.ifaceP[i].F.tke.im/EPS.im;";
    codeStr ~= "iface.dFdU[blk.OMEGA][" ~ posInArray ~ "] = blk.ifaceP[i].F.omega.im/EPS.im;";        
    codeStr ~= "}";
    codeStr ~= "}";
    codeStr ~= "foreach (i, cell; pcell.jacobian_cell_stencil) {";
    //codeStr ~= "writeln(cell.Q);";
    codeStr ~= "cell.dQdU[blk.MASS][" ~ posInArray ~ "] = cell.Q.mass.im/EPS.im;";         
    codeStr ~= "cell.dQdU[blk.X_MOM][" ~ posInArray ~ "] = cell.Q.momentum.x.im/EPS.im;";
    codeStr ~= "cell.dQdU[blk.Y_MOM][" ~ posInArray ~ "] = cell.Q.momentum.y.im/EPS.im;";
    codeStr ~= "if (blk.myConfig.dimensions == 3) {";
    codeStr ~= "cell.dQdU[blk.Z_MOM][" ~ posInArray ~ "] = cell.Q.momentum.z.im/EPS.im;";
    codeStr ~= "}";
    codeStr ~= "cell.dQdU[blk.TOT_ENERGY][" ~ posInArray ~ "] = cell.Q.total_energy.im/EPS.im;";
    codeStr ~= "if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {";
    codeStr ~= "cell.dQdU[blk.TKE][" ~ posInArray ~ "] = cell.Q.tke.im/EPS.im;";
    codeStr ~= "cell.dQdU[blk.OMEGA][" ~ posInArray ~ "] = cell.Q.omega.im/EPS.im;";        
    codeStr ~= "}";
    codeStr ~= "}";
    codeStr ~= "pcell.copy_values_from(blk.cellSave, CopyDataOption.all);";
    /*
    codeStr ~= "foreach(cell; pcell.jacobian_cell_stencil) {";
    codeStr ~= "cell.fs.mu_t = cell.fs.mu_t.re;";
    codeStr ~= "cell.fs.k_t = cell.fs.k_t.re;";
    codeStr ~= "cell.gradients.compute_lsq_values(cell.cell_cloud, cell.ws, blk.myConfig);";
    codeStr ~= "}";
    codeStr ~= "foreach(i, iface; pcell.jacobian_face_stencil) {";
    codeStr ~= "iface.fs.copy_average_values_from(iface.left_cell.fs, iface.right_cell.fs);";
    codeStr ~= "blk.ifaceP[i].fs.copy_average_values_from(iface.left_cell.fs, iface.right_cell.fs);";
    codeStr ~= "}";
    codeStr ~= "foreach(face; pcell.jacobian_face_stencil) {";
    codeStr ~= "face.grad.gradients_leastsq(face.cloud_fs, face.cloud_pos, face.ws_grad);";
    codeStr ~= "}";
    */
    codeStr ~= "compute_flux(pcell, blk, orderOfJacobian, pcell.jacobian_cell_stencil, pcell.jacobian_face_stencil, blk.ifaceP);"; 
    //codeStr ~= "steadystate_core.evalRHS(0.0, 0);";
    return codeStr;
}

void compute_flux(FVCell pcell, FluidBlock blk, size_t orderOfJacobian, ref FVCell[] cell_list, FVInterface[] iface_list, FVInterface[] ifaceP_list) {
    //writeln("COMPUTE FLUX");
    foreach(iface; iface_list) iface.F.clear();
    foreach(iface; ifaceP_list) iface.F.clear();
    foreach(cell; cell_list) cell.clear_source_vector();
    if (orderOfJacobian > 1) {
        // TODO: add in missing MLP code. 
        // compute gradients for reconstruction
        foreach(c; cell_list) {
	    //writeln(pcell.global_id, ", ", c.global_id, ", ", c.cell_cloud.length, ", ", blk.id);
            c.gradients.compute_lsq_values(c.cell_cloud, c.ws, blk.myConfig);
        }

	if (GlobalConfig.frozen_limiter == false) {
	    foreach(c; cell_list) {
		// It is more efficient to determine limiting factor here for some usg limiters.
		final switch (blk.myConfig.unstructured_limiter) {
		case UnstructuredLimiter.van_albada:
		    // do nothing now
		    break;
		case UnstructuredLimiter.min_mod:
		    // do nothing now
		    break;
		case UnstructuredLimiter.mlp:
		    c.gradients.mlp_limit(c.cell_cloud, c.ws, blk.myConfig);
		    break;
		case UnstructuredLimiter.barth:
		    c.gradients.barth_limit(c.cell_cloud, c.ws, blk.myConfig);
		    break;
		case UnstructuredLimiter.heuristic_van_albada:
		    c.gradients.heuristic_van_albada_limit(c.cell_cloud, c.ws, blk.myConfig, 0);
		    break;
		case UnstructuredLimiter.venkat:
		    c.gradients.venkat_limit(c.cell_cloud, c.ws, blk.myConfig, 0);
		    break;
		} // end switch
	    } // end foreach c
	}

        // Fill in gradients for ghost cells so that left- and right- cells at all faces,
        // including those along block boundaries, have the latest gradient values.
        // Note that we DO NOT copy gradients from neighbour blocks even if the boundary
        // has a mapped cell, this is due to efficiency, and block-parallel formation reasons.
        // So we just use the fall back for all boundary conditions.
        foreach (bcond; blk.bc) {
	if (bcond.type != "exchange_using_mapped_cells") {
	    // There are no other cells backing the ghost cells on this boundary.
		// Fill in ghost-cell gradients from the other side of the face.
		foreach (i, f; bcond.faces) {
		    if (bcond.outsigns[i] == 1) {
			f.right_cell.gradients.copy_values_from(f.left_cell.gradients);
		    } else {
			f.left_cell.gradients.copy_values_from(f.right_cell.gradients);
		    }
		} // end foreach f
	    } // end foreach bcond
	}

	foreach(c; cell_list) {
            c.gradients.compute_lsq_values(c.cell_cloud, c.ws, blk.myConfig);
	}

    } // end if interpolation_order > 1
    // Convective flux update
    //writeln("BEGIN INVISCID");
    foreach(iface; iface_list) {
        auto ublk = cast(UFluidBlock) blk;
        size_t gtl = 0;
        bool allow_high_order_interpolation = true;
        if (iface.left_cell && iface.right_cell) {
            ublk.lsq.interp_both(iface, gtl, ublk.Lft, ublk.Rght, allow_high_order_interpolation);
	    compute_interface_flux(ublk.Lft, ublk.Rght, iface, ublk.myConfig, ublk.omegaz);
        } else if (iface.right_cell) {
            ublk.lsq.interp_right(iface, gtl, ublk.Rght, allow_high_order_interpolation);
	    compute_flux_at_left_wall(ublk.Rght, iface, ublk.myConfig, ublk.omegaz);
        } else if (iface.left_cell) {
            ublk.lsq.interp_left(iface, gtl, ublk.Lft, allow_high_order_interpolation);
	    compute_flux_at_right_wall(ublk.Lft, iface, ublk.myConfig, ublk.omegaz);
        } else {
                assert(0, "oops, a face without attached cells");
        }
        
    }
    //writeln("INVISCID CHECKS OUT");
    /*
    foreach (f; iface_list) {
	if(pcell.global_id == 4 || pcell.global_id == 12) writeln("---- face check: ", pcell.global_id, ", ", f.global_id, ", ", f.pos.x.re, ", ", f.pos.y.re, ", ", blk.id, ", ", f.fs.vel.x, ", ", f.fs.vel.y, ", ", f.fs.gas.rho, ", ", f.fs.gas.p);
	//f.average_cell_deriv_values(0);
    }
    foreach(c; cell_list) {
        c.grad.gradients_leastsq(c.cloud_fs, c.cloud_pos, c.ws_grad); // blk.flow_property_spatial_derivatives(0); 
        if(pcell.global_id == 4 || pcell.global_id == 12) {
            writeln("fs check 0 : ", blk.id);
            foreach(i; 0..c.cloud_fs.length) {
                writeln(pcell.global_id, ", ", c.global_id, ", ", c.cloud_pos[i].x.re, ", ", c.cloud_pos[i].y.re, ", ", c.cloud_fs[i].vel.x, ", ", c.cloud_fs[i].vel.y, ", ", c.cloud_fs[i].gas.rho, ", ", c.cloud_fs[i].gas.p, ", ", c.cloud_fs[i].gas.T);
            }
            writef("grad check:  %d    %d    %.12e   %.12e   %.12e   %.12e   %.12e   %.12e\n", pcell.global_id, c.global_id, c.grad.vel[0][0], c.grad.vel[0][1], c.grad.vel[1][0], c.grad.vel[1][1], c.grad.T[0], c.grad.T[1]);
	    writeln(c.ws_grad.wx);
	    writeln(c.ws_grad.wy);
	}
    }
    */
    
    //writeln("END INVISCID");
    blk.applyPostConvFluxAction(0.0, 0, 0);
    // Viscous flux update
    //writeln("BEGIN VISCOUS");
    if (blk.myConfig.viscous) {
        // currently only for least-squares at faces
        blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
	//writeln("STAGE 1");
        foreach(c; cell_list) {
	    //writeln(c.global_id, " in");
	    //writeln("fs: ", c.cloud_fs);
	    //writeln("pos: ", c.cloud_pos);
	    //writeln("ws_grad: ", c.ws_grad);
	    c.grad.gradients_leastsq(c.cloud_fs, c.cloud_pos, c.ws_grad); // blk.flow_property_spatial_derivatives(0); 
	    //writeln(c.global_id, " out");
	}
	//writeln("SPATIAL GRAD CHECKS OUT");
	//writeln("STAGE 2");
        foreach (bcond; blk.bc) {
	    if (bcond.type != "exchange_using_mapped_cells") {
		// There are no other cells backing the ghost cells on this boundary.
		// Fill in ghost-cell gradients from the other side of the face.
		foreach (i, f; bcond.faces) {
		    if (bcond.outsigns[i] == 1) {
			f.right_cell.grad.copy_values_from(f.left_cell.grad);
		    } else {
			f.left_cell.grad.copy_values_from(f.right_cell.grad);
		    }
		} // end foreach f
	    }
	}

	//writeln("SPATIAL GRAD COPY CHECKS OUT");
        if(pcell.global_id == 17) //writef("pcell info: %d   %d    %d \n", blk.id, pcell.jacobian_cell_stencil.length, pcell.jacobian_face_stencil.length);
	foreach(c; cell_list) {
            c.grad.gradients_leastsq(c.cloud_fs, c.cloud_pos, c.ws_grad); // blk.flow_property_spatial_derivatives(0); 
	    if(pcell.global_id == 17) //writef("%d   %.12e   %.12e   %.12e   %.12e   %.12e   %.12e\n", c.global_id, c.grad.vel[0][0], c.grad.vel[0][1], c.grad.vel[1][0], c.grad.vel[1][1], c.grad.T[0], c.grad.T[1]);
            if(pcell.global_id == 17 && c.global_id == 17) {
                foreach(i; 0..c.cloud_fs.length) {
                    //writeln("c_cloud: ", c.cloud_pos[i].x.re, ", ", c.cloud_pos[i].y.re, ", ", c.cloud_fs[i].gas.rho, ", ", c.cloud_fs[i].vel.x, ", ", c.cloud_fs[i].vel.y, ", ", c.cloud_fs[i].gas.p);
                }
	    }
	}

	if(pcell.global_id == 17) {
            foreach(face; pcell.iface) {
                //writeln("face check 0: ", face.left_cell.global_id, ", ", face.left_cell.fs.gas.rho, ", ", face.right_cell.global_id, ", ", face.right_cell.fs.gas.rho);
            }
	}

	//writeln("SPATIAL GRAD 2 CHECKS OUT");
	//writeln("STAGE 3");
        foreach (f; iface_list) {
	    //writeln("in......blk: ", blk.id, ", f.id: ", f.global_id, ", ", f.left_cell, ", ", f.right_cell);
	    //if(pcell.global_id == 14) writeln("check: ", f.id, ", ", blk.id, ", ", f.F.mass);
            f.average_cell_deriv_values(0);
	    //writeln("out......blk: ", blk.id, ", f.id: ", f.global_id, ", ", f.left_cell, ", ", f.right_cell);
        }

        if(pcell.global_id == 17) {
            foreach(f; pcell.jacobian_face_stencil) {
                //writeln("face check 1: ", f.global_id, ", ", f.grad.vel[0][0], " ,", f.grad.vel[0][1], ", ", f.grad.vel[1][0], ", ", f.grad.vel[1][1], ", ", f.grad.T[0], ", ", f.grad.T[1]);
            }
        }

	//writeln("SPATIAL GRAD AVG CHECKS OUT");
	//writeln("STAGE 4");
        final switch (blk.myConfig.turbulence_model) {
        case TurbulenceModel.none:
            foreach (cell; cell_list) cell.turbulence_viscosity_zero();
            break;
        case TurbulenceModel.baldwin_lomax:
            throw new FlowSolverException("need to port baldwin_lomax_turbulence_model");
        case TurbulenceModel.spalart_allmaras:
            throw new FlowSolverException("Should implement Spalart-Allmaras some day.");
        case TurbulenceModel.k_omega:
            foreach (cell; cell_list) cell.turbulence_viscosity_k_omega();
            break;
        }

        foreach (cell; cell_list) {
            cell.turbulence_viscosity_factor(blk.myConfig.transient_mu_t_factor);
            cell.turbulence_viscosity_limit(blk.myConfig.max_mu_t_factor);
            cell.turbulence_viscosity_zero_if_not_in_zone();
        }

        if(pcell.global_id == 17) {
            foreach(face; pcell.jacobian_face_stencil) {
		//writeln("face check 2: ", face.global_id, ", ", face.fs.gas.rho, ", ", face.fs.vel.x, ", ", face.fs.vel.y, ", ", face.fs.gas.p, ", ", face.fs.gas.T);

            }
        }

        foreach(iface; iface_list) {
            iface.viscous_flux_calc();
        }

        if(pcell.global_id == 17) {
            foreach(face; pcell.jacobian_face_stencil) {
                //writeln("face check 3: ", face.global_id, ", ", face.F.mass, ", ", face.F.momentum.x, ", ", face.F.momentum.y, ", ", face.F.total_energy);
            }
        }
        
	//writeln("VISCOUS CALC CHECKS OUT");
        blk.applyPostDiffFluxAction(0.0, 0, 0);
    }
    bool local_with_k_omega = (blk.myConfig.turbulence_model == TurbulenceModel.k_omega);
    foreach (i, cell; cell_list) {
        cell.add_inviscid_source_vector(0, 0.0);
        if (blk.myConfig.viscous) {
            cell.add_viscous_source_vector(local_with_k_omega);
        }
        if (blk.myConfig.udf_source_terms) {
            size_t i_cell = cell.id;
            size_t j_cell = 0;
            size_t k_cell = 0;
            if (blk.grid_type == Grid_t.structured_grid) {
                auto sblk = cast(SFluidBlock) blk;
                assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                i_cell = ijk_indices[0];
                j_cell = ijk_indices[1];
                k_cell = ijk_indices[2];
            }
            addUDFSourceTermsToCell(blk.myL, cell, 0, 
                                    0.0, blk.myConfig.gmodel,
                                    blk.id, i_cell, j_cell, k_cell);
        }
    }

    foreach (f; iface_list) {
	//if(pcell.global_id == 14 && f.id == 9) writeln("check 1: ", f.id, ", ", f.global_id, ", ", f.pos.x.re, ", ", f.pos.y.re, ", ", blk.id, ", ", f.F.mass);
	//f.average_cell_deriv_values(0);
    }

    // copy perturbed flux
    foreach(i, iface; iface_list) {
        ifaceP_list[i].copy_values_from(iface, CopyDataOption.all);
    }
}


void compute_design_variable_partial_derivatives(Vector3[] design_variables, ref number[] g, size_t nPrimitive, bool with_k_omega, number EPS, string jobName) {
    size_t nDesignVars = design_variables.length;
    int gtl; int ftl; number objFcnEvalP; number objFcnEvalM; string varID; number dP; number P0;

    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;

    foreach (i; 0..nDesignVars) {
        foreach (myblk; localFluidBlocks) {
            ensure_directory_is_present(make_path_name!"grid"(0));
            string gridFileName = make_file_name!"grid"(jobName, myblk.id, 0, GlobalConfig.gridFileExt = "gz");
            myblk.read_new_underlying_grid(gridFileName);
            myblk.sync_vertices_from_underlying_grid(0);
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }

        /*
        foreach (myblk; localFluidBlocks) {
            myblk.read_solution(make_file_name!"flow"("ramp", myblk.id, 0, GlobalConfig.flowFileExt), false);
            
            // We can apply a special initialisation to the flow field, if requested.
            //if (GlobalConfig.diffuseWallBCsOnInit) {
            //writeln("Applying special initialisation to blocks: wall BCs being diffused into domain.");
            //writefln("%d passes of the near-wall flow averaging operation will be performed.", GlobalConfig.nInitPasses);
            //foreach (blk; parallel(localFluidBlocks,1)) {
            //    diffuseWallBCsIntoBlock(blk, GlobalConfig.nInitPasses, GlobalConfig.initTWall);
            //}
            //}
            
            foreach (cell; myblk.cells) {
                cell.encode_conserved(0, 0, myblk.omegaz);
                // Even though the following call appears redundant at this point,
                // fills in some gas properties such as Prandtl number that is
                // needed for both the cfd_check and the BaldwinLomax turbulence model.
                cell.decode_conserved(0, 0, myblk.omegaz);
            }
        }
        */
        //steadystate_core.evalRHS(0.0, 0);
        
        // perturb design variable +ve
        gtl = 1; ftl = 1;
        
        // perturb design variable in complex plan
        P0 = design_variables[i].y; 
        design_variables[i].refy = P0 + EPS;
        
        // perturb grid
        gridUpdate(design_variables, 1);

        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, vtx; myblk.vertices) {
                vtx.pos[0].refx = vtx.pos[1].x;
                vtx.pos[0].refy = vtx.pos[1].y;
            }
        }

        //exchange_ghost_cell_boundary_data(0.0, 0, 0);
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
            
            //foreach ( face; myblk.faces )
            //    foreach ( j; 0..face.cloud_pos.length) writef("%d    %.16f    %.16f \n", face.id, face.ws_grad.wx[j], face.ws_grad.wy[j]); 
        }

        exchange_ghost_cell_boundary_data(0.0, 0, 0);
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            //myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
            
            //foreach ( face; myblk.faces )
            //    foreach ( j; 0..face.cloud_pos.length) writef("%d    %.16f    %.16f \n", face.id, face.ws_grad.wx[j], face.ws_grad.wy[j]); 
        }
        
        //evalRHS(0.0, ftl, 0, with_k_omega);
        steadystate_core.evalRHS(0.0, ftl);
        
        objFcnEvalP = objective_function_evaluation(0);
        
        // compute cost function sensitivity
        g[i] = (objFcnEvalP.im)/(EPS.im);

        // compute residual sensitivity
        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, cell; myblk.cells) {
                myblk.rT[i, j*nPrimitive+myblk.MASS] = to!number((-cell.dUdt[ftl].mass.im)/(EPS.im));
                myblk.rT[i, j*nPrimitive+myblk.X_MOM] = to!number((-cell.dUdt[ftl].momentum.x.im)/(EPS.im));
                myblk.rT[i, j*nPrimitive+myblk.Y_MOM] = to!number((-cell.dUdt[ftl].momentum.y.im)/(EPS.im));
                if (myblk.myConfig.dimensions == 3) 
                    myblk.rT[i, j*nPrimitive+myblk.Z_MOM] = to!number((-cell.dUdt[ftl].momentum.z.im)/(EPS.im));
                myblk.rT[i, j*nPrimitive+myblk.TOT_ENERGY] = to!number((-cell.dUdt[ftl].total_energy.im)/(EPS.im));                
                
                if (myblk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    myblk.rT[i, j*nPrimitive+myblk.TKE] = to!number((-cell.dUdt[ftl].tke.im)/(EPS.im));
                    myblk.rT[i, j*nPrimitive+myblk.OMEGA] = to!number((-cell.dUdt[ftl].omega.im)/(EPS.im));
                }
            }
        }
        
        // restore design variable
        design_variables[i].refy = P0;
    }
}

/**************************/
/*  OBJECTIVE FUNCTIONS   */
/**************************/
number objective_function_evaluation(int gtl=0, string bndaryForSurfaceIntergral = "objective_function_surface") {

    number ObjFcn = 0.0;    
    foreach (myblk; parallel(localFluidBlocks,1)) {
        myblk.locObjFcn = 0.0;
        foreach (bndary; myblk.bc) {
            if (bndary.group == bndaryForSurfaceIntergral) {
                foreach (i, f; bndary.faces) {
                    FVCell cell;
                    // collect interior cell
                    if (bndary.outsigns[i] == 1) {
                        cell = f.left_cell;
                    } else {
                        cell = f.right_cell;
                    }
                    myblk.locObjFcn += cell.fs.gas.p*f.area[gtl];//*f.n.x;
                }
            }
        }
    }
    foreach ( myblk; localFluidBlocks) ObjFcn += myblk.locObjFcn;
    // writef("drag: %.16f    constraint: %.16f    scale: %.16f \n", 2.0*PI*abs(ObjFcn), ConsFcn, scale);
    return fabs(ObjFcn);
}

void form_objective_function_sensitivity(FluidBlock blk, size_t np, number EPS, string bndaryForSurfaceIntergral = "objective_function_surface") {

    // for now we have hard coded the pressure drag in the x-direction as the objective function
    size_t nLocalCells = blk.cells.length;
    blk.f.length = nLocalCells * np;

    foreach(cell; blk.cells) {
        for ( size_t ip = 0; ip < np; ++ip ) {
            blk.f[cell.id*np + ip] = 0.0;
        }
    }
    
    foreach (bndary; blk.bc) {
        if (bndary.group == bndaryForSurfaceIntergral) {            
	    foreach (i, f; bndary.faces) {
		FVCell cell; number origValue; number ObjFcnM; number ObjFcnP; number h;
                // collect interior cell
                if (bndary.outsigns[i] == 1) {
		    cell = f.left_cell;
		} else {
		    cell = f.right_cell;
		}
                // for current objective function only perturbations in pressure have any effect
                origValue = cell.fs.gas.p;
                cell.fs.gas.p = origValue + EPS;
                ObjFcnP = objective_function_evaluation();
                //writeln(ObjFcnP);
                blk.f[cell.id*np + blk.TOT_ENERGY] = (ObjFcnP.im)/(EPS.im);
                cell.fs.gas.p = origValue;
            }
        }
    }
}

/**********************************/
/*  GRID PERTURBATION FUNCTIONs   */
/**********************************/
void fit_design_parameters_to_surface(ref Vector3[] designVars)
{
    // fitting tolerances
    double tol = GlobalConfig.sscOptions.tolBezierCurveFit;
    int maxSteps = GlobalConfig.sscOptions.maxStepsBezierCurveFit;
    int nCntrlPts;
    
    // collect vertices along design surface (may cross multiple blocks)
    Vector3[] orderedList; Vector3[] unorderedList;
    size_t[string] origPosId; // used to identify where the point is in the unordered list
    number[] xPosition;
    foreach ( blk; localFluidBlocks ) {
        size_t[] idList;
        foreach ( bndary; blk.bc ) {
            if (bndary.is_design_surface) {
                nCntrlPts = bndary.num_cntrl_pts;
                foreach ( face; bndary.faces) {
                    foreach ( vtx; face.vtx) {
                        // check x-position uniqueness
                        bool uniqueXPos = true ;
                        foreach ( i; 0..unorderedList.length) {
                            number diff = abs(vtx.pos[0].x - xPosition[i]);
                            if ( diff < ESSENTIALLY_ZERO) uniqueXPos = false;
                        }
                        if (uniqueXPos) {                   
                            unorderedList ~= Vector3(vtx.pos[0].x, vtx.pos[0].y, vtx.pos[0].z);
                            xPosition ~= vtx.pos[0].x;
                            string xPosIdx = to!string(vtx.pos[0].x);
                            origPosId[xPosIdx] = unorderedList.length-1;
                        }
                    }
                }
            }
        }
    }

    // order points in ascending x-position (WARNING: it is assumed that, for a given design surface boundary, each x-coordinate is unique).
    xPosition.sort();
    foreach(x; xPosition) orderedList ~= unorderedList[origPosId[to!string(x)]];

    double[] ts;
    Bezier bezier = optimiseBezierPoints(orderedList, nCntrlPts, ts, tol, maxSteps);
    // first and last control points are not design variables
    foreach ( i; 1..bezier.B.length-1) {
        designVars ~= bezier.B[i];
    }

    // transmit global bezier to all relevant blocks
    // copy bezier curve object to each design surface boundary, along with relevant portion of the ts array
    foreach (myblk; parallel(localFluidBlocks,1)) {
        foreach (bndary; myblk.bc) {
            if (bndary.is_design_surface) {
                bndary.bezier = bezier;
                number xi = bndary.faces[0].vtx[0].pos[0].x;
                number xf = bndary.faces[0].vtx[0].pos[0].x;
                foreach ( face; bndary.faces ) {
                    foreach ( vtx; face.vtx ) {
                        if (vtx.pos[0].x < xi) xi = vtx.pos[0].x;
                        else if (vtx.pos[0].x > xf) xf = vtx.pos[0].x;
                        else continue;
                    }
                }
                size_t idxi; size_t idxf;
                foreach (i, t; ts ) {
                    if ( fabs(orderedList[i].x - xi) < ESSENTIALLY_ZERO) idxi = i;
                    else if ( fabs(orderedList[i].x - xf) < ESSENTIALLY_ZERO) idxf = i;
                    else {} //do nothing
                } // foreach (ts)
                bndary.ts.length = ts[idxi..idxf+1].length;
                bndary.ts[] = ts[idxi..idxf+1];
            } // end if
        } // end foreach blk.bc
    } // end foreach blk
    writeBezierDataToFile();
} // end parameterise_design_surfaces

void gridUpdate(Vector3[] designVars, size_t gtl, bool gridUpdate = false, string jobName = "") {
    size_t nDesignVars = designVars.length;
    
    foreach (myblk; parallel(localFluidBlocks,1)) {
        foreach(j, vtx; myblk.vertices) {
            vtx.pos[gtl].refx = vtx.pos[0].x;
            vtx.pos[gtl].refy = vtx.pos[0].y;
        }
    }

    Vector3[] bndaryVtxInitPos;
    foreach (myblk; localFluidBlocks) {
        size_t[] idList;
        foreach(bndary; myblk.bc) {
            foreach( face; bndary.faces) {
                foreach ( vtx; face.vtx) {
                    if (idList.canFind(vtx.id) == false && bndary.type != "exchange_using_mapped_cells") {            
                        myblk.boundaryVtxIndexList ~= vtx.id;
                        //
                        bool update = true;
                        foreach(pos; bndaryVtxInitPos) {
                            if ( fabs(pos.x.re - vtx.pos[0].x.re) < 1.0e-50 && fabs(pos.y.re - vtx.pos[0].y.re) < 1.0e-50 ) {
                                update = false;
                            }
                        }
                        //
                        if(update) {
                            bndaryVtxInitPos ~= vtx.pos[0];
                            //myblk.boundaryVtxIndexList ~= vtx.id;
                            idList ~= vtx.id;
                        }
                    }
                }
            }
        }
    }
    //writeln("length: ", bndaryVtxInitPos.length);
    foreach (myblk; localFluidBlocks) {
        foreach( bndary; myblk.bc ) {
                if (bndary.is_design_surface) {
                    foreach ( i; 1..bndary.bezier.B.length-1) {
                        // y-variable
                        bndary.bezier.B[i].refy = designVars[i-1].y;
                    }
                    
                    foreach(j, vtx; bndary.vertices) {
                        if (gridUpdate) {
                            vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                            vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                        } else {
                            version(complex_numbers) vtx.pos[gtl].refx = complex(vtx.pos[gtl].x.re, bndary.bezier(bndary.ts[j]).x.im);
                            version(complex_numbers) vtx.pos[gtl].refy = complex(vtx.pos[gtl].y.re, bndary.bezier(bndary.ts[j]).y.im);
                        }
                    }
                }
        }
    }
    
    Vector3[] bndaryVtxNewPos;
    foreach (myblk; localFluidBlocks) {
        size_t[] idList;
        foreach(bndary; myblk.bc) {
            foreach( face; bndary.faces) {
                foreach ( vtx; face.vtx) {
                    if (idList.canFind(vtx.id) == false  && bndary.type != "exchange_using_mapped_cells") {

                        //
                        bool update = true;
                        foreach(pos; bndaryVtxNewPos) {
                            if ( fabs(pos.x.re - vtx.pos[gtl].x.re) < 1.0e-50 && fabs(pos.y.re - vtx.pos[gtl].y.re) < 1.0e-50 ) {
                                update = false;
                            }
                        }
                        //
                        if(update) {
                            bndaryVtxNewPos ~= vtx.pos[gtl];
                            idList ~= vtx.id;
                        }
                    }
                }
            }
        }
    }
    
    foreach (myblk; localFluidBlocks) {
        inverse_distance_weighting(myblk, bndaryVtxInitPos, bndaryVtxNewPos, gtl);
    }

    foreach (myblk; localFluidBlocks) {
        myblk.boundaryVtxIndexList = [];
    }

    if (gridUpdate) {
        foreach (myblk; localFluidBlocks) {
            foreach(j, vtx; myblk.vertices) {
                vtx.pos[0].refx = vtx.pos[gtl].x;
                vtx.pos[0].refy = vtx.pos[gtl].y;
            }
            
            // write out grid
            // save mesh
            myblk.sync_vertices_to_underlying_grid(0);
            ensure_directory_is_present(make_path_name!"grid"(0));
            auto fileName = make_file_name!"grid"(jobName, myblk.id, 0, GlobalConfig.gridFileExt = "gz");
            myblk.write_underlying_grid(fileName);
        }
    }
}

void collect_boundary_vertices(FluidBlock blk)
{
    // make a block local collection of the vertices along domain boundaries
    foreach (bndary; blk.bc) {
        if ( bndary.type != "exchange_using_mapped_cells") {
            if (bndary.is_design_surface) { // we need to order the vertices by x-position
                FVVertex[] vtxOrdered; FVVertex[] vtxUnordered;
                size_t[string] origPosId; // used to identify where the point is in the unordered list
                number[] xPosition; size_t[] listOfAddedVtxIds;
                foreach(face; bndary.faces) {
                    foreach(vtx; face.vtx) {
                        if (!listOfAddedVtxIds.canFind(vtx.id)) {
                            blk.boundaryVtxIndexList ~= vtx.id;
                            vtxUnordered ~= vtx;
                            xPosition ~= vtx.pos[0].x;
                            listOfAddedVtxIds ~= vtx.id;
                            string xPosIdx = to!string(vtx.pos[0].x);
                            origPosId[xPosIdx] = vtxUnordered.length-1;
                        } // end if
                    } // foreach vtx
                } // end foreach face
                
                // order points in ascending x-position (WARNING: it is assumed that, for a given design surface boundary, each x-coordinate is unique).
                xPosition.sort();
                foreach(x; xPosition) vtxOrdered ~= vtxUnordered[origPosId[to!string(x)]];
                bndary.vertices ~= vtxOrdered;
            } else { // we don't need the vertices in any particular order
                size_t[] listOfAddedVtxIds;
                foreach(face; bndary.faces) {
                    foreach(vtx; face.vtx) {
                        if (!listOfAddedVtxIds.canFind(vtx.id)) {
                            bndary.vertices ~= vtx;
                            listOfAddedVtxIds ~= vtx.id;
                            blk.boundaryVtxIndexList ~= vtx.id;
                        } // end if
                    } // end foreach vtx
                } // end foreach face
            } // end else
        } // end if
    } // end foreach bndary
} // end collect_boundary_vertices


/*************************/
/*  EVALUATE RHS @ gtl   */
/*************************/
/*
void evalRHS(double pseudoSimTime, int ftl, int gtl, bool with_k_omega)
{
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
    }
    
    exchange_ghost_cell_boundary_data(pseudoSimTime, gtl, ftl);
    
    foreach (blk; localFluidBlocks) {
        blk.applyPreReconAction(pseudoSimTime, gtl, ftl);
    }
    
    // We don't want to switch between flux calculator application while
    // doing the Frechet derivative, so we'll only search for shock points
    // at ftl = 0, which is when the F(U) evaluation is made.
    if ( ftl == 0 && (GlobalConfig.flux_calculator == FluxCalculator.adaptive_hanel_ausmdv ||
		      GlobalConfig.flux_calculator == FluxCalculator.adaptive_hlle_roe) ||
		      GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.detect_shock_points();
        }
    }

    bool allow_high_order_interpolation = true;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase0(allow_high_order_interpolation, gtl);
    }
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase1(allow_high_order_interpolation, gtl);
    }
    foreach (blk; localFluidBlocks) {
        blk.applyPostConvFluxAction(pseudoSimTime, gtl, ftl);
    }
    if (GlobalConfig.viscous) {
        foreach (blk; localFluidBlocks) {
            blk.applyPreSpatialDerivActionAtBndryFaces(pseudoSimTime, gtl, ftl);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.flow_property_spatial_derivatives(gtl); 
            blk.estimate_turbulence_viscosity();
            blk.viscous_flux();
        }
        foreach (blk; localFluidBlocks) {
            blk.applyPostDiffFluxAction(pseudoSimTime, gtl, ftl);
        }
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        foreach (i, cell; blk.cells) {
            cell.add_inviscid_source_vector(gtl, 0.0);
            if (blk.myConfig.viscous) {
                cell.add_viscous_source_vector(local_with_k_omega);
            }
            if (blk.myConfig.udf_source_terms) {
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                if (blk.grid_type == Grid_t.structured_grid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                    i_cell = ijk_indices[0];
                    j_cell = ijk_indices[1];
                    k_cell = ijk_indices[2];
                }
                addUDFSourceTermsToCell(blk.myL, cell, gtl, 
                                        pseudoSimTime, blk.myConfig.gmodel,
                                        blk.id, i_cell, j_cell, k_cell);
            }
            cell.time_derivatives(gtl, ftl, local_with_k_omega);
        }
    }
}
*/

/**********************/
/*  GMRES FUNCTIONS   */
/**********************/
string dot_over_blocks(string dot, string A, string B)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.dotAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.dotAcc += blk.`~A~`[k].re*blk.`~B~`[k].re;
   }
}
`~dot~` = 0.0;
foreach (blk; localFluidBlocks) `~dot~` += blk.dotAcc;`;

}

string norm2_over_blocks(string norm2, string blkMember)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.normAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.normAcc += blk.`~blkMember~`[k].re*blk.`~blkMember~`[k].re;
   }
}
`~norm2~` = 0.0;
foreach (blk; localFluidBlocks) `~norm2~` += blk.normAcc;
`~norm2~` = sqrt(`~norm2~`);`;
}

void rpcGMRES_solve0(size_t nPrimitive) {    

    auto fileName = "e4ssc.diagnostics.dat";
    auto outFile = File(fileName, "w");
    outFile.close();

    // restarted-GMRES settings
    size_t maxIters = GlobalConfig.sscOptions.maxOuterIterations; // maxOuterIters
    size_t m = maxIters;
    number outerTol = GlobalConfig.sscOptions.stopOnRelativeGlobalResidual;
    size_t maxRestarts = GlobalConfig.sscOptions.maxRestarts;
    size_t iterCount;
    number resid;
    size_t nRestarts;
    size_t r;
    // allocate GMRES arrays attached to the block objectcs
    foreach (blk; localFluidBlocks) {
        size_t n = nPrimitive*blk.cells.length;
        blk.nvars = n;
        // Now allocate arrays and matrices
        blk.psi.length = n;
        blk.r0.length = n;
        blk.x0.length = n;
        blk.v.length = n;
        blk.w.length = n;
        blk.wext.length = n;
        blk.z.length = n;
        blk.V = new Matrix!number(n, m+1);
        blk.Q1 = new Matrix!number(m+1, m+1);
        blk.g0.length = m+1;
        blk.g1.length = m+1;
    }    

    // allocate global GMRES arrays
    g0.length = m+1;
    g1.length = m+1;
    h.length = m+1;
    hR.length = m+1;
    H0 = new Matrix!number(m+1, m);
    H1 = new Matrix!number(m+1, m);
    Gamma = new Matrix!number(m+1, m+1);
    Q0 = new Matrix!number(m+1, m+1);
    Q1 = new Matrix!number(m+1, m+1);
    
    // Initialise some global arrays and matrices that have already been allocated
    g0[] = to!number(0.0);
    g1[] = to!number(0.0);
    H0.zeros();
    H1.zeros();

    number[] Z; // global array used in the matrix-vector product
    
    // 1. Evaluate r0, beta, v1
    // r0 = b - A*x0
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.x0[] = to!number(1.0); 
    }

    // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
    // 1. exchange interface data
    // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
    Z = [];                        
    foreach (blk; localFluidBlocks) {
        Z ~= blk.x0[];
    }
    foreach (blk; localFluidBlocks) {
        blk.Z.length = Z.length;
        blk.Z[] = Z[];
    }
    // 2. local product
    foreach (blk; parallel(localFluidBlocks,1)) {
        multiply(blk.JlocT, blk.x0, blk.r0);
    }
    // 3. external product
    foreach (blk; parallel(localFluidBlocks,1)) {
         multiply(blk.JextT, blk.Z, blk.wext);
    }
    // 4. sum the two contributions
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.r0[] += blk.wext[];
    }
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.f[k] - blk.r0[k];}
    }
    
    // Then compute v = r0/||r0||
    number betaTmp;
    mixin(norm2_over_blocks("betaTmp", "r0"));
    number beta = betaTmp;
    number residRef = beta;
    g0[0] = beta;
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) {
            blk.v[k] = blk.r0[k]/beta;
            blk.V[k,0] = blk.v[k];
        }
    }
    
    // Compute tolerance
    //auto outerTol = eta*beta;
    
    // 2. Start outer-loop of restarted GMRES
    int step = 0;
    for ( r = 0; r < maxRestarts; r++ ) {
        // 2a. Begin iterations
        foreach (j; 0 .. m) {
            iterCount = j+1;

            // compute z (preconditioning step);
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.z[] = blk.v[];
                solve(blk.P, blk.z);
            }

            // compute w
            // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
            // 1. exchange interface data
            // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
            Z = [];                        
            foreach (blk; localFluidBlocks) {
                Z ~= blk.z[];
            }
            foreach (blk; localFluidBlocks) {
                blk.Z.length = Z.length;
                blk.Z[] = Z[];
            }
            // 2. local product
            foreach (blk; parallel(localFluidBlocks,1)) {
                multiply(blk.JlocT, blk.z, blk.w);
            }
            // 3. external product
            foreach (blk; parallel(localFluidBlocks,1)) {
                multiply(blk.JextT, blk.Z, blk.wext);
            }
            // 4. sum the two contributions
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.w[] += blk.wext[];
            }
            
            // The remainder of the algorithm looks a lot like any standard
            // GMRES implementation (for example, see smla.d)
            foreach (i; 0 .. j+1) {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    // Extract column 'i'
                    foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.V[k,i]; 
                }
                number H0_ij_tmp;
                mixin(dot_over_blocks("H0_ij_tmp", "w", "v"));
                number H0_ij = H0_ij_tmp;
                H0[i,j] = H0_ij;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k]; 
                }
            }
            number H0_jp1j_tmp;
            mixin(norm2_over_blocks("H0_jp1j_tmp", "w"));
            number H0_jp1j = H0_jp1j_tmp;
            H0[j+1,j] = H0_jp1j;
        
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) {
                    blk.v[k] = blk.w[k]/H0_jp1j;
                    blk.V[k,j+1] = blk.v[k];
                }
            }

            // Build rotated Hessenberg progressively
            if ( j != 0 ) {
                // Extract final column in H
                foreach (i; 0 .. j+1) h[i] = H0[i,j];
                // Rotate column by previous rotations (stored in Q0)
                nm.bbla.dot!number(Q0, j+1, j+1, h, hR);
                // Place column back in H
                foreach (i; 0 .. j+1) H0[i,j] = hR[i];
            }
            // Now form new Gamma
            Gamma.eye();
            auto denom = sqrt(H0[j,j]*H0[j,j] + H0[j+1,j]*H0[j+1,j]);
            auto s_j = H0[j+1,j]/denom; 
            auto c_j = H0[j,j]/denom;
            Gamma[j,j] = c_j; Gamma[j,j+1] = s_j;
            Gamma[j+1,j] = -s_j; Gamma[j+1,j+1] = c_j;
            // Apply rotations
            nm.bbla.dot!number(Gamma, j+2, j+2, H0, j+1, H1);
            nm.bbla.dot!number(Gamma, j+2, j+2, g0, g1);
            // Accumulate Gamma rotations in Q.
            if ( j == 0 ) {
                copy(Gamma, Q1);
            }
            else {
                nm.bbla.dot!number(Gamma, j+2, j+2, Q0, j+2, Q1);
            }
            // Prepare for next step
            copy(H1, H0);
            g0[] = g1[];
            copy(Q1, Q0);
            // Get residual
            resid = fabs(g1[j+1]);
            // DEBUG:
            //      writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
            if ( resid <= outerTol*residRef ) {
                m = j+1;
                // DEBUG:
                //      writefln("OUTER: TOL ACHIEVED restart-count= %d iteration-count= %d, resid= %e", r, m, resid);
                break;
            }
        }
        
        if (iterCount == maxIters)
            m = maxIters;
        // At end H := R up to row m
        //        g := gm up to row m
        nm.bbla.upperSolve!number(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot!number(blk.V, blk.nvars, m, blk.g1, blk.psi);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.smla.solve(blk.P, blk.psi);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) blk.psi[k] += blk.x0[k];
        }
        writef("step: %d    relative residual: %.4e    absolute residual: %.4e \n", step, (resid/residRef).re, (resid).re);
        outFile = File(fileName, "a");
        outFile.writef("%d %.16e %.16e \n", step, resid/residRef, resid);
        outFile.close();
        step += 1;
        if ( resid <= outerTol*residRef || r+1 == maxRestarts ) {
            // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
            // DEBUG:  writefln("Breaking restart loop.");
            break;
        }
        // Else, we prepare for restart by setting x0 and computing r0
        // Computation of r0 as per Fraysee etal (2005)
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.x0[] = blk.psi[];
        }
        /*
        foreach (blk; localFluidBlocks) copy(Q1, blk.Q1);
        // Set all values in g0 to 0.0 except for final (m+1) value
        foreach (i; 0 .. m) g0[i] = 0.0;
        foreach (blk; localFluidBlocks) blk.g0[] = g0[];
        
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.Q1, m, m+1, blk.g0, blk.g1);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.V, blk.nvars, m+1, blk.g1, blk.r0);
        }
        */

        // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
        // 1. exchange interface data
        // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
        Z = [];                        
        foreach (blk; localFluidBlocks) {
            Z ~= blk.x0[];
        }
        foreach (blk; localFluidBlocks) {
            blk.Z.length = Z.length;
            blk.Z[] = Z[];
        }
        // 2. local product
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.smla.multiply(blk.JlocT, blk.x0, blk.r0);
        }
        // 3. external product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.JextT, blk.Z, blk.wext);
        }
        // 4. sum the two contributions
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.r0[] += blk.wext[];
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.f[k] - blk.r0[k];}
        }
        mixin(norm2_over_blocks("betaTmp", "r0"));
        beta = betaTmp;
        // DEBUG: writefln("OUTER: ON RESTART beta= %e", beta);
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.r0[k]/beta;
                blk.V[k,0] = blk.v[k];
            }
        }
        // Re-initialise some vectors and matrices for restart
        g0[] = to!number(0.0);
        g1[] = to!number(0.0);
        H0.zeros();
        H1.zeros();
        // And set first residual entry
        g0[0] = beta;

    }
    nRestarts = to!int(r);
    writeln(nRestarts, " restarts.");
}

void rpcGMRES_solve1(size_t nPrimitive) {    

    auto fileName = "e4ssc.diagnostics.dat";
    auto outFile = File(fileName, "w");
    outFile.close();
    
    // restarted-GMRES settings
    size_t maxIters = GlobalConfig.sscOptions.maxOuterIterations;
    size_t m = maxIters;
    number outerTol = GlobalConfig.sscOptions.stopOnRelativeGlobalResidual;
    size_t maxRestarts = GlobalConfig.sscOptions.maxRestarts;
    double eta = GlobalConfig.sscOptions.eta;
    size_t iterCount;
    number resid;
    size_t nRestarts;
    size_t r;
    double CFL = GlobalConfig.sscOptions.cfl0;
    double dt = steadystate_core.determine_initial_dt(CFL);
    
    // allocate GMRES arrays attached to the block objectcs
    foreach (blk; localFluidBlocks) {
        size_t n = nPrimitive*blk.cells.length;
        blk.nvars = n;
        // Now allocate arrays and matrices
        blk.psi.length = n;
        blk.r0.length = n;
        blk.x0.length = n;
        blk.b.length = n;
        blk.delpsi.length = n;
        blk.v.length = n;
        blk.w.length = n;
        blk.wext.length = n;
        blk.z.length = n;
        blk.V = new Matrix!number(n, m+1);
        blk.Q1 = new Matrix!number(m+1, m+1);
        blk.g0.length = m+1;
        blk.g1.length = m+1;
    }    

    // allocate global GMRES arrays
    g0.length = m+1;
    g1.length = m+1;
    h.length = m+1;
    hR.length = m+1;
    H0 = new Matrix!number(m+1, m);
    H1 = new Matrix!number(m+1, m);
    Gamma = new Matrix!number(m+1, m+1);
    Q0 = new Matrix!number(m+1, m+1);
    Q1 = new Matrix!number(m+1, m+1);
    
    // Initialise some global arrays and matrices that have already been allocated
    g0[] = to!number(0.0);
    g1[] = to!number(0.0);
    H0.zeros();
    H1.zeros();

    number[] Z; // global array used in the matrix-vector product

    // set initial guess for the adjoint variables (either 0 or 1 is popular in the literature)
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.x0[] = to!number(0.0); // this is delpsi guess
        blk.delpsi[] = to!number(0.0);
        blk.psi[] = to!number(1.0); // this is psi guess
    }

    // compute original adjoint system residual (resid = b - Ax0)
    number residualRef = 0.0;
    number residual = 0.0;
    // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
    // 1. exchange interface data
    // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
    Z = [];                        
    foreach (blk; localFluidBlocks) {
        Z ~= blk.psi[];
    }
    foreach (blk; localFluidBlocks) {
        blk.Z.length = Z.length;
        blk.Z[] = Z[];
    }
    // 2. local product
    foreach (blk; parallel(localFluidBlocks,1)) {
        multiply(blk.JlocT, blk.psi, blk.r0);
    }
    // 3. external product
    foreach (blk; parallel(localFluidBlocks,1)) {
         multiply(blk.JextT, blk.Z, blk.wext);
    }
    // 4. sum the two contributions
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.r0[] += blk.wext[];
    }
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.f[k] - blk.r0[k];}
    }

    mixin(norm2_over_blocks("residualRef", "r0"));

    residual = residualRef;
    writef("reference residual: %.13e \n", residualRef.re);
    //start outer loop
    int step = 0;
    while (residual/residualRef > outerTol) {
        
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.x0[] = to!number(0.0); // this is delpsi guess
        }

        // 1. Evaluate r0, beta, v1
        // where r0 = b* - A*x0 (* refers to pseudo-time augumented system)

        // compute b = -dRdQ*psi_n + f
        // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
        // 1. exchange interface data
        // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
        Z = [];                        
        foreach (blk; localFluidBlocks) {
            Z ~= blk.psi[];
        }
        foreach (blk; localFluidBlocks) {
            blk.Z.length = Z.length;
            blk.Z[] = Z[];
        }
        // 2. local product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.JlocT, blk.psi, blk.r0);
        }
        // 3. external product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.JextT, blk.Z, blk.wext);
        }
        // 4. sum the two contributions
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.r0[] += blk.wext[];
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) { blk.b[k] = blk.f[k] - blk.r0[k];}
        }

        // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
        // 1. exchange interface data
        // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
        Z = [];                        
        foreach (blk; localFluidBlocks) {
            Z ~= blk.x0[];
        }
        foreach (blk; localFluidBlocks) {
            blk.Z.length = Z.length;
            blk.Z[] = Z[];
        }
        // 2. local product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.A, blk.x0, blk.r0);
        }
        // 3. external product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.Aext, blk.Z, blk.wext);
        }
        // 4. sum the two contributions
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.r0[] += blk.wext[];
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.b[k] - blk.r0[k];}
        }
        
        // Then compute v = r0/||r0||
        number betaTmp;
        mixin(norm2_over_blocks("betaTmp", "r0"));
        number beta = betaTmp;
        //writef("beta: %.13e \n", beta.re);
        g0[0] = beta;
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.r0[k]/beta;
                blk.V[k,0] = blk.v[k];
            }
        }

        // Compute tolerance
        //auto outerTol = eta*beta;
        
        // 2. Start inner-loop of restarted GMRES
        for ( r = 0; r < maxRestarts; r++ ) {
            // 2a. Begin iterations
            foreach (j; 0 .. m) {
                iterCount = j+1;
                
                // compute z (preconditioning step);
                foreach (blk; parallel(localFluidBlocks,1)) {
                    blk.z[] = blk.v[];
                    solve(blk.P, blk.z);
                }

                // compute w
                // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
                // 1. exchange interface data
                // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
                Z = [];                        
                foreach (blk; localFluidBlocks) {
                    Z ~= blk.z[];
                }
                foreach (blk; localFluidBlocks) {
                    blk.Z.length = Z.length;
                    blk.Z[] = Z[];
                }
                // 2. local product
                foreach (blk; parallel(localFluidBlocks,1)) {
                    multiply(blk.A, blk.z, blk.w);
                }
                // 3. external product
                foreach (blk; parallel(localFluidBlocks,1)) {
                    multiply(blk.Aext, blk.Z, blk.wext);
                }
                // 4. sum the two contributions
                foreach (blk; parallel(localFluidBlocks,1)) {
                    blk.w[] += blk.wext[];
                }
                //number wtmp;
                //mixin(norm2_over_blocks("wtmp", "w"));
                //writef("w: %.13e \n", wtmp.re);

                // The remainder of the algorithm looks a lot like any standard
                // GMRES implementation (for example, see smla.d)
                foreach (i; 0 .. j+1) {
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        // Extract column 'i'
                        foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.V[k,i]; 
                    }
                    number H0_ij_tmp;
                    mixin(dot_over_blocks("H0_ij_tmp", "w", "v"));
                    number H0_ij = H0_ij_tmp;
                    H0[i,j] = H0_ij;
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k]; 
                    }
                }
                number H0_jp1j_tmp;
                mixin(norm2_over_blocks("H0_jp1j_tmp", "w"));
                number H0_jp1j = H0_jp1j_tmp;
                H0[j+1,j] = H0_jp1j;
                
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (k; 0 .. blk.nvars) {
                        blk.v[k] = blk.w[k]/H0_jp1j;
                        blk.V[k,j+1] = blk.v[k];
                    }
                }
                
                // Build rotated Hessenberg progressively
                if ( j != 0 ) {
                    // Extract final column in H
                    foreach (i; 0 .. j+1) h[i] = H0[i,j];
                    // Rotate column by previous rotations (stored in Q0)
                    nm.bbla.dot!number(Q0, j+1, j+1, h, hR);
                    // Place column back in H
                    foreach (i; 0 .. j+1) H0[i,j] = hR[i];
                }
                // Now form new Gamma
                Gamma.eye();
                auto denom = sqrt(H0[j,j]*H0[j,j] + H0[j+1,j]*H0[j+1,j]);
                auto s_j = H0[j+1,j]/denom; 
                auto c_j = H0[j,j]/denom;
                Gamma[j,j] = c_j; Gamma[j,j+1] = s_j;
                Gamma[j+1,j] = -s_j; Gamma[j+1,j+1] = c_j;
                // Apply rotations
                nm.bbla.dot!number(Gamma, j+2, j+2, H0, j+1, H1);
                nm.bbla.dot!number(Gamma, j+2, j+2, g0, g1);
                // Accumulate Gamma rotations in Q.
                if ( j == 0 ) {
                    copy(Gamma, Q1);
                }
                else {
                    nm.bbla.dot!number(Gamma, j+2, j+2, Q0, j+2, Q1);
                }
                // Prepare for next step
                copy(H1, H0);
                g0[] = g1[];
                copy(Q1, Q0);
                // Get residual
                resid = fabs(g1[j+1]);
                // DEBUG:
                //      writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
                if ( resid <= beta*eta ) { // inner tolerance of 0.01
                    m = j+1;
                    // DEBUG:
                    //      writefln("OUTER: TOL ACHIEVED restart-count= %d iteration-count= %d, resid= %e", r, m, resid);
                    break;
                }
            }
            
            if (iterCount == maxIters)
                m = maxIters;
            // At end H := R up to row m
            //        g := gm up to row m
            nm.bbla.upperSolve!number(H1, to!int(m), g1);
            // In serial, distribute a copy of g1 to each block
            foreach (blk; localFluidBlocks) blk.g1[] = g1[];
            foreach (blk; parallel(localFluidBlocks,1)) {
                nm.bbla.dot!number(blk.V, blk.nvars, m, blk.g1, blk.delpsi);
            }
            foreach (blk; parallel(localFluidBlocks,1)) {
                nm.smla.solve(blk.P, blk.delpsi);
            }
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) blk.delpsi[k] += blk.x0[k];
            }
            //writef("global residual: %.16e \n",  resid);
            if ( resid <= beta*eta || r+1 == maxRestarts ) {
                // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
                // DEBUG:  writefln("Breaking restart loop.");
                break;
            }
            //writeln("RESTARTING");
            // Else, we prepare for restart by setting x0 and computing r0
            // Computation of r0 as per Fraysee etal (2005)
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.x0[] = blk.delpsi[];
            }
            /*
              foreach (blk; localFluidBlocks) copy(Q1, blk.Q1);
              // Set all values in g0 to 0.0 except for final (m+1) value
              foreach (i; 0 .. m) g0[i] = 0.0;
              foreach (blk; localFluidBlocks) blk.g0[] = g0[];
              
              foreach (blk; parallel(localFluidBlocks,1)) {
              nm.bbla.dot(blk.Q1, m, m+1, blk.g0, blk.g1);
              }
              foreach (blk; parallel(localFluidBlocks,1)) {
              nm.bbla.dot(blk.V, blk.nvars, m+1, blk.g1, blk.r0);
              }

            */

            // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
            // 1. exchange interface data
            // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
            Z = [];                        
            foreach (blk; localFluidBlocks) {
                Z ~= blk.psi[];
            }
            foreach (blk; localFluidBlocks) {
                blk.Z.length = Z.length;
                blk.Z[] = Z[];
            }
            // 2. local product
            foreach (blk; parallel(localFluidBlocks,1)) {
                multiply(blk.JlocT, blk.psi, blk.r0);
            }
            // 3. external product
            foreach (blk; parallel(localFluidBlocks,1)) {
                multiply(blk.JextT, blk.Z, blk.wext);
            }
            // 4. sum the two contributions
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.r0[] += blk.wext[];
            }
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) { blk.b[k] = blk.f[k] - blk.r0[k];}
            }
            
            // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
            // 1. exchange interface data
            // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
            Z = [];                        
            foreach (blk; localFluidBlocks) {
                Z ~= blk.x0[];
            }
            foreach (blk; localFluidBlocks) {
                blk.Z.length = Z.length;
                blk.Z[] = Z[];
            }
            // 2. local product
            foreach (blk; parallel(localFluidBlocks,1)) {
                multiply(blk.A, blk.x0, blk.r0);
            }
            // 3. external product
            foreach (blk; parallel(localFluidBlocks,1)) {
                multiply(blk.Aext, blk.Z, blk.wext);
            }
            // 4. sum the two contributions
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.r0[] += blk.wext[];
            }
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.b[k] - blk.r0[k];}
            }

            mixin(norm2_over_blocks("betaTmp", "r0"));
            beta = betaTmp;
            // DEBUG: writefln("OUTER: ON RESTART beta= %e", beta);
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) {
                    blk.v[k] = blk.r0[k]/beta;
                    blk.V[k,0] = blk.v[k];
                }
            }
            // Re-initialise some vectors and matrices for restart
            g0[] = to!number(0.0);
            g1[] = to!number(0.0);
            H0.zeros();
            H1.zeros();
            // And set first residual entry
            g0[0] = beta;
            
        }
        // should have a decent estimate for delpsi at this point, so update psi
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) { blk.psi[k] = blk.psi[k] + blk.delpsi[k];}
        }
        
        // parallel matrix-vector product; Saad, Krylov Subspace Methods in Distributed Computing Environments
        // 1. exchange interface data
        // let's take a shortcut for now by grabbing the global z array, this won't work for MPI
        Z = [];                        
        foreach (blk; localFluidBlocks) {
            Z ~= blk.psi[];
        }
        foreach (blk; localFluidBlocks) {
            blk.Z.length = Z.length;
            blk.Z[] = Z[];
        }
        // 2. local product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.JlocT, blk.psi, blk.r0);
        }
        // 3. external product
        foreach (blk; parallel(localFluidBlocks,1)) {
            multiply(blk.JextT, blk.Z, blk.wext);
        }
        // 4. sum the two contributions
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.r0[] += blk.wext[];
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.f[k] - blk.r0[k];}
        }
        number residual_tmp = 0.0;
        mixin(norm2_over_blocks("residual_tmp", "r0"));
        writef("step: %d    CFL: %.4e    relative residual: %.4e    absolute residual: %.4e \n", step, CFL, (residual_tmp/residualRef).re, (residual_tmp).re);
        outFile = File(fileName, "a");
        outFile.writef("%d %.16e %.16e $.4e\n", step, residual_tmp/residualRef, residual, CFL);
        outFile.close();
        step += 1;

        CFL = CFL;
        residual = residual_tmp;
    }

    nRestarts = to!int(r);
    writeln("COMPLETE");
}

/**********************/
/*    IO FUNCTIONS    */
/**********************/
void writeBezierDataToFile()
{
    foreach ( myblk; localFluidBlocks) {
        foreach (bndary; myblk.bc) {
            if (bndary.is_design_surface) {
                string fileName = "blk" ~ to!string(myblk.id) ~ ".bezier";
                if (exists(fileName)) {
                    string error_msg = format(".bezier files already exist. Please remove before proceeding.");
                    throw new FlowSolverException(error_msg);
                }
                auto outFile = File(fileName, "a");
                foreach ( point; bndary.bezier.B) {
                    outFile.writef("%.16e %.16e %.16e \n", point.x.re, point.y.re, point.z.re);
                }
                foreach( t; bndary.ts ) {
                    outFile.writef("%.16e \n", t);
                }
            } 
        } 
    } 
} // end writeBexierDataToFile

void readBezierDataFromFile(ref Vector3[] designVars)
{
    foreach (myblk; parallel(localFluidBlocks,1)) {
        collect_boundary_vertices(myblk);
    }

    foreach ( myblk; localFluidBlocks) {
        foreach (bndary; myblk.bc) {
            if (bndary.is_design_surface) {
                string fileName = "blk" ~ to!string(myblk.id) ~ ".bezier";
                if (!exists(fileName)) {
                    string error_msg = format(".bezier file does not exist");
                    throw new FlowSolverException(error_msg);
                }
                auto fR = File(fileName, "r");
                //while (!fR.eof) {
                Vector3[] bezPts;
                foreach ( i; 0..bndary.num_cntrl_pts) {
                    auto line = fR.readln().strip();
                    auto tokens = line.split();
                    Vector3 pt;
                    pt.refx = to!number(tokens[0]); pt.refy = to!number(tokens[1]); pt.refz = to!number(tokens[2]);
                    bezPts ~= pt;
                }
                bndary.bezier = new Bezier(bezPts);
                //writeln(myblk.id, ", ", bndary.vertices.length, ", ", bndary.which_boundary, ", ", bndary.group, ", ", bndary.type);
                foreach ( i; 0..bndary.vertices.length) {
                    auto line = fR.readln().strip();
                    auto tokens = line.split();
                    bndary.ts ~= to!double(tokens[0]);
                    //}
                }
                if (designVars.length < 1) {
                    foreach ( i; 1..bndary.bezier.B.length-1) {
                        Vector3 dvar;
                        dvar.refx = bndary.bezier.B[i].x;
                        dvar.refy = bndary.bezier.B[i].y;
                        designVars ~= dvar;
                    }
                }
            } // end if
        } // end foreach bndary
    } // end foreach myblk
} // end readBezierDataFromFile

void writeDesignVarsToDakotaFile(Vector3[] design_variables, string jobName) {
    size_t ndvars = design_variables.length;
    // write subset of bezier control points out to DAKOTA input script
    string dakota_input_filename = jobName ~ ".in";
    string pathToDakotaInputFile = "../" ~ dakota_input_filename;
    if (exists(pathToDakotaInputFile)) {
        string error_msg = format("DAKOTA input file already exists. Please remove the .in file before proceeding");
        throw new FlowSolverException(error_msg);
    }
    string dakota_input_tmplt_filename = jobName ~ ".tmplt";
    string pathToDakotaInputTmpltFile = "../" ~ dakota_input_tmplt_filename;
    if (!exists(pathToDakotaInputTmpltFile)) {
        string error_msg = format("DAKOTA input template file does not exist.");
        throw new FlowSolverException(error_msg);
    }
    
    // generate design variable information string
    
    string designVariableInputContent;
    designVariableInputContent ~= "    continuous_design = " ~ to!string(ndvars) ~ "\n";
    designVariableInputContent ~= "    initial_point    ";
    foreach ( i; 0..design_variables.length) designVariableInputContent ~= format!"%.16e    %.16e    "(design_variables[i].x.re, design_variables[i].y.re);
    designVariableInputContent ~= " \n"; 
    designVariableInputContent ~= "    descriptors    ";
    foreach ( i; 0..design_variables.length) {
                string descriptor;
                // y-variable
                descriptor = "'" ~ "design_var" ~ "_x" ~ to!string(i) ~ "'" ~ "    " ~ "'" ~ "design_var" ~ "_y" ~ to!string(i) ~ "'";
                designVariableInputContent ~= descriptor ~ "    ";
    }

    designVariableInputContent ~= " \n"; 
    
    auto fR = File(pathToDakotaInputTmpltFile, "r");
    auto fW = File(pathToDakotaInputFile, "w");
    
    while (!fR.eof) {
        auto line = fR.readln().strip();
        if (line == "variables") {
            fW.writeln(line);
            fW.writef(designVariableInputContent);
        }
        else fW.writeln(line);
    }
}

void readDesignVarsFromDakotaFile(ref Vector3[] design_variables)
{
    // read in new control points (note we should have pre-filled the design variables array with the original points)
    auto f = File("params.in", "r");
    auto line = f.readln().strip;// read first line & do nothing 
    auto tokens = line.split();
    foreach ( i; 0..design_variables.length) {
	// x-variable
        line = f.readln().strip;
        tokens = line.split();
        design_variables[i].refx = to!double(tokens[0]);
        // y-variable
        line = f.readln().strip;
        tokens = line.split();
        design_variables[i].refy = to!double(tokens[0]);
    }
    // assign design variables to bezier curve
    foreach ( myblk; localFluidBlocks) {
	foreach ( bndary; myblk.bc) {
	    if ( bndary.is_design_surface) {
		foreach ( i; 1..bndary.bezier.B.length-1) {
		    bndary.bezier.B[i].refx = design_variables[i-1].x;
		    bndary.bezier.B[i].refy = design_variables[i-1].y;
		} // end foreach i
	    } // end if
	} // end foreach bndary
    } // end foreach myblk
} // end readDesignVarsFromDakotaFile

void write_adjoint_variables_to_file(FluidBlock blk, size_t np, string jobName) {
    // Make a stack-local copy of conserved quantities info
    //size_t nConserved = nConservedQuantities;
    //size_t MASS = massIdx;
    //size_t X_MOM = xMomIdx;
    //size_t Y_MOM = yMomIdx;
    //size_t Z_MOM = zMomIdx;
    //size_t TOT_ENERGY = totEnergyIdx;
    //size_t TKE = tkeIdx;
    //size_t OMEGA = omegaIdx;

    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;
    // write out adjoint variables in VTK-format
    if (blk.grid_type == Grid_t.structured_grid) throw new FlowSolverException("Shape Sensitivity Calculator not implemented for structured grids, yet."); 
    if (blk.grid_type == Grid_t.unstructured_grid) {
        auto fileName = "adjointVars" ~ to!string(blk.id) ~ ".vtk";
        auto ublk = cast(UFluidBlock) blk; 
        auto outFile = File(fileName, "w");
        outFile.writef("# vtk DataFile Version 3.0 \n");
        outFile.writef("%s \n", jobName);
        outFile.writef("ASCII \n");
        outFile.writef("DATASET UNSTRUCTURED_GRID \n");
        outFile.writef("POINTS %d double \n", nvertices);
        // write grid data
        foreach(i, vtx; blk.vertices) {
            outFile.writef("%.16f %.16f %.16f \n", vtx.pos[0].x.re, vtx.pos[0].y.re, vtx.pos[0].z.re); 
        }
        // write cell connectivity
        size_t connections = 0;
        foreach ( cell; blk.cells) {
            connections += cell.iface.length;
        }
        size_t size = ncells + connections; // TODO: only for quads, need to generalise

        outFile.writef("CELLS %d %d \n", ncells, size);
        foreach(i, cell; ublk.grid.cells) {
            outFile.writef("%d ", cell.vtx_id_list.length);
            foreach(vid; cell.vtx_id_list) {
                outFile.writef("%d ", vid);
            }
            outFile.writef("\n");
        }
        outFile.writef("CELL_TYPES %d \n", ncells);
        foreach(i, cell; ublk.grid.cells) {
            outFile.writef("%d \n", ublk.grid.vtk_element_type_for_cell(i)); //cell.cell_type);
        }
        
        // write cell data
        outFile.writef("CELL_DATA %d \n", ncells);

        outFile.writef("SCALARS adjoint_density double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+blk.MASS].re);
        }

        outFile.writef("SCALARS adjoint_velx double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+blk.X_MOM].re);
        }

        outFile.writef("SCALARS adjoint_vely double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+blk.Y_MOM].re);
        }

        if (blk.myConfig.dimensions == 3) {
            outFile.writef("SCALARS adjoint_vely double \n");
            outFile.writef("LOOKUP_TABLE default \n");
            foreach(i; 0..ncells) {
                outFile.writef("%.16f \n", blk.psi[np*i+blk.Z_MOM].re);
            }
        }
        
        outFile.writef("SCALARS adjoint_pressure double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) { 
            outFile.writef("%.16f \n", blk.psi[np*i+blk.TOT_ENERGY].re);
        }
 
        if ( blk.myConfig.turbulence_model == TurbulenceModel.k_omega ) {
            outFile.writef("SCALARS adjoint_pressure double \n");
            outFile.writef("LOOKUP_TABLE default \n");
            foreach(i; 0..ncells) { 
                outFile.writef("%.16f \n", blk.psi[np*i+blk.TKE].re);
            }
            
            outFile.writef("SCALARS adjoint_pressure double \n");
            outFile.writef("LOOKUP_TABLE default \n");
            foreach(i; 0..ncells) { 
                outFile.writef("%.16f \n", blk.psi[np*i+blk.OMEGA].re);
            }
        }
    }
}

void write_objective_fn_to_file(string fileName, number objFnEval) {
    auto outFile = File(fileName, "w");
    outFile.writef("%.16e f\n", objFnEval.re); 
}

void write_gradients_to_file(string fileName, number[] grad) {
    auto outFile = File(fileName, "a");
    outFile.writef("[ ");
    foreach( i; 0..grad.length ) {
        outFile.writef("%.16e %.16e ", 0.0, grad[i].re);
    }
    outFile.writef(" ]\n"); 
}


/*****************************/
/*  DIRECT GRADIENT METHOD   */
/*****************************/
void compute_direct_complex_step_derivatives(string jobName, int last_tindx, int maxCPUs, Vector3[] design_variables, number EPS) {
    writeln(" ");
    writeln("------------------------------------------------------");
    writeln("----EVALUATING DERIVATIVES VIA DIRECT COMPLEX STEP----");
    writeln("------------------------------------------------------");
    writeln(" ");
    
    size_t nDesignVars = design_variables.length;
    double[] gradients; number P0; number objFcnP; number objFcnM; 

    foreach ( i; 0..nDesignVars) {
        writeln("----- Computing Gradient for variable: ", i);
        ensure_directory_is_present(make_path_name!"grid"(0));
        foreach (myblk; localFluidBlocks) {
            string gridFileName = make_file_name!"grid"(jobName, myblk.id, 0, GlobalConfig.gridFileExt = "gz");
            myblk.read_new_underlying_grid(gridFileName);
            myblk.sync_vertices_from_underlying_grid(0);
        }

        exchange_ghost_cell_boundary_data(0.0, 0, 0);
        
        foreach (myblk; localFluidBlocks) {
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }
        
        
        foreach (myblk; localFluidBlocks) {
            myblk.read_solution(make_file_name!"flow"(jobName, myblk.id, 0, GlobalConfig.flowFileExt), false);
            foreach (cell; myblk.cells) {
                cell.encode_conserved(0, 0, myblk.omegaz);
                // Even though the following call appears redundant at this point,
                // fills in some gas properties such as Prandtl number that is
                // needed for both the cfd_check and the BaldwinLomax turbulence model.
                cell.decode_conserved(0, 0, myblk.omegaz);
            }
        }

        // perturb design variable in complex plane
        P0 = design_variables[i].y; 
        design_variables[i].refy = P0 + EPS;
        
        // perturb grid
        gridUpdate(design_variables, 1); // gtl = 1

        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, vtx; myblk.vertices) {
                vtx.pos[0].refx = vtx.pos[1].x;
                vtx.pos[0].refy = vtx.pos[1].y;
            }
        }

        exchange_ghost_cell_boundary_data(0.0, 0, 0);
                
        foreach (myblk; localFluidBlocks) {
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }

        exchange_ghost_cell_boundary_data(0.0, 0, 0);

        foreach (myblk; localFluidBlocks) {
            // save mesh
            myblk.sync_vertices_to_underlying_grid(0);
            ensure_directory_is_present(make_path_name!"grid-p"(0));
            auto fileName = make_file_name!"grid-p"(jobName, myblk.id, 0, GlobalConfig.gridFileExt = "gz");
            myblk.write_underlying_grid(fileName);
        }

	if (GlobalConfig.diffuseWallBCsOnInit) {
	    //        if (GlobalConfig.viscous) {
            // We can apply a special initialisation to the flow field, if requested.
            //if (GlobalConfig.diffuseWallBCsOnInit) {
            //writeln("Applying special initialisation to blocks: wall BCs being diffused into domain.");
            //writefln("%d passes of the near-wall flow averaging operation will be performed.", GlobalConfig.nInitPasses);
            foreach (blk; parallel(localFluidBlocks,1)) {
                diffuseWallBCsIntoBlock(blk, GlobalConfig.nInitPasses, GlobalConfig.initTWall);
            }
        }
        
        // Additional memory allocation specific to steady-state solver
        allocate_global_workspace();
        foreach (myblk; localFluidBlocks) {
            myblk.allocate_GMRES_workspace();
        }

        foreach (myblk; localFluidBlocks) {
            //myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }
	
	if (GlobalConfig.sscOptions.read_frozen_limiter_values_from_file) {
	    steadystate_core.evalRHS(0.0, 0); // fill in some values
	    auto fileName = "frozen_limiter_values.dat";
	    auto outFile = File(fileName, "r");
	    foreach (blk; localFluidBlocks) {
		foreach (cell; blk.cells) {
		    auto line = outFile.readln().strip();
		    auto token = line.split();
		    cell.gradients.rhoPhi = to!double(token[0]);
		    
		    line = outFile.readln().strip();
		    token = line.split();
		    cell.gradients.velxPhi = to!double(token[0]);
		    
		    line = outFile.readln().strip();
		    token = line.split();
		    cell.gradients.velyPhi = to!double(token[0]);
		    
		    if (blk.myConfig.dimensions == 3) {
			line = outFile.readln().strip();
			token = line.split();
			cell.gradients.velzPhi = to!double(token[0]);
		    }
		    
		    line = outFile.readln().strip();
		    token = line.split();
		    cell.gradients.pPhi = to!double(token[0]);
		    
		    if (blk.myConfig.turbulence_model == TurbulenceModel.k_omega) {
			line = outFile.readln().strip();
			token = line.split();
			cell.gradients.tkePhi = to!double(token[0]);
			
			line = outFile.readln().strip();
			token = line.split();
			cell.gradients.omegaPhi = to!double(token[0]);
		    }
		}
	    }
	    outFile.close();
	    GlobalConfig.frozen_limiter = true;
	}
	

        // run steady-state solver
        iterate_to_steady_state(0, maxCPUs); // snapshotStart = 0
        //GlobalConfig.report_residuals = true;
        //sim_time = 0.0;
        //integrate_in_time(GlobalConfig.max_time);
        
        // compute objective function gradient
        objFcnP = objective_function_evaluation();
        
        // return value to original state
        design_variables[i].refy = P0;
        
        // compute objective function gradient
        objFcnM = objective_function_evaluation();
        gradients ~= (objFcnP.im)/(EPS.im);
        
        // return value to original state
        design_variables[i].refy = P0;
    }
    foreach ( i; 0..nDesignVars) {
        writef("gradient for variable %d: %.16e \n", i, gradients[i]);
    }
    writeln("simulation complete.");
}

/*
void compute_direct_complex_step_derivatives(string jobName, int last_tindx, int maxCPUs, Vector3[] design_variables, number EPS) {
    writeln(" ");
    writeln("------------------------------------------------------");
    writeln("----EVALUATING DERIVATIVES VIA DIRECT COMPLEX STEP----");
    writeln("------------------------------------------------------");
    writeln(" ");
        
    size_t nDesignVars = design_variables.length;
    double[] gradients; number P0; number objFcnP; number objFcnM; 

    foreach ( i; 0..nDesignVars) {
        writeln("----- Computing Gradient for variable: ", i);
        foreach (myblk; localFluidBlocks) {
            ensure_directory_is_present(make_path_name!"grid"(0));
            string gridFileName = make_file_name!"grid"(jobName, myblk.id, 0, gridFileExt = "gz");
            myblk.read_new_underlying_grid(gridFileName);
            myblk.sync_vertices_from_underlying_grid(0);
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }
    
        foreach (myblk; localFluidBlocks) {
            myblk.read_solution(make_file_name!"flow"(jobName, myblk.id, 0, flowFileExt), false);

            foreach (cell; myblk.cells) {
                cell.encode_conserved(0, 0, myblk.omegaz);
                // Even though the following call appears redundant at this point,
                // fills in some gas properties such as Prandtl number that is
                // needed for both the cfd_check and the BaldwinLomax turbulence model.
                cell.decode_conserved(0, 0, myblk.omegaz);
            }
        }
        
        // perturb design variable in complex plane
        P0 = design_variables[i].y; 
        design_variables[i].refy = P0 + EPS.im;
        
        // perturb grid
        gridUpdate(design_variables, 1); // gtl = 1

        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, vtx; myblk.vertices) {
                vtx.pos[0].refx = vtx.pos[1].x;
                vtx.pos[0].refy = vtx.pos[1].y;
            }
        }
        
        foreach (myblk; localFluidBlocks) {
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }

        foreach (myblk; localFluidBlocks) {
            // save mesh
            myblk.sync_vertices_to_underlying_grid(0);
            ensure_directory_is_present(make_path_name!"grid-p"(0));
            auto fileName = make_file_name!"grid-p"(jobName, myblk.id, 0, gridFileExt = "gz");
            myblk.write_underlying_grid(fileName);
        }
            
        // Additional memory allocation specific to steady-state solver
        allocate_global_workspace();
        foreach (myblk; localFluidBlocks) {
            myblk.allocate_GMRES_workspace();
        }
        
        // run steady-state solver
        iterate_to_steady_state(0, maxCPUs); // snapshotStart = 0
        //GlobalConfig.report_residuals = true;
        //sim_time = 0.0;
        //integrate_in_time(GlobalConfig.max_time);
        
        // compute objective function gradient
        objFcnP = objective_function_evaluation();
                
        // return value to original state
        design_variables[i].refy = P0;


        foreach (myblk; localFluidBlocks) {
            ensure_directory_is_present(make_path_name!"grid"(0));
            string gridFileName = make_file_name!"grid"(jobName, myblk.id, 0, gridFileExt = "gz");
            myblk.read_new_underlying_grid(gridFileName);
            myblk.sync_vertices_from_underlying_grid(0);
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }
    
        foreach (myblk; localFluidBlocks) {
            myblk.read_solution(make_file_name!"flow"(jobName, myblk.id, 0, flowFileExt), false);

            foreach (cell; myblk.cells) {
                cell.encode_conserved(0, 0, myblk.omegaz);
                // Even though the following call appears redundant at this point,
                // fills in some gas properties such as Prandtl number that is
                // needed for both the cfd_check and the BaldwinLomax turbulence model.
                cell.decode_conserved(0, 0, myblk.omegaz);
            }
        }
        
        // perturb design variable in complex plane
        P0 = design_variables[i].y; 
        design_variables[i].refy = P0 - EPS.im;
        
        // perturb grid
        gridUpdate(design_variables, 1); // gtl = 1

        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, vtx; myblk.vertices) {
                vtx.pos[0].refx = vtx.pos[1].x;
                vtx.pos[0].refy = vtx.pos[1].y;
            }
        }
        
        foreach (myblk; localFluidBlocks) {
            myblk.compute_primary_cell_geometric_data(0);
            myblk.compute_least_squares_setup(0);
        }

        foreach (myblk; localFluidBlocks) {
            // save mesh
            myblk.sync_vertices_to_underlying_grid(0);
            ensure_directory_is_present(make_path_name!"grid-p"(0));
            auto fileName = make_file_name!"grid-p"(jobName, myblk.id, 0, gridFileExt = "gz");
            myblk.write_underlying_grid(fileName);
        }
        
        // run steady-state solver
        iterate_to_steady_state(0, maxCPUs); // snapshotStart = 0
        //GlobalConfig.report_residuals = true;
        //sim_time = 0.0;
        //integrate_in_time(GlobalConfig.max_time);
        
        // compute objective function gradient
        objFcnM = objective_function_evaluation();
        gradients ~= (objFcnP.re-objFcnM.re)/(2.0*EPS.im);
        
        // return value to original state
        design_variables[i].refy = P0;
    }
    foreach ( i; 0..nDesignVars) {
        writef("gradient for variable %d: %.16e \n", i, gradients[i]);
    }
    writeln("simulation complete.");
}
*/
/*****************************************/
/*  STEADY-STATE SOLVER PRECONDITIONER   */
/*****************************************/

// NB. ONLY FOR UNSTRUCTURED SOLVER

void sss_preconditioner_initialisation(ref FluidBlock blk, size_t nConservative) {
    final switch (blk.myConfig.sssOptions.preconditionMatrixType) {
    case PreconditionMatrixType.block_diagonal:
        // block-diagonal preconditioner has a special residual stencil
        foreach (cell; blk.cells) {
            cell.jacobian_cell_stencil ~= cell;
            foreach ( face; cell.iface) cell.jacobian_face_stencil ~= face;
        }
        // initialise objects
        blk.Minv = new Matrix!number(nConservative, nConservative);
        foreach (cell; blk.cells) {
            cell.dPrimitive = new Matrix!number(nConservative, nConservative);
            cell.dConservative = new Matrix!number(nConservative, nConservative);
        }
        blk.cellSave = new FVCell(blk.myConfig);
        foreach (i; 0..blk.MAX_PERTURBED_INTERFACES) {
            blk.ifaceP[i] = new FVInterface(blk.myConfig, false);
        }
        break;
    case PreconditionMatrixType.ilu:
        //initialise objects
        blk.Minv = new Matrix!number(nConservative, nConservative);
	blk.P = new SMatrix!number();
        blk.cellSave = new FVCell(blk.myConfig);
        foreach(i; 0..blk.MAX_PERTURBED_INTERFACES) blk.ifaceP[i] = new FVInterface(blk.myConfig, false);
        break;
    } // end switch
}

void sss_preconditioner(ref FluidBlock blk, size_t np, double dt, size_t orderOfJacobian=1) {
    final switch (blk.myConfig.sssOptions.preconditionMatrixType) {
    case PreconditionMatrixType.block_diagonal:
        assert(blk.myConfig.sssOptions.preconditionMatrixType != PreconditionMatrixType.block_diagonal, "Error: Block Diagonal precondition matrix currently not available");
        block_diagonal_preconditioner(blk, np, dt, orderOfJacobian);
        break;
    case PreconditionMatrixType.ilu:
        ilu_preconditioner(blk, np, dt, orderOfJacobian);
        break;
    } // end switch
}

void ilu_preconditioner(ref FluidBlock blk, size_t np, double dt, size_t orderOfJacobian=1) {
    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    version(complex_numbers) Complex!double EPS = complex(0.0, 1.0e-30);
    else double EPS;
    blk.myConfig.interpolation_order = 1;
    /* form conservative Jacobian transpose */
    blk.P.aa = [];
    blk.P.ja = [];
    blk.P.ia = [];
    local_flow_jacobian_transpose(blk.P, blk, np, 1, EPS, true, true);
    number dtInv = 1.0/dt;
    foreach (i; 0 .. np*blk.cells.length) {
        blk.P[i,i] = blk.P[i,i] + dtInv;
    }
    
    /* perform ILU0 decomposition */
    int level_of_fill_in = blk.myConfig.sssOptions.iluFill;
    if (level_of_fill_in == 0) decompILU0(blk.P);
    else decompILUp(blk.P, level_of_fill_in);
    scaleLU(blk.P);
    // reset interpolation order to the global setting
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
}

void block_diagonal_preconditioner(FluidBlock blk, size_t np, double dt, size_t orderOfJacobian=1) {
    /*
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;
    size_t OMEGA = omegaIdx;
    
    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    bool transformToConserved = true;
    version(complex_numbers) Complex!double EPS = complex(0.0, 1.0e-30);
    else double EPS;
    
    blk.myConfig.interpolation_order = 1;
    
    // compute diagonal of 1st order Jacobian
    foreach(pcell; blk.cells) {
        // 0th perturbation: rho
        mixin(computeFluxDerivativesAroundCell("gas.rho", "0", true));
        // 1st perturbation: u
        mixin(computeFluxDerivativesAroundCell("vel.refx", "1", false));
        // 2nd perturbation: v
        mixin(computeFluxDerivativesAroundCell("vel.refy", "2", false));
        // 3rd perturbation: P
        mixin(computeFluxDerivativesAroundCell("gas.p", "3", true));

        // transform face flux Jacobians from primitive to conservative form
        if (transformToConserved) {
            auto gmodel = blk.myConfig.gmodel;
            foreach (f; pcell.jacobian_face_stencil) { 
                // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
                number gamma = gmodel.gamma(pcell.fs.gas);
                // form inverse transformation matrix
                blk.Minv[0,0] = to!number(1.0);
                blk.Minv[0,1] = to!number(0.0);
                blk.Minv[0,2] = to!number(0.0);
                blk.Minv[0,3] = to!number(0.0);
                // second row
                blk.Minv[1,0] = -pcell.fs.vel.x/pcell.fs.gas.rho;
                blk.Minv[1,1] = 1.0/pcell.fs.gas.rho;
                blk.Minv[1,2] = to!number(0.0);
                blk.Minv[1,3] = to!number(0.0);
                // third row
                blk.Minv[2,0] = -pcell.fs.vel.y/pcell.fs.gas.rho;
                blk.Minv[2,1] = to!number(0.0);
                blk.Minv[2,2] = 1.0/pcell.fs.gas.rho;
                blk.Minv[2,3] = to!number(0.0);
                // fourth row
                blk.Minv[3,0] = 0.5*(gamma-1.0)*(pcell.fs.vel.x*pcell.fs.vel.x+pcell.fs.vel.y*pcell.fs.vel.y);
                blk.Minv[3,1] = -pcell.fs.vel.x*(gamma-1);
                blk.Minv[3,2] = -pcell.fs.vel.y*(gamma-1);
                blk.Minv[3,3] = gamma-1.0;

                number[][] tmp;
                tmp.length = np;
                foreach (ref a; tmp) a.length = np;                        
                for (size_t i = 0; i < np; i++) {
                    for (size_t j = 0; j < np; j++) {
                        tmp[i][j] = to!number(0.0);
                        for (size_t k = 0; k < np; k++) {
                            tmp[i][j] += f.dFdU[i][k]*blk.Minv[k,j];
                        }
                    }
                }
                foreach (i; 0..np) {
                    foreach (j; 0..np) {
                        f.dFdU[i][j] = tmp[i][j];
                    }
                }
            }
        }
        
        number integral;
        number volInv = 1.0 / pcell.volume[0];
        for ( size_t ip = 0; ip < np; ++ip ) {
            for ( size_t jp = 0; jp < np; ++jp ) {
                integral = 0.0;
                foreach(fi, iface; pcell.iface) {
                    integral -= pcell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                }
                number entry = volInv * integral;                    
                pcell.dConservative[ip,jp] = -entry;
            }
        }
        // clear the interface flux Jacobian entries
        foreach (iface; pcell.jacobian_face_stencil) {
            foreach (i; 0..iface.dFdU.length) {
                foreach (j; 0..iface.dFdU[i].length) {
                    iface.dFdU[i][j] = 0.0;
                }
            }
        }
    }

    // boundary correction
    apply_boundary_conditions_for_sss_preconditioner(blk, np, orderOfJacobian, EPS, true);

    foreach (cell; blk.cells) {
        number dtInv = 1.0/dt;
        foreach (i; 0 .. np) {
            cell.dConservative[i,i] += dtInv;
        }

        foreach ( i; 0..np) {
            foreach ( j; 0..np) {
                cell.dPrimitive[i,j] = cell.dConservative[i,j];
            }
        }
        // Get an inverse ready for repeated solves.
        Matrix!number tmp = nm.bbla.inverse(cell.dConservative);
        cell.dConservative = tmp;
    }
    // reset interpolation order to the global setting
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
    */
}

void apply_boundary_conditions_for_sss_preconditioner(FluidBlock blk, size_t np, size_t orderOfJacobian, number EPS, bool transformToConserved) {
    /*
    // initialise some re-used data objects here
    number[][] dRdq; number[][] dqdQ; number[][] Aext; number[] qP;
    qP.length = np; dRdq.length = np; dqdQ.length = np; Aext.length = np; 
    foreach (ref a; dRdq) a.length = np;
    foreach (ref a; dqdQ) a.length = np;
    foreach (ref a; Aext) a.length = np;

    foreach ( bndary; blk.bc ) {
        if (bndary.type != "exchange_using_mapped_cells") {
            foreach ( bi, bface; bndary.faces) {                
                // collect interior boundary cells (bcells) and exterior ghost cell (pcell)
                FVCell[] bcells; FVCell pcell;
                if (bndary.outsigns[bi] == 1) {
                    bcells ~= bface.left_cell;
                    pcell = bface.right_cell;
                } else {
                    bcells ~= bface.right_cell;
                    pcell = bface.left_cell;
                }
                
                // form dqdQ - ghost cell derivatives 

                // 0th perturbation: rho
                mixin(computeGhostCellDerivatives("gas.rho", "0", true));
                // 1st perturbation: u
                mixin(computeGhostCellDerivatives("vel.refx", "1", false));
                // 2nd perturbation: v
                mixin(computeGhostCellDerivatives("vel.refy", "2", false));
                // 3rd perturbation: P
                mixin(computeGhostCellDerivatives("gas.p", "3", true));

                // form dRdq //       
                pcell.jacobian_cell_stencil ~= bcells;
                size_t[] idList;
                foreach ( bcell; bcells) {
                    foreach ( face; bcell.iface) {
                        if ( idList.canFind(face.id) == false ) {
                            pcell.jacobian_face_stencil ~= face;
                            idList ~= face.id;
                        }
                    }
                }
                
                // 0th perturbation: rho
                mixin(computeFluxDerivativesAroundCell("gas.rho", "0", true));
                // 1st perturbation: u
                mixin(computeFluxDerivativesAroundCell("vel.refx", "1", false));
                // 2nd perturbation: v
                mixin(computeFluxDerivativesAroundCell("vel.refy", "2", false));
                // 3rd perturbation: P
                mixin(computeFluxDerivativesAroundCell("gas.p", "3", true));
                
                foreach(bcell; pcell.jacobian_cell_stencil) {
                    number integral;
                    number volInv = 1.0 / bcell.volume[0];
                    for ( size_t ip = 0; ip < np; ++ip ) {
                        for ( size_t jp = 0; jp < np; ++jp ) {
                            integral = 0.0;
                            foreach(fi, iface; bcell.iface) {
                                integral -= bcell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                            }
                            number entry = volInv * integral;                    
                            dRdq[ip][jp] = -entry;
                        }
                    }

                    //writeln(dRdq);

                    // perform matrix-matrix multiplication
                    for (size_t i = 0; i < np; i++) {
                        for (size_t j = 0; j < np; j++) {
                            Aext[i][j] = 0;
                            for (size_t k = 0; k < np; k++) {
                                Aext[i][j] += dRdq[i][k]*dqdQ[k][j];
                            }
                        }
                    }
                    
                    // transform
                    if (transformToConserved) {
                        auto gmodel = blk.myConfig.gmodel;
                        // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
                        number gamma = gmodel.gamma(bcells[0].fs.gas);
                        // form inverse transformation matrix
                        blk.Minv[0,0] = to!number(1.0);
                        blk.Minv[0,1] = to!number(0.0);
                        blk.Minv[0,2] = to!number(0.0);
                        blk.Minv[0,3] = to!number(0.0);
                        // second row
                        blk.Minv[1,0] = -bcells[0].fs.vel.x/bcells[0].fs.gas.rho;
                        blk.Minv[1,1] = 1.0/pcell.fs.gas.rho;
                        blk.Minv[1,2] = to!number(0.0);
                        blk.Minv[1,3] = to!number(0.0);
                        // third row
                        blk.Minv[2,0] = -bcells[0].fs.vel.y/bcells[0].fs.gas.rho;
                        blk.Minv[2,1] = to!number(0.0);
                        blk.Minv[2,2] = 1.0/bcells[0].fs.gas.rho;
                        blk.Minv[2,3] = to!number(0.0);
                        // fourth row
                        blk.Minv[3,0] = 0.5*(gamma-1.0)*(bcells[0].fs.vel.x*bcells[0].fs.vel.x+bcells[0].fs.vel.y*bcells[0].fs.vel.y);
                        blk.Minv[3,1] = -bcells[0].fs.vel.x*(gamma-1);
                        blk.Minv[3,2] = -bcells[0].fs.vel.y*(gamma-1);
                        blk.Minv[3,3] = gamma-1.0;

                        number[][] tmp;
                        tmp.length = np;
                        foreach (ref a; tmp) a.length = np;
                        
                        for (size_t i = 0; i < np; i++) {
                            for (size_t j = 0; j < np; j++) {
                                tmp[i][j] = to!number(0.0);
                                for (size_t k = 0; k < np; k++) {
                                    tmp[i][j] += Aext[i][k]*blk.Minv[k,j];
                                }
                            }
                        }
                        foreach (i; 0..np) {
                            foreach (j; 0..np) {
                                Aext[i][j] = tmp[i][j];
                            }
                        }
                    }
                    
                    // add correction to boundary entry in Jacobian
                    size_t I, J;
                    for ( size_t ip = 0; ip < np; ++ip ) {
                        I = bcell.id*np + ip; // column index
                        for ( size_t jp = 0; jp < np; ++jp ) {
                            J = bcells[0].id*np + jp; // row index
                            bcells[0].dConservative[ip,jp] = bcells[0].dConservative[ip,jp] + Aext[ip][jp];
                        }
                    }
                }
                    
                // clear the interface flux Jacobian entries
                foreach (iface; pcell.jacobian_face_stencil) {
                    foreach (i; 0..iface.dFdU.length) {
                        foreach (j; 0..iface.dFdU[i].length) {
                            iface.dFdU[i][j] = 0.0;
                        }
                    }
                }
                
                pcell.jacobian_cell_stencil = [];
                pcell.jacobian_face_stencil = [];
            }
        }
    }
*/
}
