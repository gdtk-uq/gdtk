/** shape_sensitivity_core.d
 * 
 * Eilmer4 shape sensitivity calculator core functions.
 *
 * Author: Kyle D.
**/

module shape_sensitivity;

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

enum ghost_cell_start_id = 1_000_000_000;
immutable double ESSENTIALLY_ZERO = 1.0e-50;
// some data objects used in forming the Jacobian
immutable size_t MAX_PERTURBED_INTERFACES = 40;
FVCell cellOrig;
FVInterface[MAX_PERTURBED_INTERFACES] ifaceOrig;
FVInterface[MAX_PERTURBED_INTERFACES] ifacePp;
FVInterface[MAX_PERTURBED_INTERFACES] ifacePm;

// Module-local, global memory arrays and matrices for GMRES
double[] g0;
double[] g1;
double[] h;
double[] hR;
Matrix H0;
Matrix H1;
Matrix Gamma;
Matrix Q0;
Matrix Q1;

private lua_State* L; // module-local Lua interpreter


//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                   FLOW JACOBIAN FUNCTIONS
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
void form_external_flow_jacobian_block_phase0(FluidBlock blk, size_t np, int orderOfJacobian, double EPSILON, double MU) {
    // In this phase we compute the external effects for the neighbouring blocks from perturbing cells in those blocks. The
    // effects are attached to boundary faces rather than the boundary cells (mapped cells). This is because when partitioning
    // an unstructured mesh, it is possible that a ghost cell maybe shared by more than one neighbour cell. However, the interfaces
    // are a unique connection between two cells.

    size_t nDim = GlobalConfig.dimensions;
    if (orderOfJacobian == 0) return; // no external effects
    foreach (bc; blk.bc) {
        if (bc.type == "exchange_using_mapped_cells") {
            foreach (i, bf; bc.faces) {
                FVCell intcell;
                FVCell ghostcell;
                // construct stencil of effected cells
                if (bc.outsigns[i] == 1) {
                    intcell = bf.left_cell;
                    ghostcell = bf.right_cell;
                } else {
                    intcell = bf.right_cell;
                    ghostcell = bf.left_cell;
                }
                if (orderOfJacobian == 1) {
                    ghostcell.jacobian_cell_stencil ~= intcell;
                    foreach (f; intcell.iface) {
                        ghostcell.jacobian_face_stencil ~= f;
                    }
                } 
                else { // higher-order
                    size_t[] face_id_list;
                    foreach(cell; intcell.cell_cloud) {
                        if (cell.id < ghost_cell_start_id) ghostcell.jacobian_cell_stencil ~= cell;
                        foreach (f; cell.iface) {
                            if (face_id_list.canFind(f.id) == false) {
                                ghostcell.jacobian_face_stencil ~= f;
                                face_id_list ~= f.id;
                            } // end if
                        } // end foreach cell.iface
                    } // end foreach intcell.cell_cloud
                } // end else
                construct_flow_jacobian(ghostcell, bf, blk, nDim, np, orderOfJacobian, EPSILON, MU);
            } // end foreach bc.faces
        } // end if
    } // end foreach blk.bc
} // end form_external_flow_jacobian_block_phase0()

void form_external_flow_jacobian_block_phase1(ref SMatrix A, FluidBlock blk, size_t np, int orderOfJacobian, double EPSILON, double MU) {
    // In this phase we first reach across and grab the external effects from the neighbour blocks. We then make a record of what cells in the
    // block are on the boundary, and then finally loop through all cells in the block filling our the external effects matrix used in the
    // domain decomposition Jacobian. Note that this matrix has rows of 0 for non-boundary cells, these are stored in CRS by storing a 0
    // as the initial row value.
    
    size_t nDim = GlobalConfig.dimensions;
    // retrieve neighbour data (not MPI compatible)
    foreach (bc; blk.bc) {
        if (bc.type == "exchange_using_mapped_cells") {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellMappedCellCopy) gce;
                foreach (i, bf; bc.faces) {
                    FVCell mapped_cell = mygce.mapped_cells[i];
                    foreach ( f; mapped_cell.iface) {
                        if (f.is_on_boundary) {
                            if (abs(bf.pos.x-f.pos.x) < ESSENTIALLY_ZERO && abs(bf.pos.y-f.pos.y) < ESSENTIALLY_ZERO) {
                                bc.ghostcells[i].idList = f.idList;
                                bc.ghostcells[i].aa = f.aa;
                            } // end if
                        } // end if
                    } // end foreach mapped_cell.iface
                } // end for each bc.faces
            } // end foreach bc.preReconAction
        } // end if
    } // end foreach blk.bc

    // collect & sort boundary ids in ascending order
    size_t[] bndary_cell_id_list;
    foreach (bc; blk.bc) {
        if (bc.type == "exchange_using_mapped_cells") {
            foreach (i, bf; bc.faces) {
                FVCell cell;
                FVCell ghostcell;
                if (bc.outsigns[i] == 1) {
                    cell = bf.left_cell;
                    ghostcell = bf.right_cell;
                } else {
                    cell = bf.right_cell;
                    ghostcell = bf.left_cell;
                }
                if (bndary_cell_id_list.canFind(cell.id) == false) bndary_cell_id_list ~= cell.id;
            } // end ofreach bc.faces
        } // end if
    } // end foreach blk.bc
    bndary_cell_id_list.sort();
    
    // construct flow Jacobian external matrix
    A.ia ~= 0;
    foreach ( cell; blk.cells) {
        if (bndary_cell_id_list.canFind(cell.id)) {
            size_t[size_t] pos_array; // position of unsorted array
            size_t[] idList;
            double[] entries;
            foreach (i, f; cell.iface) {
                if (f.is_on_boundary) {
                    FVCell ghostcell;
                    if (cell.outsign[i] == 1) {
                        ghostcell = f.right_cell;
                    } else {
                        ghostcell = f.left_cell;
                    }
                    idList ~= ghostcell.idList;
                    entries ~= ghostcell.aa;
                }
            }
            foreach (i, cid; idList) pos_array[cid] = i;
            idList.sort(); // we want these in global id order
            for ( size_t ip = 0; ip < np; ++ip ) {
                size_t aaLen = A.aa.length;
                foreach (cid; idList) {
                    for ( size_t jp = 0; jp < np; ++jp ) {
                        size_t idx = pos_array[cid]*np*np+np*jp+ip;
                        size_t J = cid*np + jp;
                        if (abs(entries[idx]) > ESSENTIALLY_ZERO) {
                            A.aa ~= entries[idx];
                            A.ja ~= J;
                        }
                    }
                }
                if (aaLen == A.aa.length) {
                    // no entries have been added for this row, for CRS we must store a dummy 0 entry
                    A.aa ~= 0.0;
                    A.ja ~= 0;
                }
                A.ia ~= A.aa.length;
            }
        }
        else { // internal cell (no entries are required this row, for CRS we must store a dummy 0 entry)
            for ( size_t ip = 0; ip < np; ++ip ) {
                A.aa ~= 0.0;
                A.ja ~= 0;
                A.ia ~= A.aa.length;
            }
        } // end else
    } // end for each blk.cells
} // end form_external_flow_jacobian_block_phase1()

void form_local_flow_jacobian_block(ref SMatrix A, FluidBlock blk, size_t nPrimitive, int orderOfJacobian, double EPSILON, double MU) {
    // construct internal flow Jacobian matrix used in the domain decomposition Jacobian.

    size_t nDim = GlobalConfig.dimensions;
    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    blk.myConfig.interpolation_order = orderOfJacobian;
    
    // build jacobian stencils
    if (GlobalConfig.viscous) construct_viscous_flow_jacobian_stencils(blk, orderOfJacobian);
    else construct_inviscid_flow_jacobian_stencils(blk, orderOfJacobian);
    
    // compute transpose Jacobian entries
    construct_flow_jacobian(blk, nDim, nPrimitive, orderOfJacobian, EPSILON, MU);
    
    // build Sparse transposed Jacobian
    size_t ia = 0;
    foreach(i; 0 .. blk.cells.length*nPrimitive) { // 0..nrows
        A.aa ~= blk.aa[i];
        A.ja ~= blk.ja[i];
        A.ia ~= ia;
        ia += blk.aa[i].length;
    }
    // clear local Jacobian memory
    blk.aa = [][];
    blk.ja = [][];
    A.ia ~= A.aa.length;

    // clear jacobian stencils
    foreach ( cell; blk.cells) {
        cell.jacobian_cell_stencil = [];
        cell.jacobian_face_stencil = [];
    }

    // reset interpolation order
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
}

void construct_flow_jacobian(FVCell cell, FVInterface face, FluidBlock blk, size_t ndim, size_t np, size_t orderOfJacobian, double EPSILON, double MU) {
    /++
     + computes and stores a block local transpose Jacobian in  Compressed
     Row Storage (CSR) format.     
     + initial method for efficient computation of the Jacobian -- predecessor
     to colourings.

     TODO: turbulence, 3D
     ++/
    // initialise some variables used in the finite difference perturbation
    double h; double diff;

    // initialise objects
    cellOrig = new FVCell(blk.myConfig);
    foreach(i; 0..MAX_PERTURBED_INTERFACES) {
        ifaceOrig[i] = new FVInterface(blk.myConfig, false);
        ifacePp[i] = new FVInterface(blk.myConfig, false);
        ifacePm[i] = new FVInterface(blk.myConfig, false);
    }

    // 0th perturbation: rho
    mixin(computeFluxFlowVariableDerivativesAroundCell("gas.rho", "0", true));
    // 1st perturbation: u
    mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refx", "1", false));
    // 2nd perturbation: v
    mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refy", "2", false));
    // 3rd perturbation: P
    mixin(computeFluxFlowVariableDerivativesAroundCell("gas.p", "3", true));
    // -----------------------------------------------------
    // loop through influenced cells and fill out Jacobian 
    // -----------------------------------------------------
    // at this point we can use the cell counter ci to access the correct stencil
    // because the stencil is not necessarily in cell id order, we need to
    // do some shuffling
    int count = 0;
    //writef("perturbed cell: %d effected cells: ", cell.id);
    foreach(c; cell.jacobian_cell_stencil) {
        face.idList ~= blk.globalCellId(c.id);
        double integral;
        double volInv = 1.0 / c.volume[0];
        for ( size_t ip = 0; ip < np; ++ip ) {
            for ( size_t jp = 0; jp < np; ++jp ) {
                integral = 0.0;
                foreach(fi, iface; c.iface) {
                    integral -= c.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                }
                double JacEntry = volInv * integral;
                face.aa ~= JacEntry;
            }
        }
    }
    // clear the interface flux Jacobian entries
    foreach (iface; cell.jacobian_face_stencil) {
        foreach (i; 0..iface.dFdU.length) {
            foreach (j; 0..iface.dFdU[i].length) {
                iface.dFdU[i][j] = 0.0;
            }
        }
    }
}

void construct_flow_jacobian(FluidBlock blk, size_t ndim, size_t np, size_t orderOfJacobian, double EPSILON, double MU) {
    /++
     + computes and stores a block local transpose Jacobian in  Compressed
     Row Storage (CSR) format.     
     + initial method for efficient computation of the Jacobian -- predecessor
     to colourings.

     TODO: turbulence, 3D
     ++/
    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;
    // initialise some variables used in the finite difference perturbation
    double h; double diff;

    // initialise objects
    cellOrig = new FVCell(blk.myConfig);
    foreach(i; 0..MAX_PERTURBED_INTERFACES) {
        ifaceOrig[i] = new FVInterface(blk.myConfig, false);
        ifacePp[i] = new FVInterface(blk.myConfig, false);
        ifacePm[i] = new FVInterface(blk.myConfig, false);
    }
        
    foreach(i; 0..np*ncells) {
        blk.aa ~= [null];
        blk.ja ~= [null];
    }
    foreach(ci, cell; blk.cells) {
        // 0th perturbation: rho
        mixin(computeFluxFlowVariableDerivativesAroundCell("gas.rho", "0", true));
        // 1st perturbation: u
        mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refx", "1", false));
        // 2nd perturbation: v
        mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refy", "2", false));
        // 3rd perturbation: P
        mixin(computeFluxFlowVariableDerivativesAroundCell("gas.p", "3", true));
        // -----------------------------------------------------
        // loop through influenced cells and fill out Jacobian 
        // -----------------------------------------------------
        // at this point we can use the cell counter ci to access the correct stencil
        // because the stencil is not necessarily in cell id order, we need to
        // do some shuffling
        int count = 0;
        //writef("perturbed cell: %d effected cells: ", cell.id);
        foreach(c; cell.jacobian_cell_stencil) {
            size_t I, J; // indices in Jacobian matrix
            double integral;
            double volInv = 1.0 / c.volume[0];
            for ( size_t ip = 0; ip < np; ++ip ) {
                I = c.id*np + ip; // row index
                for ( size_t jp = 0; jp < np; ++jp ) {
                    integral = 0.0;
                    J = cell.id*np + jp; // column index
                    foreach(fi, iface; c.iface) {
                        integral -= c.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                    }
                    double JacEntry = volInv * integral;
                    
                    if (abs(JacEntry) > ESSENTIALLY_ZERO) {
                        blk.aa[J] ~= JacEntry;
                        blk.ja[J] ~= I;
                    }
                }
            }
        }
        // clear the interface flux Jacobian entries
        foreach (iface; cell.jacobian_face_stencil) {
            foreach (i; 0..iface.dFdU.length) {
                foreach (j; 0..iface.dFdU[i].length) {
                    iface.dFdU[i][j] = 0.0;
                }
            }
        }
    } // end foreach cell
}

void construct_inviscid_flow_jacobian_stencils(FluidBlock blk, size_t orderOfJacobian) {
    /++
     This is the stencil of surrounding cells (and faces) that are effected by a perturbation in a 
     cells flowstate variables, for inviscid flow.
     
     We must have the cells stored in id order, for efficient storage of the Jacobian.
     
     The general formula is to first collect the faces that will have a perturbed flux from the perturbation of the
     parent cell, then collect the left and right neighbouring cells of faces in the face stencil. For 2nd order structured
     grids though it is more efficient to first collect collect cells first.

     pcell = perturbed cell
     ++/

    // for first-order simulations the structured, unstructured stencils are identical
    if (orderOfJacobian == 0) {
        foreach(pcell; blk.cells) {
            pcell.jacobian_cell_stencil ~= pcell; 
            foreach(face; pcell.iface) {
                pcell.jacobian_face_stencil ~= face;
            }
        }
    }
    else if (orderOfJacobian == 1) { 
        foreach(pcell; blk.cells) {
            FVCell[] refs_ordered;
            FVCell[] refs_unordered;
            size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
            size_t[] cell_ids;
            
            // collect faces
            foreach(face; pcell.iface) {
                pcell.jacobian_face_stencil ~= face;
            }
            
            // for each effected face, add the neighbouring cells
            foreach(face; pcell.jacobian_face_stencil) {
                // collect (non-ghost) neighbour cells
                if (cell_ids.canFind(face.left_cell.id) == false && face.left_cell.id < ghost_cell_start_id) {
                    refs_unordered ~= face.left_cell;
                    pos_array[face.left_cell.id] = refs_unordered.length-1;
                    cell_ids ~= face.left_cell.id;
                }
                if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.id < ghost_cell_start_id) {
                    refs_unordered ~= face.right_cell;
                    pos_array[face.right_cell.id] = refs_unordered.length-1;
                    cell_ids ~= face.right_cell.id;
                }
                else continue;
            }
            
            // sort ids, and store sorted cell references
            cell_ids.sort();
            foreach(id; cell_ids) refs_ordered ~= refs_unordered[pos_array[id]];
            pcell.jacobian_cell_stencil ~= refs_ordered;
            
        } // end foreach cell
    } // end if interpolation order < 2
    
    else { // higher-order
        foreach(pcell; blk.cells) {
            FVCell[] refs_ordered;
            FVCell[] refs_unordered;
            size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
            size_t[] cell_ids;
            size_t[] face_ids;
            
            if (blk.grid_type == Grid_t.structured_grid) throw new FlowSolverException("Shape Sensitivity Calculator not implemented for structured grids, yet.");
            else { // unstructured grid
                foreach(cell; pcell.cell_cloud) {
                    // collect faces
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
                    if (cell_ids.canFind(face.left_cell.id) == false && face.left_cell.id < ghost_cell_start_id) {
                        refs_unordered ~= face.left_cell;
                        pos_array[face.left_cell.id] = refs_unordered.length-1;
                        cell_ids ~= face.left_cell.id;
                    }
                    if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.id < ghost_cell_start_id) {
                        refs_unordered ~= face.right_cell;
                        pos_array[face.right_cell.id] = refs_unordered.length-1;
                        cell_ids ~= face.right_cell.id;
                    }
                    else continue;
                }
            }
            
            // finally sort ids, and store sorted cell references
            cell_ids.sort();
            foreach(id; cell_ids) {
                refs_ordered ~= refs_unordered[pos_array[id]];
            }
            pcell.jacobian_cell_stencil ~= refs_ordered;            
        }
    }
}

void construct_viscous_flow_jacobian_stencils(FluidBlock blk, size_t orderOfJacobian) {
    /++

     This is the stencil of surrounding cells (and faces) that are effected by a perturbation in a 
     cells flowstate variables, for viscous flow.
     
     We must have the cells stored in id order, for efficient storage of the Jacobian.
     
     The general formula is to first collect the faces that will have a perturbed flux from the perturbation of the
     parent cell, then collect the left and right neighbouring cells of faces in the face stencil. For 2nd order structured
     grids though it is more efficient to first collect collect cells first.
     
     pcell = perturbed cell
     
     ++/

    // for first-order simulations the structured, unstructured stencils are identical
    foreach(pcell; blk.cells) {
        FVCell[] refs_ordered;
        FVCell[] refs_unordered;
        size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
        size_t[] cell_ids;
        size_t[] face_ids;

        if (blk.grid_type == Grid_t.structured_grid) throw new FlowSolverException("Shape Sensitivity Calculator not implemented for structured grids, yet.");
        else { // unstructured grid
            foreach(cell; pcell.cell_cloud) {
                // collect faces
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
                if (cell_ids.canFind(face.left_cell.id) == false && face.left_cell.id < ghost_cell_start_id) {
                    refs_unordered ~= face.left_cell;
                    pos_array[face.left_cell.id] = refs_unordered.length-1;
                    cell_ids ~= face.left_cell.id;
                }
                if (cell_ids.canFind(face.right_cell.id) == false && face.right_cell.id < ghost_cell_start_id) {
                    refs_unordered ~= face.right_cell;
                    pos_array[face.right_cell.id] = refs_unordered.length-1;
                    cell_ids ~= face.right_cell.id;
                }
                else continue;
            }
        }
        
        // finally sort ids, and store sorted cell references
        cell_ids.sort();
        foreach(id; cell_ids) {
            refs_ordered ~= refs_unordered[pos_array[id]];
        }
        pcell.jacobian_cell_stencil ~= refs_ordered;
    } // end foreach cell
}


void compute_perturbed_flux(FVCell pcell, FluidBlock blk, size_t orderOfJacobian, FVCell[] cell_list, FVInterface[] iface_list, FVInterface[] ifaceP_list) {
    /++
     Computes the fluxes for a subset of fvinterfaces, provided in a list.

     Note: 
     This method is used to construct the transposed flow Jacobians used in the adjoint solver, as well as to form the 1st order (diagonal) Jacobian
     for the steady-state solver preconditioner. Currently the the adjoint solver is only functional for unstructured grids, however the preconditioner
     functions for both structured, and unstructured grids.
     
     WARNING: this method should resemble the convective, and viscous flux routines found in the flow solver.
     
     ++/
    foreach(iface; iface_list) iface.F.clear_values();

    // Applies BCs for entire block -- could be more efficient
    blk.applyPreReconAction(0.0, 0, 0);  // assume sim_time = 0.0, gtl = 0, ftl = 0
    
    // Convective flux update
    if (blk.grid_type == Grid_t.structured_grid) {
        // WARNING: this method only works for forming the steady-state solver precondition matrix

        // we need to cast an SFluidBlock here to reach some methods and data
        auto sblk = cast(SFluidBlock) blk;
        size_t imin = sblk.imin; size_t imax = sblk.imax; size_t jmin = sblk.jmin; size_t jmax = sblk.jmax;

        FVCell cell = cell_list[0]; // for the preconditioner we know we only have one cell in the list (the perturbed cell)
        size_t[3] ijk = sblk.cell_id_to_ijk_indices(cell.id);
            
        // assume list of faces is in order: [i, i+1, j, j+1, k, k+1]
        size_t i; size_t j; size_t k; FVInterface iface;
        FVCell cL0; FVCell cL1; FVCell cR0; FVCell cR1; 
        // East-facing interfaces
        iface = iface_list[0]; 
        i = ijk[0]; j = ijk[1]; k = ijk[2];
        cL0 = sblk.get_cell(i-1,j,k); cL1 = sblk.get_cell(i-2,j,k);
        cR0 = sblk.get_cell(i,j,k); cR1 = sblk.get_cell(i+1,j,k);
        if ((i == imin) && (sblk.bc[Face.west].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cR0.fs); sblk.Rght.copy_values_from(cR0.fs);
        } else if ((i == imin+1) && (sblk.bc[Face.west].ghost_cell_data_available == false)) {
            sblk.one_d.interp_right(iface, cL0, cR0, cR1, cL0.iLength, cR0.iLength, cR1.iLength, sblk.Lft, sblk.Rght);
        } else if ((i == imax) && (sblk.bc[Face.east].ghost_cell_data_available == false)) {
            sblk.one_d.interp_left(iface, cL1, cL0, cR0, cL1.iLength, cL0.iLength, cR0.iLength, sblk.Lft, sblk.Rght);
        } else if ((i == imax+1) && (sblk.bc[Face.east].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cL0.fs); sblk.Rght.copy_values_from(cL0.fs);
        } else { // General symmetric reconstruction.
            sblk.one_d.interp_both(iface, cL1, cL0, cR0, cR1, cL1.iLength, cL0.iLength, cR0.iLength, cR1.iLength, sblk.Lft, sblk.Rght);
        }
        iface.fs.copy_average_values_from(sblk.Lft, sblk.Rght);
        if ((i == imin) && (sblk.bc[Face.west].convective_flux_computed_in_bc == true)) {} // do nothing
        else if ((i == imax+1) && (sblk.bc[Face.east].convective_flux_computed_in_bc == true)) {} // do nothing
        else compute_interface_flux(sblk.Lft, sblk.Rght, iface, sblk.myConfig, sblk.omegaz);
        
        iface = iface_list[1]; 
        i = ijk[0] + 1; j = ijk[1]; k = ijk[2];
        cL0 = sblk.get_cell(i-1,j,k); cL1 = sblk.get_cell(i-2,j,k);
        cR0 = sblk.get_cell(i,j,k); cR1 = sblk.get_cell(i+1,j,k);
        if ((i == imin) && (sblk.bc[Face.west].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cR0.fs); sblk.Rght.copy_values_from(cR0.fs);
        } else if ((i == imin+1) && (sblk.bc[Face.west].ghost_cell_data_available == false)) {
            sblk.one_d.interp_right(iface, cL0, cR0, cR1, cL0.iLength, cR0.iLength, cR1.iLength, sblk.Lft, sblk.Rght);
        } else if ((i == imax) && (sblk.bc[Face.east].ghost_cell_data_available == false)) {
            sblk.one_d.interp_left(iface, cL1, cL0, cR0, cL1.iLength, cL0.iLength, cR0.iLength, sblk.Lft, sblk.Rght);
        } else if ((i == imax+1) && (sblk.bc[Face.east].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cL0.fs); sblk.Rght.copy_values_from(cL0.fs);
        } else { // General symmetric reconstruction.
            sblk.one_d.interp_both(iface, cL1, cL0, cR0, cR1, cL1.iLength, cL0.iLength, cR0.iLength, cR1.iLength, sblk.Lft, sblk.Rght);
        }
        iface.fs.copy_average_values_from(sblk.Lft, sblk.Rght);
        if ((i == imin) && (sblk.bc[Face.west].convective_flux_computed_in_bc == true)) {} // do nothing
        else if ((i == imax+1) && (sblk.bc[Face.east].convective_flux_computed_in_bc == true)) {} // do nothing
        else compute_interface_flux(sblk.Lft, sblk.Rght, iface, sblk.myConfig, sblk.omegaz);
        
        // North-facing interfaces
        iface = iface_list[2]; 
        i = ijk[0]; j = ijk[1]; k = ijk[2];
        cL0 = sblk.get_cell(i,j-1,k); cL1 = sblk.get_cell(i,j-2,k);
        cR0 = sblk.get_cell(i,j,k); cR1 = sblk.get_cell(i,j+1,k);
        if ((j == jmin) && (sblk.bc[Face.south].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cR0.fs); sblk.Rght.copy_values_from(cR0.fs);
        } else if ((j == jmin+1) && (sblk.bc[Face.south].ghost_cell_data_available == false)) {
            sblk.one_d.interp_right(iface, cL0, cR0, cR1, cL0.jLength, cR0.jLength, cR1.jLength, sblk.Lft, sblk.Rght);
        } else if ((j == jmax) && (sblk.bc[Face.north].ghost_cell_data_available == false)) {
            sblk.one_d.interp_left(iface, cL1, cL0, cR0, cL1.jLength, cL0.jLength, cR0.jLength, sblk.Lft, sblk.Rght);
        } else if ((j == jmax+1) && (sblk.bc[Face.north].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cL0.fs); sblk.Rght.copy_values_from(cL0.fs);
        } else { // General symmetric reconstruction.
            sblk.one_d.interp_both(iface, cL1, cL0, cR0, cR1, cL1.jLength, cL0.jLength, cR0.jLength, cR1.jLength, sblk.Lft, sblk.Rght);
        }
        iface.fs.copy_average_values_from(sblk.Lft, sblk.Rght);
        if ((j == jmin) && (sblk.bc[Face.south].convective_flux_computed_in_bc == true)) {} // do nothing
        else if ((j == jmax+1) && (sblk.bc[Face.north].convective_flux_computed_in_bc == true)) {} // do nothing
        else compute_interface_flux(sblk.Lft, sblk.Rght, iface, sblk.myConfig, sblk.omegaz);
        
        iface = iface_list[3]; 
        i = ijk[0]; j = ijk[1]; k = ijk[2];
        cL0 = sblk.get_cell(i,j-1,k); cL1 = sblk.get_cell(i,j-2,k);
        cR0 = sblk.get_cell(i,j,k); cR1 = sblk.get_cell(i,j+1,k);
        if ((j == jmin) && (sblk.bc[Face.south].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cR0.fs); sblk.Rght.copy_values_from(cR0.fs);
        } else if ((j == jmin+1) && (sblk.bc[Face.south].ghost_cell_data_available == false)) {
            sblk.one_d.interp_right(iface, cL0, cR0, cR1, cL0.jLength, cR0.jLength, cR1.jLength, sblk.Lft, sblk.Rght);
        } else if ((j == jmax) && (sblk.bc[Face.north].ghost_cell_data_available == false)) {
            sblk.one_d.interp_left(iface, cL1, cL0, cR0, cL1.jLength, cL0.jLength, cR0.jLength, sblk.Lft, sblk.Rght);
        } else if ((j == jmax+1) && (sblk.bc[Face.north].ghost_cell_data_available == false)) {
            sblk.Lft.copy_values_from(cL0.fs); sblk.Rght.copy_values_from(cL0.fs);
        } else { // General symmetric reconstruction.
            sblk.one_d.interp_both(iface, cL1, cL0, cR0, cR1, cL1.jLength, cL0.jLength, cR0.jLength, cR1.jLength, sblk.Lft, sblk.Rght);
        }
        iface.fs.copy_average_values_from(sblk.Lft, sblk.Rght);
        if ((j == jmin) && (sblk.bc[Face.south].convective_flux_computed_in_bc == true)) {} // do nothing
        else if ((j == jmax+1) && (sblk.bc[Face.north].convective_flux_computed_in_bc == true)) {} // do nothing
        else compute_interface_flux(sblk.Lft, sblk.Rght, iface, sblk.myConfig, sblk.omegaz);
    }
    else { // unstructured grid
        if (orderOfJacobian > 1) {
            // for the MLP limiter we need to first loop over the vertices
            if (blk.myConfig.unstructured_limiter == UnstructuredLimiter.mlp) {
                FVVertex[] vtx_list;
                size_t[] vtx_id_list; 
                foreach ( cell; cell_list) {
                    foreach ( vtx; cell.vtx){
                        if (vtx_id_list.canFind(vtx.id) == false) {
                            vtx_list ~= vtx;
                            vtx_id_list ~= vtx.id;
                        }
                    }
                }
                foreach (vtx; vtx_list) {
                    FVCell[] cell_cloud;
                    foreach(cid; blk.cellIndexListPerVertex[vtx.id]) cell_cloud ~= blk.cells[cid];
                    vtx.gradients.store_max_min_values_for_mlp_limiter(cell_cloud, blk.myConfig);
                }
            }
            // compute gradients for reconstruction
            foreach(c; cell_list) {
                c.gradients.compute_lsq_values(c.cell_cloud, c.ws, blk.myConfig);
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
                case UnstructuredLimiter.venkat:
                    c.gradients.venkat_limit(c.cell_cloud, c.ws, blk.myConfig, 0);
                    break;
                } // end switch
            } // end foreach c
            
            // Fill in gradients for ghost cells so that left- and right- cells at all faces,
            // including those along block boundaries, have the latest gradient values.
            foreach (bcond; blk.bc) {
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
                if (!found_mapped_cell_bc) {
                    // There are no other cells backing the ghost cells on this boundary.
                    // Fill in ghost-cell gradients from the other side of the face.
                    foreach (i, f; bcond.faces) {
                        if (bcond.outsigns[i] == 1) {
                            f.right_cell.gradients.copy_values_from(f.left_cell.gradients);
                        } else {
                            f.left_cell.gradients.copy_values_from(f.right_cell.gradients);
                        }
                    } // end foreach f
                } // end if !found_mapped_cell_bc
            } // end foreach bcond
        } // end if interpolation_order > 1
        // compute flux
        foreach(iface; iface_list) {
            auto ublk = cast(UFluidBlock) blk;
            ublk.lsq.interp_both(iface, 0, ublk.Lft, ublk.Rght); // gtl assumed 0
            iface.fs.copy_average_values_from(ublk.Lft, ublk.Rght);
            compute_interface_flux(ublk.Lft, ublk.Rght, iface, ublk.myConfig, ublk.omegaz);
        }
    }
    // Applies BCs for entire block -- could be more efficient
    blk.applyPostConvFluxAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0

    // Viscous flux update
    if (GlobalConfig.viscous) {
        // Applies BCs for entire block -- could be more efficient
        blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);

        // currently only for least-squares at faces
        foreach(iface; iface_list) {
            iface.grad.gradients_leastsq(iface.cloud_fs, iface.cloud_pos, iface.ws_grad); // blk.flow_property_spatial_derivatives(0); 
        }

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
        foreach(iface; iface_list) {
            iface.viscous_flux_calc();
        }
        
        // Applies BCs for entire block -- could be more efficient
        blk.applyPostDiffFluxAction(0.0, 0, 0);
    }
    
    // copy perturbed flux
    foreach(i, iface; iface_list) {
        ifaceP_list[i].copy_values_from(iface, CopyDataOption.all);
    }
}

string computeFluxFlowVariableDerivativesAroundCell(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "h = (abs(cell.fs."~varName~") + MU) * EPSILON;";
    codeStr ~= "cellOrig.copy_values_from(cell, CopyDataOption.all);";
    // ------------------ negative perturbation ------------------
    codeStr ~= "cell.fs."~varName~" -= h;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(cell, blk, orderOfJacobian, cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePm);"; 
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(cell, blk, orderOfJacobian, cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePp);"; 
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "foreach (i, iface; cell.jacobian_face_stencil) {";
    codeStr ~= "diff = ifacePp[i].F.mass - ifacePm[i].F.mass;";
    codeStr ~= "iface.dFdU[0][" ~ posInArray ~ "] = diff/(2.0*h);";         
    codeStr ~= "diff = ifacePp[i].F.momentum.x - ifacePm[i].F.momentum.x;";
    codeStr ~= "iface.dFdU[1][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp[i].F.momentum.y - ifacePm[i].F.momentum.y;";
    codeStr ~= "iface.dFdU[2][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp[i].F.total_energy - ifacePm[i].F.total_energy;";
    codeStr ~= "iface.dFdU[3][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "}";
    return codeStr;
}

void compute_design_variable_partial_derivatives(Vector3[] design_variables, ref double[] g, size_t nDesignVars, size_t nPrimitive, bool with_k_omega, double ETA) {
    foreach ( i; 0..nDesignVars) {
        //evalRHS(0.0, 0, 0, with_k_omega, myblk);
        int gtl; int ftl; double objFcnEvalP; double objFcnEvalM; string varID; double P0; double dP;

        // store origianl value, and compute perturbation
        P0 = design_variables[i].y;
        dP = ETA;
        
        // perturb design variable +ve
        gtl = 1; ftl = 1;
        design_variables[i].refy = P0 + dP;
        
        // perturb grid
        gridUpdate(true, false, nDesignVars, design_variables, gtl);
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
            if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                (myblk.myConfig.interpolation_order > 1)) { 
                auto myUBlock = cast(UFluidBlock) myblk;
                myUBlock.compute_least_squares_setup(gtl);
            }
        }

        evalRHS(0.0, ftl, gtl, with_k_omega);
                
        objFcnEvalP = objective_function_evaluation(gtl);
        
        // perturb design variable -ve
        gtl = 2; ftl = 2;
        design_variables[i].refy = P0 - dP;
        
        // perturb grid
        gridUpdate(true, false, nDesignVars, design_variables, gtl);
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
            if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                (myblk.myConfig.interpolation_order > 1)) { 
                auto myUBlock = cast(UFluidBlock) myblk;
                myUBlock.compute_least_squares_setup(gtl);
            }
        }
        
        evalRHS(0.0, ftl, gtl, with_k_omega);
        
        objFcnEvalM = objective_function_evaluation(gtl);
        
        // compute cost function sensitivity
        g[i] = (objFcnEvalP-objFcnEvalM)/(2.0*dP);
	
        // compute residual sensitivity
        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, cell; myblk.cells) {
                myblk.rT[i, j*nPrimitive] = (cell.dUdt[1].mass - cell.dUdt[2].mass)/(2.0*dP);
                myblk.rT[i, j*nPrimitive+1] = (cell.dUdt[1].momentum.x - cell.dUdt[2].momentum.x)/(2.0*dP);
                myblk.rT[i, j*nPrimitive+2] = (cell.dUdt[1].momentum.y - cell.dUdt[2].momentum.y)/(2.0*dP);
                myblk.rT[i, j*nPrimitive+3] = (cell.dUdt[1].total_energy - cell.dUdt[2].total_energy)/(2.0*dP);
            }
        }
        
        // restore design variable
        design_variables[i].refy = P0;
    }
}


//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                   GMRES FUNCTIONS
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
string dot_over_blocks(string dot, string A, string B)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.dotAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.dotAcc += blk.`~A~`[k]*blk.`~B~`[k];
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
      blk.normAcc += blk.`~blkMember~`[k]*blk.`~blkMember~`[k];
   }
}
`~norm2~` = 0.0;
foreach (blk; localFluidBlocks) `~norm2~` += blk.normAcc;
`~norm2~` = sqrt(`~norm2~`);`;

}



void rpcGMRES_solve(size_t nPrimitive) {    
    // restarted-GMRES settings
    size_t maxIters = GlobalConfig.sscOptions.gmresRestartInterval; // maxOuterIters
    size_t m = maxIters;
    double outerTol = GlobalConfig.sscOptions.stopOnRelativeGlobalResidual;
    size_t maxRestarts = 1000;
    size_t iterCount;
    double resid;
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
        blk.V = new Matrix(n, m+1);
        blk.Q1 = new Matrix(m+1, m+1);
        blk.g0.length = m+1;
        blk.g1.length = m+1;
    }    

    // allocate global GMRES arrays
    g0.length = m+1;
    g1.length = m+1;
    h.length = m+1;
    hR.length = m+1;
    H0 = new Matrix(m+1, m);
    H1 = new Matrix(m+1, m);
    Gamma = new Matrix(m+1, m+1);
    Q0 = new Matrix(m+1, m+1);
    Q1 = new Matrix(m+1, m+1);
    
    // Initialise some global arrays and matrices that have already been allocated
    g0[] = 0.0;
    g1[] = 0.0;
    H0.zeros();
    H1.zeros();

    double[] Z; // global array used in the matrix-vector product
    
    // 1. Evaluate r0, beta, v1
    // r0 = b - A*x0
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.x0[] = 1.0; 
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
    double betaTmp;
    mixin(norm2_over_blocks("betaTmp", "r0"));
    double beta = betaTmp;
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
                double H0_ij_tmp;
                mixin(dot_over_blocks("H0_ij_tmp", "w", "v"));
                double H0_ij = H0_ij_tmp;
                H0[i,j] = H0_ij;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k]; 
                }
            }
            double H0_jp1j_tmp;
            mixin(norm2_over_blocks("H0_jp1j_tmp", "w"));
            double H0_jp1j = H0_jp1j_tmp;
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
                nm.bbla.dot(Q0, j+1, j+1, h, hR);
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
            nm.bbla.dot(Gamma, j+2, j+2, H0, j+1, H1);
            nm.bbla.dot(Gamma, j+2, j+2, g0, g1);
            // Accumulate Gamma rotations in Q.
            if ( j == 0 ) {
                copy(Gamma, Q1);
            }
            else {
                nm.bbla.dot(Gamma, j+2, j+2, Q0, j+2, Q1);
            }
            // Prepare for next step
            copy(H1, H0);
            g0[] = g1[];
            copy(Q1, Q0);
            // Get residual
            resid = fabs(g1[j+1]);
            writef("global residual: %.16e \n",  resid);
            // DEBUG:
            //      writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
            if ( resid <= outerTol ) {
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
        upperSolve(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.V, blk.nvars, m, blk.g1, blk.psi);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            solve(blk.P, blk.psi);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) blk.psi[k] += blk.x0[k];
        }
        if ( resid <= outerTol || r+1 == maxRestarts ) {
            // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
            // DEBUG:  writefln("Breaking restart loop.");
            break;
        }
        writeln("RESTARTING");
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
        g0[] = 0.0;
        g1[] = 0.0;
        H0.zeros();
        H1.zeros();
        // And set first residual entry
        g0[0] = beta;

    }
    nRestarts = to!int(r);
    writeln(nRestarts, " restarts.");
}

//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                   WRITE TO FILE FUNCTIONS
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
void write_adjoint_variables_to_file(FluidBlock blk, size_t np, string jobName) {
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
            outFile.writef("%.16f %.16f %.16f \n", vtx.pos[0].x, vtx.pos[0].y, vtx.pos[0].z); 
        }
        // write cell connectivity
        size_t size = ncells + 4*ncells; // TODO: only for quads, need to generalise
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
            outFile.writef("%.16f \n", blk.psi[np*i]);
        }

        outFile.writef("SCALARS adjoint_velx double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+1]);
        }

        outFile.writef("SCALARS adjoint_vely double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+2]);
        }

        outFile.writef("SCALARS adjoint_pressure double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) { 
            outFile.writef("%.16f \n", blk.psi[np*i+3]);
        }

    }
}

void writeBezierCntrlPtsToDakotaFile(Vector3[] design_variables, size_t ndvars, string jobName) {
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
    foreach ( i; 0..design_variables.length) designVariableInputContent ~= format!"%.16e    "(design_variables[i].y);
    designVariableInputContent ~= " \n"; 
    designVariableInputContent ~= "    descriptors    ";
    foreach ( i; 0..design_variables.length) {
                string descriptor;
                // y-variable
                descriptor = "'" ~ "design_var" ~ "_y" ~ to!string(i) ~ "'";
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

void readBezierCntrlPtsFromDakotaFile(size_t nDesignVars, ref Vector3[] design_variables) {
    // read in new control points
    auto f = File("params.in", "r");
    auto line = f.readln().strip;// read first line & do nothing 
    auto tokens = line.split();
    foreach ( i; 0..nDesignVars) {
        Vector3 pt;
        // x-variable
        line = f.readln().strip;
        tokens = line.split();
        pt.refx = to!double(tokens[0]);
        // y-variable
        line = f.readln().strip;
        tokens = line.split();
        pt.refy = to!double(tokens[0]);
        // z-variable
        pt.refz = 0.0;
        design_variables ~= pt;
    }
}

double finite_difference_grad(string jobName, int last_tindx, string varID) {    
    // run simulation
    string command = "bash run-flow-solver.sh";
    auto output = executeShell(command);
    
    // store simulation diagnostics file
    string commandCpy0 = "cp e4sss.diagnostics.dat adjoint/" ~ varID;
    auto outputCpy0 = executeShell(commandCpy0);

    string commandCpy1 = "cp -r plot adjoint/plot-" ~ varID;
    auto outputCpy1 = executeShell(commandCpy1);

    // read perturbed solution
    foreach (blk; localFluidBlocks) {
        blk.read_solution(make_file_name!"flow"(jobName, blk.id, last_tindx, flowFileExt), false);
    }
    
    // compute cost function
    double J;
    J = objective_function_evaluation();

    // read original grid in
    foreach (blk; localFluidBlocks) {
        ensure_directory_is_present(make_path_name!"grid-original"(0));
        string gridFileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
        blk.read_new_underlying_grid(gridFileName);
        blk.sync_vertices_from_underlying_grid(0);
    }
    // clear old simulation files
    string command0 = "bash clear.sh";
    output = executeShell(command0);

    return J;
}

void compute_design_sensitivities_with_finite_differences(string jobName, int last_tindx, Vector3[] design_variables, size_t nDesignVars, double[] adjointGradients, double DELTA) {

    double[] finiteDiffGradients;
    
    // save a copy of the original mesh for later use
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.sync_vertices_to_underlying_grid(0);
        ensure_directory_is_present(make_path_name!"grid-original"(0));
        auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
        blk.write_underlying_grid(fileName);
    }
    
    foreach ( i; 0..nDesignVars) {
        
        double objFcnEvalP; double objFcnEvalM; string varID; double P0; double dP; double err;
        
        // store origianl value, and compute perturbation
        P0 = design_variables[i].y;
        dP = DELTA;
        
        // perturb design variable +ve
        design_variables[i].refy = P0 + dP;

        // perturb grid
        gridUpdate(false, false, nDesignVars, design_variables, 1, jobName);

        // run simulation
        varID = "p-y-" ~ to!string(i);
        objFcnEvalP = finite_difference_grad(jobName, last_tindx, varID);
                
        // perturb design variable -ve
        design_variables[i].refy = P0 - dP;

        // perturb grid
        gridUpdate(false, false, nDesignVars, design_variables, 1, jobName);
        
        // run simulation
        varID = "m-y-" ~ to!string(i);
        objFcnEvalM = finite_difference_grad(jobName, last_tindx, varID);

        finiteDiffGradients ~= (objFcnEvalP - objFcnEvalM)/(2.0*dP); 
        
        err = ((adjointGradients[i] - finiteDiffGradients[i])/finiteDiffGradients[i]) * 100.0;
        writeln("finite-difference gradient for variable ", i+1, ": ", finiteDiffGradients[i], ", % error: ", abs(err));
    }
}

/*
void initLuaStateForUserDefinedObjFunc()
{
    L = init_lua_State();
    registerVector3(L);
    registerBBLA(L);
    // Load some Lua modules using 'require'.
    // There is no convenient C API expression to do the equivalent of "require"
    luaL_dostring(L, "require 'lua_helper'");
    // Set some globally available constants for the Lua state.
    lua_pushnumber(L, GlobalConfig.nFluidBlocks);
    lua_setglobal(L, "nFluidBlocks");
    lua_pushnumber(L, n_ghost_cell_layers);
    lua_setglobal(L, "nGhostCellLayers");
    // Give the user a table that holds information about
    // all of the blocks in the full simulation.
    // Note that not all of these blocks may be fully present
    // in an MPI-parallel simulation.
    lua_newtable(L);
    foreach (int i, blk; globalFluidBlocks) {
        lua_newtable(L);
        lua_pushnumber(L, blk.cells.length);
        lua_setfield(L, -2, "nCells");
        lua_pushnumber(L, blk.vertices.length);
        lua_setfield(L, -2, "nVertices");
        if ( blk.grid_type == Grid_t.structured_grid ) {
            auto sblk = cast(SFluidBlock) blk;
            assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
            lua_pushnumber(L, sblk.nicell);
            lua_setfield(L, -2, "niCells");
            lua_pushnumber(L, sblk.njcell);
            lua_setfield(L, -2, "njCells");
            lua_pushnumber(L, sblk.nkcell);
            lua_setfield(L, -2, "nkCells");
            lua_pushnumber(L, sblk.imin);
            lua_setfield(L, -2, "vtxImin");
            lua_pushnumber(L, sblk.imax+1);
            lua_setfield(L, -2, "vtxImax");
            lua_pushnumber(L, sblk.jmin);
            lua_setfield(L, -2, "vtxJmin");
            lua_pushnumber(L, sblk.jmax+1);
            lua_setfield(L, -2, "vtxJmax");
            lua_pushnumber(L, sblk.kmin);
            lua_setfield(L, -2, "vtxKmin");
            lua_pushnumber(L, sblk.kmax+1);
            lua_setfield(L, -2, "vtxKmax");
        }
        lua_rawseti(L, -2, i);
    }
    lua_setglobal(L, "blockData");
    //
    setSampleHelperFunctions(L);
    doLuaFile(L, GlobalConfig.sscOptions.userDefinedObjectiveFile);
}

double objectiveFunctionEvaluation()
{
    lua_getglobal(L, "objective_function");
    int numberArgs = 0;
    int numberResults = 1;
    if (lua_pcall(L, numberArgs, numberResults, 0) != 0) {
        luaL_error(L, "error running user-defined 'objective_function' in file: %s\n", GlobalConfig.sscOptions.userDefinedObjectiveFile.toStringz);
    }
    // Now assume a double is sitting at top of stack
    if (!lua_isnumber(L, -1)) {
         luaL_error(L, "error in return from user-defined 'objective_function' in file: %s\n A single float value was expected.", GlobalConfig.sscOptions.userDefinedObjectiveFile.toStringz);
    }
    double val = luaL_checknumber(L, -1);
    lua_settop(L, 0);
    return val;
}
*/

double objective_function_evaluation(int gtl=0, string bndaryForSurfaceIntergral = "objective_function_surface",) {

    double ObjFcn = 0.0;    
    foreach (myblk; parallel(localFluidBlocks,1)) {
        double locObjFcn = 0.0;
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
                    locObjFcn += cell.fs.gas.p*f.area[gtl]*f.n.x;
                }
            }
        }
        ObjFcn += locObjFcn;
    }
    return abs(ObjFcn);
}

void form_objective_function_sensitivity(FluidBlock blk, size_t np, double EPSILON, double MU, string bndaryForSurfaceIntergral = "objective_function_surface") {

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
		FVCell cell; double origValue; double ObjFcnM; double ObjFcnP; double h;
                // collect interior cell
                if (bndary.outsigns[i] == 1) {
		    cell = f.left_cell;
		} else {
		    cell = f.right_cell;
		}
                // for current objectrive function only perturbations in pressure have any effect
                origValue = cell.fs.gas.p;
                h = (abs(cell.fs.gas.p) + MU) * EPSILON;
                cell.fs.gas.p = origValue + h;
                ObjFcnP = objective_function_evaluation();
                cell.fs.gas.p = origValue - h;
                ObjFcnM = objective_function_evaluation();
                blk.f[cell.id*np + 3] = (ObjFcnP-ObjFcnM)/(2.0*h);
                cell.fs.gas.p = origValue;
            }
        }
    }
}





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
    if ( ftl == 0 && (GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv ||
		      GlobalConfig.flux_calculator == FluxCalculator.adaptive_hlle_ausmdv)) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.detect_shock_points();
        }
    }

     foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase0(gtl);
    }
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase1(gtl);
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
                addUDFSourceTermsToCell(blk.myL, cell, gtl, 
                                        pseudoSimTime, blk.myConfig.gmodel);
            }
            cell.time_derivatives(gtl, ftl, local_with_k_omega);
        }
    }
}


void parameteriseSurfaces(bool parameteriseSurfacesFlag, ref Vector3[] design_variables, ref size_t ndvars, double tol, int max_steps) {
    // collect boundary vertices for each design surface
    foreach (myblk; parallel(localFluidBlocks,1)) {
        foreach (bndary; myblk.bc) {
            if (bndary.is_design_surface) {
                // some arrays used in ordering the boundary vertice points
                size_t[] vtx_id_list;
                foreach(face; bndary.faces) {
                    foreach(i, vtx; face.vtx) {
                        if (vtx_id_list.canFind(vtx.id) == false) {
                            bndary.surfacePoints ~= Vector3(vtx.pos[0].x, vtx.pos[0].y, vtx.pos[0].z);
                            vtx_id_list ~= vtx.id;
                        }
                    }
                }
            }
        }
    }

    // collect all surface points for each design surface in a global array
    Vector3[] globalSurfacePoints;
    double[] ts;
    Vector3[] pts_ordered;
    Vector3[] pts_unordered;
    size_t[string] pos_array; // used to identify where the point is in the unordered list
    double[] x_pos_array;
    string[] idx_list;
    Bezier bezier;
    int num_cntrl_pts;
    foreach (myblk; localFluidBlocks) {
        foreach (bndary; myblk.bc) {
            if (bndary.is_design_surface) {
                num_cntrl_pts = bndary.num_cntrl_pts;
                foreach (pt; bndary.surfacePoints) {
                    string idx = to!string(pt.x);
                    if (idx_list.canFind(idx)) {} // do nothing
                    else {
                        pts_unordered ~= pt;
                        x_pos_array ~= pt.x;
                        pos_array[idx] = pts_unordered.length-1;
                        idx_list ~= idx;
                    } // end else
                } // end foreach bndary.surfacePoints
            } // end if
        } // end foreach blk.bc
    } // end foreach localFluidBlocks
    x_pos_array.sort();
    foreach(i; x_pos_array) pts_ordered ~= pts_unordered[pos_array[to!string(i)]];
    globalSurfacePoints.length = pts_ordered.length;
    globalSurfacePoints[] = pts_ordered[];
    // perform bezier curve fit
    bezier = optimiseBezierPoints(globalSurfacePoints, num_cntrl_pts, ts, tol, max_steps);
    foreach ( i; 1..bezier.B.length-1) {
        design_variables ~= bezier.B[i];
    }
    ndvars = design_variables.length;
    if (parameteriseSurfacesFlag) return;
    foreach (myblk; parallel(localFluidBlocks,1)) {
        prepForMeshPerturbation(myblk, bezier, globalSurfacePoints, ts);
    }
 } // end parameteriseSurfaces()



void prepForMeshPerturbation(FluidBlock blk, Bezier bezier, Vector3[] globalSurfacePoints, double[] ts) {
    // prepare arrays for mesh perturbation
    foreach (bndary; blk.bc) {
        if ( bndary.type == "exchange_using_mapped_cells") {} // do nothing
        else if (bndary.is_design_surface) {
            // some arrays used in ordering the boundary vertice points
            FVVertex[] vtx_ordered;
            FVVertex[] vtx_unordered;
            size_t[double] pos_array; // used to identify where the point is in the unordered list
            size_t[] vtx_id_list;
            double[] x_pos_array;
            foreach(face; bndary.faces) {
                foreach(i, vtx; face.vtx) {
                    if (vtx_id_list.canFind(vtx.id) == false) {
                        blk.boundaryVtxIndexList ~= vtx.id;
                        vtx_unordered ~= vtx;
                        x_pos_array ~= vtx.pos[0].refx;
                        vtx_id_list ~= vtx.id;
                        pos_array[vtx.pos[0].x] = vtx_unordered.length-1;
                    }
                }
            }
            
            // order points in ascending x-position (WARNING: it is assumed that, for a given boundary,  each x-coordinate is unique).
            x_pos_array.sort();
            foreach(i; x_pos_array) vtx_ordered ~= vtx_unordered[pos_array[i]];
            bndary.vertices ~= vtx_ordered;
        }
        else {
            size_t[] vtx_id_list;
            foreach(face; bndary.faces) {
                foreach(vtx; face.vtx) {
                    if (vtx_id_list.canFind(vtx.id) == false) {
                        bndary.vertices ~= vtx;
                        vtx_id_list ~= vtx.id;
                        blk.boundaryVtxIndexList ~= vtx.id;
                    }
                }
            }
        }
    }
    
    // copy bezier curve object to each design surface boundary, along with relevant portion of the ts array
    foreach (bndary; blk.bc) {
        if (bndary.is_design_surface) {
            bndary.bezier = bezier;
            double xi = bndary.vertices[0].pos[0].x;
            double xf = bndary.vertices[$-1].pos[0].x;
            size_t idxi; size_t idxf;
            foreach (i, t; ts ) {
                if ( abs(globalSurfacePoints[i].x - xi) < ESSENTIALLY_ZERO) idxi = i;
                else if ( abs(globalSurfacePoints[i].x - xf) < ESSENTIALLY_ZERO) idxf = i;
                else {} //do nothing
            } // foreach (ts)
            bndary.ts.length = ts[idxi..idxf+1].length;
            bndary.ts[] = ts[idxi..idxf+1];
        } // end if
    } // end foreach blk.bc
} // end prepMeshForPerturbation



void gridUpdate(bool doNotWriteGridToFile, bool readDesignVarsFromFile, size_t nDesignVars, ref Vector3[] design_variables, size_t gtl, string jobName = "") {
    if (readDesignVarsFromFile) readBezierCntrlPtsFromDakotaFile(nDesignVars, design_variables);
    // assign new control points
    // initialise all positions
    Vector3[] bndaryVtxInitPos;
    Vector3[] bndaryVtxNewPos;
    foreach (myblk; localFluidBlocks) {
        foreach(bndary; myblk.bc) {
            if ( bndary.type != "exchange_using_mapped_cells") {
                foreach(j, vtx; bndary.vertices) {
                    bndaryVtxInitPos ~= vtx.pos[0];
                    vtx.pos[1].refx = vtx.pos[0].x;
                    vtx.pos[1].refy = vtx.pos[0].y;
                    vtx.pos[2].refx = vtx.pos[0].x;
                    vtx.pos[2].refy = vtx.pos[0].y;
                }
            }
        }
            
        foreach( bndary; myblk.bc ) {
            if (bndary.is_design_surface) {
                foreach ( i; 1..bndary.bezier.B.length-1) {
                    // y-variable
                    bndary.bezier.B[i].refy = design_variables[i-1].y;
                }
                
                foreach(j, vtx; bndary.vertices) {
                    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                }
            }
        }
        foreach(bndary; myblk.bc) {
            if ( bndary.type != "exchange_using_mapped_cells") {
                foreach(j, vtx; bndary.vertices) {
                    bndaryVtxNewPos ~= vtx.pos[gtl];
                }
            }
        }
    }
    /*
    Vector3[] tempInit;
    Vector3[] tempNew;
    tempInit.length = bndaryVtxInitPos.length;
    tempNew.length = bndaryVtxNewPos.length;
    tempInit[] = bndaryVtxInitPos[];
    tempNew[] = bndaryVtxNewPos[];
    bndaryVtxInitPos = [];
    bndaryVtxNewPos = [];
    
    foreach ( i; 0..tempInit.length) {
        foreach ( j; 0..tempInit.length) {
            if ( i!=j && abs(tempInit[i].x - tempInit[j].x) < 1.0e-50 && abs(tempInit[i].y - tempInit[j].y) < 1.0e-50) {
                tempInit.remove(j);
                tempNew.remove(j);
            }
        }
    }
    foreach( i; 0..tempInit.length) {
        if (!isNaN(tempInit[i].x)) {
            bndaryVtxInitPos ~= tempInit[i];
            bndaryVtxNewPos ~= tempNew[i];
        }
    }
    */
    foreach (myblk; parallel(localFluidBlocks,1)) {
        inverse_distance_weighting(myblk, bndaryVtxInitPos, bndaryVtxNewPos, gtl);
    }

    if (doNotWriteGridToFile) return;
    foreach (myblk; parallel(localFluidBlocks,1)) {
        foreach(j, vtx; myblk.vertices) {
            vtx.pos[0].refx = vtx.pos[gtl].x;
            vtx.pos[0].refy = vtx.pos[gtl].y;
        }
        
        // save mesh
        myblk.sync_vertices_to_underlying_grid(0);
        ensure_directory_is_present(make_path_name!"grid"(0));
        auto fileName = make_file_name!"grid"(jobName, myblk.id, 0, gridFileExt = "gz");
        myblk.write_underlying_grid(fileName);
    }
}

void write_objective_fn_to_file(string fileName, double objFnEval) {
    auto outFile = File(fileName, "w");
    outFile.writef("%.16e f\n", objFnEval); 
}

void write_gradients_to_file(string fileName, double[] grad) {
    auto outFile = File(fileName, "a");
    outFile.writef("[ ");
    foreach( i; 0..grad.length ) {
        outFile.writef("%.16e ", grad[i]);
    }
    outFile.writef(" ]\n"); 
}



//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                   STEADY-STATE SOLVER PRECONDITIONER
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
void sss_preconditioner_initialisation(FluidBlock blk, size_t nConservative) {
    if (blk.grid_type == Grid_t.structured_grid) {
        auto sblk = cast(SFluidBlock) blk;
        foreach ( cell; blk.cells) {
            cell.jacobian_cell_stencil ~= cell;
            size_t[3] ijk = sblk.cell_id_to_ijk_indices(cell.id);
            size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2]; 
            // assume list of faces is in order: [i, i+1, j, j+1]
            cell.jacobian_face_stencil ~= sblk.get_ifi(i, j, k);
            cell.jacobian_face_stencil ~= sblk.get_ifi(i+1, j, k);
            cell.jacobian_face_stencil ~= sblk.get_ifj(i, j, k);
            cell.jacobian_face_stencil ~= sblk.get_ifj(i, j+1, k);
        }
    }
    else { // unstructured_grid
        foreach( cell; blk.cells) {
            cell.jacobian_cell_stencil ~= cell;
            foreach ( face; cell.iface) cell.jacobian_face_stencil ~= face;
        }
    }
    // initialise objects
    blk.transform = new Matrix(nConservative, nConservative);
    foreach( cell; blk.cells) {
        cell.dPrimitive = new Matrix(nConservative, nConservative);
        cell.dConservative = new Matrix(nConservative, nConservative);
    }
    cellOrig = new FVCell(blk.myConfig);
    foreach(i; 0..MAX_PERTURBED_INTERFACES) {
        ifaceOrig[i] = new FVInterface(blk.myConfig, false);
        ifacePp[i] = new FVInterface(blk.myConfig, false);
        ifacePm[i] = new FVInterface(blk.myConfig, false);
    }
}

void sss_preconditioner(FluidBlock blk, size_t np, double EPSILON, double MU, int orderOfJacobian=1) {

    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    blk.myConfig.interpolation_order = orderOfJacobian;

    // initialise some variables used in the finite difference perturbation
    double h; double diff;
    
    // compute diagonal of 1st order Jacobian (w.r.t. primitive variables)
    foreach(cell; blk.cells) {
        // 0th perturbation: rho
        mixin(computeFluxFlowVariableDerivativesAroundCell("gas.rho", "0", true));
        // 1st perturbation: u
        mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refx", "1", false));
        // 2nd perturbation: v
        mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refy", "2", false));
        // 3rd perturbation: P
        mixin(computeFluxFlowVariableDerivativesAroundCell("gas.p", "3", true));
        
        double integral;
        double volInv = 1.0 / cell.volume[0];
        for ( size_t ip = 0; ip < np; ++ip ) {
            for ( size_t jp = 0; jp < np; ++jp ) {
                integral = 0.0;
                foreach(fi, iface; cell.iface) {
                    integral -= cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                }
                double entry = volInv * integral;                    
                cell.dPrimitive[ip,jp] = entry;
            }
        }
        // clear the interface flux Jacobian entries
        foreach (iface; cell.jacobian_face_stencil) {
            foreach (i; 0..iface.dFdU.length) {
                foreach (j; 0..iface.dFdU[i].length) {
                    iface.dFdU[i][j] = 0.0;
                }
            }
        }
    }
    
    // multiply by transform matrix diagonal (transforming primitive to conservative form)
    foreach ( cell; blk.cells) {
        // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
        double gamma = cell.fs.gas.p/(cell.fs.gas.rho * cell.fs.gas.u); // ratio of specific heats minus 1

        // first row
        blk.transform[0,0] = 1.0;
        blk.transform[0,1] = 0.0;
        blk.transform[0,2] = 0.0;
        blk.transform[0,3] = 0.0;
        // second row
        blk.transform[1,0] = -cell.fs.vel.x/cell.fs.gas.rho;
        blk.transform[1,1] = 1.0/cell.fs.gas.rho;
        blk.transform[1,2] = 0.0;
        blk.transform[1,3] = 0.0;
        // third row
        blk.transform[2,0] = -cell.fs.vel.y/cell.fs.gas.rho;
        blk.transform[2,1] = 0.0;
        blk.transform[2,2] = 1.0/cell.fs.gas.rho;
        blk.transform[2,3] = 0.0;
        // fourth row
        blk.transform[3,0] = 0.5*(gamma-1.0)*(cell.fs.vel.x*cell.fs.vel.x+cell.fs.vel.y*cell.fs.vel.y);
        blk.transform[3,1] = -cell.fs.vel.x*(gamma-1);
        blk.transform[3,2] = -cell.fs.vel.y*(gamma-1);
        blk.transform[3,3] = gamma-1.0;

        dot(cell.dPrimitive, blk.transform, cell.dConservative);
    }
    
    // reset interpolation order to the global setting
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
}
