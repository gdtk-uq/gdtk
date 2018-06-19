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

enum ghost_cell_start_id = 1_000_000_000;
immutable double ESSENTIALLY_ZERO = 1.0e-50;
// some data objects used in forming the Jacobian
immutable size_t MAX_PERTURBED_INTERFACES = 40;
FVCell cellOrig;
FVInterface[MAX_PERTURBED_INTERFACES] ifaceOrig;
FVInterface[MAX_PERTURBED_INTERFACES] ifacePp;
FVInterface[MAX_PERTURBED_INTERFACES] ifacePm;

version(complex_numbers)
{
    Complex!double EPS = complex(0.0, 1.0e-20); //  0 + iEPSILON
} else {
    double EPS = 1.0; //  0 + iEPSILON
}

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


//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                   FLOW JACOBIAN FUNCTIONS
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
string computeBndaryVariableDerivatives(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "cellOrig.copy_values_from(int_cell, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "int_cell.fs."~varName~" += EPS;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(int_cell.fs.gas);";
    }
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "qp[0] = cell.fs.gas.rho;";
    codeStr ~= "qp[1] = cell.fs.vel.x;";
    codeStr ~= "qp[2] = cell.fs.vel.y;";
    codeStr ~= "qp[3] = cell.fs.gas.p;";
    codeStr ~= "int_cell.copy_values_from(cellOrig, CopyDataOption.all);";
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "dqdQ[0][" ~ posInArray ~ "] = qp[0].im/(EPS.im);";
    codeStr ~= "dqdQ[1][" ~ posInArray ~ "] = qp[1].im/(EPS.im);";
    codeStr ~= "dqdQ[2][" ~ posInArray ~ "] = qp[2].im/(EPS.im);";
    codeStr ~= "dqdQ[3][" ~ posInArray ~ "] = qp[3].im/(EPS.im);";         
    return codeStr;
}

void jacobian_bndary_correction(FluidBlock blk, ref SMatrix!number L, size_t np,
                                int orderOfJacobian) {

    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    blk.myConfig.interpolation_order = orderOfJacobian;
    
    // define, and initialise data structures
    number h; number diff;
    number[][] dRdq; number[][] dqdQ; number[][] Lext; number[] qp; number[] qm;
    qp.length = np;
    qm.length = np;
    dRdq.length = np; // number of conserved variables
    dqdQ.length = np; // number of conserved variables
    Lext.length = np; // number of conserved variables
    foreach (ref a; dRdq) a.length = np;
    foreach (ref a; dqdQ) a.length = np;
    foreach (ref a; Lext) a.length = np;
    foreach (i; 0..np) {
        foreach (j; 0..np) {
            dRdq[i][j] = 0.0;
            dqdQ[i][j] = 0.0;
            Lext[i][j] = 0.0;
        }
    }

    foreach ( bndary; blk.bc) {
        if (bndary.type != "exchange_using_mapped_cells") {
            foreach ( fj,f; bndary.faces) {
                
                // collect int and ext cells
                FVCell int_cell;
                FVCell cell;
                if (bndary.outsigns[fj] == 1) {
                    int_cell = f.left_cell;
                    cell = f.right_cell;
                } else {
                    int_cell = f.right_cell;
                    cell = f.left_cell;
                }
                
                // form dqdQ
                // 0th perturbation: rho
                mixin(computeBndaryVariableDerivatives("gas.rho", "0", true));
                // 1st perturbation: u
                mixin(computeBndaryVariableDerivatives("vel.refx", "1", false));
                // 2nd perturbation: v
                mixin(computeBndaryVariableDerivatives("vel.refy", "2", false));
                // 3rd perturbation: P
                mixin(computeBndaryVariableDerivatives("gas.p", "3", true));
                
                // form dRdq
                cell.jacobian_cell_stencil ~= int_cell;
                foreach ( face; int_cell.iface) cell.jacobian_face_stencil ~= face;
                // 0th perturbation: rho
                mixin(computeFluxFlowVariableDerivativesAroundCell("gas.rho", "0", true));
                // 1st perturbation: u
                mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refx", "1", false));
                // 2nd perturbation: v
                mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refy", "2", false));
                // 3rd perturbation: P
                mixin(computeFluxFlowVariableDerivativesAroundCell("gas.p", "3", true));
                
                number integral;
                number volInv = 1.0 / int_cell.volume[0];
                for ( size_t ip = 0; ip < np; ++ip ) {
                    for ( size_t jp = 0; jp < np; ++jp ) {
                    integral = 0.0;
                    foreach(fi, iface; int_cell.iface) {
                        integral -= int_cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                    }
                    number entry = volInv * integral;                    
                    dRdq[ip][jp] = entry;
                    }
                }
                
                // perform matrix-matrix multiplication
                for (size_t i = 0; i < np; i++) {
                    for (size_t j = 0; j < np; j++) {
                        Lext[i][j] = 0;
                        for (size_t k = 0; k < np; k++) {
                            Lext[i][j] += dRdq[i][k]*dqdQ[k][j];
                        }
                    }
                }
                
                // add correction to boundary entry in Jacobian
                size_t I, J;
                for ( size_t ip = 0; ip < np; ++ip ) {
                    I = int_cell.id*np + ip; // column index
                    for ( size_t jp = 0; jp < np; ++jp ) {
                        J = int_cell.id*np + jp; // column index
                        L[J,I] = L[J,I] + Lext[ip][jp];
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
                
                // clear working matrices
                foreach (i; 0..np) {
                    foreach (j; 0..np) {
                        dRdq[i][j] = 0.0;
                        dqdQ[i][j] = 0.0;
                        Lext[i][j] = 0.0;
                    }
                }
                cell.jacobian_cell_stencil = [];
                cell.jacobian_face_stencil = [];
            }
        }
    }
    // reset interpolation order to the global setting
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
}

void local_flow_jacobian_transpose(ref SMatrix!number A, FluidBlock blk, size_t nPrimitive,
                                    int orderOfJacobian) {
    // construct internal flow Jacobian matrix used in the domain decomposition Jacobian.

    size_t nDim = GlobalConfig.dimensions;
    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    blk.myConfig.interpolation_order = orderOfJacobian;
    
    // build jacobian stencils
    construct_inviscid_flow_jacobian_stencils(blk, orderOfJacobian);
    
    // compute transpose Jacobian entries
    construct_flow_jacobian(blk, nDim, nPrimitive, orderOfJacobian);
    
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
    jacobian_bndary_correction(blk, A, nPrimitive, blk.myConfig.interpolation_order);
}

void construct_flow_jacobian(FluidBlock blk, size_t ndim, size_t np, size_t orderOfJacobian) {
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
    number h; number diff;

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
            number integral;
            number volInv = 1.0 / c.volume[0];
            for ( size_t ip = 0; ip < np; ++ip ) {
                I = c.id*np + ip; // row index
                for ( size_t jp = 0; jp < np; ++jp ) {
                    integral = 0.0;
                    J = cell.id*np + jp; // column index
                    foreach(fi, iface; c.iface) {
                        integral -= c.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                    }
                    number JacEntry = volInv * integral;

                    //if (abs(JacEntry) > ESSENTIALLY_ZERO) {
                    blk.aa[J] ~= JacEntry;
                    blk.ja[J] ~= I;
                    //}
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
    
    // for first-order simulations the structured, unstructured stencils are identical
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
}

void compute_perturbed_flux(FVCell pcell, FluidBlock blk, size_t orderOfJacobian, FVCell[] cell_list, FVInterface[] iface_list, FVInterface[] ifaceP_list) {
    
    foreach(iface; iface_list) iface.F.clear_values();
    foreach(iface; ifaceP_list) iface.F.clear_values();
    
    // Convective flux update
    foreach(iface; iface_list) {
        auto ublk = cast(UFluidBlock) blk;
        ublk.lsq.interp_both(iface, 0, ublk.Lft, ublk.Rght); // gtl assumed 0
        iface.fs.copy_average_values_from(ublk.Lft, ublk.Rght);
        compute_interface_flux(ublk.Lft, ublk.Rght, iface, ublk.myConfig, ublk.omegaz);
    }
    
    // copy perturbed flux
    foreach(i, iface; iface_list) {
        ifaceP_list[i].copy_values_from(iface, CopyDataOption.all);
    }
}

string computeFluxFlowVariableDerivativesAroundCell(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "cellOrig.copy_values_from(cell, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "cell.fs."~varName~" += EPS;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(cell, blk, orderOfJacobian, cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePp);"; 
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "foreach (i, iface; cell.jacobian_face_stencil) {";
    codeStr ~= "iface.dFdU[0][" ~ posInArray ~ "] = ifacePp[i].F.mass.im/(EPS.im);";         
    codeStr ~= "iface.dFdU[1][" ~ posInArray ~ "] = ifacePp[i].F.momentum.x.im/(EPS.im);";
    codeStr ~= "iface.dFdU[2][" ~ posInArray ~ "] = ifacePp[i].F.momentum.y.im/(EPS.im);";
    codeStr ~= "iface.dFdU[3][" ~ posInArray ~ "] = ifacePp[i].F.total_energy.im/(EPS.im);";
    codeStr ~= "}";
    return codeStr;
}

void compute_design_variable_partial_derivatives(number[] design_variables, ref number[] g, size_t nPrimitive, bool with_k_omega) {
    size_t nDesignVars = design_variables.length;
    int gtl; int ftl; number objFcnEvalP; number objFcnEvalM; string varID; number dP; number P0;

    foreach (i; 0..nDesignVars) {
        
        // perturb design variable +ve
        gtl = 1; ftl = 1;
        
        // perturb design variable in complex plan
        P0 = design_variables[i]; 
        design_variables[i] = P0 + EPS;
        
        // perturb grid
        gridUpdate(design_variables, gtl);
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.sync_vertices_to_underlying_grid(gtl);
            myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
            myblk.compute_least_squares_setup(gtl);
        }

        /*
        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach (cell; myblk.cells) {
                writeln(cell.volume);
                foreach (face; cell.iface) {
                    writeln(face.n, ", ", face.t1, ", ", face.t2, ", ", face.area);
                }
            }
        }
        */
        
        evalRHS(0.0, ftl, gtl, with_k_omega);
        
        objFcnEvalP = objective_function_evaluation(gtl);
        
        // compute cost function sensitivity
        g[i] = (objFcnEvalP.im)/(EPS.im);
        
        // compute residual sensitivity
        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(j, cell; myblk.cells) {
                myblk.rT[i, j*nPrimitive] = to!number((cell.dUdt[1].mass.im)/(EPS.im));
                myblk.rT[i, j*nPrimitive+1] = to!number((cell.dUdt[1].momentum.x.im)/(EPS.im));
                myblk.rT[i, j*nPrimitive+2] = to!number((cell.dUdt[1].momentum.y.im)/(EPS.im));
                myblk.rT[i, j*nPrimitive+3] = to!number((cell.dUdt[1].total_energy.im)/(EPS.im));
            }
        }
        
        // restore design variable
        design_variables[i] = P0;
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



void rpcGMRES_solve(size_t nPrimitive) {    
    // restarted-GMRES settings
    size_t maxIters = GlobalConfig.sscOptions.gmresRestartInterval; // maxOuterIters
    size_t m = maxIters;
    number outerTol = GlobalConfig.sscOptions.stopOnRelativeGlobalResidual;
    size_t maxRestarts = 1000;
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
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) { blk.r0[k] = blk.f[k] - blk.r0[k];}
    }
    
    // Then compute v = r0/||r0||
    number betaTmp;
    mixin(norm2_over_blocks("betaTmp", "r0"));
    number beta = betaTmp;
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
        writef("global residual: %.16e \n",  resid);
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
            nm.smla.multiply(blk.JlocT, blk.x0, blk.r0);
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
            outFile.writef("%.16f \n", blk.psi[np*i].re);
        }

        outFile.writef("SCALARS adjoint_velx double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+1].re);
        }

        outFile.writef("SCALARS adjoint_vely double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", blk.psi[np*i+2].re);
        }

        outFile.writef("SCALARS adjoint_pressure double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) { 
            outFile.writef("%.16f \n", blk.psi[np*i+3].re);
        }

    }
}

void compute_direct_complex_step_derivatives(string jobName, int last_tindx, int maxCPUs, number[] design_variables) {

    writeln(" ");
    writeln("------------------------------------------------------");
    writeln("----EVALUATING DERIVATIVES VIA DIRECT COMPLEX STEP----");
    writeln("------------------------------------------------------");
    writeln(" ");
        
    size_t nDesignVars = design_variables.length;
    double[] gradients; number P0; number objFcn; 

    foreach ( i; 0..nDesignVars) {
        writeln("----- Computing Gradient for variable: ", i);
        foreach (myblk; localFluidBlocks) {
            ensure_directory_is_present(make_path_name!"grid"(0));
            string gridFileName = make_file_name!"grid"(jobName, myblk.id, 0, gridFileExt = "gz");
            myblk.read_new_underlying_grid(gridFileName);
            myblk.sync_vertices_from_underlying_grid(0);
            myblk.compute_primary_cell_geometric_data(0);
        }
    
        foreach (myblk; localFluidBlocks) {
            myblk.read_solution(make_file_name!"flow"(jobName, myblk.id, 0, flowFileExt), false);
        }
        
        // perturb design variable in complex plane
        P0 = design_variables[i]; 
        design_variables[i] = P0 + EPS;
    
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
        
        // compute objective function gradient
        objFcn = objective_function_evaluation();
        gradients ~= objFcn.im/EPS.im;
        
        // return value to original state
        design_variables[i] = P0;
    }
    
    writef("gradient for variable %d: %.16e \n", 1, gradients[0]);
    writef("gradient for variable %d: %.16e \n", 2, gradients[1]);
    writef("gradient for variable %d: %.16e \n", 3, gradients[2]);
    writeln("simulation complete.");
}


number objective_function_evaluation(int gtl=0, string bndaryForSurfaceIntergral = "objective_function_surface") {

    
    double[] p_target;
    foreach (myblk; parallel(localFluidBlocks,1)) {
        // target pressure distribution saved in file target.dat
        auto file = File("target.dat", "r");
        foreach(line; 0..myblk.cells.length) {
            auto lineContent = file.readln().strip();
            auto tokens = lineContent.split();
            p_target ~= 1000.0; //myblk.cells[line].fs.gas.p.re*line; //to!double(tokens[8]);
        }
    }
    
    number ObjFcn = 0.0;    
    foreach (myblk; parallel(localFluidBlocks,1)) {
        myblk.locObjFcn = 0.0;
        foreach (i, cell; myblk.cells) myblk.locObjFcn += 0.5*(cell.fs.gas.p-to!number(p_target[i]))*(cell.fs.gas.p-to!number(p_target[i])); 
    }
    foreach ( myblk; localFluidBlocks) ObjFcn += myblk.locObjFcn;
    return abs(ObjFcn);
}

void form_objective_function_sensitivity(FluidBlock blk, size_t np, string bndaryForSurfaceIntergral = "objective_function_surface") {

    // for now we have hard coded the pressure drag in the x-direction as the objective function
    size_t nLocalCells = blk.cells.length;
    blk.f.length = nLocalCells * np;

    foreach(cell; blk.cells) {
        for ( size_t ip = 0; ip < np; ++ip ) {
            blk.f[cell.id*np + ip] = 0.0;
        }
    }
    
    foreach (cell; blk.cells) {
        number origValue; number ObjFcnM; number ObjFcnP; number h;
        // for current objective function only perturbations in pressure have any effect
        origValue = cell.fs.gas.p;
        cell.fs.gas.p = origValue + EPS;
        ObjFcnP = objective_function_evaluation();
        blk.f[cell.id*np + 3] = (ObjFcnP.im)/(EPS.im);
        cell.fs.gas.p = origValue;
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


void gridUpdate(number[] designVars, size_t gtl) {
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
                    if (idList.canFind(vtx.id) == false) {
                        bndaryVtxInitPos ~= vtx.pos[0];
                        myblk.boundaryVtxIndexList ~= vtx.id;
                    }
                }
            }
        }
    }
                      
    number y0 = 0.105;
    number b = designVars[0];
    number c = designVars[1];
    number d = designVars[2];
    number scale = 1.0;
    number a = y0 - b*tanh(-d/scale);
    
    foreach (myblk; localFluidBlocks) {
        size_t[] idList;
        foreach(bndary; myblk.bc) {
            if (bndary.is_design_surface) {
                foreach( face; bndary.faces) {
                    foreach ( vtx; face.vtx) {
                        if (idList.canFind(vtx.id) == false) {
                            vtx.pos[gtl].refx = vtx.pos[0].x;
                            vtx.pos[gtl].refy = complex(vtx.pos[0].y.re, (a + b*tanh((c*vtx.pos[0].x-d)/scale)).im);
                            //vtx.pos[gtl].refy = a + b*tanh((c*vtx.pos[0].x-d)/scale);
                            
                            /*
                            number xp = sqrt(vtx.pos[0].x^^2+(vtx.pos[0].y-y0)^^2);
                            number yp = a + b*tanh((c*xp-d)/scale) - y0;
                            number theta = to!number(PI)/4.0;

                            vtx.pos[gtl].refx = xp*cos(theta) - yp*sin(theta);
                            vtx.pos[gtl].refy = xp*sin(theta) + yp*cos(theta) + y0;
                            */

                            idList ~= vtx.id;
                        }
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
                    if (idList.canFind(vtx.id) == false) {
                        bndaryVtxNewPos ~= vtx.pos[gtl];
                    }
                }
            }
        }
    }
    
    foreach (myblk; localFluidBlocks) {
        inverse_distance_weighting(myblk, bndaryVtxInitPos, bndaryVtxNewPos, gtl);
    }

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
        foreach (cell; blk.cells) {
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
        foreach (cell; blk.cells) {
            cell.jacobian_cell_stencil ~= cell;
            foreach ( face; cell.iface) cell.jacobian_face_stencil ~= face;
        }
    }
    // initialise objects
    blk.transform = new Matrix!number(nConservative, nConservative);
    foreach (cell; blk.cells) {
        cell.dPrimitive = new Matrix!number(nConservative, nConservative);
        cell.dConservative = new Matrix!number(nConservative, nConservative);
        cell.pivot.length = nConservative;
    }
    cellOrig = new FVCell(blk.myConfig);
    foreach (i; 0..MAX_PERTURBED_INTERFACES) {
        ifaceOrig[i] = new FVInterface(blk.myConfig, false);
        ifacePp[i] = new FVInterface(blk.myConfig, false);
        ifacePm[i] = new FVInterface(blk.myConfig, false);
    }
}

void jacobian_bndary_correction_for_sss_preconditioner(FluidBlock blk, size_t np, int orderOfJacobian) {

    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    blk.myConfig.interpolation_order = orderOfJacobian;
    
    // define, and initialise data structures
    number h; number diff;
    number[][] dRdq; number[][] dqdQ; number[][] Lext; number[] qp; number[] qm;
    qp.length = np;
    qm.length = np;
    dRdq.length = np; // number of conserved variables
    dqdQ.length = np; // number of conserved variables
    Lext.length = np; // number of conserved variables
    foreach (ref a; dRdq) a.length = np;
    foreach (ref a; dqdQ) a.length = np;
    foreach (ref a; Lext) a.length = np;
    foreach (i; 0..np) {
        foreach (j; 0..np) {
            dRdq[i][j] = 0.0;
            dqdQ[i][j] = 0.0;
            Lext[i][j] = 0.0;
        }
    }

    foreach ( bndary; blk.bc) {
        if (bndary.type != "exchange_using_mapped_cells") {
            foreach ( fj,f; bndary.faces) {
                
                // collect int and ext cells
                FVCell int_cell;
                FVCell cell;
                if (bndary.outsigns[fj] == 1) {
                    int_cell = f.left_cell;
                    cell = f.right_cell;
                } else {
                    int_cell = f.right_cell;
                    cell = f.left_cell;
                }
                
                // form dqdQ
                // 0th perturbation: rho
                mixin(computeBndaryVariableDerivatives("gas.rho", "0", true));
                // 1st perturbation: u
                mixin(computeBndaryVariableDerivatives("vel.refx", "1", false));
                // 2nd perturbation: v
                mixin(computeBndaryVariableDerivatives("vel.refy", "2", false));
                // 3rd perturbation: P
                mixin(computeBndaryVariableDerivatives("gas.p", "3", true));
                
                // form dRdq
                cell.jacobian_cell_stencil ~= int_cell;
                foreach ( face; int_cell.iface) cell.jacobian_face_stencil ~= face;
                // 0th perturbation: rho
                mixin(computeFluxFlowVariableDerivativesAroundCell("gas.rho", "0", true));
                // 1st perturbation: u
                mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refx", "1", false));
                // 2nd perturbation: v
                mixin(computeFluxFlowVariableDerivativesAroundCell("vel.refy", "2", false));
                // 3rd perturbation: P
                mixin(computeFluxFlowVariableDerivativesAroundCell("gas.p", "3", true));
                
                number integral;
                number volInv = 1.0 / int_cell.volume[0];
                for ( size_t ip = 0; ip < np; ++ip ) {
                    for ( size_t jp = 0; jp < np; ++jp ) {
                    integral = 0.0;
                    foreach(fi, iface; int_cell.iface) {
                        integral -= int_cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                    }
                    number entry = volInv * integral;                    
                    dRdq[ip][jp] = entry;
                    }
                }
                
                // perform matrix-matrix multiplication
                for (size_t i = 0; i < np; i++) {
                    for (size_t j = 0; j < np; j++) {
                        Lext[i][j] = 0;
                        for (size_t k = 0; k < np; k++) {
                            Lext[i][j] += dRdq[i][k]*dqdQ[k][j];
                        }
                    }
                }

                // add correction to boundary entry in Jacobian
                size_t I, J;
                for ( size_t ip = 0; ip < np; ++ip ) {
                    I = int_cell.id*np + ip; // column index
                    for ( size_t jp = 0; jp < np; ++jp ) {
                        J = int_cell.id*np + jp; // column index
                        int_cell.dPrimitive[ip,jp] = int_cell.dPrimitive[ip,jp] + Lext[ip][jp];
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
                
                // clear working matrices
                foreach (i; 0..np) {
                    foreach (j; 0..np) {
                        dRdq[i][j] = 0.0;
                        dqdQ[i][j] = 0.0;
                        Lext[i][j] = 0.0;
                    }
                }
                cell.jacobian_cell_stencil = [];
                cell.jacobian_face_stencil = [];
            }
        }
    }
    // reset interpolation order to the global setting
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
}

void sss_preconditioner(FluidBlock blk, size_t np, double dt, double EPSILON, double MU, int orderOfJacobian=1) {
    // temporarily switch the interpolation order of the config object to that of the Jacobian 
    blk.myConfig.interpolation_order = orderOfJacobian;

    // initialise some variables used in the finite difference perturbation
    number h; number diff;
    
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
        
        number integral;
        number volInv = 1.0 / cell.volume[0];
        for ( size_t ip = 0; ip < np; ++ip ) {
            for ( size_t jp = 0; jp < np; ++jp ) {
                integral = 0.0;
                foreach(fi, iface; cell.iface) {
                    integral -= cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                }
                number entry = volInv * integral;                    
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

    // boundary correction
    jacobian_bndary_correction_for_sss_preconditioner(blk, np, orderOfJacobian);
    auto gmodel = blk.myConfig.gmodel;
    // multiply by transform matrix (transforming primitive to conservative form)
    foreach (cell; blk.cells) {
        // form transformation matrix (TODO: genearlise, currently only for 2D Euler/Laminar Navier-Stokes).
        number gamma = gmodel.gamma(cell.fs.gas);
        // first row
        blk.transform[0,0] = to!number(1.0);
        blk.transform[0,1] = to!number(0.0);
        blk.transform[0,2] = to!number(0.0);
        blk.transform[0,3] = to!number(0.0);
        // second row
        blk.transform[1,0] = -cell.fs.vel.x/cell.fs.gas.rho;
        blk.transform[1,1] = 1.0/cell.fs.gas.rho;
        blk.transform[1,2] = to!number(0.0);
        blk.transform[1,3] = to!number(0.0);
        // third row
        blk.transform[2,0] = -cell.fs.vel.y/cell.fs.gas.rho;
        blk.transform[2,1] = to!number(0.0);
        blk.transform[2,2] = 1.0/cell.fs.gas.rho;
        blk.transform[2,3] = to!number(0.0);
        // fourth row
        blk.transform[3,0] = 0.5*(gamma-1.0)*(cell.fs.vel.x*cell.fs.vel.x+cell.fs.vel.y*cell.fs.vel.y);
        blk.transform[3,1] = -cell.fs.vel.x*(gamma-1);
        blk.transform[3,2] = -cell.fs.vel.y*(gamma-1);
        blk.transform[3,3] = gamma-1.0;

        
        nm.bbla.dot!number(cell.dPrimitive, blk.transform, cell.dConservative);

        number dtInv = 1.0/dt;
        foreach (i; 0 .. np) {
            cell.dConservative[i,i] += dtInv;
        }

        // Get an LU decomposition ready for repeated solves.
        nm.bbla.LUDecomp!number(cell.dConservative, cell.pivot);
    }
    
    // reset interpolation order to the global setting
    blk.myConfig.interpolation_order = GlobalConfig.interpolation_order;
}
