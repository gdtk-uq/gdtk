/** adjoint_solver.d
 * Code to construct, and solve the adjoint equations.
 *
 * Author: Kyle D.
 * Date: 2017-09-18
 *
 *
 * History:
 *   2017-09-18 : started -- 2D, 1st order, single-block, structured grid, Euler solver.
**/

import core.stdc.stdlib : exit;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.array;
import std.math;
import std.datetime;
import std.process;

import nm.smla;
import nm.bbla;
import nm.rsla;

import fluxcalc;
import grid;
import sgrid;
import usgrid;
import block;
import sblock;
import ublock;
import fvcell;
import fvinterface;
import fvvertex;
import globaldata;
import globalconfig;
import simcore;
import fvcore;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess;
import loads;
import gzip;
import geom;
import gas;
import flowsolution;

// EPSILON parameter for numerical differentiation of flux jacobian
// Value used based on Vanden and Orkwis (1996), AIAA J. 34:6 pp. 1125-1129
immutable double EPSILON = 1.0e-05;
immutable double ESSENTIALLY_ZERO = 1.0e-15;

string adjointDir = "adjoint";

void init_adjoint_dir()
{
    ensure_directory_is_present(adjointDir);
}

void main(string[] args) {
    
    init_adjoint_dir();
    
    writeln("Eilmer compressible-flow simulation code -- adjoint solver:");
    
    // -----------------------------------------------------
    // 1. Read in flow solution
    // -----------------------------------------------------

    string msg = "Usage:                              Comment:\n";
    msg       ~= "e4adjoint  [--job=<string>]     name of job\n";
    
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    try {
	getopt(args,
	       "job", &jobName,
	       );
    } catch (Exception e) {
	writeln("Problem parsing command-line options.");
	writeln("Arguments not processed:");
	args = args[1 .. $]; // Dispose of program in first arg
	foreach (arg; args) writeln("   arg: ", arg);
	write(msg);
	exit(1);
    }
    GlobalConfig.base_file_name = jobName;
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    auto last_tindx = tindx_list[$-1];

    // TODO: currently assume we only want 1 CPU
    int maxCPUs = 1;
    int maxWallClock = 5*24*3600; // 5 days default
    init_simulation(last_tindx, 0, maxCPUs, maxWallClock);
    
    writeln("simulation initialised");

    // save a copy of the original mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }

    // -----------------------------------------------------
    // 2. store the stencil of effected cells for each cell
    // -----------------------------------------------------

    // TODO: high order interpolation, currently only for 1st order interpolation    
    FVCell[][] cellStencil;
    foreach (blk; gasBlocks) {
	foreach(i, cell; blk.cells) {
	    FVCell[] cell_refs_ordered;
	    FVCell[] cell_refs_unordered;
	    size_t[size_t] array_pos; // a dictionary which uses the cellid
	                        // as an index to retrieve the array pos in cell_refs
	    size_t[] ids;
	    cell_refs_unordered ~= cell; // add the parent cell as the first reference
	    array_pos[cell.id] = cell_refs_unordered.length-1;
	    ids ~= cell.id;
	    foreach(f; cell.iface) {
		if (f.left_cell.id != cell.id &&
		    f.left_cell.id < ghost_cell_start_id) {
		    cell_refs_unordered ~= f.left_cell;
		    array_pos[f.left_cell.id] = cell_refs_unordered.length-1;
		    ids ~= f.left_cell.id;
		}
		if (f.right_cell.id != cell.id &&
		    f.right_cell.id < ghost_cell_start_id) {
		    cell_refs_unordered ~= f.right_cell;
		    array_pos[f.right_cell.id] = cell_refs_unordered.length-1;
		    ids ~= f.right_cell.id;
		}
	    }
	    ids.sort(); // sort ids
	    foreach(id; ids) {
		cell_refs_ordered ~= cell_refs_unordered[array_pos[id]];
	    }
	    cellStencil ~= cell_refs_ordered;
	}
    }

    // ------------------------------------------------------------
    // 3. Compute and store perturbed flux (form residual Jacobian)
    // ------------------------------------------------------------
    FVCell cellPp; FVCell cellPm; FVCell cellR; FVCell cellL;
    FVInterface ifacePp; FVInterface ifacePm;
    double h; double diff;
    
    // TODO: make this automatic, hard-coded for Inviscid solver
    size_t nc = 4; // number of primitive variables
    
    // TODO: need to fix this code to be block agnostic
    size_t ncells = gasBlocks[0].cells.length;
    size_t nvertices = gasBlocks[0].vertices.length;
    size_t ndim = gasBlocks[0].myConfig.dimensions;

    // TODO: sparse matrix storage, currently stores the entire Jacobian
    //SMatrix JacT;
    //JacT = new SMatrix();
    double[][] aa;
    size_t[][] ja;
    bool[] first_nonzero_filled;
    foreach(i; 0..nc*ncells) {
	aa ~= [null];
	ja ~= [null];
    }
    
    foreach (blk; gasBlocks) {
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
	    foreach(c; cellStencil[ci]) {
		size_t I, J; // indices in Jacobian matrix
		double integral;
		double volInv = 1.0 / c.volume[0];
		for ( size_t ic = 0; ic < nc; ++ic ) {
		    I = c.id*nc + ic; // row index
		    for ( size_t jc = 0; jc < nc; ++jc ) {
			integral = 0.0;
			J = cell.id*nc + jc; // column index
			foreach(fi, iface; c.iface) {
			    integral -= c.outsign[fi] * iface.dFdU[ic][jc] * iface.area[0]; // gtl=0
			}
			double JacEntry = volInv * integral;
			if (JacEntry != 0.0) {
			    aa[J] ~= volInv * integral;
			    ja[J] ~= I; //J;
			}
		    }
		}
	    }
	    // clear the interface flux Jacobian entries
	    foreach (iface; cell.iface) {
		foreach (i; 0..iface.dFdU.length) {
		    foreach (j; 0..iface.dFdU[i].length) {
			iface.dFdU[i][j] = 0.0;
		    }
		}
	    }
	} // end foreach cell
    } // end foreach block
    SMatrix JacT = new SMatrix();
    size_t ia = 0;
    foreach(i; 0 .. nc*ncells) {
	JacT.aa ~= aa[i];
	JacT.ja ~= ja[i];
	JacT.ia ~= ia;
	ia += aa[i].length;
    }
    JacT.ia ~= JacT.aa.length;
    //--------------------------------------------------------
    // Transpose Jac
    //--------------------------------------------------------
    /*
    Matrix JacT;
    JacT = transpose(Jac);

    double[] aa;
    size_t[] ja;
    size_t[] ia;

    foreach(i; 0..nc*ncells) {
	bool first_nonzero_val_in_row = false;
	foreach(j; 0..nc*ncells) {
	    if (JacT[i,j] != 0.0) {
		aa ~= JacT[i,j]; 
		ja ~= j;
		if (first_nonzero_val_in_row == false) {
		    ia ~= aa.length-1;
		    first_nonzero_val_in_row = true;
		}
	    }
	}
	if (i == nc*ncells-1) ia ~= aa.length;
    }
    */
    
    //auto JacT_sparse = new SMatrix(aa, ja, ia);
    
    // -----------------------------------------------------
    //  Form cost function sensitvity
    // -----------------------------------------------------

    // TODO: how to handle this for multiple cost functions...
    
    // Analytically form dJdV by hand differentiation
    // cost function is defined as: J(Q) = 0.5*integral[0->l] (p-p*)^2
    double[] dJdV;
    double[] p_target;

    // target pressure distribution saved in file target.dat
    auto file = File("target.dat", "r");
    foreach(i; 0 .. ncells) {
	auto lineContent = file.readln().strip();
	auto tokens = lineContent.split();
	p_target ~= to!double(tokens[8]);
    }
    writeln("target pressure imported");
    foreach (blk; gasBlocks) {
	foreach(i, cell; blk.cells) {
	    dJdV ~= 0.0;
	    dJdV ~= 0.0;
	    dJdV ~= 0.0;
	    dJdV ~= 0.5*(2.0*cell.fs.gas.p - 2.0*p_target[i]);
	}
    }
    // -----------------------------------------------------
    // 4. Solve adjoint equations
    // -----------------------------------------------------

    /*
    // form augmented matrix aug = [A|B] = [Jac|dJdQ]
    size_t ncols = nc*ncells+1;
    size_t nrows = nc*ncells;
    Matrix aug;
    aug = new Matrix(nrows, ncols);
    foreach (i; 0 .. nrows) {
	foreach (j; 0 .. ncols ) {
	    if (j < nrows) aug[i,j] =  JacT[i,j];
	    else aug[i,j] = -dJdV[i];
	}
    }
   
    // solve for adjoint variables
    gaussJordanElimination(aug);
    
    double[] psi;
    foreach (i; 0 .. nrows) {
	psi ~= aug[i,ncols-1];
    }
    */

    int maxInnerIters = 85;
    int maxOuterIters = 1000;
    int nIter = 0;
    double normRef = 0.0;
    double normNew = 0.0;
    double residTol = 1.0e-10;
    double[] psi0;
    double[] psiN;
    foreach(i; 0..nc*ncells) psi0 ~= 1.0;
    foreach(i; 0..nc*ncells) psiN ~= 1.0;
    double[] residVec;
    residVec.length = dJdV.length;
    foreach(i; 0..dJdV.length) dJdV[i] = -1.0 * dJdV[i];

    // compute ILU[0] for preconditioning
    SMatrix m = new SMatrix(JacT);
    auto M = decompILUp(m, 6);

    // compute reference norm
    multiply(JacT, psi0, residVec);
    foreach (i; 0..nc*ncells) residVec[i] = dJdV[i] - residVec[i];
    foreach (i; 0..nc*ncells) normRef += residVec[i]*residVec[i];
    normRef = sqrt(fabs(normRef));
    auto gws = GMRESWorkSpace(psi0.length, maxInnerIters);
    
    while (nIter < maxOuterIters) {
	// compute psi
	//psiN = gmres2(JacT, dJdV, psi0, maxInnerIters, residTol);
	rpcGMRES(JacT, M, dJdV, psi0, psiN, maxInnerIters, residTol, gws);
	    
	// compute new norm
	normNew = 0.0;
	multiply(JacT, psiN, residVec);
	foreach (i; 0..nc*ncells) residVec[i] = dJdV[i] - residVec[i];
	foreach (i; 0..nc*ncells) normNew += residVec[i]*residVec[i];
	normNew = sqrt(fabs(normNew));
	
	writeln("iter = ", nIter, ", resid = ", normNew/normRef,
		", adjoint: rho = ", psiN[0], ", velx = ", psiN[1], ", vely = ", psiN[2], ", p = ", psiN[3]);
	nIter += 1;
	// tolerance check
	if (normNew/normRef < residTol) {
	    writeln("final residual: ", normNew/normRef);
	    break;
	}
	foreach(i; 0..nc*ncells) psi0[i] = psiN[i];
    }
    
    //writeln(psiN[0]);
    double[] psi;
    psi.length = psiN.length;
    foreach(i; 0..nc*ncells) psi[i] = psiN[i];
    
    //writeln(Jac);
    //writeln(psiN);

    // store adjoint variables
    // TODO: formalise this
    
    foreach(i; 0 .. ncells) {
	FVCell cell = gasBlocks[0].cells[i];
	auto writer = format("%f %f %f %f %f \n", cell.pos[0].x, psi[i*nc], psi[i*nc+1], psi[i*nc+2], psi[i*nc+3]);
	append("e4_adjoint_vars.dat", writer);
    }
    
    // ---------------------------------------------------------------------
    // 5. form dR/dX -- sensitivity of residual to perturbations of the mesh 
    // ---------------------------------------------------------------------
    
    // store the stencil of effected cells for each cell
    // TODO: higher order reconsturction, currently only for 1st order
    FVCell[][] vtxStencil;
    foreach (blk; gasBlocks) {
	foreach(i, vtx; blk.vertices) {
	    FVCell[] cell_refs;
	    foreach (cid; blk.cellIndexListPerVertex[vtx.id]) {
		cell_refs ~= blk.cells[cid];
	    }
	    vtxStencil ~= cell_refs;
	}
    }

    // TODO: sparse matrix storage, currently stores the entire Jacobian
    Matrix dRdX;
    dRdX = new Matrix(ncells*nc, nvertices*ndim);
    dRdX.zeros();

    FVInterface ifaceOrig;
    foreach (blk; gasBlocks) {
	// we need to compute the converged flow solution fluxes
	// TODO: handle viscous fluxes
	blk.applyPreReconAction(0.0, 0, 0);
	blk.convective_flux_phase0();
	blk.convective_flux_phase1();
	blk.applyPostConvFluxAction(0.0, 0, 0);
	foreach(vi, vtx; blk.vertices) {
	    // 0th perturbation: x
	    mixin(computeFluxMeshPointDerivativesAroundCell("pos[0].refx", "0"));
	    // 1st perturbation: y
	    mixin(computeFluxMeshPointDerivativesAroundCell("pos[0].refy", "1"));
	    // -----------------------------------------------------
	    // loop through influenced cells and fill out Jacobian 
	    // -----------------------------------------------------
	    // at this point we can use the vertex counter vi to access
	    // the correct stencil
	    foreach(c; vtxStencil[vi]) {
		size_t I, J; // indices in Jacobian matrix
		double integral;
		double volInv = 1.0 / c.volume[0];
		for ( size_t ic = 0; ic < nc; ++ic ) {
		    I = c.id*nc + ic; // row index
		    for ( size_t jc = 0; jc < ndim; ++jc ) {
			// there are three contributions, dRdX = dRdF * dFdX + dRdA * dAdX + dRdV * dVdX
			// 1. dRdF * dFdX ---------------------------
			integral = 0.0;
			J = vtx.id*ndim + jc; //vtx.id*nc + jc; // column index
			foreach(fi, iface; c.iface) {
			    integral -= c.outsign[fi] * iface.dFdU[ic][jc]*iface.area[0]; // gtl=0
			}
			dRdX[I,J] = volInv * integral;
			// 2. dRdA * dAdX ---------------------------
			integral = 0.0;
			double dAdX; double A0; double A1;
			if (ic == 0 ) {
			    foreach(fi, iface; c.iface) {
				if (jc == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
				else if (jc == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
				integral -= c.outsign[fi]*iface.F.mass*dAdX;
			    }
			}
			else if (ic == 1) {
			    foreach(fi, iface; c.iface) {
				if (jc == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
				else if (jc == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
				integral -= c.outsign[fi]*iface.F.momentum.x*dAdX;
			    }
			}
			else if (ic == 2) {
			    foreach(fi, iface; c.iface){
				if (jc == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
				else if (jc == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
				integral -= c.outsign[fi]*iface.F.momentum.y*dAdX;
			    }
			}
			else if (ic == 3) {
			    foreach(fi, iface; c.iface){
				if (jc == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
				else if (jc == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
				integral -= c.outsign[fi]*iface.F.total_energy*dAdX;
			    }
			}
			dRdX[I,J] += volInv * integral;

			// 3. dRdV * dVdX ---------------------------
			double dVdX; double V0; double V1;
			integral = 0.0;
			if (ic == 0 ) {
			    foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.mass*iface.area[0]; }
			}
			else if (ic == 1) {
			    foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.momentum.x*iface.area[0]; }
			}
			else if (ic == 2) {
			    foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.momentum.y*iface.area[0]; }
			}
			else if (ic == 3) {
			    foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.total_energy*iface.area[0]; }
			}
			if (jc == 0) mixin(computeCellVolumeSensitivity("pos[0].refx")); // x-dimension
			else if (jc == 1) mixin(computeCellVolumeSensitivity("pos[0].refy")); // y-dimension
			//dVdX = (V1-V0)/(2*h);
			dRdX[I,J] -= volInv*volInv*integral * dVdX;
		    }
		}
	    }
	    // clear the flux Jacobian entries
	    foreach (iface; blk.faceIndexListPerVertex[vtx.id]) {
		foreach (i; 0..blk.faces[iface].dFdU.length) {
		    foreach (j; 0..blk.faces[iface].dFdU[i].length) {
			blk.faces[iface].dFdU[i][j] = 0.0;
		    }
		}
	    }
	} // end foreach cell
    } // end foreach block
    
    // -----------------------------------------------------
    // form dX/dD -- mesh perturbation specific code
    // -----------------------------------------------------
    size_t nvar = 3; // number of design variables
    size_t nsurfnodes = 101; // number of surface nodes
    // sensitivity of mesh points to movements of the surface mesh points
    Matrix dXdXb;
    dXdXb = new Matrix(nvertices*ndim, nsurfnodes*ndim);
    dXdXb.zeros();
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	foreach(vi, vtx; sblk.vertices) {
	    ulong i; ulong j;
	    j = vtx.id/(sblk.imax);
	    i = vtx.id - j*(sblk.imax);
	    foreach(vbi; 0..nsurfnodes) {
		if (i != 0 && i == vbi) dXdXb[2*i+1,2*vbi+1] = -1.0*(1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1));
	    }
	}
    } 
    // sensitivity of surface mesh points to movements of design variables
    Matrix dXbdD;
    dXbdD = new Matrix(ndim*nsurfnodes, nvar);
    dXbdD.zeros();
    // shape parameters
    double scale = 1.5;
    double b = 0.07;
    double c = 0.8;
    double d = 3.8;
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	foreach(vi; 0..nsurfnodes) {
	    FVVertex vtx = sblk.get_vtx(vi+2,sblk.jmin,sblk.kmin);
	    dXbdD[vi*ndim+0,0] = 0.0;
	    dXbdD[vi*ndim+0,1] = 0.0;
	    dXbdD[vi*ndim+0,2] = 0.0;
	    dXbdD[vi*ndim+1,0] = tanh(d/scale) + tanh((c*vtx.pos[0].x - d)/scale);
	    dXbdD[vi*ndim+1,1] = b*vtx.pos[0].x*(-pow(tanh((c*vtx.pos[0].x - d)/scale),2) + 1)/scale;
	    dXbdD[vi*ndim+1,2] = b*(-pow(tanh(d/scale),2) + 1)/scale - b*(-pow(tanh((c*vtx.pos[0].x - d)/scale),2) + 1)/scale;
	}
    } 
    // sensitivity of mesh points to movements of the design variables
    Matrix dXdD;
    dXdD = new Matrix(ndim*nvertices, nvar);
    dot(dXdXb, dXbdD, dXdD);
    // compute transposes
    Matrix dXdD_T; Matrix dRdX_T;
    dXdD_T = transpose(dXdD);
    dRdX_T = transpose(dRdX);
    // temp matrix multiplication
    Matrix tempMatrix;
    tempMatrix = new Matrix(nvar, ncells*nc);
    dot(dXdD_T, dRdX_T, tempMatrix);
    // compute gradient
    double[3] grad;
    dot(tempMatrix, psi, grad);

    writeln("adjoint gradients [dLdB, dLdC, dLdD] = ", grad);
    
    // -----------------------------------------------------
    // Finite difference verification
    // -----------------------------------------------------
    /*
    // save original mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }
    */
    // read original grid in
    foreach (blk; gasBlocks) { 
    	//blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	blk.read_grid(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"), 0);
    }

    
    //double J0 = 0.0;
    //foreach (blk; gasBlocks) {
    // 	foreach (i, cell; blk.cells) {
    //	    J0 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
    //	}
    //}
    
    // perturb b +ve --------------------------------------------------------------------------------
    double del_b = b*EPSILON + EPSILON;
    b += del_b;
    // perturb mesh
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	    double delta;
	    double y_old;
	    for ( size_t j = sblk.jmax+1; j >= sblk.jmin; --j ) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		    //writeln(j, ", ", sblk.jmax, ", ", delta, ", ", y_old, ", ", y_new);
		} else {
		    vtx.pos[0].refy += (1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1)) * delta;
		    //writeln(j, ", ", sblk.jmax, ", ", delta, ", ", vtx.pos[0].refy);
		}
	    }
	}
    }
    // save mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }

    // run simulation
    string command = "bash run-perturb.sh";
    auto output = executeShell(command);
    //writeln(output[1]);
    //int status = output[0];

    foreach (blk; gasBlocks) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }
    
    double J1 = 0.0;
    foreach (blk; gasBlocks) {
	foreach (i, cell; blk.cells) {
	    J1 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }

    
    // read original grid in
    foreach (blk; gasBlocks) { 
    	//blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	blk.read_grid(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"), 0);
    }

    // clear old simulation files
    string command0 = "bash clear.sh";
    output = executeShell(command0);
    
    // perturb b -ve --------------------------------------------------------------------------------
    b -= 2.0*del_b;
    // perturb mesh
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	    double delta;
	    double y_old;
	    for ( size_t j = sblk.jmax+1; j >= sblk.jmin; --j ) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		} else {
		    vtx.pos[0].refy += (1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1)) * delta;
		}
	    }
	}
    }
    
    // save mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }
    
    // run simulation
    output = executeShell(command);
    //writeln(output[1]);
    //int status = output[0];

    foreach (blk; gasBlocks) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }
    
    double J0 = 0.0;
    foreach (blk; gasBlocks) {
	foreach (i, cell; blk.cells) {
	    J0 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }
    //writef("%.16f, %.16f \n", J1, J0);
    double grad_b = (J1 -J0)/(2.0*del_b);
    writeln("FD dLdB = ", grad_b, ", % error = ", abs((grad_b - grad[0])/grad_b * 100));

    // clear old simulation files
    output = executeShell(command0);

    // perturb c +ve --------------------------------------------------------------------------------
    // read original grid in
    foreach (blk; gasBlocks) { 
    	//blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	blk.read_grid(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"), 0);
    }

    b += del_b; // return b to original state
    double del_c = c*EPSILON + EPSILON;
    c += del_c;
    // perturb mesh
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	    double delta;
	    double y_old;
	    for ( size_t j = sblk.jmax+1; j >= sblk.jmin; --j ) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		} else {
		    vtx.pos[0].refy += (1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1)) * delta;
		}
	    }
	}
    }
    // save mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }

    // run simulation
    command = "bash run-perturb.sh";
    output = executeShell(command);
    //writeln(output[1]);
    //int status = output[0];

    foreach (blk; gasBlocks) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }
    
    J1 = 0.0;
    foreach (blk; gasBlocks) {
	foreach (i, cell; blk.cells) {
	    J1 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }

    
    // read original grid in
    foreach (blk; gasBlocks) { 
    	//blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	blk.read_grid(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"), 0);
    }

    // clear old simulation files
    command0 = "bash clear.sh";
    output = executeShell(command0);
    
    // perturb c -ve --------------------------------------------------------------------------------
    c -= 2.0*del_c;
    // perturb mesh
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	    double delta;
	    double y_old;
	    for ( size_t j = sblk.jmax+1; j >= sblk.jmin; --j ) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		} else {
		    vtx.pos[0].refy += (1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1)) * delta;
		}
	    }
	}
    }
    
    // save mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }
    
    // run simulation
    output = executeShell(command);
    //writeln(output[1]);
    //int status = output[0];

    foreach (blk; gasBlocks) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }
    
    J0 = 0.0;
    foreach (blk; gasBlocks) {
	foreach (i, cell; blk.cells) {
	    J0 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }
    //writef("%.16f, %.16f \n", J1, J0);
    double grad_c = (J1-J0)/(2.0*del_c);
    writeln("FD dLdC = ", grad_c, ", % error = ", abs((grad_c - grad[1])/grad_c * 100));

    output = executeShell(command0);
    
    // perturb d +ve --------------------------------------------------------------------------------
    // read original grid in
    foreach (blk; gasBlocks) { 
    	//blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	blk.read_grid(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"), 0);
    }

    c += del_c; // return c to original state
    double del_d = d*EPSILON + EPSILON;
    d += del_d;
    // perturb mesh
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	    double delta;
	    double y_old;
	    for ( size_t j = sblk.jmax+1; j >= sblk.jmin; --j ) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		} else {
		    vtx.pos[0].refy += (1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1)) * delta;
		}
	    }
	}
    }
    // save mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }

    // run simulation
    command = "bash run-perturb.sh";
    output = executeShell(command);
    //writeln(output[1]);
    //int status = output[0];

    foreach (blk; gasBlocks) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }
    
    J1 = 0.0;
    foreach (blk; gasBlocks) {
	foreach (i, cell; blk.cells) {
	    J1 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }

    
    // read original grid in
    foreach (blk; gasBlocks) { 
    	//blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	blk.read_grid(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"), 0);
    }

    // clear old simulation files
    command0 = "bash clear.sh";
    output = executeShell(command0);
    
    // perturb c -ve --------------------------------------------------------------------------------
    d -= 2.0*del_d;
        // perturb mesh
    foreach (blk; gasBlocks) {
	SBlock sblk = cast(SBlock) blk;
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	    double delta;
	    double y_old;
	    for ( size_t j = sblk.jmax+1; j >= sblk.jmin; --j ) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		} else {
		    vtx.pos[0].refy += (1.0 - ((sblk.jmax-1) - (j-2))/(sblk.jmax-1)) * delta;
		}
	    }
	}
    }
    // save mesh
    foreach (blk; gasBlocks) {
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
    }
    
    // run simulation
    output = executeShell(command);
    //writeln(output[1]);
    //int status = output[0];

    foreach (blk; gasBlocks) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }
    
    J0 = 0.0;
    foreach (blk; gasBlocks) {
	foreach (i, cell; blk.cells) {
	    J0 += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }
    //writef("%.16f, %.16f \n", J1, J0);
    double grad_d = (J1-J0)/(2.0*del_d);
    writeln("FD dLdD = ", grad_d, ", % error = ", abs((grad_d - grad[2])/grad_d * 100));
    
    writeln("Done simulation.");
}

string computeCellVolumeSensitivity(string varName)
{
    string codeStr;
    codeStr ~= "h = vtx."~varName~" * EPSILON + EPSILON;";
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "c.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "V1 = c.volume[0];";
    codeStr ~= "vtx."~varName~" -= 2*h;";
    codeStr ~= "c.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "V0 = c.volume[0];";
    codeStr ~= "dVdX = (V1-V0)/(2*h);";
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "c.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    return codeStr;
}

string computeInterfaceAreaSensitivity(string varName)
{
    string codeStr;
    codeStr ~= "h = vtx."~varName~" * EPSILON + EPSILON;";
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "A1 = iface.area[0];";
    codeStr ~= "vtx."~varName~" -= 2*h;";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "A0 = iface.area[0];";
    codeStr ~= "dAdX = (A1-A0)/(2*h);";
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    return codeStr;
}

string computeFluxFlowVariableDerivativesAroundCell(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "cellPp = new FVCell(dedicatedConfig[blk.id]);";
    codeStr ~= "ifacePp = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "cellPm = new FVCell(dedicatedConfig[blk.id]);";
    codeStr ~= "ifacePm = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "h = cell.fs."~varName~" * EPSILON + EPSILON;";
    codeStr ~= "cellPm.copy_values_from(cell, CopyDataOption.all);";
    codeStr ~= "cellPm.fs."~varName~" -= h;";
    codeStr ~= "cellPp.copy_values_from(cell, CopyDataOption.all);";
    codeStr ~= "cellPp.fs."~varName~" += h;";
    codeStr ~= "foreach(iface; cell.iface) {";
    // ------------------ negative perturbation ------------------
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cellPm.fs.gas);";
    }
    // ------------------ apply cell effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    codeStr ~= "cell.fs."~varName~" -= h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);";  // assume sim_time = 0.0, gtl = 0, ftl = 0
    //codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "}";
    // ------------------ compute interflux ------------------
    codeStr ~= "if(iface.left_cell.id == cellPm.id) {";
    codeStr ~= "cellR = iface.right_cell;";
    codeStr ~= "cellL = cellPm;";
    codeStr ~= "}";
    codeStr ~= "else {";
    codeStr ~= "cellR = cellPm;";
    codeStr ~= "cellL = iface.left_cell;";
    codeStr ~= "}";
    codeStr ~= "blk.Lft.copy_values_from(cellL.fs);";
    codeStr ~= "blk.Rght.copy_values_from(cellR.fs);";
    codeStr ~= "compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);";
    // ------------------ apply interface effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    codeStr ~= "cell.fs."~varName~" -= h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    //codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);"; // assume sim_time = 0.0, gtl = 0, ftl = 0
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "}";
    codeStr ~= "ifacePm.copy_values_from(iface, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    // update thermo
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cellPp.fs.gas);";
    }
    // ------------------ apply cell effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);"; // assume sim_time = 0.0, gtl = 0, ftl = 0
    //codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "cell.fs."~varName~" -= h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "}";
    // ------------------ compute interface flux ------------------
    codeStr ~= "if(iface.left_cell.id == cellPp.id) {";
    codeStr ~= "cellR = iface.right_cell;";
    codeStr ~= "cellL = cellPp;";
    codeStr ~= "}";
    codeStr ~= "else {";
    codeStr ~= "cellR = cellPp;";
    codeStr ~= "cellL = iface.left_cell;";
    codeStr ~= "}";
    codeStr ~= "blk.Lft.copy_values_from(cellL.fs);";
    codeStr ~= "blk.Rght.copy_values_from(cellR.fs);";
    codeStr ~= "compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);";
    // ------------------ apply interface effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    //codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);"; // assume sim_time = 0.0, gtl = 0, ftl = 0
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "cell.fs."~varName~" -= h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "}";
    codeStr ~= "ifacePp.copy_values_from(iface, CopyDataOption.all);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "diff = ifacePp.F.mass - ifacePm.F.mass;";
    codeStr ~= "iface.dFdU[0][" ~ posInArray ~ "] = diff/(2.0*h);";	    
    codeStr ~= "diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;";
    codeStr ~= "iface.dFdU[1][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;";
    codeStr ~= "iface.dFdU[2][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp.F.total_energy - ifacePm.F.total_energy;";
    codeStr ~= "iface.dFdU[3][" ~ posInArray ~ "] = diff/(2.0*h);";
    //
    codeStr ~= "}";

    return codeStr;
}

string computeFluxMeshPointDerivativesAroundCell(string varName, string posInArray)
{
    string codeStr;
    codeStr ~= "ifacePp = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "ifacePm = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "ifaceOrig = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "h = vtx."~varName~" * EPSILON + EPSILON;";
    codeStr ~= "foreach (faceid; blk.faceIndexListPerVertex[vtx.id]) { ";
    codeStr ~= "FVInterface iface = blk.faces[faceid];";
    codeStr ~= "ifaceOrig.copy_values_from(iface, CopyDataOption.all);";
    // ------------------ negative perturbation ------------------
    codeStr ~= "vtx."~varName~" -= h;";
    // ------------------ apply grid metrics ------------------
    //codeStr ~= "foreach (cid; blk.cellIndexListPerVertex[vtx.id]) { blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric); }";
    //codeStr ~= "foreach (fid; blk.faceIndexListPerVertex[vtx.id]) { blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric); }";
    codeStr ~= "iface.right_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.left_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";

    //codeStr ~= "blk.compute_primary_cell_geometric_data(0);";
    //codeStr ~= "if (GlobalConfig.do_compute_distance_to_nearest_wall) {";
    //codeStr ~= "blk.compute_distance_to_nearest_wall_for_all_cells(0);";
    //codeStr ~= "}";
    //codeStr ~= "if ((blk.grid_type == Grid_t.unstructured_grid) &&";
    //codeStr ~= "(blk.myConfig.interpolation_order > 1)) {"; 
    //codeStr ~= "blk.compute_least_squares_setup_for_reconstruction(0);";
    //codeStr ~= "}";
    // ------------------ apply cell effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);";  // assume sim_time = 0.0, gtl = 0, ftl = 0
    //codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "}";
    // ------------------ compute interface flux ------------------
    codeStr ~= "cellR = iface.right_cell;";
    codeStr ~= "cellL = iface.left_cell;";
    codeStr ~= "blk.Lft.copy_values_from(cellL.fs);";
    codeStr ~= "blk.Rght.copy_values_from(cellR.fs);";
    codeStr ~= "compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);";
    // ------------------ apply interface effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    //codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);"; // assume sim_time = 0.0, gtl = 0, ftl = 0
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "}";
    codeStr ~= "ifacePm.copy_values_from(iface, CopyDataOption.all);";
    codeStr ~= "vtx."~varName~" += h;";
    // ------------------ positive perturbation ------------------
    codeStr ~= "vtx."~varName~" += h;";
    // ------------------ apply grid metrics ------------------
    //codeStr ~= "foreach (cid; blk.cellIndexListPerVertex[vtx.id]) { blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric); }";
    //codeStr ~= "foreach (fid; blk.faceIndexListPerVertex[vtx.id]) { blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric); }";

    codeStr ~= "iface.right_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.left_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    
    //codeStr ~= "if (GlobalConfig.do_compute_distance_to_nearest_wall) {";
    //codeStr ~= "blk.compute_distance_to_nearest_wall_for_all_cells(0);";
    //codeStr ~= "}";
    //codeStr ~= "if ((blk.grid_type == Grid_t.unstructured_grid) &&";
     //codeStr ~= "(blk.myConfig.interpolation_order > 1)) {"; 
    //codeStr ~= "blk.compute_least_squares_setup_for_reconstruction(0);";
    //codeStr ~= "}";
    // ------------------ apply cell effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);"; // assume sim_time = 0.0, gtl = 0, ftl = 0
    //codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "}";
    // ------------------ compute interface flux ------------------
    codeStr ~= "cellR = iface.right_cell;";
    codeStr ~= "cellL = iface.left_cell;";
    codeStr ~= "blk.Lft.copy_values_from(cellL.fs);";
    codeStr ~= "blk.Rght.copy_values_from(cellR.fs);";
    codeStr ~= "compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);";
    // ------------------ apply interface effect bcs ------------------
    codeStr ~= "if (iface.is_on_boundary) {";
    //codeStr ~= "blk.applyPreReconAction(0.0, 0, 0);"; // assume sim_time = 0.0, gtl = 0, ftl = 0
    codeStr ~= "blk.applyPostConvFluxAction(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);";
    //codeStr ~= "blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);";
    //codeStr ~= "blk.applyPostDiffFluxAction(0.0, 0, 0);";
    codeStr ~= "}";
    codeStr ~= "ifacePp.copy_values_from(iface, CopyDataOption.all);";
    codeStr ~= "vtx."~varName~" -= h;";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "diff = ifacePp.F.mass - ifacePm.F.mass;";
    codeStr ~= "iface.dFdU[0][" ~ posInArray ~ "] = diff/(2.0*h);";	    
    codeStr ~= "diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;";
    codeStr ~= "iface.dFdU[1][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;";
    codeStr ~= "iface.dFdU[2][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp.F.total_energy - ifacePm.F.total_energy;";
    codeStr ~= "iface.dFdU[3][" ~ posInArray ~ "] = diff/(2.0*h);";
    // ------------------ restore original geometry ------------------
    //codeStr ~= "foreach (cid; blk.cellIndexListPerVertex[vtx.id]) { blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric); }";
    //codeStr ~= "foreach (fid; blk.faceIndexListPerVertex[vtx.id]) { blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric); }";
    codeStr ~= "iface.right_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.left_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    //
    codeStr ~= "iface.copy_values_from(ifaceOrig, CopyDataOption.all);";
    codeStr ~= "}";

    return codeStr;
}

