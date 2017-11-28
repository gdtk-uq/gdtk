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
import core.memory;
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
import geom;
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
import gas;
import flowsolution;
import onedinterp;

// EPSILON parameter for numerical differentiation of flux jacobian
// Value used based on Vanden and Orkwis (1996), AIAA J. 34:6 pp. 1125-1129
immutable double EPSILON = 1.0e-4;
immutable double MU = 1.0e-4;
immutable double ESSENTIALLY_ZERO = 1.0e-15;
immutable bool withFiniteDiffVerification = true;

string adjointDir = "adjoint";

void init_adjoint_dir()
{
    ensure_directory_is_present(adjointDir);
}

void main(string[] args) {

    // -------------------------------------------
    // 0. Some house-keeping, and initialisations
    // -------------------------------------------	 

    init_adjoint_dir();
    writeln("Eilmer compressible-flow simulation code -- adjoint solver:");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                              Comment:\n";
    msg       ~= "e4adjoint  [--job=<string>]            name of job\n";
    msg       ~= "           [--max-cpus=<int>]          defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine\n";
    msg       ~= "           [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "           [--help]                    writes this message\n";
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;
    try {
	getopt(args,
	       "job", &jobName,
	       "max-cpus", &maxCPUs,
	       "max-wall-clock", &maxWallClock,
	       "help", &helpWanted
	       );
    } catch (Exception e) {
	writeln("Problem parsing command-line options.");
	writeln("Arguments not processed:");
	args = args[1 .. $]; // Dispose of program in first arg
	foreach (arg; args) writeln("   arg: ", arg);
	write(msg);
	exit(1);
    }
    if (helpWanted) {
	write(msg);
	exit(0);
    }
    if (jobName.length == 0) {
	writeln("Need to specify a job name.");
	write(msg);
	exit(1);
    }
    GlobalConfig.base_file_name = jobName;
    maxCPUs = min(max(maxCPUs, 1), totalCPUs); // don't ask for more than available
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    auto last_tindx = tindx_list[$-1];

    writefln("Initialising simulation from tindx: %d", last_tindx);
    init_simulation(last_tindx, 0, maxCPUs, maxWallClock);
    
    // save a copy of the original mesh
    foreach (blk; parallel(gasBlocks,1)) {
	blk.sync_vertices_to_underlying_grid(0);
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_underlying_grid(fileName);
    }
    
    // ----------------------------------------------------
    // 1.a store stencils used in forming the flow Jacobian
    //-----------------------------------------------------

    foreach (myblk; parallel(gasBlocks,1)) {
	construct_flow_jacobian_stencils(myblk);
    }
    
    // ------------------------------------------------------------
    // 1.b construct local flow Jacobian via Frechet derivatives
    // ------------------------------------------------------------
    size_t ndim = GlobalConfig.dimensions;
    // number of primitive variables
    bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
    size_t np;
    if (GlobalConfig.dimensions == 2) {
	if (with_k_omega) np = 6; // rho, momentum(x,y), total energy, tke, omega
	else np = 4; // rho, momentum(x,y), total energy
    }
    else if (GlobalConfig.dimensions == 3) {
	if (with_k_omega) np = 7; // rho, momentum(x,y,z), total energy, tke, omega
	else np = 5; // rho, momentum(x,y,z), total energy, tke, omega
    }

    foreach (myblk; parallel(gasBlocks,1)) {
	construct_flow_jacobian(myblk, ndim, np);
    }

    // ------------------------------------------------------------
    // 1.c construct global flow Jacobian transpose
    // ------------------------------------------------------------
    SMatrix globalJacobianT = new SMatrix();    
    // we must do this in serial
    size_t ia = 0;
    foreach (myblk; gasBlocks) {
	size_t nrows = myblk.cells.length*np;
	foreach(i; 0 .. nrows) {
	    globalJacobianT.aa ~= myblk.aa[i];
	    globalJacobianT.ja ~= myblk.ja[i];
	    globalJacobianT.ia ~= ia;
	    ia += myblk.aa[i].length;
	}
	// clear local Jacobian memory
	myblk.aa = [];
	myblk.ja = [];
    }
    globalJacobianT.ia ~= globalJacobianT.aa.length;
    // -----------------------------------------------------
    // 2. Form cost function sensitvity
    // -----------------------------------------------------
     // Analytically form dJdV by hand differentiation
    // cost function is defined as: J(Q) = 0.5*integral[0->l] (p-p*)^2
    double[] dJdV;
    double[] p_target;
    size_t ncells = gasBlocks[0].cells.length;
    size_t nvertices = gasBlocks[0].vertices.length;

    
    // target pressure distribution saved in file target.dat
    auto file = File("target.dat", "r");
    foreach(line; 0..ncells) {
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
    // 3. Solve for the adjoint variables 
    // -----------------------------------------------------
    // restarted-GMRES settings
    int maxInnerIters = 85;
    int maxOuterIters = 1000;
    int nIter = 0;
    double normRef = 0.0;
    double normNew = 0.0;
    double residTol = 1.0e-10;
    double[] psi0;
    double[] psiN;
    foreach(i; 0..np*ncells) psi0 ~= 1.0;
    foreach(i; 0..np*ncells) psiN ~= 1.0;
    double[] residVec;
    residVec.length = dJdV.length;
    foreach(i; 0..dJdV.length) dJdV[i] = -1.0 * dJdV[i];

    // compute ILU[0] for preconditioning
    SMatrix m = new SMatrix(globalJacobianT);
    auto M = decompILUp(m, 6);

    // compute reference norm
    multiply(globalJacobianT, psi0, residVec);
    foreach (i; 0..np*ncells) residVec[i] = dJdV[i] - residVec[i];
    foreach (i; 0..np*ncells) normRef += residVec[i]*residVec[i];
    normRef = sqrt(fabs(normRef));
    auto gws = GMRESWorkSpace(psi0.length, maxInnerIters);
    
    while (nIter < maxOuterIters) {
	// compute psi
	rpcGMRES(globalJacobianT, M, dJdV, psi0, psiN, maxInnerIters, residTol, gws);
	    
	// compute new norm
	normNew = 0.0;
	multiply(globalJacobianT, psiN, residVec);
	foreach (i; 0..np*ncells) residVec[i] = dJdV[i] - residVec[i];
	foreach (i; 0..np*ncells) normNew += residVec[i]*residVec[i];
	normNew = sqrt(fabs(normNew));
	
	writeln("iter = ", nIter, ", resid = ", normNew/normRef,
		", adjoint: rho = ", psiN[0], ", velx = ", psiN[1], ", vely = ", psiN[2], ", p = ", psiN[3]);
	nIter += 1;
	// tolerance check
	if (normNew/normRef < residTol) {
	    writeln("final residual: ", normNew/normRef);
	    break;
	}
	foreach(i; 0..np*ncells) psi0[i] = psiN[i];
    }
    
    double[] psi;
    psi.length = psiN.length;
    foreach(i; 0..np*ncells) psi[i] = psiN[i];
    
    foreach(i; 0 .. ncells) {
	FVCell cell = gasBlocks[0].cells[i];
	auto writer = format("%f %f %f %f %f \n", cell.pos[0].x, psi[i*np], psi[i*np+1], psi[i*np+2], psi[i*np+3]);
	append("e4_adjoint_vars.dat", writer);
    }

    // clear the global Jacobian from memory
    destroy(globalJacobianT);
    GC.minimize();
    // ---------------------------------------------------------------------
    // 4.a store stencils used in forming the mesh Jacobian 
    // ---------------------------------------------------------------------

    foreach (myblk; parallel(gasBlocks,1)) {
	construct_mesh_jacobian_stencils(myblk);
    }
    // ------------------------------------------------------------
    // 4.b construct local mesh Jacobian via Frechet derivatives
    // ------------------------------------------------------------

    foreach (myblk; parallel(gasBlocks,1)) {
	construct_mesh_jacobian(myblk, ndim, np);
    }
    // ------------------------------------------------------------
    // 4.c construct global mesh Jacobian transpose
    // ------------------------------------------------------------
    SMatrix globaldRdXT = new SMatrix();    
    // we must do this in serial
    ia = 0;
    foreach (myblk; gasBlocks) {
	size_t nrows = myblk.cells.length*np;
	foreach(i; 0 .. nrows) {
	    globaldRdXT.aa ~= myblk.aa[i];
	    globaldRdXT.ja ~= myblk.ja[i];
	    globaldRdXT.ia ~= ia;
	    ia += myblk.aa[i].length;
	}
	// clear local Jacobian memory
	myblk.aa = [];
	myblk.ja = [];
    }
    globaldRdXT.ia ~= globaldRdXT.aa.length;
    // -----------------------------------------------------
    // 5. form dX/dD via finite-differences
    // -----------------------------------------------------
    double[3] D = [0.07, 0.8, 3.8]; // array of design variables 
    auto D0 = D;
    size_t nvar = D.length; 
    Matrix dXdD_T;
    dXdD_T = new Matrix(nvar, nvertices*ndim);
    foreach (i; 0..D.length) {
	foreach (myblk; parallel(gasBlocks,1)) {
	    double dD = (abs(D[i]) + MU)*EPSILON;
	    D[i] = D0[i] + dD; 
	    auto meshPp = perturbMesh(myblk, D);
	    D[i] = D0[i] - dD; 
	    auto meshPm = perturbMesh(myblk, D);
	    foreach(j; 0..meshPm.length) {
		dXdD_T[i, j*ndim] = (meshPp[j][0] - meshPm[j][0])/(2.0*dD);
		dXdD_T[i, j*ndim+1] = (meshPp[j][1] - meshPm[j][1])/(2.0*dD);
	    }
	    D[i] = D0[i];
	}
    }

    // -----------------------------------------------------
    // 7. Compute adjoint gradients
    // -----------------------------------------------------

    // temp matrix multiplication
    Matrix tempMatrix;
    tempMatrix = new Matrix(nvar, ncells*np);
    for (int i = 0; i < nvar; i++) {
        for (int j = 0; j < np*ncells; j++) {
            tempMatrix[i,j] = 0;
            for (int k = 0; k < nvertices*ndim; k++) {
                tempMatrix[i,j] += dXdD_T[i,k]*globaldRdXT[k,j];
	    }
	}
    }
    double[3] grad;
    dot(tempMatrix, psi, grad);    
    writeln("adjoint gradients [dLdB, dLdC, dLdD] = ", grad);

    // clear the sensitivity matrices from memory
    destroy(globaldRdXT);
    destroy(dXdD_T);
    destroy(tempMatrix);
    GC.minimize();
    
    // -----------------------------------------------------
    // Finite difference verification
    // -----------------------------------------------------

    if (withFiniteDiffVerification) {
    
	double J1; double J0; double EPSILON0 = 1.0e-07;
	
	// read original grid in
	foreach (blk; gasBlocks) { 
	    //blk.init_grid_and_flow_arrays(make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz"));
	    ensure_directory_is_present(make_path_name!"grid-original"(0));
	    string gridFileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
	    blk.read_new_underlying_grid(gridFileName);
	    blk.sync_vertices_from_underlying_grid(0);
	}
	// perturb b +ve --------------------------------------------------------------------------------
	double del_b = EPSILON0;
	D[0] = D0[0] + del_b;
	J1 = finite_difference_grad(jobName, last_tindx, gasBlocks, p_target, D);
	// perturb b -ve --------------------------------------------------------------------------------
	D[0] = D0[0] - del_b;
	J0 = finite_difference_grad(jobName, last_tindx, gasBlocks, p_target, D);
	D[0] = D0[0];
	double grad_b = (J1 -J0)/(2.0*del_b);
	writeln("FD dLdB = ", grad_b, ", % error = ", abs((grad_b - grad[0])/grad_b * 100));
	
	// perturb c +ve --------------------------------------------------------------------------------
	double del_c = EPSILON0;
	D[1] = D0[1] + del_c;
	J1 = finite_difference_grad(jobName, last_tindx, gasBlocks, p_target, D);
	// perturb c -ve --------------------------------------------------------------------------------
	D[1] = D0[1] - del_c;
	J0 = finite_difference_grad(jobName, last_tindx, gasBlocks, p_target, D);
	D[1] = D0[1];
	double grad_c = (J1 -J0)/(2.0*del_c);
	writeln("FD dLdC = ", grad_c, ", % error = ", abs((grad_c - grad[1])/grad_c * 100));
	
	// perturb d +ve --------------------------------------------------------------------------------
	double del_d = EPSILON0;
	D[2] = D0[2] + del_d;
	J1 = finite_difference_grad(jobName, last_tindx, gasBlocks, p_target, D);
	// perturb d -ve --------------------------------------------------------------------------------
	D[2] = D0[2] - del_d;
	J0 = finite_difference_grad(jobName, last_tindx, gasBlocks, p_target, D);
	D[2] = D0[2];
	double grad_d = (J1 -J0)/(2.0*del_d);
	writeln("FD dLdD = ", grad_d, ", % error = ", abs((grad_d - grad[2])/grad_d * 100));
    }
    writeln("Simulation complete.");
}

double[][] perturbMesh(Block blk, double[] D) {
    // compute new nozzle surface, and delta
    SBlock sblk = cast(SBlock) blk;
    double y0; double y1;
    double[size_t] delta;
    double scale = 1.5;
    double b = D[0]; double c = D[1]; double d = D[2];
    double yo = 0.105;
    double a = yo - b*tanh(-d/scale);
    for (size_t i = sblk.imin; i <= sblk.imax+1; ++i ) {
	FVVertex vtx = sblk.get_vtx(i,sblk.jmax+1,0);
	y0 = vtx.pos[0].refy;
	y1 = a + b*tanh((c*vtx.pos[0].x-d)/scale);
	delta[i] = y1-y0;
    }
    double[][] meshP;
    meshP.length = sblk.vertices.length;
    // perturb mesh by linearly scaling the nozzle surface movement
    size_t vtxID = 0;
    for ( size_t j = sblk.jmin; j <= sblk.jmax+1; ++j ) {
	for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i) {
	    meshP[vtxID].length = 2;
	    double[2] point;
	    FVVertex vtx = sblk.get_vtx(i,j,0);
	    double jmax = to!double(sblk.jmax);
	    double jd = to!double(j)-2.0;
	    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
	    point = [vtx.pos[0].refx, vtx.pos[0].refy+p*delta[i]];
	    meshP[vtxID][0] = point[0];
	    meshP[vtxID][1] = point[1];
	    vtxID += 1;
	}
    }
    return meshP;
}

double finite_difference_grad(string jobName, int last_tindx, Block[] gasBlocks, double[] p_target, double[] D) {

    foreach (blk; parallel(gasBlocks,1)) {
	auto meshP = perturbMesh(blk, D);
	SBlock sblk = cast(SBlock) blk;
	size_t vtxID = 0;
	for ( size_t j = sblk.jmin; j <= sblk.jmax+1; ++j ) {
	    for ( size_t i = sblk.imin; i <= sblk.imax+1; ++i) {
		FVVertex vtx = sblk.get_vtx(i,j,0);
		vtx.pos[0].refy = meshP[vtxID][1];
		vtxID += 1;
	    }
	}
    }
    
    // save mesh
    foreach (blk; parallel(gasBlocks,1)) {
	blk.sync_vertices_to_underlying_grid(0);
	ensure_directory_is_present(make_path_name!"grid-perturb"(0));
	auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_underlying_grid(fileName);
    }
    
    // run simulation
    string command = "bash run-perturb.sh";
    auto output = executeShell(command);
    
    // read perturbed solution
    foreach (blk; parallel(gasBlocks,1)) {
	blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    }

    // compute cost function
    double J = 0.0;
    foreach (blk; parallel(gasBlocks,1)) {
	foreach (i, cell; blk.cells) {
	    J += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
	}
    }

    // read original grid in
    foreach (blk; parallel(gasBlocks,1)) { 
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

void construct_flow_jacobian(Block blk, size_t ndim, size_t np) {
    /++
     + computes and stores a block local transpose Jacobian in  Compressed
     Row Storage (CSR) format.     
     + initial method for efficient computation of the Jacobian -- predecessor
     to colourings.

     TODO: turbulence, 3D
     ++/

    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;
    
    // some data objects used in forming the Jacobian
    FVCell cellOrig; FVCell cellL; FVCell cellR;
    FVInterface ifaceOrig; FVInterface ifacePp; FVInterface ifacePm;
    double h; double diff;

    cellOrig = new FVCell(dedicatedConfig[blk.id]);
    ifaceOrig = new FVInterface(dedicatedConfig[blk.id], false);
    ifacePp = new FVInterface(dedicatedConfig[blk.id], false);
    ifacePm = new FVInterface(dedicatedConfig[blk.id], false);

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
	foreach(c; cell.jacobian_stencil) {
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
		    if (JacEntry != 0.0) {
			blk.aa[J] ~= JacEntry;
			blk.ja[J] ~= I; //J;
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
}

void construct_flow_jacobian_stencils(Block blk) {
    /++
     
     This stencil holds references to the cells effected by a 
     perturbation in the parent cell.
     
     For inviscid simulations, stencils are made up of the cells used in the
     reconstruction step for each of the cells interfaces.
     
     For viscous simulations, stencils are made up of the inviscid stencil, plus
     any cells that are additionally used in the viscous flux computations.
     
     NB. we need the stencils in cell id order, so that we can sequentially fill
     a row in the transposed Jacobian in Compressed Row Storage format.
     
     ++/
    foreach(c; blk.cells) {
	FVCell[] refs_ordered;
	FVCell[] refs_unordered;
	size_t[size_t] pos_array; // this is a dictionary that uses a cell id to reference the position of that cell in the unordered reference array
	size_t[] cell_ids;

	if (blk.myConfig.interpolation_order < 2) { 
	    // add the parent cell as the first reference
	    refs_unordered ~= c;
	    pos_array[c.id] = refs_unordered.length-1;
	    cell_ids ~= c.id;
	    foreach(f; c.iface) {
		if (f.left_cell.id != c.id && f.left_cell.id < ghost_cell_start_id) {
		    refs_unordered ~= f.left_cell;
		    pos_array[f.left_cell.id] = refs_unordered.length-1;
		    cell_ids ~= f.left_cell.id;
		}
		if (f.right_cell.id != c.id && f.right_cell.id < ghost_cell_start_id) {
		    refs_unordered ~= f.right_cell;
		    pos_array[f.right_cell.id] = refs_unordered.length-1;
		    cell_ids ~= f.right_cell.id;
		} else continue;
	    } // end foreach
	} // end if interpolation order < 2
	else { // higher-order interpolation
	    if (blk.grid_type == Grid_t.structured_grid) {
		// add the parent cell as the first reference
		refs_unordered ~= c;
		pos_array[c.id] = refs_unordered.length-1;
		cell_ids ~= c.id;
		//throw new Error("adjoint: 2nd order interpolation for structured_grid() not yet implemented");
		SBlock sblk = cast(SBlock) blk;
		FVCell c0;
		size_t[3] ijk = sblk.cell_id_to_ijk_indices(c.id);
		size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2]; 
		FVCell[] cells;
		cells ~= sblk.get_cell(i-2, j, k);
		cells ~= sblk.get_cell(i-1, j, k);
		cells ~= sblk.get_cell(i+1, j, k);
		cells ~= sblk.get_cell(i+2, j, k);
		cells ~= sblk.get_cell(i, j-2, k);
		cells ~= sblk.get_cell(i, j-1, k);
		cells ~= sblk.get_cell(i, j+1, k);
		cells ~= sblk.get_cell(i, j+2, k);
		foreach(cell; cells) {
		    if (cell.id < ghost_cell_start_id) {
			refs_unordered ~= cell;
			pos_array[cell.id] = refs_unordered.length-1;
			cell_ids ~= cell.id;
		    } else continue;
		}
	    }
	    else { // unstructured grid
		throw new Error("adjoint: 2nd order interpolation for unstructured_grid() not yet implemented");
	    }
	}
	// finally sort ids, and store sorted cell references
	cell_ids.sort();
	foreach(id; cell_ids) {
	    refs_ordered ~= refs_unordered[pos_array[id]];
	}
	c.jacobian_stencil ~= refs_ordered;
    }
}
 
void construct_mesh_jacobian_stencils(Block blk) {
    /++
     
     This stencil holds references to the cells effected by a 
     perturbation in the parent vertex.
     
     NB. we need the stencils in cell id order, so that we can sequentially fill
     a row in the transposed Jacobian in Compressed Row Storage format.
     
     ++/
    
    if (blk.myConfig.interpolation_order < 2) { // first-order
	foreach(i, vtx; blk.vertices) {
	    FVCell[] cell_refs;
	    foreach (cid; blk.cellIndexListPerVertex[vtx.id]) {
		cell_refs ~= blk.cells[cid];
	    }
	    vtx.jacobian_stencil ~= cell_refs;
	}
    } else { // higher-order interpolation
	if (blk.grid_type == Grid_t.structured_grid) {
	    SBlock sblk = cast(SBlock) blk;
	    foreach(vtx; blk.vertices) {
		FVCell[] refs_ordered;
		FVCell[] refs_unordered;
		size_t[size_t] pos_array; // this is a dictionary that uses a cell id to reference the position of that cell in the unordered reference array
		size_t[] cell_ids;
		foreach (cid; sblk.cellIndexListPerVertex[vtx.id]) {
		    size_t[3] ijk = sblk.cell_id_to_ijk_indices(cid);
		    size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2];  
		    FVCell[] cells;
		    cells ~= sblk.get_cell(i, j, k);
		    cells ~= sblk.get_cell(i-1, j, k);
		    cells ~= sblk.get_cell(i-2, j, k);
		    cells ~= sblk.get_cell(i+1, j, k);
		    cells ~= sblk.get_cell(i+2, j, k);
		    cells ~= sblk.get_cell(i, j-1, k);
		    cells ~= sblk.get_cell(i, j-2, k);
		    cells ~= sblk.get_cell(i, j+1, k);
		    cells ~= sblk.get_cell(i, j+2, k);
		    foreach(c; cells) {
			if (cell_ids.canFind(c.id) == false && c.id < ghost_cell_start_id) {
			    refs_unordered ~= c;
			    pos_array[c.id] = refs_unordered.length-1;
			    cell_ids ~= c.id;
			} else continue;
		    }
		}
		cell_ids.sort();
		foreach(id; cell_ids) {
		    refs_ordered ~= refs_unordered[pos_array[id]];
		}
		vtx.jacobian_stencil ~= refs_ordered;
	    }
	}
	else { // unstructured grid
	    throw new Error("adjoint: 2nd order interpolation for unstructured_grid() not yet implemented");
	}
    }
}

void construct_mesh_jacobian(Block blk, size_t ndim, size_t np) {
    /++
     + computes and stores a block local transpose Jacobian in  Compressed
     Row Storage (CSR) format.     
     + initial method for efficient computation of the Jacobian -- predecessor
     to colourings.

     TODO: turbulence, 3D
     ++/
    
    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;
    
    // some data objects used in forming the Jacobian
    FVCell cellOrig; FVCell cellL; FVCell cellR;
    FVInterface ifaceOrig; FVInterface ifacePp; FVInterface ifacePm;
    double h; double diff;

    cellOrig = new FVCell(dedicatedConfig[blk.id]);
    ifaceOrig = new FVInterface(dedicatedConfig[blk.id], false);
    ifacePp = new FVInterface(dedicatedConfig[blk.id], false);
    ifacePm = new FVInterface(dedicatedConfig[blk.id], false);

    foreach(i; 0..np*ncells) {
	blk.aa ~= [null];
	blk.ja ~= [null];
    }

    blk.clear_fluxes_of_conserved_quantities();
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
	foreach(ci, c; vtx.jacobian_stencil) {
	    size_t I, J; // indices in Jacobian matrix
	    double integral;
	    double volInv = 1.0 / c.volume[0];
	    for ( size_t ip = 0; ip < np; ++ip ) {
		I = c.id*np + ip; // row index
		for ( size_t jp = 0; jp < ndim; ++jp ) {
		    double JacEntry = 0.0;
		    // there are three contributions, dRdX = dRdF * dFdX + dRdA * dAdX + dRdV * dVdX
		    // 1. dRdF * dFdX ---------------------------
		    integral = 0.0;
		    J = vtx.id*ndim + jp; //vtx.id*nc + jc; // column index
		    foreach(fi, iface; c.iface) {
			integral -= c.outsign[fi] * iface.dFdU[ip][jp]*iface.area[0]; // gtl=0
		    }
		    //dRdX[I,J] = volInv * integral;
		    JacEntry += volInv * integral;
		    // 2. dRdA * dAdX ---------------------------
		    integral = 0.0;
		    double dAdX; double A0; double A1;
		    if (ip == 0 ) { // mass
			foreach(fi, iface; c.iface) {
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.mass*dAdX;
			}
		    }
		    else if (ip == 1) { // x-momentum
			foreach(fi, iface; c.iface) {
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.momentum.x*dAdX;
			}
		    }
		    else if (ip == 2) { // y-momentum
			foreach(fi, iface; c.iface){
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.momentum.y*dAdX;
			}
		    }
		    else if (ip == 3) { // total energy
			foreach(fi, iface; c.iface){
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.total_energy*dAdX;
			}
		    }
		    //dRdX[I,J] += volInv * integral;
		    JacEntry += volInv * integral;
		    // 3. dRdV * dVdX ---------------------------
		    double dVdX; double V0; double V1;
		    integral = 0.0;
		    if (ip == 0 ) { // mass
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.mass*iface.area[0]; }
		    }
		    else if (ip == 1) { // x-momentum
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.momentum.x*iface.area[0]; }
		    }
		    else if (ip == 2) { // y-momentum
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.momentum.y*iface.area[0]; }
		    }
		    else if (ip == 3) { // total energy
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.total_energy*iface.area[0]; }
		    }
		    if (jp == 0) mixin(computeCellVolumeSensitivity("pos[0].refx")); // x-dimension
		    else if (jp == 1) mixin(computeCellVolumeSensitivity("pos[0].refy")); // y-dimension
		    JacEntry -= volInv*volInv*integral * dVdX;
		    if (JacEntry != 0.0) {
			blk.aa[J] ~= JacEntry;
			blk.ja[J] ~= I; //J;
		    }
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
    } // end foreach vtx
}

void compute_perturbed_flux(Block blk, FVInterface iface, FVInterface ifaceP) {
    /++
     This method computes a perturbed flux at a given interface. 

     NB. The accuracy of the adjoint solver is largely dependent on how well this
     routine replicates the flow solver flux update procedures.
     ++/
    
    // pre-reconstrucion stage
    //if (iface.is_on_boundary) { // only apply bc's for cells on boundary
    blk.applyPreReconAction(0.0, 0, 0);  // assume sim_time = 0.0, gtl = 0, ftl = 0
    //}

    // Convective flux update
    if (blk.grid_type == Grid_t.structured_grid) {
	// we need to cast an SBlock here to reach some methods and data
	SBlock sblk = cast(SBlock) blk;
	size_t imin = sblk.imin; size_t imax = sblk.imax; size_t jmin = sblk.jmin; size_t jmax = sblk.jmax;
	foreach(faceIdentity, face; iface.left_cell.iface) { // use of left_cell is arbitrary -- could use right_cell
	    if (face == iface) {
		if (faceIdentity == 1 || faceIdentity == 3) { // east-facing faces 
		    // get iface ijk indices
		    size_t[3] ijk;
		    size_t i; size_t j; size_t k; 
		    if (iface.left_cell.id < ghost_cell_start_id) {
			ijk = sblk.cell_id_to_ijk_indices(iface.left_cell.id);
			if (faceIdentity == 1) i = ijk[0]+1;
			if (faceIdentity == 3) i = ijk[0];
			j = ijk[1]; k = ijk[2];
		    }
		    else {
			ijk = sblk.cell_id_to_ijk_indices(iface.right_cell.id);
			if (faceIdentity == 1) i = ijk[0]+1-1;
			if (faceIdentity == 3) i = ijk[0]-1;
			j = ijk[1]; k = ijk[2];
		    }
		    auto cL0 = sblk.get_cell(i-1,j,k); auto cL1 = sblk.get_cell(i-2,j,k);
		    auto cR0 = sblk.get_cell(i,j,k); auto cR1 = sblk.get_cell(i+1,j,k);
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
		    if ((i == imin) && (sblk.bc[Face.west].convective_flux_computed_in_bc == true)) continue;
		    if ((i == imax+1) && (sblk.bc[Face.east].convective_flux_computed_in_bc == true)) continue;
		    compute_interface_flux(sblk.Lft, sblk.Rght, iface, sblk.myConfig, sblk.omegaz);
		}
		else if (faceIdentity == 0 || faceIdentity == 2) { // north-facing faces 
		    // get iface ijk indices
		    size_t[3] ijk;
		    size_t i; size_t j; size_t k; 
		    if (iface.left_cell.id < ghost_cell_start_id) {
			ijk = sblk.cell_id_to_ijk_indices(iface.left_cell.id);
			if (faceIdentity == 0) j = ijk[1]+1;
			if (faceIdentity == 2) j = ijk[1];
			i = ijk[0]; k = ijk[2];
		    }
		    else {
			ijk = sblk.cell_id_to_ijk_indices(iface.right_cell.id);
			if (faceIdentity == 0) j = ijk[1]+1-1;
			if (faceIdentity == 2) j = ijk[1]-1;
			i = ijk[0]; k = ijk[2];
		    }
		    auto cL0 = sblk.get_cell(i,j-1,k); auto cL1 = sblk.get_cell(i,j-2,k);
		    auto cR0 = sblk.get_cell(i,j,k); auto cR1 = sblk.get_cell(i,j+1,k);
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
		    if ((j == jmin) && (sblk.bc[Face.south].convective_flux_computed_in_bc == true)) continue;
		    if ((j == jmax+1) && (sblk.bc[Face.north].convective_flux_computed_in_bc == true)) continue;
		    compute_interface_flux(sblk.Lft, sblk.Rght, iface, sblk.myConfig, sblk.omegaz);
		}
	    }
	}
    }
    else { // unstructured grid
	throw new Error("adjoint: 2nd order interpolation for unstructured_grid() not yet implemented");
    }
    
    //if (iface.is_on_boundary) {
    blk.applyPostConvFluxAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
	//}

    // copy perturbed flux
    ifaceP.copy_values_from(iface, CopyDataOption.all);
}

string computeFluxFlowVariableDerivativesAroundCell(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "h = (abs(cell.fs."~varName~") + MU) * EPSILON;";
    codeStr ~= "cellOrig.copy_values_from(cell, CopyDataOption.all);";
    codeStr ~= "foreach(stencilCell; cell.jacobian_stencil) {";
    codeStr ~= "foreach(iface; stencilCell.iface) {";
    codeStr ~= "ifaceOrig.copy_values_from(iface, CopyDataOption.all);";
    // ------------------ negative perturbation ------------------
    codeStr ~= "cell.fs."~varName~" -= h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(blk, iface, ifacePm);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
	codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(blk, iface, ifacePp);";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "diff = ifacePp.F.mass - ifacePm.F.mass;";
    codeStr ~= "iface.dFdU[0][" ~ posInArray ~ "] = diff/(2.0*h);";	    
    codeStr ~= "diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;";
    codeStr ~= "iface.dFdU[1][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;";
    codeStr ~= "iface.dFdU[2][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp.F.total_energy - ifacePm.F.total_energy;";
    codeStr ~= "iface.dFdU[3][" ~ posInArray ~ "] = diff/(2.0*h);";
    // ------------------ restore original values ------------------
    codeStr ~= "iface.copy_values_from(ifaceOrig, CopyDataOption.all);";
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    codeStr ~= "}";
    codeStr ~= "}";
    return codeStr;
}

string computeCellVolumeSensitivity(string varName)
{
    string codeStr;
    codeStr ~= "h = (abs(vtx."~varName~") + MU) * EPSILON;";
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
    codeStr ~= "h = (abs(vtx."~varName~") + MU) * EPSILON;";
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

string computeFluxMeshPointDerivativesAroundCell(string varName, string posInArray)
{
    string codeStr;
    codeStr ~= "h = (abs(vtx."~varName~") + MU) * EPSILON;";
    codeStr ~= "foreach (stencilCell; vtx.jacobian_stencil) {";
    codeStr ~= "foreach (iface; stencilCell.iface) { ";
    codeStr ~= "ifaceOrig.copy_values_from(iface, CopyDataOption.all);";
    // ------------------ negative perturbation ------------------
    codeStr ~= "vtx."~varName~" -= h;";
    codeStr ~= "foreach (cid; blk.cellIndexListPerVertex[vtx.id]) blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "foreach (fid; blk.faceIndexListPerVertex[vtx.id]) blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "compute_perturbed_flux(blk, iface, ifacePm);";
    codeStr ~= "vtx."~varName~" += h;";
    // ------------------ positive perturbation ------------------
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "foreach (cid; blk.cellIndexListPerVertex[vtx.id]) blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "foreach (fid; blk.faceIndexListPerVertex[vtx.id]) blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "compute_perturbed_flux(blk, iface, ifacePp);";
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
    codeStr ~= "iface.copy_values_from(ifaceOrig, CopyDataOption.all);";
    // ------------------ restore original values ------------------
    codeStr ~= "foreach (cid; blk.cellIndexListPerVertex[vtx.id]) blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "foreach (fid; blk.faceIndexListPerVertex[vtx.id]) blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "}";
    codeStr ~= "}";
    return codeStr;
}

