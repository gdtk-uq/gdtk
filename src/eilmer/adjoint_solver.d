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
immutable double EPSILON = 1.0e-06;
immutable double ESSENTIALLY_ZERO = 1.0e-15;

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
	ensure_directory_is_present(make_path_name!"grid-original"(0));
	auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
	blk.write_grid(fileName, 0.0, 0);
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
	construct_residual_jacobian(myblk, ndim, np);
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
		if (i == vbi) {
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j);
		    double entry = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0)); 
		    //writeln(i, ", ", j, ", ", entry, ", ", 2*vi+1, ", ", 2*vbi+1);
		    dXdXb[vi*ndim+1,vbi*ndim+1] = entry;
		}
	    }
	}
    }
    // sensitivity of surface mesh points to movements of design variables
    Matrix dXbdD;
    dXbdD = new Matrix(nsurfnodes*ndim, nvar);
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
	    dXbdD[vi*ndim,0] = 0.0;
	    dXbdD[vi*ndim,1] = 0.0;
	    dXbdD[vi*ndim,2] = 0.0;
	    dXbdD[vi*ndim+1,0] = tanh(d/scale) + tanh((c*vtx.pos[0].x - d)/scale);
	    dXbdD[vi*ndim+1,1] = b*vtx.pos[0].x*(-pow(tanh((c*vtx.pos[0].x - d)/scale),2) + 1)/scale;
	    dXbdD[vi*ndim+1,2] = b*(-pow(tanh(d/scale),2) + 1)/scale - b*(-pow(tanh((c*vtx.pos[0].x - d)/scale),2) + 1)/scale;
	}
    } 
    // sensitivity of mesh points to movements of the design variables
    Matrix dXdD;
    dXdD = new Matrix(nvertices*ndim, nvar);
    //writeln(dXbdD);
    //writeln("-------");
    //writeln(dXdXb);
    dot(dXdXb, dXbdD, dXdD);
    //writeln("-------");
    //writeln(dXdD);
    // compute transposes
    Matrix dXdD_T; //Matrix dRdX_T;
    dXdD_T = transpose(dXdD);
    //dRdX_T = transpose(dRdX);
    // temp matrix multiplication
    Matrix tempMatrix;
    tempMatrix = new Matrix(nvar, ncells*np);
    //matrixMultiply(dXdD_T, dRdX_T, tempMatrix);
    for (int i = 0; i < nvar; i++) {
        for (int j = 0; j < np*ncells; j++) {
            tempMatrix[i,j] = 0;
            for (int k = 0; k < nvertices*ndim; k++) {
                tempMatrix[i,j] += dXdD_T[i,k]*globaldRdXT[k,j];
	    }
	}
    }
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
		//writeln(i, ", ", j, ", ", vtx.pos[0].x, ", ", vtx.pos[0].y);
		y_old = vtx.pos[0].refy;
		if (j == sblk.jmax+1) {
		    double yo = 0.105;
		    double a = yo - b*tanh(-d/scale);
		    double y_new = a + b*tanh((c*vtx.pos[0].x-d)/scale);
		    vtx.pos[0].refy = y_new;
		    delta = y_new - y_old;
		    //writeln(j, ", ", sblk.jmax, ", ", delta, ", ", y_old, ", ", y_new);
		} else {
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j)-2.0;
		    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
		    //writeln(i, ", ", j, ", ", jmax, ", ", jd, ", ", p);
		    vtx.pos[0].refy += p*delta;
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
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j)-2.0;
		    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
		    vtx.pos[0].refy += p*delta;
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
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j)-2.0;
		    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
		    vtx.pos[0].refy += p*delta;
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
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j)-2.0;
		    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
		    vtx.pos[0].refy += p*delta;
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
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j)-2.0;
		    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
		    vtx.pos[0].refy += p*delta;
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
		    double jmax = to!double(sblk.jmax);
		    double jd = to!double(j)-2.0;
		    double p = 1.0*(1.0 - ((jmax-1.0) - (jd))/(jmax-1.0));
		    vtx.pos[0].refy += p*delta;
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

void construct_residual_jacobian(Block blk, size_t ndim, size_t np) {
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
     + This stencil holds references to the cells effected by a 
     perturbation in the parent cell.
     + For inviscid simulations, stencils are made up of the cells used in the
     reconstruction step for each of the cells interfaces.
     + For viscous simulations, stencils are made up of the inviscid stencil, plus
     any cells that are additionally used in the viscous flux computations.
     + NB. we need the stencils in cell id order, so that we can sequentially fill
     a row in the transposed Jacobian in Compressed Row Storage format.

     TODO: high order interpolation
     ++/
    foreach(c; blk.cells) {
	FVCell[] refs_ordered;
	FVCell[] refs_unordered;
	size_t[size_t] pos_array;
	// pos_array is a dictionary which uses the cell id
	// as an index to retrieve the array position in refs_unordered
	size_t[] cell_ids;

	// add the parent cell as the first reference
	refs_unordered ~= c;
	pos_array[c.id] = refs_unordered.length-1;
	cell_ids ~= c.id;

	// now fill out the rest of rest stencil
	if (blk.myConfig.interpolation_order < 2) { 
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
		}
	    }
	}
	else { // higher-order interpolation
	    if (blk.grid_type == Grid_t.structured_grid) {
		//throw new Error("adjoint: 2nd order interpolation for structured_grid() not yet implemented");
		SBlock sblk = cast(SBlock) blk;
		FVCell c0;
		size_t[3] ijk = sblk.cell_id_to_ijk_indices(c.id);
		size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2]; 
		c0 = sblk.get_cell(i-2, j, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i-1, j, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i+1, j, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i+2, j, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i, j-2, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i, j-1, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i, j+1, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
		}
		c0 = sblk.get_cell(i, j+2, k);
		if (c0.id < ghost_cell_start_id) {
		    refs_unordered ~= c0;
		    pos_array[c0.id] = refs_unordered.length-1;
		    cell_ids ~= c0.id;
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
     + This stencil holds references to the cells effected by a 
     perturbation in the parent vertice.
     + NB. we need the stencils in cell id order, so that we can sequentially fill
     a row in the transposed Jacobian in Compressed Row Storage format -- fortuitosly
     the cellIndexListPerVertex is already in the correct order.
     ++/
    
    foreach(i, vtx; blk.vertices) {
	FVCell[] cell_refs;
	foreach (cid; blk.cellIndexListPerVertex[vtx.id]) {
	    cell_refs ~= blk.cells[cid];
	}
	vtx.jacobian_stencil ~= cell_refs;
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
    foreach(i; 0..np*ncells) {
	blk.aa ~= [null];
	blk.ja ~= [null];
    }

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
	foreach(c; vtx.jacobian_stencil) {
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
		    if (ip == 0 ) {
			foreach(fi, iface; c.iface) {
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.mass*dAdX;
			}
		    }
		    else if (ip == 1) {
			foreach(fi, iface; c.iface) {
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.momentum.x*dAdX;
			}
		    }
		    else if (ip == 2) {
			foreach(fi, iface; c.iface){
			    if (jp == 0) mixin(computeInterfaceAreaSensitivity("pos[0].refx"));  // x-dimension
			    else if (jp == 1) mixin(computeInterfaceAreaSensitivity("pos[0].refy")); // y-dimension
			    integral -= c.outsign[fi]*iface.F.momentum.y*dAdX;
			}
		    }
		    else if (ip == 3) {
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
		    if (ip == 0 ) {
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.mass*iface.area[0]; }
		    }
		    else if (ip == 1) {
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.momentum.x*iface.area[0]; }
		    }
		    else if (ip == 2) {
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.momentum.y*iface.area[0]; }
		    }
		    else if (ip == 3) {
			foreach(fi, iface; c.iface) { integral -= c.outsign[fi]*iface.F.total_energy*iface.area[0]; }
		    }
		    if (jp == 0) mixin(computeCellVolumeSensitivity("pos[0].refx")); // x-dimension
		    else if (jp == 1) mixin(computeCellVolumeSensitivity("pos[0].refy")); // y-dimension
		    //dVdX = (V1-V0)/(2*h);
		    //dRdX[I,J] -= volInv*volInv*integral * dVdX;
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
    } // end foreach cell
}

void compute_perturbed_flux(Block blk, FVInterface iface, FVInterface ifaceP) {
    FVCell cellL;
    FVCell cellR;
    // apply cell effect bcs ------------------
    if (iface.is_on_boundary) {
	blk.applyPreReconAction(0.0, 0, 0);  // assume sim_time = 0.0, gtl = 0, ftl = 0
    }
    // compute interflux ------------------
    cellR = iface.right_cell;
    cellL = iface.left_cell;
    blk.Lft.copy_values_from(cellL.fs);
    blk.Rght.copy_values_from(cellR.fs);
    compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);
    // apply interface effect bcs ------------------
    if (iface.is_on_boundary) {
	blk.applyPostConvFluxAction(0.0, 0, 0);
    }
    ifaceP.copy_values_from(iface, CopyDataOption.all);
}

string computeFluxFlowVariableDerivativesAroundCell(string varName, string posInArray, bool includeThermoUpdate)
{
    string codeStr;
    codeStr ~= "cellOrig = new FVCell(dedicatedConfig[blk.id]);";
    codeStr ~= "ifaceOrig = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "ifacePp = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "ifacePm = new FVInterface(dedicatedConfig[blk.id], false);";
    codeStr ~= "h = cell.fs."~varName~" * EPSILON + EPSILON;";
    codeStr ~= "cellOrig.copy_values_from(cell, CopyDataOption.all);";
    codeStr ~= "foreach(iface; cell.iface) {";
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
    return codeStr;
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
    codeStr ~= "iface.right_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.left_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "compute_perturbed_flux(blk, iface, ifacePm);";
    codeStr ~= "vtx."~varName~" += h;";
    // ------------------ positive perturbation ------------------
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "iface.right_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.left_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
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
    codeStr ~= "iface.right_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.left_cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "iface.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);";
    codeStr ~= "}";

    return codeStr;
}

