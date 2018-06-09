/** shape_sensitivity_calc.d
 * 
 * Eilmer4 shape sensitivity calculator top-level function.
 *
 * Note: This is the first attempt at 'production' code.
 * Some test implementations began on 2017-09-18.
 *
 * Date: 2018-10-02 -- implemented shared memory parallelism.
 *
 * Author: Kyle D.
**/

import core.stdc.stdlib : exit;
import core.memory;
import std.stdio;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;

import nm.smla;
import nm.bbla;
import nm.complex;
import nm.number;

import shape_sensitivity_core;
import fileutil;
import globalconfig;
import postprocess;
import simcore;
import geom;
import globaldata;


string adjointDir = "adjoint";
void init_adjoint_dir()
{
    ensure_directory_is_present(adjointDir);
}

void main(string[] args) {

    // -------------------------------------------
    // Simulation Initialisation
    // -------------------------------------------       

    init_adjoint_dir();

    writeln("Eilmer shape sensitivity calculator:");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    version(mpi_parallel) {
        assert(0, "Adjoint solver is not MPI parallel, yet.");
    }

    string msg = "Usage:                                       Comment:\n";
    msg       ~= "e4ssc      [--job=<string>]               name of job\n";
    msg       ~= "           [--return-objective-function]             \n";
    msg       ~= "           [--parameterise-surfaces                  \n";
    msg       ~= "           [--grid-update]                           \n";
    msg       ~= "           [--direct-method                          \n";
    msg       ~= "           [--verification                           \n";
    msg       ~= "           [--max-cpus=<int>]              defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine               \n";
    msg       ~= "           [--max-wall-clock=<int>]        in seconds\n";
    msg       ~= "           [--help]               writes this message\n";

    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        write(msg);
        exit(1);
    }

    string jobName = "";
    bool returnObjFcnFlag = false;
    bool parameteriseSurfacesFlag = false;
    bool gridUpdateFlag = false;
    bool directMethodFlag = false;
    bool verificationFlag = false;
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;

    try {
        getopt(args,
               "job", &jobName,
               "return-objective-function", &returnObjFcnFlag,
               "parameterise-surfaces", &parameteriseSurfacesFlag,
               "grid-update", &gridUpdateFlag,
               "direct-method", &directMethodFlag,
               "verification", &verificationFlag,
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

    // read simulation details to initialise stored simulation
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    auto last_tindx = tindx_list[$-1];
    writefln("Initialising simulation from tindx: %d", last_tindx);
    // TODO: MPI implementation (we can't assume threadsPerMPITask = 1)
    if (directMethodFlag) init_simulation(0, 0, maxCPUs, 1, maxWallClock);
    else init_simulation(last_tindx, 0, maxCPUs, 1, maxWallClock);
    
    // set some global config values specified in user input lua script
    GradientMethod gradientMethod = GlobalConfig.sscOptions.gradientMethod;
    // finite-difference perturbation parameters
    immutable double EPSILON = GlobalConfig.sscOptions.epsilon; // flow Jacobian
    immutable double MU = GlobalConfig.sscOptions.mu;           // flow Jacobian
    immutable double ETA = GlobalConfig.sscOptions.eta;         // residual sensitivity
    immutable double DELTA = GlobalConfig.sscOptions.delta;     // finite-difference gradient
    // bezier curve parameters
    double bezierCurveFitTol = GlobalConfig.sscOptions.tolBezierCurveFit;
    int bezierCurveFitMaxSteps = GlobalConfig.sscOptions.maxStepsBezierCurveFit;
    // order of flow Jacobian used in preconditioning
    int orderOfPreconditioningMatrix = 0;
    bool inviscidPreconditioningMatrix = true;
    bool viscousConfigSave = GlobalConfig.viscous; 
        
    // some global variables
    //size_t nDesignVars = 0;
    Vector3[] designVars;
    bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
    
    // Initialise Lua state for calling user-defined objective function.
    // TODO: complete user-defined objective functions
    //initLuaStateForUserDefinedObjFunc();
    
    // -----------------------------
    // -----------------------------
    // GEOMETRY PARAMETERISATION
    // -----------------------------
    // -----------------------------
    if (parameteriseSurfacesFlag) {
        parameterise_design_surfaces(designVars, bezierCurveFitTol, bezierCurveFitMaxSteps);
        writeDesignVarsToDakotaFile(designVars, jobName);
        writeBezierDataToFile();
        return; // --parameterise-surfaces complete
    }
    
    
    // -----------------------------
    // -----------------------------
    // Grid Update
    // -----------------------------
    // -----------------------------
    if (gridUpdateFlag) {
	// fill-in bezier curve with previous bezier points
        readBezierDataFromFile(designVars);
	// update bezier curve with new design variables
	readDesignVarsFromDakotaFile(designVars);
	gridUpdate(false, false, designVars, 1, jobName); // gtl=1
        return; // --grid-update complete
    }
    
    // -----------------------------
    // -----------------------------
    // OBJECTIVE FUNCTION EVALUATION
    //------------------------------
    // -----------------------------
    number objFnEval;
    if (returnObjFcnFlag) {
        objFnEval = objective_function_evaluation();
        writeln("objective fn evaluation: ", objFnEval);
        write_objective_fn_to_file("results.out", objFnEval);
        
        return; // --return-objective-function
    }
    

    // -------------------------------
    // -------------------------------
    // DIRECT COMPLEX STEP DERIVATIVES
    // -------------------------------
    // -------------------------------

    if (directMethodFlag) {
        parameterise_design_surfaces(designVars, bezierCurveFitTol, bezierCurveFitMaxSteps);
        compute_direct_complex_step_derivatives(jobName, last_tindx, maxCPUs, designVars, EPSILON);
        return; //        
    }

    // ----------------------------
    // ----------------------------
    // SHAPE SENSITIVITY CALCULATOR
    // -----------------------------
    // ----------------------------

    // --------------
    // ADJOINT SOLVER
    //---------------
    size_t nPrimitive; 
    if (GlobalConfig.dimensions == 2) nPrimitive = 4; // density, velocity(x,y), pressure
    else nPrimitive = 5;                               // density, velocity(x,y,z), pressure
    
    
    // TODO: MPI implementation (can't rely on this method for computing total number of cells)
    size_t nGlobalCells = 0;
    foreach ( myblk; localFluidBlocks) {
        nGlobalCells += myblk.cells.length;
    }

    // we make sure that the ghost cells are filled here (we don't want to update ghost cells during the flow Jacobian formation)
    // TODO: currently the code works for second order interpolation for internal interfaces, with first order interpolation at boundaries.
    // Should think about how to extend this to second order at block boundaries (i.e. manage the copying of gradients correctly).
    exchange_ghost_cell_boundary_data(0.0, 0, 0); // pseudoSimTime = 0.0; gtl = 0; ftl = 0
    
    foreach (myblk; parallel(localFluidBlocks,1)) {
        myblk.applyPreReconAction(0.0, 0, 0);
        myblk.applyPostConvFluxAction(0.0, 0, 0);
        myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
        myblk.applyPostDiffFluxAction(0.0, 0, 0);
        
        myblk.JlocT = new SMatrix!number();
        form_local_flow_jacobian_block(myblk.JlocT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPSILON, MU);
        jacobian_bndary_correction(myblk, myblk.JlocT, nPrimitive, EPSILON, MU, myblk.myConfig.interpolation_order);        
    }
    
    foreach (myblk; parallel(localFluidBlocks,1)) {
        form_external_flow_jacobian_block_phase0(myblk, nPrimitive, myblk.myConfig.interpolation_order, EPSILON, MU); // orderOfJacobian=interpolation_order
    }

    foreach (myblk; parallel(localFluidBlocks,1)) {
        myblk.JextT = new SMatrix!number();
        form_external_flow_jacobian_block_phase1(myblk.JextT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPSILON, MU); // orderOfJacobian=interpolation_order
    }
    
    if (inviscidPreconditioningMatrix)
        { GlobalConfig.viscous = false; }

    foreach (myblk; parallel(localFluidBlocks,1)) {
        myblk.P = new SMatrix!number();
        if (orderOfPreconditioningMatrix == 0) { // block-diagonal of a first order flow Jacobian (TODO: something isn't quite right for the new parallel implementation)
            form_local_flow_jacobian_block(myblk.P, myblk, nPrimitive, 0, EPSILON, MU); // orderOfJacobian=0
            jacobian_bndary_correction(myblk, myblk.P, nPrimitive, EPSILON, MU, 0);
        }
        else { // first order flow Jacobian {
            form_local_flow_jacobian_block(myblk.P, myblk, nPrimitive, 1, EPSILON, MU); // orderOfJacobian=1
            jacobian_bndary_correction(myblk, myblk.P, nPrimitive, EPSILON, MU, 1);
        }
        decompILU0(myblk.P);
    }
    
    if (inviscidPreconditioningMatrix)
        { GlobalConfig.viscous = viscousConfigSave; }
    
    // Surface intergal objective functions can be computed in parallel with a reduction process across blocks to gather the final value,
    // however the sensitivity w.r.t to primitive variables cannot be computed in parallel, since we are only computing a single objective calue (for example drag).
    // TODO: think about how this will effect user-defined objective functions
    foreach (myblk; localFluidBlocks) {
        form_objective_function_sensitivity(myblk, nPrimitive, EPSILON, MU);
    }

    // solve the adjoint system
    rpcGMRES_solve(nPrimitive);

    // Write out adjoint variables for visualisation 
    foreach (myblk; localFluidBlocks) {
        write_adjoint_variables_to_file(myblk, nPrimitive, jobName);
    }

    // clear some expensive data structures from memory
    foreach (myblk; parallel(localFluidBlocks,1)) {
        destroy(myblk.JlocT);
        destroy(myblk.JextT);
        destroy(myblk.P);
        GC.minimize();
    }

    
    // ------------------------------------------------------
    // RESIDUAL/OBJECTIVE SENSITIVITY W.R.T. DESIGN VARIABLES
    // ------------------------------------------------------
    if (verificationFlag) parameterise_design_surfaces(designVars, bezierCurveFitTol, bezierCurveFitMaxSteps);    
    else { 
	readBezierDataFromFile(designVars);
        readDesignVarsFromDakotaFile(designVars);
    }
    
    // objective function sensitivity w.r.t. design variables
    number[] g;
    g.length = designVars.length;
    foreach (myblk; localFluidBlocks) {
        size_t nLocalCells = myblk.cells.length;
        size_t nDesignVars = designVars.length;
        myblk.rT = new Matrix!number(nDesignVars, nLocalCells*nPrimitive);
    }

    compute_design_variable_partial_derivatives(designVars, g, nPrimitive, with_k_omega, ETA);

    // ---------------------
    // COMPUTE SENSITIVITIES
    // ---------------------
    foreach (myblk; parallel(localFluidBlocks,1)) {
        size_t nDesignVars = designVars.length;
        myblk.rTdotPsi.length = nDesignVars;
        dot(myblk.rT, myblk.psi, myblk.rTdotPsi);
    }
    number[] adjointGradients;
    adjointGradients.length = designVars.length;
    adjointGradients[] = to!number(0.0);
    foreach (myblk; localFluidBlocks) adjointGradients[] += myblk.rTdotPsi[];
    adjointGradients[] = g[] - adjointGradients[];
    foreach ( i; 0..adjointGradients.length) writef("gradient for variable %d: %.16f \n", i+1, adjointGradients[i]);
    write_gradients_to_file("results.out", adjointGradients);
    // clear some expensive data structures from memory
    foreach (myblk; parallel(localFluidBlocks,1)) {
        destroy(myblk.rT);
        GC.minimize();
    }
    
    writeln("Simulation complete.");
}
