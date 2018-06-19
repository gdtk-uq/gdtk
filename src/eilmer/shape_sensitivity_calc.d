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

import steadystate_core;
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

    /* Simulation Initialisation */
    
    writeln("Eilmer shape sensitivity calculator:");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    version(mpi_parallel) {
        assert(0, "Adjoint solver is not MPI parallel, yet.");
    }

    string msg = "Usage:                                       Comment:\n";
    msg       ~= "e4ssc      [--job=<string>]                  name of job                                              \n";
    msg       ~= "           [--return-objective-function]     evaluates objective function                             \n";
    msg       ~= "           [--parameterise-surfaces]         parameterises all design surfaces using Bezier curves    \n";
    msg       ~= "           [--update-grid]                   perturbs grid according to updated Bezier control points \n";
    msg       ~= "           [--direct-method]                 evaluates sensitivities via direct complex step method   \n";
    msg       ~= "           [--adjoint-method]                evalautes sensitivities via adjoint method               \n";
    msg       ~= "           [--max-cpus=<int>]                defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine                                                                \n";
    msg       ~= "           [--max-wall-clock=<int>]          in seconds                                               \n";
    msg       ~= "           [--help]                          writes this message                                      \n";

    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        write(msg);
        exit(1);
    }

    string jobName = "";
    bool returnObjFcnFlag = false;
    bool parameteriseSurfacesFlag = false;
    bool updateGridFlag = false;
    bool directMethodFlag = false;
    bool adjointMethodFlag = true;
    bool verifyFlowJacobianFlag = false;
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;

    try {
        getopt(args,
               "job", &jobName,
               "return-objective-function", &returnObjFcnFlag,
               "parameterise-surfaces", &parameteriseSurfacesFlag,
               "update-grid", &updateGridFlag,
               "direct-method", &directMethodFlag,
               "adjoint-method", &adjointMethodFlag,
               "verify-flow-jacobian", &verifyFlowJacobianFlag,
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
    sort(tindx_list);
    auto last_tindx = tindx_list[$-1];
    writefln("Initialising simulation from tindx: %d", last_tindx);
    // TODO: MPI implementation (we can't assume threadsPerMPITask = 1) here.
    if (directMethodFlag) init_simulation(0, 0, maxCPUs, 1, maxWallClock); // initialise tindx=0
    else init_simulation(last_tindx, 0, maxCPUs, 1, maxWallClock);
    writefln("Initialising simulation from tindx: %d", last_tindx);

    // check some flag option compatibilities
    assert(directMethodFlag == adjointMethodFlag, "Error: Incompatible command line flags: direct-method & adjoint-method");
    assert(verifyFlowJacobianFlag != adjointMethodFlag, "Error: Incompatible command line flags: verify-flow-jacobian & adjoint-method");
    
    /* some global variables */    
    bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
    number[] designVars;
    designVars ~= to!number(0.347);
    designVars ~= to!number(0.6);
    designVars ~= to!number(3.8);
    size_t nDesignVars = designVars.length;

    /* Evaluate Sensitivities via Direct Complex Step */
    if (directMethodFlag) {
        compute_direct_complex_step_derivatives(jobName, last_tindx, maxCPUs, designVars);
        return;        
    }

    /* Evaluate Sensitivities via Discrete Adjoint Method */
    if (adjointMethodFlag) {

        size_t nPrimitive; 
        if (GlobalConfig.dimensions == 2) nPrimitive = 4;  // density, velocity(x,y), pressure
        else nPrimitive = 5;                               // density, velocity(x,y,z), pressure
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            // make sure ghost cells are filled before proceeding...
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
            
            myblk.JlocT = new SMatrix!number();
            local_flow_jacobian_transpose(myblk.JlocT, myblk, nPrimitive, myblk.myConfig.interpolation_order);
        }
              
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.P = new SMatrix!number();
            local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 0); // orderOfJacobian=0
            decompILU0(myblk.P);
        }
        
        // Surface intergal objective functions can be computed in parallel with a reduction process across blocks to gather the final value,
        // however the sensitivity w.r.t to primitive variables cannot be computed in parallel, since we are only computing a single objective calue (for example drag).
        // TODO: think about how this will effect user-defined objective functions
        foreach (myblk; localFluidBlocks) {
            form_objective_function_sensitivity(myblk, nPrimitive);
        }
        
        // solve the adjoint system
        rpcGMRES_solve(nPrimitive);
        
        // Write out adjoint variables for visualisation 
        foreach (myblk; localFluidBlocks) {
            //writeln(myblk.psi);
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
        
        // objective function sensitivity w.r.t. design variables
        number[] g;
        g.length = nDesignVars;
        foreach (myblk; localFluidBlocks) {
            size_t nLocalCells = myblk.cells.length;
            myblk.rT = new Matrix!number(nDesignVars, nLocalCells*nPrimitive);
        }
        
        compute_design_variable_partial_derivatives(designVars, g, nPrimitive, with_k_omega);
        
        // ---------------------
        // COMPUTE SENSITIVITIES
        // ---------------------
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.rTdotPsi.length = nDesignVars;
            dot(myblk.rT, myblk.psi, myblk.rTdotPsi);
            //writeln(myblk.rT);
        }
        
        number[] adjointGradients;
        adjointGradients.length = nDesignVars;
        foreach (myblk; localFluidBlocks) adjointGradients[] = myblk.rTdotPsi[];
        adjointGradients[] = g[] - adjointGradients[];
        
        foreach ( i; 0..adjointGradients.length) writef("gradient for variable %d: %.16e \n", i+1, adjointGradients[i]);
        
        // clear some expensive data structures from memory
        foreach (myblk; parallel(localFluidBlocks,1)) {
            destroy(myblk.rT);
            GC.minimize();
        }
        
        writeln("Simulation complete.");
    }

    /* Verify Flow Jacobian via Frechet Derivative */
    if (verifyFlowJacobianFlag) {

        /* Construct Flow Jacobian */
        size_t nPrimitive; 
        if (GlobalConfig.dimensions == 2) nPrimitive = 4;  // density, velocity(x,y), pressure
        else nPrimitive = 5;                               // density, velocity(x,y,z), pressure
        
        foreach (myblk; parallel(localFluidBlocks,1)) {
            // make sure ghost cells are filled before proceeding...
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
            
            myblk.JlocT = new SMatrix!number();
            local_flow_jacobian_transpose(myblk.JlocT, myblk, nPrimitive, myblk.myConfig.interpolation_order);
        }
        
        foreach (blk; localFluidBlocks) {

            // compute directional vector
            number[] v;
            v.length = blk.cells.length*nPrimitive;
            
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                v[i] = i+1; // theoretically can be any vector.
            }
            
            // normalise directional vector
            number norm = 0.0;
            foreach( i; 0..v.length) {
                norm += v[i]*v[i];
            }
            norm = sqrt(norm);
            foreach( i; 0..v.length) {
                v[i] = v[i]/norm;
            }

            // result vectors
            number[] p1;
            p1.length = blk.cells.length*nPrimitive;
            number[] p2;
            p2.length = blk.cells.length*nPrimitive;
            
            // explicit multiplication of Jv (note we want to nultiply the Jacobian NOT transpose Jacobian by v)
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                p1[i] = 0.0;
                foreach ( j; 0..blk.cells.length*nPrimitive) {
                    p1[i] += blk.JlocT[j,i]*v[j];
                }
            }
            
            // Frechet derivative of Jv
            evalJacobianVecProd(blk, nPrimitive, v, p2);

            string fileName = "flow_jacobian_test.output";
            auto outFile = File(fileName, "w");
            foreach( i; 0..v.length ) {
                outFile.writef("%.16e    %.16f    %.16f \n", (p1[i]-p2[i])/p2[i], p1[i], p2[i]);
            }
                        
            assert(approxEqualNumbers(p2, p1, 1.0e-14), "Flow Jacobian Test: FAILED.");
            
        }

        writeln("Flow Jacobian Test: PASSED.");
        return;
    }

    /* Program should terminate before here */
    assert(0, "Oh dear. The Eilmer Shape Sensitivity Calculator executed correctly, however, nothing meaningful has been computed. Please check command line flags.");
}
