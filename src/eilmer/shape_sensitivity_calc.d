 /** shape_sensitivity_calc.d
 * 
 * Eilmer4 shape sensitivity calculator top-level function.
 *
 * This file coordinates all the functions required for the Eilmer4 sensitivity analysis code, including:
 *     + complex-variable discrete adjoint solver
 *     + direct differentiation via complex-variable finitie difference solver (used for verification)
 *     + surface parameterisation
 *     + grid perturbation
 *     + objective function evaluation
 *
 * Author: Kyle D.
**/

// standard library
import core.stdc.stdlib : exit;
import core.memory;
import std.stdio;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.file;

// numerical methods
import nm.smla;
import nm.bbla;
import nm.complex;
import nm.number;

// eilmer
import fvcell;
import bc;
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
    msg       ~= "           [--adjoint-verification           flag for verifying the adjoint sensitivities             \n";
    msg       ~= "           [--verify-primitive-jacobian]     test the formation of the primtive Jacobian              \n";
    msg       ~= "           [--verify-conservative-jacobian   test the formation of the conservative Jacobian          \n";
    msg       ~= "           [--verify-sss-preconditioner      test the formation of the steady-state preconditioner    \n";
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
    bool adjointVerificationFlag = false;
    bool verifyPrimitiveJacobianFlag = false;
    bool verifyConservativeJacobianFlag = false;
    bool verifySSSPreconditionerFlag = false;
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
               "adjoint-verification", &adjointVerificationFlag,
               "verify-primitive-jacobian", &verifyPrimitiveJacobianFlag,
               "verify-conservative-jacobian", &verifyConservativeJacobianFlag,
               "verify-sss-preconditioner", &verifySSSPreconditionerFlag,
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
    //
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    sort(tindx_list);
    auto last_tindx = tindx_list[$-1];
    writefln("Initialising simulation from tindx: %d", last_tindx);
    // TODO: MPI implementation (we can't assume threadsPerMPITask = 1) here.
    if (directMethodFlag) init_simulation(0, 0, maxCPUs, 1, maxWallClock); // initialise tindx=0
    else init_simulation(last_tindx, 0, maxCPUs, 1, maxWallClock);

    // check some flag option compatibilities
    if (verifyPrimitiveJacobianFlag || verifyConservativeJacobianFlag)
        assert(verifyPrimitiveJacobianFlag != adjointMethodFlag || verifyConservativeJacobianFlag != adjointMethodFlag, "Error: Incompatible command line flags: verify-flow-jacobian & adjoint-method");
    else if (verifySSSPreconditionerFlag)
        assert(verifySSSPreconditionerFlag != adjointMethodFlag, "Error: Incompatible command line flags: verify-flow-jacobian & adjoint-method");
    else
        assert(directMethodFlag != adjointMethodFlag, "Error: Incompatible command line flags: direct-method & adjoint-method");
    
    /* initilise some global variables */    

    Vector3[] designVars;
    version(complex_numbers) number EPS = complex(0.0, GlobalConfig.sscOptions.epsilon);
    else number EPS = -1.0;
    assert(EPS != -1.0, "Error: complex step size incorrectly set");
    writeln("/* Complex Step Size: ", EPS.im, ", i */");
    with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
    

    /* Geometry parameterisation */
    
    if (parameteriseSurfacesFlag) {
        fit_design_parameters_to_surface(designVars);
        if (!adjointVerificationFlag) { 
            writeDesignVarsToDakotaFile(designVars, jobName);
        }
        return; // --parameterise-surfaces complete
    }

    /* update grid */

    if (updateGridFlag) {
	// fill-in bezier curve with previous bezier points
        readBezierDataFromFile(designVars);
	// update bezier curve with new design variables
	readDesignVarsFromDakotaFile(designVars);
	gridUpdate(designVars, 1, updateGridFlag, jobName); // gtl=1
        return; // --grid-update complete
    }

    /* objective function evaluation */

    number objFnEval;
    if (returnObjFcnFlag) {
        objFnEval = objective_function_evaluation();
        writeln("objective fn evaluation: ", objFnEval);
        write_objective_fn_to_file("results.out", objFnEval);
        
        return; // --return-objective-function
    }
    
    /* Evaluate Sensitivities via Direct Complex Step */

    if (directMethodFlag) {
        readBezierDataFromFile(designVars);
        compute_direct_complex_step_derivatives(jobName, last_tindx, maxCPUs, designVars, EPS);
        return;        
    }

    /* Evaluate Sensitivities via Discrete Adjoint Method */

    if (adjointMethodFlag) {
	
	// store global ids
	foreach (myblk; parallel(localFluidBlocks,1)) {
	    foreach (cell; myblk.cells) {
		cell.global_id = to!int(myblk.globalCellId(to!size_t(cell.id)));
		//writeln(cell.global_id, ", ", cell.pos[0].x.re, ", ", cell.pos[0].y.re);
	    }
	}
	
        foreach (blk; localFluidBlocks) {
            foreach (bc; blk.bc) {
                if (bc.type != "exchange_using_mapped_cells") {
                    //writeln("boundary: ", blk.id, ", ", bc.which_boundary, ", ", bc.faces.length, " ------");
                    foreach (i, face0; bc.faces) {
                        FVCell ghostcell;
                        if (bc.outsigns[i] == 1) {
                            ghostcell = face0.right_cell;
                        } else {
                            ghostcell = face0.left_cell;
                        }
                        string cid = to!string(ghostcell.id);
                        string bid = to!string(blk.id);
                        string gid = cid ~ bid;
                        ghostcell.global_id = to!size_t(gid);
                    }
                }
            }
        }
	
	foreach (blk; localFluidBlocks) {
	    foreach (bcond; blk.bc) {
		foreach (gce; bcond.preReconAction) {
		    auto mygce = cast(GhostCellMappedCellCopy)gce;
		    if (mygce && !blk.myConfig.in_mpi_context) {
			foreach (i, f; bcond.faces) {
			    // Only FVCell objects in an unstructured-grid are expected to have
			    // precomputed gradients.  There will be an initialized reference
			    // in the FVCell object of a structured-grid block, so we need to
			    // test and avoid copying from such a reference.
			    auto mapped_cell = mygce.get_mapped_cell(i);
			    FVCell ghost_cell;
			    if (bcond.outsigns[i] == 1) {
				ghost_cell = f.right_cell;
			    } else {
				ghost_cell = f.left_cell;
			    }
			    // writeln(ghost_cell.global_id, ", ", mapped_cell.global_id);
			    ghost_cell.global_id = mapped_cell.global_id;
			    // writeln(ghost_cell.global_id, ", ", mapped_cell.global_id);
			}
		    }
		}
	    }
	}
	
	foreach (myblk; localFluidBlocks) {
	    foreach(face; myblk.faces) {
		string id;
		if(face.left_cell.global_id > face.right_cell.global_id) id = to!string(face.left_cell.global_id) ~ "_" ~ to!string(face.right_cell.global_id);
		else id = to!string(face.right_cell.global_id) ~ "_" ~ to!string(face.left_cell.global_id);
		face.global_id = id;
	    }
	}

        foreach (myblk; localFluidBlocks) {
            foreach(cell; myblk.cells) {
                //writeln(cell.global_id, ", ", cell.pos[0].x.re, ", ", cell.pos[0].y.re, ", ", myblk.id);
            }
            

            // Make a stack-local copy of conserved quantities info
            myblk.nConserved = nConservedQuantities;
            myblk.MASS = massIdx;
            myblk.X_MOM = xMomIdx;
            myblk.Y_MOM = yMomIdx;
            myblk.Z_MOM = zMomIdx;
            myblk.TOT_ENERGY = totEnergyIdx;
            myblk.TKE = tkeIdx;
            myblk.OMEGA = omegaIdx;
        }
        
        bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
        size_t nPrimitive; 

        nPrimitive = nConservedQuantities;  // density, velocity(x,y), pressure
        //if (GlobalConfig.dimensions == 3) nPrimitive += 1; // velocity(z)
        //if (with_k_omega) nPrimitive += 2; // tke, omega

        /* Flow Jacobian */
        
        // we make sure that the ghost cells are filled here (we don't want to update ghost cells during the flow Jacobian formation)
        // TODO: currently the code works for second order interpolation for internal interfaces, with first order interpolation at boundaries.
        // Should think about how to extend this to second order at block boundaries (i.e. manage the copying of gradients correctly).
	
	if (GlobalConfig.sscOptions.read_frozen_limiter_values_from_file) {
	    //steadystate_core.evalRHS(0.0, 0); // fill in some values
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
		    } else {
			cell.gradients.velzPhi = 1.0;
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

	exchange_ghost_cell_boundary_data(0.0, 0, 0); // pseudoSimTime = 0.0; gtl = 0; ftl = 0
	

        foreach (myblk; localFluidBlocks) {
            initialisation(myblk, nPrimitive, myblk.myConfig.interpolation_order);

            // make sure ghost cells are filled
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
        }

        steadystate_core.evalRHS(0.0, 0);
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
        }

	foreach (myblk; parallel(localFluidBlocks,1)) {
            fill_boundary_conditions(myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS, false, false);
	}	

        //if (GlobalConfig.interpolation_order > 1) {
	foreach (myblk; localFluidBlocks) {
	    ghost_cell_connectivity_for_gradients(myblk);
	}

        foreach (myblk; parallel(localFluidBlocks,1)) {
	    local_flow_jacobian_transpose(myblk.JlocT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS);
	}

        foreach (myblk; localFluidBlocks) {
	    form_external_flow_jacobian_block_phase0(myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS);
        }

        foreach (myblk; localFluidBlocks) {
	    form_external_flow_jacobian_block_phase1(nPrimitive, myblk.myConfig.interpolation_order, EPS);
        }

        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.JextT = new SMatrix!number();
            form_external_flow_jacobian_block_phase2(myblk.JextT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS);
	}

        foreach (myblk; parallel(localFluidBlocks,1)) {
            form_external_flow_jacobian_block_phase3(myblk.JextT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS);
	}

	int ngcells = 0;
	foreach (myblk; localFluidBlocks) {
	    ngcells += myblk.cells.length;
	}
	
	/*
        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach(icell; myblk.cells) {
                // local
		foreach(jcell; myblk.cells) {
                    writeln(icell.global_id, ", ", jcell.global_id, " ------------------- ", myblk.id);
                    foreach(i; 0..nPrimitive) {
                        foreach(j; 0..nPrimitive) {
                            size_t I = icell.id*(nPrimitive) + i;
                            size_t J = jcell.id*(nPrimitive) + j;
                            writef("%.12e ", myblk.JlocT[I,J].re);
                        }
                        writef("\n");
                    }
                }
		// external
		foreach(jcell; 0..ngcells) {
                    writeln(icell.global_id, ", ", jcell, " external ------------------ ", myblk.id);
                    foreach(i; 0..nPrimitive) {
                        foreach(j; 0..nPrimitive) {
                            size_t I = icell.id*(nPrimitive) + i;
                            size_t J = jcell*(nPrimitive) + j;
                            writef("%.12e ", myblk.JextT[I,J].re);
                        }
                        writef("\n");
                    }
                }
		
            } 
        }
	
	*/
        steadystate_core.evalRHS(0.0, 0);
        exchange_ghost_cell_boundary_data(0.0, 0, 0); // pseudoSimTime = 0.0; gtl = 0; ftl = 0
        
        if (GlobalConfig.sscOptions.pseudotime) {
            double CFL = GlobalConfig.sscOptions.cfl0;
            double dt = steadystate_core.determine_initial_dt(CFL);

            int interpolation_order = 2;
            final switch (GlobalConfig.sscOptions.pseudotime_lhs_jacobian_order) {
            case 1:
                interpolation_order = 1;
                break;
            case 2:
                interpolation_order = 2;
                break;
            }
            /* low order Jacobian */
            foreach (myblk; parallel(localFluidBlocks,1)) {
                local_flow_jacobian_transpose(myblk.A, myblk, nPrimitive, interpolation_order, EPS, false, false);
                foreach (i; 0 .. nPrimitive*myblk.cells.length) {
                    myblk.A[i,i] = myblk.A[i,i] + (1.0/dt);
                }
            }

            foreach (myblk; localFluidBlocks) {
                form_external_flow_jacobian_block_phase0(myblk, nPrimitive, interpolation_order, EPS);
            }
            
            foreach (myblk; localFluidBlocks) {
                form_external_flow_jacobian_block_phase1(nPrimitive, interpolation_order, EPS);
            }
            
            foreach (myblk; parallel(localFluidBlocks,1)) {
                myblk.Aext = new SMatrix!number();
                form_external_flow_jacobian_block_phase2(myblk.Aext, myblk, nPrimitive, interpolation_order, EPS);
            }
            
            foreach (myblk; parallel(localFluidBlocks,1)) {
                form_external_flow_jacobian_block_phase3(myblk.Aext, myblk, nPrimitive, interpolation_order, EPS);
            }

            /* Preconditioner */
            final switch (GlobalConfig.sscOptions.adjoint_precondition_matrix_order) {
            case 0:
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    //myblk.P = new SMatrix!number();
                    local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 1, EPS, true, false); // orderOfJacobian=0
                    foreach (i; 0 .. nPrimitive*myblk.cells.length) {
                        myblk.P[i,i] = myblk.P[i,i] + (1.0/dt);
                    }
                    decompILU0(myblk.P);
                }
                break;
            case 1:
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    //myblk.P = new SMatrix!number();
                    local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 1, EPS, false, false); // orderOfJacobian=0
                    foreach (i; 0 .. nPrimitive*myblk.cells.length) {
                        myblk.P[i,i] = myblk.P[i,i] + (1.0/dt);
                    }
                    decompILU0(myblk.P);
                }
                break;
            case 2:
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    //myblk.P = new SMatrix!number();
                    local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 2, EPS, false, false); // orderOfJacobian=0
                    foreach (i; 0 .. nPrimitive*myblk.cells.length) {
                        myblk.P[i,i] = myblk.P[i,i] + (1.0/dt);
                    }
                    decompILU0(myblk.P);
                }
                break;
            }
        } else {
            /* Preconditioner */
            final switch (GlobalConfig.sscOptions.adjoint_precondition_matrix_order) {
            case 0:
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    //myblk.P = new SMatrix!number();
                    local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 1, EPS, true, false); // orderOfJacobian=0
                    decompILU0(myblk.P);
                }
                break;
            case 1:
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    //myblk.P = new SMatrix!number();
                    local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 1, EPS, false, false); // orderOfJacobian=0
                    decompILU0(myblk.P);
                }
                break;
            case 2:
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    //myblk.P = new SMatrix!number();
                    local_flow_jacobian_transpose(myblk.P, myblk, nPrimitive, 2, EPS, false, false); // orderOfJacobian=0
                    decompILU0(myblk.P);
                }
                break;
            }
        }
        
        /* Objective Function Sensitivity */

        // Surface intergal objective functions can be computed in parallel with a reduction process across blocks to gather the final value,
        // however the sensitivity w.r.t to primitive variables cannot be computed in parallel, since we are only computing a single objective calue (for example drag).
        foreach (myblk; localFluidBlocks) {
            form_objective_function_sensitivity(myblk, nPrimitive, EPS);
        }
        
        /* solve the adjoint system */

        if (GlobalConfig.sscOptions.pseudotime)
            rpcGMRES_solve1(nPrimitive);
        else
            rpcGMRES_solve0(nPrimitive);
        
        /* Write out adjoint variables for visualisation and to eilmer native format*/ 
        ensure_directory_is_present("adjoint-solution");
        foreach (myblk; localFluidBlocks) {
            write_adjoint_variables_to_file(myblk, nPrimitive, jobName);
            auto fileName = format!"adjoint-solution/%s.adjoint_vars.b%04d.gz"(jobName, myblk.id);
            myblk.write_adjoint_variables(fileName);
        }

        /* clear some expensive data structures from memory */

        foreach (myblk; parallel(localFluidBlocks,1)) {
            destroy(myblk.JlocT);
            destroy(myblk.JextT);
            destroy(myblk.P);
            GC.minimize();
        }
        
        // ------------------------------------------------------
        // RESIDUAL/OBJECTIVE SENSITIVITY W.R.T. DESIGN VARIABLES
        // ------------------------------------------------------
        /*
        if (adjointVerificationFlag) fit_design_parameters_to_surface(designVars);    
        else { 
            readBezierDataFromFile(designVars);
            readDesignVarsFromDakotaFile(designVars);
        }
        */
        readBezierDataFromFile(designVars);
        if (!adjointVerificationFlag) { 
            readDesignVarsFromDakotaFile(designVars);
        }
        
        // objective function sensitivity w.r.t. design variables
        number[] g;
        //readBezierDataFromFile(designVars);
        g.length = designVars.length;
        foreach (myblk; localFluidBlocks) {
            size_t nLocalCells = myblk.cells.length;
            myblk.rT = new Matrix!number(designVars.length, nLocalCells*nPrimitive);
        }

        compute_design_variable_partial_derivatives(designVars, g, nPrimitive, with_k_omega, EPS, jobName);
        
        // ---------------------
        // COMPUTE SENSITIVITIES
        // ---------------------
        foreach (myblk; parallel(localFluidBlocks,1)) {
            myblk.rTdotPsi.length = designVars.length;
            dot(myblk.rT, myblk.psi, myblk.rTdotPsi);
        }
        
        number[] adjointGradients;
        adjointGradients.length = designVars.length;
        adjointGradients[] = to!number(0.0);
        foreach (myblk; localFluidBlocks) adjointGradients[] += myblk.rTdotPsi[];
        adjointGradients[] = g[] - adjointGradients[];

        foreach ( i; 0..adjointGradients.length) writef("gradient for variable %d: %.16e \n", i+1, adjointGradients[i]);

        write_gradients_to_file("results.out", adjointGradients);
        // clear some expensive data structures from memory
        foreach (myblk; parallel(localFluidBlocks,1)) {
            destroy(myblk.rT);
            GC.minimize();
        }
        
        writeln("Simulation complete.");
        return;
    }

    /* Check Accuracy of Primitive Jacobian routines via Frechet Derivative Comparison */

    if (verifyPrimitiveJacobianFlag) {
        	// store global ids
	foreach (myblk; parallel(localFluidBlocks,1)) {
	    foreach (cell; myblk.cells) {
		cell.global_id = to!int(myblk.globalCellId(to!size_t(cell.id)));
		//writeln(cell.global_id, ", ", cell.pos[0].x.re, ", ", cell.pos[0].y.re);
	    }
	}
	
        foreach (blk; localFluidBlocks) {
            foreach (bc; blk.bc) {
                if (bc.type != "exchange_using_mapped_cells") {
                    //writeln("boundary: ", blk.id, ", ", bc.which_boundary, ", ", bc.faces.length, " ------");
                    foreach (i, face0; bc.faces) {
                        FVCell ghostcell;
                        if (bc.outsigns[i] == 1) {
                            ghostcell = face0.right_cell;
                        } else {
                            ghostcell = face0.left_cell;
                        }
                        string cid = to!string(ghostcell.id);
                        string bid = to!string(blk.id);
                        string gid = cid ~ bid;
                        ghostcell.global_id = to!size_t(gid);
                    }
                }
            }
        }
	
	foreach (blk; localFluidBlocks) {
	    foreach (bcond; blk.bc) {
		foreach (gce; bcond.preReconAction) {
		    auto mygce = cast(GhostCellMappedCellCopy)gce;
		    if (mygce && !blk.myConfig.in_mpi_context) {
			foreach (i, f; bcond.faces) {
			    // Only FVCell objects in an unstructured-grid are expected to have
			    // precomputed gradients.  There will be an initialized reference
			    // in the FVCell object of a structured-grid block, so we need to
			    // test and avoid copying from such a reference.
			    auto mapped_cell = mygce.get_mapped_cell(i);
			    FVCell ghost_cell;
			    if (bcond.outsigns[i] == 1) {
				ghost_cell = f.right_cell;
			    } else {
				ghost_cell = f.left_cell;
			    }
			    // writeln(ghost_cell.global_id, ", ", mapped_cell.global_id);
			    ghost_cell.global_id = mapped_cell.global_id;
			    // writeln(ghost_cell.global_id, ", ", mapped_cell.global_id);
			}
		    }
		}
	    }
	}
	
	foreach (myblk; localFluidBlocks) {
	    foreach(face; myblk.faces) {
		string id;
		if(face.left_cell.global_id > face.right_cell.global_id) id = to!string(face.left_cell.global_id) ~ "_" ~ to!string(face.right_cell.global_id);
		else id = to!string(face.right_cell.global_id) ~ "_" ~ to!string(face.left_cell.global_id);
		face.global_id = id;
	    }
	}

        foreach (myblk; localFluidBlocks) {
            // Make a stack-local copy of conserved quantities info
            myblk.nConserved = nConservedQuantities;
            myblk.MASS = massIdx;
            myblk.X_MOM = xMomIdx;
            myblk.Y_MOM = yMomIdx;
            myblk.Z_MOM = zMomIdx;
            myblk.TOT_ENERGY = totEnergyIdx;
            myblk.TKE = tkeIdx;
            myblk.OMEGA = omegaIdx;
        }

        // set number of primitive variables
        bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
        size_t nPrimitive; 
        nPrimitive = nConservedQuantities; //4;  // density, velocity(x,y), pressure
        //if (GlobalConfig.dimensions == 3) nPrimitive += 1; // velocity(z)
        //if (with_k_omega) nPrimitive += 2; // tke, omega
        
        // construct the transposed primitive Jacobian
        foreach (myblk; localFluidBlocks) {
            initialisation(myblk, nPrimitive, myblk.myConfig.interpolation_order);

            // make sure ghost cells are filled
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);

            steadystate_core.evalRHS(0.0, 0);

            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
            
            myblk.Minv = new Matrix!number(nPrimitive, nPrimitive);
            myblk.JcT = new SMatrix!number();
            myblk.A = new SMatrix!number();
            myblk.P = new SMatrix!number();
            //myblk.JlocT = new SMatrix!number();
            //myblk.JlocT = new SMatrix!number();
            local_flow_jacobian_transpose(myblk.JlocT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS);

            // transpose
            myblk.P.ia.length = myblk.JlocT.ia.length;
            myblk.P.ja.length = myblk.JlocT.ja.length;
            myblk.P.aa.length = myblk.JlocT.aa.length;
            nm.smla.transpose(myblk.JlocT.ia, myblk.JlocT.ja, myblk.JlocT.aa, myblk.P.ia, myblk.P.ja, myblk.P.aa);

            steadystate_core.evalRHS(0.0, 0);
        }
        
        foreach (blk; localFluidBlocks) {

            // compute arbitrary vector
            number[] v;
            v.length = blk.cells.length*nPrimitive;
            
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                v[i] = i+1; // theoretically can be any vector.
            }
            
            // normalise the vector
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
            
            // explicit multiplication of Jv
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                p1[i] = 0.0;
                foreach ( j; 0..blk.cells.length*nPrimitive) {
                    p1[i] += blk.P[i,j]*v[j];
                }
            }
                        
            // Frechet derivative of Jv
            evalPrimitiveJacobianVecProd(blk, nPrimitive, v, p2, EPS);

            // write out results for the 1st order Jacobian test
            {
                string fileName = "primitive_jacobian_test.output";
                auto outFile = File(fileName, "w");
                foreach( i; 0..v.length ) {
                    size_t id = i/nPrimitive;
                    outFile.writef("%d    %d    %.16e    %.16e    %.16f    %.16f \n", i, id, fabs((p1[i]-p2[i])/p1[i]).re, fabs(p1[i]-p2[i]).re, p1[i], p2[i].re);
                }
            }
        }
        writeln("Primitive Jacobian Test: COMPLETED");
        return;
    }

    /* Check Accuracy of Conservative Jacobian routines via Frechet Derivative Comparison */
    
    if (verifyConservativeJacobianFlag) {
        foreach (myblk; localFluidBlocks) {
            // Make a stack-local copy of conserved quantities info
            myblk.nConserved = nConservedQuantities;
            myblk.MASS = massIdx;
            myblk.X_MOM = xMomIdx;
            myblk.Y_MOM = yMomIdx;
            myblk.Z_MOM = zMomIdx;
            myblk.TOT_ENERGY = totEnergyIdx;
            myblk.TKE = tkeIdx;
            myblk.OMEGA = omegaIdx;
        }
        
        bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
        size_t nPrimitive; 
        nPrimitive = nConservedQuantities; //4;  // density, velocity(x,y), pressure
        //if (GlobalConfig.dimensions == 3) nPrimitive += 1; // velocity(z)
        //if (with_k_omega) nPrimitive += 2; // tke, omega
        
        // make sure ghost cells are filled
        foreach (myblk; parallel(localFluidBlocks,1)) {
            initialisation(myblk, nPrimitive, myblk.myConfig.interpolation_order);
            
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
            
            steadystate_core.evalRHS(0.0, 0);
            steadystate_core.evalRHS(0.0, 1);

            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);

            
            myblk.Minv = new Matrix!number(nPrimitive, nPrimitive);
            myblk.JcT = new SMatrix!number();
            myblk.A = new SMatrix!number();
            myblk.P = new SMatrix!number();
            myblk.JlocT = new SMatrix!number();
            local_flow_jacobian_transpose(myblk.JlocT, myblk, nPrimitive, myblk.myConfig.interpolation_order, EPS, false, true); // false, true
            // transpose
            myblk.P.ia.length = myblk.JlocT.ia.length;
            myblk.P.ja.length = myblk.JlocT.ja.length;
            myblk.P.aa.length = myblk.JlocT.aa.length;
            nm.smla.transpose(myblk.JlocT.ia, myblk.JlocT.ja, myblk.JlocT.aa, myblk.P.ia, myblk.P.ja, myblk.P.aa);
            
        }
        
        foreach (blk; localFluidBlocks) {
            
            //////////////////
            /*
            foreach (cell; blk.cells) {
                double[][] tmp;
                tmp.length = nPrimitive;
                foreach (ref a; tmp) a.length = nPrimitive;
                with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
                
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].mass += EPS;
                cell.decode_conserved(0, 1, 0.0);
                tmp[0][0] = cell.fs.gas.rho.im/EPS.im;
                tmp[1][0] = cell.fs.vel.x.im/EPS.im;
                tmp[2][0] = cell.fs.vel.y.im/EPS.im;
                tmp[3][0] = cell.fs.gas.p.im/EPS.im;
                tmp[4][0] = cell.fs.turb[0].im/EPS.im;
                tmp[5][0] = cell.fs.turb[1].im/EPS.im;
                    
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].momentum.refx += EPS;
                cell.decode_conserved(0, 1, 0.0);
                tmp[0][1] = cell.fs.gas.rho.im/EPS.im;
                tmp[1][1] = cell.fs.vel.x.im/EPS.im;
                tmp[2][1] = cell.fs.vel.y.im/EPS.im;
                tmp[3][1] = cell.fs.gas.p.im/EPS.im;
                tmp[4][1] = cell.fs.turb[0].im/EPS.im;
                tmp[5][1] = cell.fs.turb[1].im/EPS.im;
                
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].momentum.refy += EPS;
                cell.decode_conserved(0, 1, 0.0);
                tmp[0][2] = cell.fs.gas.rho.im/EPS.im;
                tmp[1][2] = cell.fs.vel.x.im/EPS.im;
                tmp[2][2] = cell.fs.vel.y.im/EPS.im;
                tmp[3][2] = cell.fs.gas.p.im/EPS.im;
                tmp[4][2] = cell.fs.turb[0].im/EPS.im;
                tmp[5][2] = cell.fs.turb[1].im/EPS.im;
                
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].total_energy += EPS;
                cell.decode_conserved(0, 1, 0.0);
                tmp[0][3] = cell.fs.gas.rho.im/EPS.im;
                tmp[1][3] = cell.fs.vel.x.im/EPS.im;
                tmp[2][3] = cell.fs.vel.y.im/EPS.im;
                tmp[3][3] = cell.fs.gas.p.im/EPS.im;
                tmp[4][3] = cell.fs.turb[0].im/EPS.im;
                tmp[5][3] = cell.fs.turb[1].im/EPS.im;
                
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].tke += EPS;
                cell.decode_conserved(0, 1, 0.0);
                tmp[0][4] = cell.fs.gas.rho.im/EPS.im;
                tmp[1][4] = cell.fs.vel.x.im/EPS.im;
                tmp[2][4] = cell.fs.vel.y.im/EPS.im;
                tmp[3][4] = cell.fs.gas.p.im/EPS.im;
                tmp[4][4] = cell.fs.turb[0].im/EPS.im;
                tmp[5][4] = cell.fs.turb[1].im/EPS.im;
                
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].omega += EPS;
                cell.decode_conserved(0, 1, 0.0);
                tmp[0][5] = cell.fs.gas.rho.im/EPS.im;
                tmp[1][5] = cell.fs.vel.x.im/EPS.im;
                tmp[2][5] = cell.fs.vel.y.im/EPS.im;
                tmp[3][5] = cell.fs.gas.p.im/EPS.im;
                tmp[4][5] = cell.fs.turb[0].im/EPS.im;
                tmp[5][5] = cell.fs.turb[1].im/EPS.im;
                
                writeln("ref: ", cell.id, ", ", tmp);
            }
            */
            //////////////////
            
            // compute arbitrary vector
            number[] v;
            v.length = blk.cells.length*nPrimitive;
            
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                v[i] = i+1; // theoretically can be any vector.
            }
            
            // normalise the vector
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

            // form 1st order Jacobian via ILU preconditioner routines
            //GlobalConfig.sssOptions.preconditionMatrixType = PreconditionMatrixType.ilu;
            //sss_preconditioner_initialisation(blk, nPrimitive);
            //sss_preconditioner(blk, nPrimitive, 1.0e+16);
            
            // form block diagonal precondition sub-matrices
            //GlobalConfig.sssOptions.preconditionMatrixType = PreconditionMatrixType.block_diagonal;
            //sss_preconditioner_initialisation(blk, nPrimitive);
            //sss_preconditioner(blk, nPrimitive, 1.0e+16);

            // explicit multiplication of Jv
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                p1[i] = 0.0;
                foreach ( j; 0..blk.cells.length*nPrimitive) {
                    p1[i] += blk.P[i,j]*v[j];
                }
            }

            // Frechet derivative of Jv
            evalConservativeJacobianVecProd(blk, nPrimitive, v, p2, EPS);

            // write out results for the 1st order Jacobian test
            {
                string fileName = "conservative_jacobian_test.output";
                auto outFile = File(fileName, "w");
                foreach( i; 0..v.length ) {
                    size_t id = i/nPrimitive;
                    outFile.writef("%d    %d    %.16e    %.16e    %.16f    %.16f \n", i, id, fabs((p1[i]-p2[i])/p1[i]).re, fabs(p1[i]-p2[i]).re, p1[i].re, p2[i].re);
                }
            }

            /*
            // write out results for testing the block-diagonal sub-matrices
            {
                string fileName = "convervative_jacobian_test1.output";
                auto outFile = File(fileName, "w");
                foreach (myblk; parallel(localFluidBlocks,1)) {
                    foreach( cell; myblk.cells ) {
                        size_t id = cell.id;
                        foreach ( i; 0..nPrimitive) {
                            foreach ( j; 0..nPrimitive) {
                                number z1 = blk.Jc[id*nPrimitive+i, id*nPrimitive+j];
                                number z2 = cell.dPrimitive[i, j];
                                outFile.writef("%d    %.16e    %.16e    %.16f    %.16f \n", id, fabs((z1-z2)/z1), fabs(z1-z2), z1, z2);
                            }
                        }
                    }
                }
            }
            */
        }
        writeln("Conservative Jacobian Test: COMPLETED");
        return;
    }

    /* Steady-state solver preconditioner routines test */
    if (verifySSSPreconditionerFlag) {
        foreach (myblk; localFluidBlocks) {
            // Make a stack-local copy of conserved quantities info
            myblk.nConserved = nConservedQuantities;
            myblk.MASS = massIdx;
            myblk.X_MOM = xMomIdx;
            myblk.Y_MOM = yMomIdx;
            myblk.Z_MOM = zMomIdx;
            myblk.TOT_ENERGY = totEnergyIdx;
            myblk.TKE = tkeIdx;
            myblk.OMEGA = omegaIdx;
        }
        
        size_t nPrimitive; 
        bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
        nPrimitive = nConservedQuantities; //4;  // density, velocity(x,y), pressure
        //if (GlobalConfig.dimensions == 3) nPrimitive += 1; // velocity(z)
        //if (with_k_omega) nPrimitive += 2; // tke, omega

        
        // Construct 1st order Flow Jacobian Transpose
        foreach (myblk; parallel(localFluidBlocks,1)) {
            // make sure ghost cells are filled before proceeding...
            myblk.applyPreReconAction(0.0, 0, 0);
            myblk.applyPostConvFluxAction(0.0, 0, 0);
            myblk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
            myblk.applyPostDiffFluxAction(0.0, 0, 0);
            
            myblk.JlocT = new SMatrix!number();
            local_flow_jacobian_transpose(myblk.JlocT, myblk, nPrimitive, 1, EPS);
        }

        
        // form the conservative Jacobian diagonal blocks used in the preconditioner
        foreach (myblk; parallel(localFluidBlocks,1)) {
            sss_preconditioner_initialisation(myblk, nPrimitive);
            sss_preconditioner(myblk, nPrimitive, 1.0e+16);
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
            
            // take the transpose of the transpose Jacobian
            Matrix!number Ap; // primitive Jacobian
            Ap = new Matrix!number(nPrimitive*blk.cells.length, nPrimitive*blk.cells.length);
            foreach ( i; 0..nPrimitive*blk.cells.length) {
                foreach ( j; 0..nPrimitive*blk.cells.length) {
                    Ap[i,j] = blk.JlocT[j,i];
                }
            }
            
            // compute transform matrix
            Matrix!number M;
            M = zeros!number(nPrimitive*blk.cells.length, nPrimitive*blk.cells.length);
            foreach ( cell; blk.cells ) {
                auto gmodel = blk.myConfig.gmodel;
                number gamma = gmodel.gamma(cell.fs.gas);
                // form inverse transformation matrix
                blk.Minv[0,0] = to!number(1.0);
                blk.Minv[0,1] = to!number(0.0);
                blk.Minv[0,2] = to!number(0.0);
                blk.Minv[0,3] = to!number(0.0);
                // second row
                blk.Minv[1,0] = -cell.fs.vel.x/cell.fs.gas.rho;
                blk.Minv[1,1] = 1.0/cell.fs.gas.rho;
                blk.Minv[1,2] = to!number(0.0);
                blk.Minv[1,3] = to!number(0.0);
                // third row
                blk.Minv[2,0] = -cell.fs.vel.y/cell.fs.gas.rho;
                blk.Minv[2,1] = to!number(0.0);
                blk.Minv[2,2] = 1.0/cell.fs.gas.rho;
                blk.Minv[2,3] = to!number(0.0);
                // fourth row
                blk.Minv[3,0] = 0.5*(gamma-1.0)*(cell.fs.vel.x*cell.fs.vel.x+cell.fs.vel.y*cell.fs.vel.y);
                blk.Minv[3,1] = -cell.fs.vel.x*(gamma-1);
                blk.Minv[3,2] = -cell.fs.vel.y*(gamma-1);
                blk.Minv[3,3] = gamma-1.0;
                
                size_t id = cell.id;
                foreach ( i; 0..nPrimitive) {
                    foreach ( j; 0..nPrimitive) {
                        M[id*nPrimitive+i, id*nPrimitive+j] = blk.Minv[i, j];
                    }
                }
            }
            
            // transform primitive => conservative
            Matrix!number Ac; // conservative Jacobian
            Ac = new Matrix!number(nPrimitive*blk.cells.length, nPrimitive*blk.cells.length);
            for (size_t i = 0; i < nPrimitive*blk.cells.length; i++) {
                for (size_t j = 0; j < nPrimitive*blk.cells.length; j++) {
                    Ac[i,j] = to!number(0.0);
                    for (size_t k = 0; k < nPrimitive*blk.cells.length; k++) {
                        Ac[i,j] += Ap[i,k]*M[k,j];
                    }
                }
            }

            // drop the off-diagonal terms
            Matrix!number P;
            P = zeros!number(nPrimitive*blk.cells.length, nPrimitive*blk.cells.length);
            foreach ( cell; blk.cells ) {
                foreach ( i; 0..nPrimitive) {
                    size_t I = cell.id*nPrimitive + i;
                    foreach ( j; 0..nPrimitive) {
                        size_t J = cell.id*nPrimitive + j;
                        P[I,J] = Ac[I,J];
                    }
                }
            }

            // invert P
            Matrix!number Pinv = inverse(P); 

            number[] z;
            z.length = nPrimitive*blk.cells.length;
            
            // explicit multiplication of P**(-1).v 
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                z[i] = 0.0;
                foreach ( j; 0..blk.cells.length*nPrimitive) {
                    z[i] += Pinv[i,j]*v[j];
                }
            }

            
            // explicit multiplication of Ac.v
            foreach ( i; 0..blk.cells.length*nPrimitive) {
                p1[i] = 0.0;
                foreach ( j; 0..blk.cells.length*nPrimitive) {
                    p1[i] += Ac[i,j]*z[j];
                }
            }
            
            
            // perform same computation with steady-state solver mechanics (i.e. use preconditioner routines & Frechet derivative)
            z[] = to!number(0.0);
            int cellCount = 0;
            number[] tmp;
            tmp.length = nPrimitive;
            foreach (cell; blk.cells) {
                LUSolve!number(cell.dConservative, cell.pivot,
                               v[cellCount..cellCount+nPrimitive], tmp);
                z[cellCount..cellCount+nPrimitive] = tmp[];
                cellCount += nPrimitive;
            }

            // Frechet derivative of Jv
            evalConservativeJacobianVecProd(blk, nPrimitive, z, p2, EPS);
            
            string fileName = "sss_preconditioner_test.output";
            auto outFile = File(fileName, "w");
            foreach( i; 0..v.length ) {
                size_t id = i/4;
                outFile.writef("%d    %d    %.16e    %.16e    %.16f    %.16f \n", i, id, fabs((p1[i]-p2[i])/p1[i]).re, fabs(p1[i]-p2[i]).re, p1[i].re, p2[i].re);
            }
        }
        writeln("Steady-state Solver Preconditioner Test: COMPLETE.");
        return;
    }
    
    /* Program should terminate before here */
    assert(0, "Oh dear. The Eilmer Shape Sensitivity Calculator executed correctly, however, nothing meaningful has been computed. Please check command line flags.");
}
