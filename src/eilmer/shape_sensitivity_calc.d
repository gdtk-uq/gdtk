/** shape_sensitivity_calc.d
 * 
 * Eilmer4 shape sensitivity calculator code, core coordination functions.
 *
 * Note: This is the first attempt at 'production' code.
 * Some test implementations began on 2017-09-18.
 *
 * Author: Kyle D.
 * Date: 2018-03-07
 *
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
import fluidblock;
import sfluidblock;
import ufluidblock;
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
import lsqinterp;
import bc;
import grid_deform;
import shape_sensitivity;

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
    writeln("Eilmer shape sensitivity calculator code:");
    writeln("Revision: PUT_REVISION_STRING_HERE");
    version(mpi_parallel) {
        assert(0, "Adjoint solver is not MPI parallel, yet.");
    }

    string msg = "Usage:                                    Comment:\n";
    msg       ~= "e4ssc      [--job=<string>]            name of job\n";
    msg       ~= "           [--prep]                               \n";
    msg       ~= "           [--return-objective-function]          \n";
    msg       ~= "           [--parameterise-surfaces               \n";
    msg       ~= "           [--grid-update]                        \n";
    msg       ~= "           [--verification]                       \n";
    msg       ~= "           [--max-cpus=<int>]           defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine            \n";
    msg       ~= "           [--max-wall-clock=<int>]     in seconds\n";
    msg       ~= "           [--help]            writes this message\n";
    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        write(msg);
        exit(1);
    }
    string jobName = "";
    bool returnObjFcnFlag = false;
    bool parameteriseSurfacesFlag = false;
    bool gridUpdateFlag = false;
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
    init_simulation(last_tindx, 0, maxCPUs, 1, maxWallClock);

    // perform some config checks
    if ( GlobalConfig.interpolation_order > 1 &&
         GlobalConfig.suppress_reconstruction_at_boundaries == false) {
        writeln("WARNING:");
        writeln("   suppress_reconstruction_at_boundaries is set to false.");
        writeln("   This setting must be true when using the shape sensitivity calculator.");
        writeln("   Its use will likely cause errorneous values for boundary fluxes.");
        writeln("   Continuing with simulation anyway.");
        writeln("END WARNING.");
    }
    
    if (GlobalConfig.diffuseWallBCsOnInit) {
        writeln("WARNING:");
        writeln("   diffuse_wall_bcs_on_init is set to true.");
        writeln("   This setting must be false when using the shape sensitivity calculator.");
        writeln("   Its use will likely cause errorneous values for gradients.");
        writeln("   Continuing with simulation anyway.");
        writeln("END WARNING.");
    }
    
    // set some global config values
    GradientMethod gradientMethod = GlobalConfig.sscOptions.gradientMethod;
    //bool gradientVerification = GlobalConfig.sscOptions.gradientVerification;
    // finite-difference perturbation parameters
    double EPSILON = GlobalConfig.sscOptions.epsilon; // flow Jacobian
    double MU = GlobalConfig.sscOptions.mu; // flow Jacobian
    double ETA = GlobalConfig.sscOptions.eta; // residual sensitivity
    double DELTA = GlobalConfig.sscOptions.delta; // finite-difference gradient
    // bezier curve parameters
    double bezier_curve_tolerance = GlobalConfig.sscOptions.tolBezierCurveFit;
    int bezier_curve_max_steps = GlobalConfig.sscOptions.maxStepsBezierCurveFit;
    
    // Initialise Lua state for calling user-defined objective function.
    //initLuaStateForUserDefinedObjFunc();

    FluidBlock myblk = localFluidBlocks[0]; // currently only compatible with single block simulations
    // number of design variables
    size_t ndvars = 0;
    size_t ndim = GlobalConfig.dimensions;

    // save a copy of the original mesh for later use
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.sync_vertices_to_underlying_grid(0);
        ensure_directory_is_present(make_path_name!"grid-original"(0));
        auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
        blk.write_underlying_grid(fileName);
    }
    

    // -----------------------------
    // -----------------------------
    // GEOMETRY PARAMETERISATION
    // -----------------------------
    // -----------------------------
    if (parameteriseSurfacesFlag) {
        parameteriseSurfaces(myblk, ndvars, ndim, bezier_curve_tolerance, bezier_curve_max_steps);    
        writeBezierCntrlPtsToDakotaFile(myblk, ndvars, jobName);
        return; // --parameterise-surfaces complete
    }
    
    // -----------------------------
    // -----------------------------
    // Grid Update
    // -----------------------------
    // -----------------------------
    Vector3[] design_variables;
    if (gridUpdateFlag) {
        parameteriseSurfaces(myblk, ndvars, ndim, bezier_curve_tolerance, bezier_curve_max_steps);
        prepForMeshPerturbation(myblk);
        gridUpdate(myblk, design_variables, jobName);
        return; // --grid-update complete
    }

    // -----------------------------
    // -----------------------------
    // OBJECTIVE FUNCTION EVALUATION
    //------------------------------
    // -----------------------------
    double objFnEval;
    if (returnObjFcnFlag) {
        foreach (bndary; myblk.bc) {
            if (bndary.group == "design") objFnEval = cost_function(bndary, 0); // assume gtl = 0
        }
        writeln("objective fn evaluation: ", objFnEval);
        write_objective_fn_to_file("results.out", objFnEval);
        
        return; // --return-objective-function
    }
    // ----------------------------
    // ----------------------------
    // SHAPE SENSITIVITY CALCULATOR
    // -----------------------------
    // ----------------------------
    
    // --------------
    // ADJOINT SOLVER
    //---------------
    
    // set the number of primitive variables
    bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
    size_t np; size_t nb = 0; 
    if (GlobalConfig.dimensions == 2) {
        if (with_k_omega) np = 6; // density, velocity(x,y), pressure, tke, omega
        else np = 4; // density, velocity(x,y), pressure
    }
    else if (GlobalConfig.dimensions == 3) {
        if (with_k_omega) np = 7; // density, velocity(x,y,z), pressure, tke, omega
        else np = 5; // density, velocity(x,y,z), pressure, tke, omega
    }

    // -----------------------------
    // Build transpose Flow Jacobian
    //------------------------------
    SMatrix L = new SMatrix();    
    build_flow_jacobian(L, myblk, ndim, np, myblk.myConfig.interpolation_order, EPSILON, MU); // orderOfJacobian=interpolation_order

    // -----------------------------
    // Build 1st order Flow Jacobian
    //------------------------------
    SMatrix L1 = new SMatrix();    
    build_flow_jacobian(L1, myblk, ndim, np, 1, EPSILON, MU); // orderOfJacobian=1

    // -----------------------------------------
    // Build approximate 1st order Flow Jacobian
    //------------------------------------------
    SMatrix L1D = new SMatrix();
    build_flow_jacobian(L1D, myblk, ndim, np, 0, EPSILON, MU); // orderOfJacobian=0
    
    // ----------------------------------
    // Construct cost function sensitvity
    // ----------------------------------
    double[] g;
    g.length = myblk.cells.length*np;
    cost_function_sensitivity(g, myblk, np, EPSILON, MU);

    // --------------------
    // Solve adjoint system 
    // --------------------
    double[] psi;
    psi = adjoint_solver(L, g, L1, L1D, myblk, np);

    // ---------------------------------------------
    // Write out adjoint variables for visualisation 
    // ---------------------------------------------
    write_adjoint_variables_to_file(myblk, psi, np, jobName);
    
    // clear some expensive data structures from memory
    destroy(L);
    destroy(L1);
    destroy(L1D);
    GC.minimize();

    // ------------------------------
    // RESIDUAL/OBJECTIVE SENSITIVITY
    // ------------------------------
    parameteriseSurfaces(myblk, ndvars, ndim, bezier_curve_tolerance, bezier_curve_max_steps);    
    prepForMeshPerturbation(myblk);
    size_t ncells = myblk.cells.length;
    double[] dJdD;
    dJdD.length = ndvars;
    Matrix dRdD_T;
    dRdD_T = new Matrix(ndvars, myblk.cells.length*np);

    foreach (bndary; myblk.bc) {
        if (bndary.is_design_surface) {
            foreach (i; 1 .. bndary.num_cntrl_pts-1) {
                int gtl; int ftl; double Jp; double Jm; string varID; double P0; double dP;

                // x
                // initialise all positions
		foreach(otherBndary; myblk.bc) {
		    foreach(j, vtx; otherBndary.vertices) {
			    vtx.pos[1].refx = vtx.pos[0].x;
			    vtx.pos[1].refy = vtx.pos[0].y;
			    vtx.pos[2].refx = vtx.pos[0].x;
			    vtx.pos[2].refy = vtx.pos[0].y;
		    }
		}
                evalRHS(0.0, 0, 0, with_k_omega, myblk);
                // store origianl value, and compute perturbation
                P0 = bndary.bezier.B[i].x;
                dP = ETA; //max(OMEGA*abs(P0), 1.0e-10);

                // perturb design variable +ve
                gtl = 1; ftl = 1;
                bndary.bezier.B[i].refx = P0 + dP;

                // update design surface
                foreach(j, vtx; bndary.vertices) {
                    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                }

                // perturb mesh, and compute new geometry
		inverse_distance_weighting(myblk, gtl);
		myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
                if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                        (myblk.myConfig.interpolation_order > 1)) { 
                        auto myUBlock = cast(UFluidBlock) myblk;
                        myUBlock.compute_least_squares_setup(gtl);
                }
                // compute perturbed flux
                myblk.clear_fluxes_of_conserved_quantities();
                evalRHS(0.0, ftl, gtl, with_k_omega, myblk);
                foreach (otherBndary; myblk.bc) {
                    if (otherBndary.group == "design") Jp = cost_function(otherBndary, gtl);
                }
		writef("perturbed %d J(D, Q(D) = %.16f , h = %.16f \n", i, Jp, dP);

                // perturb design variable -ve
                gtl = 2; ftl = 2;
		bndary.bezier.B[i].refx = P0 - dP;

                // update design variable
                foreach(j, vtx; bndary.vertices) {
                    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                }

                // perturb mesh, and compute new geometry
		inverse_distance_weighting(myblk, gtl);
		myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
                if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                    (myblk.myConfig.interpolation_order > 1)) { 
                    auto myUBlock = cast(UFluidBlock) myblk;
                    myUBlock.compute_least_squares_setup(gtl);
                }
                
                // compute perturbed flux
                myblk.clear_fluxes_of_conserved_quantities();
                evalRHS(0.0, ftl, gtl, with_k_omega, myblk);
                foreach (otherBndary; myblk.bc) {
                    if (otherBndary.group == "design") Jm = cost_function(otherBndary, gtl);
                }
		writef("perturbed %d J(D, Q(D) = %.16f , h = %.16f \n", i, Jm, dP);

                // compute cost function sensitivity
                dJdD[ndim*(i-1)] = (Jp-Jm)/(2.0*dP);
		
                // compute residual sensitivity
                foreach(j, cell; myblk.cells) {
		    dRdD_T[ndim*(i-1), j*np] = (cell.dUdt[1].mass - cell.dUdt[2].mass)/(2.0*dP);
		    dRdD_T[ndim*(i-1), j*np+1] = (cell.dUdt[1].momentum.x - cell.dUdt[2].momentum.x)/(2.0*dP);
		    dRdD_T[ndim*(i-1), j*np+2] = (cell.dUdt[1].momentum.y - cell.dUdt[2].momentum.y)/(2.0*dP);
		    dRdD_T[ndim*(i-1), j*np+3] = (cell.dUdt[1].total_energy - cell.dUdt[2].total_energy)/(2.0*dP);
		}

                // restore design variable
                bndary.bezier.B[i].refx = P0;

                // y
                // initialise all positions
		foreach(otherBndary; myblk.bc) {
		    foreach(j, vtx; otherBndary.vertices) {
			    vtx.pos[1].refx = vtx.pos[0].x;
			    vtx.pos[1].refy = vtx.pos[0].y;
			    vtx.pos[2].refx = vtx.pos[0].x;
			    vtx.pos[2].refy = vtx.pos[0].y;
		    }
		}
                evalRHS(0.0, 0, 0, with_k_omega, myblk);
                // store origianl value, and compute perturbation
                P0 = bndary.bezier.B[i].y;
                dP = ETA; //max(OMEGA*abs(P0), 1.0e-10);

                // perturb design variable +ve
                gtl = 1; ftl = 1;
                bndary.bezier.B[i].refy = P0 + dP;

                // update design surface
                foreach(j, vtx; bndary.vertices) {
                    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                }

                // perturb mesh, and compute new geometry
		inverse_distance_weighting(myblk, gtl);
		myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
                if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                        (myblk.myConfig.interpolation_order > 1)) { 
                        auto myUBlock = cast(UFluidBlock) myblk;
                        myUBlock.compute_least_squares_setup(gtl);
                }
                // compute perturbed flux
                myblk.clear_fluxes_of_conserved_quantities();
                evalRHS(0.0, ftl, gtl, with_k_omega, myblk);
                foreach (otherBndary; myblk.bc) {
                    if (otherBndary.group == "design") Jp = cost_function(otherBndary, gtl);
                }
		writef("perturbed %d J(D, Q(D) = %.16f , h = %.16f \n", i, Jp, dP);

                // perturb design variable -ve
                gtl = 2; ftl = 2;
		bndary.bezier.B[i].refy = P0 - dP;

                // update design variable
                foreach(j, vtx; bndary.vertices) {
                    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                }

                // perturb mesh, and compute new geometry
		inverse_distance_weighting(myblk, gtl);
		myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
                if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                    (myblk.myConfig.interpolation_order > 1)) { 
                    auto myUBlock = cast(UFluidBlock) myblk;
                    myUBlock.compute_least_squares_setup(gtl);
                }
                
                // compute perturbed flux
                myblk.clear_fluxes_of_conserved_quantities();
                evalRHS(0.0, ftl, gtl, with_k_omega, myblk);
                foreach (otherBndary; myblk.bc) {
                    if (otherBndary.group == "design") Jm = cost_function(otherBndary, gtl);
                }
		writef("perturbed %d J(D, Q(D) = %.16f , h = %.16f \n", i, Jm, dP);

                // compute cost function sensitivity
                dJdD[ndim*(i-1)+1] = (Jp-Jm)/(2.0*dP);
		
                // compute residual sensitivity
                foreach(j, cell; myblk.cells) {
		    dRdD_T[ndim*(i-1)+1, j*np] = (cell.dUdt[1].mass - cell.dUdt[2].mass)/(2.0*dP);
		    dRdD_T[ndim*(i-1)+1, j*np+1] = (cell.dUdt[1].momentum.x - cell.dUdt[2].momentum.x)/(2.0*dP);
		    dRdD_T[ndim*(i-1)+1, j*np+2] = (cell.dUdt[1].momentum.y - cell.dUdt[2].momentum.y)/(2.0*dP);
		    dRdD_T[ndim*(i-1)+1, j*np+3] = (cell.dUdt[1].total_energy - cell.dUdt[2].total_energy)/(2.0*dP);
		}

                // restore design variable
                bndary.bezier.B[i].refy = P0;
	    }
	}
    }

    //writeln("dRdD_T");
    //writeln("psi = ", psi);
    //writeln("dRdD_T = ", dRdD_T);
    //writeln(dJdV);
    //writeln(dJdD);
    // ------------------------------
    // 3. Compute shape sensitivities
    // ------------------------------
    writeln("compute shape sensitivities");
    double[] adjointGradients;
    adjointGradients.length = dRdD_T.nrows;
    dot(dRdD_T, psi, adjointGradients);    
    foreach(i; 0..ndvars) adjointGradients[i] = dJdD[i] - adjointGradients[i];
    writeln("adjoint gradients = ", adjointGradients);
    write_gradients_to_file("results.out", adjointGradients);
    // clear the sensitivity matrices from memory
    //destroy(globaldRdXT);
    destroy(dRdD_T);
    //destroy(tempMatrix);
    GC.minimize();
    
    // -----------------------------------------------------
    // Finite difference verification
    // -----------------------------------------------------
    if (verificationFlag) {
        double[] finiteDiffGradients; double J1; double J0;
        
        // save a copy of the original mesh for later use
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.sync_vertices_to_underlying_grid(0);
            ensure_directory_is_present(make_path_name!"grid-original"(0));
            auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
            blk.write_underlying_grid(fileName);
        }
        
        foreach (bndary; myblk.bc) {
	    if (bndary.is_design_surface) {
		foreach (i; 1 .. bndary.num_cntrl_pts-1) {
                    double dP = DELTA; //max(ETA*abs(P0), 1.0e-10);                
                    double err; double P0;
                    // x
		    writeln("computing finite difference gradient for variable ", i*ndim-1, " out of ", ndvars, " variables");
		    // could loop over all x, y, z coordinates here -- for now just use y-coordinates
		    string varID;
                    foreach(otherBndary; myblk.bc) {
                        foreach(j, vtx; otherBndary.vertices) {
                            vtx.pos[1].refx = vtx.pos[0].x;
                            vtx.pos[1].refy = vtx.pos[0].y;
                            vtx.pos[2].refx = vtx.pos[0].x;
                            vtx.pos[2].refy = vtx.pos[0].y;
                        }
                    }
		    P0 = bndary.bezier.B[i].refx;
                    bndary.bezier.B[i].refx = P0 + dP;
                    varID = "p-x-" ~ to!string(i);
		    J0 = finite_difference_grad(jobName, last_tindx, myblk, bndary, varID);
		    writef("perturbed J(D, Q(D) = %.16f , h = %.16f \n", J0, dP);
		    bndary.bezier.B[i].refx = P0 - dP;
                    varID = "m-x-" ~ to!string(i);
		    J1 = finite_difference_grad(jobName, last_tindx, myblk, bndary, varID);
		    finiteDiffGradients ~= (J0 - J1)/(2.0*dP);
		    bndary.bezier.B[i].refx = P0;
                    writef("perturbed J(D, Q(D) = %.16f , h = %.16f \n", J1, dP);

                    err = ((adjointGradients[ndim*i-2] - finiteDiffGradients[ndim*i-2])/finiteDiffGradients[ndim*i-2]) * 100.0;
		    writeln("finite-difference gradient for variable ", ndim*i-1, ": ", finiteDiffGradients[ndim*i-2], ", % error: ", abs(err));

                    // y
		    writeln("computing finite difference gradient for variable ", i*ndim, " out of ", ndvars, " variables");
		    // could loop over all x, y, z coordinates here -- for now just use y-coordinates
                    foreach(otherBndary; myblk.bc) {
                        foreach(j, vtx; otherBndary.vertices) {
                            vtx.pos[1].refx = vtx.pos[0].x;
                            vtx.pos[1].refy = vtx.pos[0].y;
                            vtx.pos[2].refx = vtx.pos[0].x;
                            vtx.pos[2].refy = vtx.pos[0].y;
                        }
                    }
		    P0 = bndary.bezier.B[i].refy;
                    bndary.bezier.B[i].refy = P0 + dP;
                    varID = "p-y-" ~ to!string(i);
		    J0 = finite_difference_grad(jobName, last_tindx, myblk, bndary, varID);
		    writef("perturbed J(D, Q(D) = %.16f , h = %.16f \n", J0, dP);
		    bndary.bezier.B[i].refy = P0 - dP;
                    varID = "m-y-" ~ to!string(i);
		    J1 = finite_difference_grad(jobName, last_tindx, myblk, bndary, varID);
		    finiteDiffGradients ~= (J0 - J1)/(2.0*dP);
		    bndary.bezier.B[i].refy = P0;
                    writef("perturbed J(D, Q(D) = %.16f , h = %.16f \n", J1, dP);

                    err = ((adjointGradients[ndim*i-1] - finiteDiffGradients[ndim*i-1])/finiteDiffGradients[ndim*i-1]) * 100.0;
		    writeln("finite-difference gradient for variable ", ndim*i, ": ", finiteDiffGradients[ndim*i-1], ", % error: ", abs(err));
                }
	    }
	}
    }
    writeln("Simulation complete.");
}
