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

/*
enum ghost_cell_start_id = 1_000_000_000;
immutable double ESSENTIALLY_ZERO = 1.0e-50;
// some data objects used in forming the Jacobian
immutable size_t MAX_PERTURBED_INTERFACES = 40;
FVCell cellOrig;
FVInterface[MAX_PERTURBED_INTERFACES] ifaceOrig;
FVInterface[MAX_PERTURBED_INTERFACES] ifacePp;
FVInterface[MAX_PERTURBED_INTERFACES] ifacePm;
*/

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

    string msg = "Usage:                              Comment:\n";
    msg       ~= "e4ssc      [--job=<string>]            name of job\n";
    msg       ~= "           [--no-gradients=<bool>]     \n";
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
    bool noGradients = false;
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "no-gradients", &noGradients,
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
    bool gradientVerification = GlobalConfig.sscOptions.gradientVerification;
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

    writeln("----------------");
    writeln("Running With Perturbation Parameters:");
    writeln("EPSILON  = ", EPSILON);
    writeln("MU  = ", MU);
    writeln("ETA = ", ETA);
    writeln("DELTA = ", DELTA);
    writeln("----------------");
    
    // save a copy of the original mesh for later use
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.sync_vertices_to_underlying_grid(0);
        ensure_directory_is_present(make_path_name!"grid-original"(0));
        auto fileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
        blk.write_underlying_grid(fileName);
    }

    FluidBlock myblk = localFluidBlocks[0]; // currently only compatible with single block simulations
    
    // -----------------------------
    // -----------------------------
    // OBJECTIVE FUNCTION EVALUATION
    //------------------------------
    // -----------------------------
    double objFnEval;
    foreach (bndary; myblk.bc) {
        if (bndary.group == "design") objFnEval = cost_function(bndary, 0); // assume gtl = 0
    }
    writeln("objective fn evaluation: ", objFnEval);
    write_objective_fn_to_file("results.out", objFnEval);
    
    if (noGradients) return;
        
    // --------------
    // --------------
    // ADJOINT SOLVER
    //---------------
    // --------------

    // set the number of primitive variables
    size_t ndim = GlobalConfig.dimensions;
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

    // ----------------------------
    // ----------------------------
    // SHAPE SENSITIVITY CALCULATOR
    // -----------------------------
    // ----------------------------

    // ----------------------------
    // Compute residual sensitivity
    // ----------------------------

    string[] NonFixedBoundaryList;
    NonFixedBoundaryList ~= "outflow";
    size_t ndvars = 0;
    foreach (bndary; myblk.bc) {
        if (bndary.is_design_surface) {
	    // gather boundary vertices
            size_t[] vtx_id_list;
            foreach(face; bndary.faces) {
                foreach(i, vtx; face.vtx) {
                    if (i == 0) {
			bndary.vertices ~= vtx;
			bndary.surfacePoints ~= Vector3(vtx.pos[0].refx, vtx.pos[0].refy, vtx.pos[0].refz);
                        vtx_id_list ~= vtx.id;
			myblk.boundaryVtxIndexList ~= vtx.id;
		    }
		    else if (i > 0 && vtx_id_list.canFind(vtx.id) == false) {
			// currently assume points in x-order
			FVVertex vtx0 = bndary.vertices[$-1];
			if (vtx.pos[0].x > vtx0.pos[0].x) {
			    bndary.vertices ~= vtx;
			    bndary.surfacePoints ~= Vector3(vtx.pos[0].refx, vtx.pos[0].refy, vtx.pos[0].refz);
                            vtx_id_list ~= vtx.id;
			    myblk.boundaryVtxIndexList ~= vtx.id;
			}
			else {
			    bndary.vertices ~= vtx0;
			    bndary.surfacePoints ~= Vector3(vtx0.pos[0].refx, vtx0.pos[0].refy, vtx0.pos[0].refz);
                            vtx_id_list ~= vtx0.id;
			    myblk.boundaryVtxIndexList ~= vtx0.id;

			    bndary.vertices[i-1] = vtx;
			    bndary.surfacePoints[i-1] = Vector3(vtx.pos[0].refx, vtx.pos[0].refy, vtx.pos[0].refz);
                            vtx_id_list[i-1] = vtx.id;
			    myblk.boundaryVtxIndexList[i-1] = vtx.id;
			}
		    }
                }
            }
            // compute bezier control points
            ndvars += bndary.num_cntrl_pts;
            bndary.bezier = optimiseBezierPoints(bndary.surfacePoints, bndary.num_cntrl_pts, bndary.ts, bezier_curve_tolerance, bezier_curve_max_steps);
        }
	else {
	    size_t[] vtx_id_list;
            foreach(face; bndary.faces) {
                foreach(vtx; face.vtx) {
                    if (vtx_id_list.canFind(vtx.id) == false) {
                        bndary.vertices ~= vtx;
			vtx_id_list ~= vtx.id;
			if (NonFixedBoundaryList.canFind(bndary.group) == false) myblk.boundaryVtxIndexList ~= vtx.id;
                    }
                }
            }
	}
    }

    size_t ncells = myblk.cells.length;
    double[] dJdD;
    dJdD.length = ndvars;
    Matrix dRdD_T;
    dRdD_T = new Matrix(ndvars, myblk.cells.length*np);

    foreach (bndary; myblk.bc) {
        if (bndary.is_design_surface) {
            foreach (i; 0 .. bndary.num_cntrl_pts) {
                int gtl; int ftl; double P0; double dP;
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
                double grad;
                foreach(j, vtx; bndary.vertices) {
                    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
                    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
                }

                // perturb mesh, and compute new geometry
		inverse_distance_weighting(myblk, NonFixedBoundaryList, gtl);
		myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
                if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                        (myblk.myConfig.interpolation_order > 1)) { 
                        auto myUBlock = cast(UFluidBlock) myblk;
                        myUBlock.compute_least_squares_setup(gtl);
                }
                // compute perturbed flux
                myblk.clear_fluxes_of_conserved_quantities();
                evalRHS(0.0, ftl, gtl, with_k_omega, myblk);
                double Jp;
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
		inverse_distance_weighting(myblk, NonFixedBoundaryList, gtl);
		myblk.compute_primary_cell_geometric_data(gtl); // need to add in 2nd order effects
                if ((myblk.grid_type == Grid_t.unstructured_grid) &&
                    (myblk.myConfig.interpolation_order > 1)) { 
                    auto myUBlock = cast(UFluidBlock) myblk;
                    myUBlock.compute_least_squares_setup(gtl);
                }
                
                // compute perturbed flux
                myblk.clear_fluxes_of_conserved_quantities();
                evalRHS(0.0, ftl, gtl, with_k_omega, myblk);
                double Jm;
                foreach (otherBndary; myblk.bc) {
                    if (otherBndary.group == "design") Jm = cost_function(otherBndary, gtl);
                }
		writef("perturbed %d J(D, Q(D) = %.16f , h = %.16f \n", i, Jm, dP);

                // compute cost function sensitivity
                dJdD[i] = (Jp-Jm)/(2.0*dP);
		
                // compute residual sensitivity
                foreach(j, cell; myblk.cells) {
		    dRdD_T[i, j*np] = (cell.dUdt[1].mass - cell.dUdt[2].mass)/(2.0*dP);
		    dRdD_T[i, j*np+1] = (cell.dUdt[1].momentum.x - cell.dUdt[2].momentum.x)/(2.0*dP);
		    dRdD_T[i, j*np+2] = (cell.dUdt[1].momentum.y - cell.dUdt[2].momentum.y)/(2.0*dP);
		    dRdD_T[i, j*np+3] = (cell.dUdt[1].total_energy - cell.dUdt[2].total_energy)/(2.0*dP);
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

    if (gradientVerification) {
        double[] finiteDiffGradients; double J1; double J0;
        
        // read original grid in
        ensure_directory_is_present(make_path_name!"grid-original"(0));
        string gridFileName = make_file_name!"grid-original"(jobName, myblk.id, 0, gridFileExt = "gz");
        myblk.read_new_underlying_grid(gridFileName);
        myblk.sync_vertices_from_underlying_grid(0);

        foreach (bndary; myblk.bc) {
	    if (bndary.is_design_surface) {
		foreach (i; 0 .. bndary.num_cntrl_pts) {
                    double err;
		    if (i == 0) {
                        double J;
                        foreach (otherBndary; myblk.bc) {
                            if (otherBndary.group == "design") J = cost_function(otherBndary, 0);
                        }
                        writef("Original J(D, Q(D) = %.16f \n", J);
                    }
		    writeln("computing finite difference gradient for variable ", i+1, " out of ", ndvars, " variables");
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
		    auto P0 = bndary.bezier.B[i].refy;
                    double dP = DELTA; //max(ETA*abs(P0), 1.0e-10);                
		    bndary.bezier.B[i].refy = P0 + dP;
                    varID = "p-y-" ~ to!string(i);
		    J0 = finite_difference_grad(jobName, last_tindx, myblk, bndary, NonFixedBoundaryList, varID);
		    writef("perturbed J(D, Q(D) = %.16f , h = %.16f \n", J0, dP);
		    bndary.bezier.B[i].refy = P0 - dP;
                    varID = "m-y-" ~ to!string(i);
		    J1 = finite_difference_grad(jobName, last_tindx, myblk, bndary, NonFixedBoundaryList, varID);
		    finiteDiffGradients ~= (J0 - J1)/(2.0*dP);
		    bndary.bezier.B[i].refy = P0;
                    writef("perturbed J(D, Q(D) = %.16f , h = %.16f \n", J1, dP);

                    err = ((adjointGradients[i] - finiteDiffGradients[i])/finiteDiffGradients[i]) * 100.0;
		    writeln("finite-difference gradient for variable ", i+1, ": ", finiteDiffGradients[i], ", % error: ", abs(err));
                }
	    }
	}
    }
    writeln("Simulation complete.");
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

