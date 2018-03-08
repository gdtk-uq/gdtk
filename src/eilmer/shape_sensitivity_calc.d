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

enum ghost_cell_start_id = 1_000_000_000;
immutable double ESSENTIALLY_ZERO = 1.0e-50;

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

    // --------------
    // --------------
    // ADJOINT SOLVER
    //---------------
    // --------------

    FluidBlock myblk = localFluidBlocks[0]; // currently only compatible with single block simulations
    size_t orderOfJacobian = myblk.myConfig.interpolation_order;
    
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
    build_flow_jacobian(L, myblk, ndim, np, orderOfJacobian, EPSILON, MU);

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
            bndary.bezier = optimiseBezierPoints(bndary.surfacePoints, bndary.num_cntrl_pts, bndary.ts);
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

double finite_difference_grad(string jobName, int last_tindx, FluidBlock blk, BoundaryCondition bndary, string[] NonFixedBoundaryList, string varID) {
    size_t gtl = 1;
    double grad;
    grad = (bndary.vertices[$-1].pos[1].y - bndary.vertices[0].pos[0].y)/(bndary.vertices[$-1].pos[0].x - bndary.vertices[0].pos[0].x);
    foreach(j, vtx; bndary.vertices) {
	vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
	vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
    }
    inverse_distance_weighting(blk, NonFixedBoundaryList, gtl);

    foreach(j, vtx; blk.vertices) {
        //writef("%d %.16f %.16f \n", vtx.id, vtx.pos[0].x, vtx.pos[0].y);
	vtx.pos[0].refx = vtx.pos[gtl].x;
	vtx.pos[0].refy = vtx.pos[gtl].y;
    }

    // save mesh
    blk.sync_vertices_to_underlying_grid(0);
    ensure_directory_is_present(make_path_name!"grid-perturb"(0));
    auto fileName = make_file_name!"grid-perturb"(jobName, blk.id, 0, gridFileExt = "gz");
    blk.write_underlying_grid(fileName);
    
    // run simulation
    string command = "bash run-flow-solver-on-perturbed-mesh.sh";
    auto output = executeShell(command);
    
    // store simulation diagnostics file
    string commandCpy0 = "cp e4sss.diagnostics.dat adjoint/" ~ varID;
    auto outputCpy0 = executeShell(commandCpy0);

    string commandCpy1 = "cp -r plot adjoint/plot-" ~ varID;
    auto outputCpy1 = executeShell(commandCpy1);

    // read perturbed solution
    blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    
    // compute cost function
    double J;
    foreach (otherBndary; blk.bc) {
        if (otherBndary.group == "design") J = cost_function(otherBndary, 0);
    }
        
    // read original grid in
    ensure_directory_is_present(make_path_name!"grid-original"(0));
    string gridFileName = make_file_name!"grid-original"(jobName, blk.id, 0, gridFileExt = "gz");
    blk.read_new_underlying_grid(gridFileName);
    blk.sync_vertices_from_underlying_grid(0);
        
    // clear old simulation files
    string command0 = "bash clear.sh";
    output = executeShell(command0);

    return J;
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
    
    // some data objects used in forming the Jacobian
    immutable size_t MAX_PERTURBED_INTERFACES = 40; 

    FVCell cellOrig; FVCell cellL; FVCell cellR;
    FVInterface[MAX_PERTURBED_INTERFACES] ifaceOrig;
    FVInterface[MAX_PERTURBED_INTERFACES] ifacePp;
    FVInterface[MAX_PERTURBED_INTERFACES] ifacePm;
    double h; double diff;

    cellOrig = new FVCell(dedicatedConfig[blk.id]);
    foreach(i; 0..MAX_PERTURBED_INTERFACES) {
        ifaceOrig[i] = new FVInterface(dedicatedConfig[blk.id], false);
        ifacePp[i] = new FVInterface(dedicatedConfig[blk.id], false);
        ifacePm[i] = new FVInterface(dedicatedConfig[blk.id], false);
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
                    
                    if (JacEntry != 0.0) {
                        if (ip == 0 && jp == 0) {
                            //writef("%d ", c.id);
                            count += 1;
                        }
                        //writeln("col: ", I, ", row: ", J, ", p-cell: ", cell.id, ", e-cell: ", c.id, ", ip: ", ip, ", jp: ", jp, ", val: ", JacEntry);
                        blk.aa[J] ~= JacEntry;
                        blk.ja[J] ~= I; //J;
                    }
                }
            }
        }
        //writef(" total effected cells: %d cell in stencil: %d \n", count, cell.jacobian_stencil.length);
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


void compute_perturbed_flux(FluidBlock blk, size_t orderOfJacobian, FVCell[] cell_list, FVInterface[] iface_list, FVInterface[] ifaceP_list) {
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
                    c.gradients.venkat_limit(c.cell_cloud, c.ws, blk.myConfig);
                    break;
                } // end switch
            } // end foreach c
        } // end if interpolation_order > 1

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
    codeStr ~= "compute_perturbed_flux(blk, orderOfJacobian, cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePm);"; 
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(blk, orderOfJacobian, cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePp);"; 
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

double cost_function(BoundaryCondition bndary, int gtl) {
    // computes normal force acting on design surface
    double J = 0.0;
    foreach (i, f; bndary.faces) {
        FVCell cell;
        // pick-up internal cell
        if (bndary.outsigns[i] == 1) {
            cell = f.left_cell;
        } else {
            cell = f.right_cell;
        }
        J += cell.fs.gas.p*f.area[gtl];
    }
    return J;
}

void cost_function_sensitivity(ref double[] dJdV, FluidBlock blk, size_t np, double EPSILON, double MU) {
    // for the moment we just have a hard-coded cost function    
    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;

    foreach(cell; blk.cells) {
        foreach(i; 0..np) {
            dJdV[cell.id*np + i] = 0.0;
        }
    }
    
    foreach (bndary; blk.bc) {
        if (bndary.group == "design") {            
	    foreach (i, f; bndary.faces) {
		FVCell cell;
                double origP;
                if (bndary.outsigns[i] == 1) {
		    cell = f.left_cell;
		} else {
		    cell = f.right_cell;
		}

                double Jm; double Jp; double h;
                // pressure
                h = (abs(cell.fs.gas.p) + MU) * EPSILON;
                origP = cell.fs.gas.p;
                cell.fs.gas.p += h;
                Jp = cost_function(bndary, 0);
                cell.fs.gas.p = origP;
                cell.fs.gas.p -= h;
                Jm = cost_function(bndary, 0);
                cell.fs.gas.p = origP;
                dJdV[cell.id*np + 3] = (Jp-Jm)/(2.0*h);
            }
        }
    }
}

double[] adjoint_solver(SMatrix globalJacobianT, double[] dJdV, SMatrix foJac, SMatrix m, FluidBlock blk, size_t np) {
    size_t ncells = blk.cells.length;
    // restarted-GMRES settings
    int maxInnerIters = GlobalConfig.sscOptions.gmresRestartInterval;
    int maxOuterIters = 1000;
    int nIter = 0;
    double normRef = 0.0;
    double normNew = 0.0;
    double residTol = GlobalConfig.sscOptions.stopOnRelativeGlobalResidual;
    double[] psi0;
    double[] psiN;
    foreach(i; 0..np*ncells) psi0 ~= 1.0;
    foreach(i; 0..np*ncells) psiN ~= 1.0;
    double[] residVec;
    residVec.length = dJdV.length;

    // compute ILU[p] for preconditioning
    //decompILUp(m, 1);
    decompILU0(m);

    // compute reference norm
    multiply(globalJacobianT, psi0, residVec);
    foreach (i; 0..np*ncells) residVec[i] = dJdV[i] - residVec[i];
    foreach (i; 0..np*ncells) normRef += residVec[i]*residVec[i];
    normRef = sqrt(fabs(normRef));
    auto gws = GMRESWorkSpace(psi0.length, maxInnerIters);
    
    while (nIter < maxOuterIters) {
        // compute psi
        rpcGMRES(globalJacobianT, m, dJdV, psi0, psiN, maxInnerIters, residTol, gws);
            
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

    destroy(m);
    GC.minimize();
    
    return psiN;
}

void write_adjoint_variables_to_file(FluidBlock blk, double[] psi, size_t np, string jobName) {
    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;
    // write out adjoint variables in VTK-format
    if (blk.grid_type == Grid_t.structured_grid) throw new FlowSolverException("Shape Sensitivity Calculator not implemented for structured grids, yet."); 
    if (blk.grid_type == Grid_t.unstructured_grid) {
        auto ublk = cast(UFluidBlock) blk; 
        auto outFile = File("adjointVars.vtk", "w");
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
            outFile.writef("%.16f \n", psi[np*i]);
        }

        outFile.writef("SCALARS adjoint_velx double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", psi[np*i+1]);
        }

        outFile.writef("SCALARS adjoint_vely double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) {
            outFile.writef("%.16f \n", psi[np*i+2]);
        }

        outFile.writef("SCALARS adjoint_pressure double \n");
        outFile.writef("LOOKUP_TABLE default \n");
        foreach(i; 0..ncells) { 
            outFile.writef("%.16f \n", psi[np*i+3]);
        }

    }
}

void evalRHS(double pseudoSimTime, int ftl, int gtl, bool with_k_omega, FluidBlock blk)
{
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    
    exchange_ghost_cell_boundary_data(pseudoSimTime, gtl, ftl);
    blk.applyPreReconAction(pseudoSimTime, gtl, ftl);
    
    // We don't want to switch between flux calculator application while
    // doing the Frechet derivative, so we'll only search for shock points
    // at ftl = 0, which is when the F(U) evaluation is made.
    if ( ftl == 0 && (GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv ||
		      GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv)) {
        blk.detect_shock_points();
    }
    
    blk.convective_flux_phase0();
    blk.convective_flux_phase1();
    blk.applyPostConvFluxAction(pseudoSimTime, gtl, ftl);
    if (GlobalConfig.viscous) {
        blk.applyPreSpatialDerivActionAtBndryFaces(pseudoSimTime, gtl, ftl);
        blk.flow_property_spatial_derivatives(gtl); 
        blk.estimate_turbulence_viscosity();
        blk.viscous_flux();
        blk.applyPostDiffFluxAction(pseudoSimTime, gtl, ftl);
    }
    
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

   
void build_flow_jacobian(SMatrix L, FluidBlock blk, size_t ndim, size_t np, size_t orderOfJacobian, double EPSILON, double MU) {

    // build jacobian stencils
    if (GlobalConfig.viscous) construct_viscous_flow_jacobian_stencils(blk, orderOfJacobian);
    else construct_inviscid_flow_jacobian_stencils(blk, orderOfJacobian);
    
    // compute transpose Jacobian entries
    construct_flow_jacobian(blk, ndim, np, orderOfJacobian, EPSILON, MU);
    
    // build Sparse transposed Jacobian
    size_t ia = 0;
    foreach(i; 0 .. blk.cells.length*np) { // 0..nrows
        L.aa ~= blk.aa[i];
        L.ja ~= blk.ja[i];
        L.ia ~= ia;
        ia += blk.aa[i].length;
    }
    // clear local Jacobian memory
    blk.aa = [];
    blk.ja = [];
    L.ia ~= L.aa.length;

    // clear jacobian stencils
    foreach ( cell; blk.cells) {
        cell.jacobian_cell_stencil = [];
        cell.jacobian_face_stencil = [];
    }

}


//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                   STEADY-STATE SOLVER PRECONDITIONER
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------

void sss_preconditioner(FluidBlock blk, size_t np, double[] aa, double EPSILON, double MU, size_t orderOfJacobian=1) {

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
    
    // some data objects used in forming the preconditioner
    immutable size_t MAX_PERTURBED_INTERFACES = 10; FVCell cellOrig;
    FVInterface[MAX_PERTURBED_INTERFACES] ifaceOrig;
    FVInterface[MAX_PERTURBED_INTERFACES] ifacePp;
    FVInterface[MAX_PERTURBED_INTERFACES] ifacePm;
    double h; double diff;

    // initialise objects
    cellOrig = new FVCell(dedicatedConfig[blk.id]);
    foreach(i; 0..MAX_PERTURBED_INTERFACES) {
        ifaceOrig[i] = new FVInterface(dedicatedConfig[blk.id], false);
        ifacePp[i] = new FVInterface(dedicatedConfig[blk.id], false);
        ifacePm[i] = new FVInterface(dedicatedConfig[blk.id], false);
    }

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
                if ( ip == jp) {
                    integral = 0.0;
                    foreach(fi, iface; cell.iface) {
                        integral -= cell.outsign[fi] * iface.dFdU[ip][jp] * iface.area[0]; // gtl=0
                    }
                    double entry = volInv * integral;                    
                    aa[cell.id*np+jp] = entry;
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

    // multiply by transform matrix diagonal (transforming primitive to conservative form)
    foreach ( cell; blk.cells) {
        aa[cell.id*np] *= 1.0;
        aa[cell.id*np+1] *= 1.0/cell.fs.gas.rho;
        aa[cell.id*np+2] *= 1.0/cell.fs.gas.rho;
        aa[cell.id*np+3] *= cell.fs.gas.p/(cell.fs.gas.rho * cell.fs.gas.u); // ratio of specific heats minus 1
    }
}
