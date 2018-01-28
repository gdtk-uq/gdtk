/** adjoint_gradient.d
 * 
 * Code to compute shape sensitivities via a discrete adjoint method.
 * Solvers/Calculators included:
 * 1. adjoint variable solver.
 * 2. mesh sensitivity calculators:
 * 2.1 sensitivity of cell residuals w.r.t. mesh point perturbations,
 * 2.2 sensitivity of mesh points w.r.t. design variable perturbations.
 * 3. gradient calculator to compute the shape sensitivities.
 *
 * Author: Kyle D.
 * Date: 2017-09-18
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
    version(mpi_parallel) {
        assert(0, "Adjoint solver is not MPI parallel, yet.");
    }

    string msg = "Usage:                              Comment:\n";
    msg       ~= "e4adjoint  [--job=<string>]            name of job\n";
    msg       ~= "           [--max-cpus=<int>]          defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine\n";
    msg       ~= "           [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "           [--epsilon=<double>]    \n";
    msg       ~= "           [--mu=<double>]    \n";
    msg       ~= "           [--epsilon-finite-diff=<double>]    \n";
    msg       ~= "           [--with-finite-diff-verification=<bool>]    \n";
    msg       ~= "           [--help]                    writes this message\n";
    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        write(msg);
        exit(1);
    }
    string jobName = "";
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    double EPSILON = 1.0e-04;
    double MU = 1.0e-04;
    double ETA = 1.0e-04;
    double OMEGA = 1.0e-04;
    bool withFiniteDiffVerification = true;
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "max-cpus", &maxCPUs,
               "max-wall-clock", &maxWallClock,
               "epsilon", &EPSILON,
               "mu", &MU,
	       "omega", &OMEGA, 
               "eta", &ETA,
               "verification", &withFiniteDiffVerification,
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
    init_simulation(last_tindx, 0, maxCPUs, 1, maxWallClock);

    writeln("Finite Difference Parameters:");
    writeln("EPSILON  = ", EPSILON);
    writeln("MU  = ", MU);
    writeln("OMEGA = ", OMEGA);
    writeln("ETA = ", ETA);

    // perform some config checks
    if ( GlobalConfig.interpolation_order > 1 &&
         GlobalConfig.suppress_reconstruction_at_boundaries == false) {
        writeln("WARNING:");
        writeln("   suppress_reconstruction_at_boundaries is set to false.");
        writeln("   This setting must be true when using the adjoint solver.");
        writeln("   Its use will likely cause errorneous values for boundary fluxes.");
        writeln("   Continuing with simulation anyway.");
        writeln("END WARNING.");
    }
    
    if (GlobalConfig.diffuseWallBCsOnInit) {
        writeln("WARNING:");
        writeln("   diffuse_wall_bcs_on_init is set to true.");
        writeln("   This setting must be false when using the adjoint solver.");
        writeln("   Its use will likely cause errorneous values for gradients.");
        writeln("   Continuing with simulation anyway.");
        writeln("END WARNING.");
    }
    
    // save a copy of the original mesh
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
    
    // ------------------------------------
    // 1.a construct flow Jacobian stencils
    //-------------------------------------
    if (GlobalConfig.viscous) construct_viscous_flow_jacobian_stencils(myblk);
    else construct_inviscid_flow_jacobian_stencils(myblk);
        
    // ------------------------------------------------------------
    // 1.b construct flow Jacobian transpose via finite differences
    // ------------------------------------------------------------
    size_t ndim = GlobalConfig.dimensions;
    bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);
    size_t np; // number of primitive variables
    if (GlobalConfig.dimensions == 2) {
        if (with_k_omega) np = 6; // rho, momentum(x,y), total energy, tke, omega
        else np = 4; // rho, momentum(x,y), total energy
    }
    else if (GlobalConfig.dimensions == 3) {
        if (with_k_omega) np = 7; // rho, momentum(x,y,z), total energy, tke, omega
        else np = 5; // rho, momentum(x,y,z), total energy, tke, omega
    }

    construct_flow_jacobian(myblk, ndim, np, EPSILON, MU);
    
    // --------------------------------------------
    // 1.c construct global flow Jacobian transpose
    // --------------------------------------------
    SMatrix globalJacobianT = new SMatrix();    
    size_t ia = 0;
    foreach(i; 0 .. myblk.cells.length*np) { // 0..nrows
        globalJacobianT.aa ~= myblk.aa[i];
        globalJacobianT.ja ~= myblk.ja[i];
        globalJacobianT.ia ~= ia;
        ia += myblk.aa[i].length;
    }
    // clear local Jacobian memory
    myblk.aa = [];
    myblk.ja = [];
    globalJacobianT.ia ~= globalJacobianT.aa.length;
    
    // --------------------------------------
    // 1.d construct cost function sensitvity
    // --------------------------------------
    double[] dJdV; // const function sensitivity
    double[] p_target;
    construct_cost_function_sensitivity(myblk, dJdV, p_target);
    
    // ------------------------
    // 1.e solve adjoint system 
    // ------------------------
    double[] psi;
    psi = adjoint_solver(globalJacobianT, dJdV, myblk, np);
    
    // -------------------------------
    // 1.f write out adjoint variables 
    // -------------------------------
    write_adjoint_variables_to_file(myblk, psi, np, jobName);
    
    // clear the global Jacobian from memory
    destroy(globalJacobianT);
    GC.minimize();

    // -------------------
    // -------------------
    // GRADIENT CALCULATOR
    //--------------------
    // -------------------

    // ------------------------------------
    // 2.a construct flow Jacobian stencils 
    // ------------------------------------
    construct_mesh_jacobian_stencils(myblk);
    
    // ------------------------------------------------------------
    // 2.b construct mesh Jacobian transpose via finite differences
    // ------------------------------------------------------------
    construct_mesh_jacobian(myblk, ndim, np, OMEGA, MU);
    
    // --------------------------------------------
    // 2.c construct global mesh Jacobian transpose
    // --------------------------------------------
    SMatrix globaldRdXT = new SMatrix();    
    // we must do this in serial
    ia = 0;
    foreach(i; 0 .. myblk.cells.length*np) { // 0..nrows
        globaldRdXT.aa ~= myblk.aa[i];
        globaldRdXT.ja ~= myblk.ja[i];
        globaldRdXT.ia ~= ia;
        ia += myblk.aa[i].length;
    }
    // clear local Jacobian memory
    myblk.aa = [];
    myblk.ja = [];

    globaldRdXT.ia ~= globaldRdXT.aa.length;

    // ------------------------------------------------------------------------------------------------
    // 2.d construct mesh point sensitivity w.r.t. design variable perturbations via finite-differences
    // ------------------------------------------------------------------------------------------------
    
    /++ 
     + To compute dXdD, we need to determine the current design variables, since we will perturb these to find the effect
     + on the mesh via finite-dfferences. Currently we parameterise surfaces to be optimised via bezier curves, 
     + that have been optimised to best fit the user specified design surfaces/boundaries. The mesh points along the 
     + design surfaces are used as the supplied data points for the bezier curve parameterisation code. 
    ++/

    // parameterise the user defined design surfaces (for now let's just tag the boundaries here)
    string[] designSurfaces;
    designSurfaces ~= "design";
    int[string] nCntrlPtsList;
    nCntrlPtsList["design"] = 5; 
    size_t ndvars = 0; // number of design variables
        
    foreach (bndary; myblk.bc) {
        if (designSurfaces.canFind(bndary.group)) {
	    // gather boundary vertices
            size_t[] vtx_id_list;
            foreach(face; bndary.faces) {
                foreach(vtx; face.vtx) {
                    if (vtx_id_list.canFind(vtx.id) == false) {
                        bndary.vertices ~= vtx;
                        bndary.surfacePoints ~= Vector3(vtx.pos[0].refx, vtx.pos[0].refy, vtx.pos[0].refz);
			vtx_id_list ~= vtx.id;
			myblk.boundaryVtxIndexList ~= vtx.id;
                    }
                }
            }
	    
	    
            // compute bezier control points
            bndary.nCntrlPts = nCntrlPtsList[bndary.group];
	    bndary.bezier = optimiseBezierPoints(bndary.surfacePoints, bndary.nCntrlPts, bndary.ts);
            ndvars += bndary.nCntrlPts*2;
        }
	else {
	    size_t[] vtx_id_list;
            foreach(face; bndary.faces) {
                foreach(vtx; face.vtx) {
                    if (vtx_id_list.canFind(vtx.id) == false) {
                        bndary.vertices ~= vtx;
			vtx_id_list ~= vtx.id;
			myblk.boundaryVtxIndexList ~= vtx.id;
                    }
                }
            }
	}
    }

    size_t nvertices = myblk.vertices.length;
    Matrix dXdD_T;
    dXdD_T = new Matrix(ndvars, nvertices*ndim);
    foreach (bndary; myblk.bc) {
        if (designSurfaces.canFind(bndary.group)) {
	    foreach (i; 0 .. bndary.bezier.B.length) {
		// could loop over all x, y, z coordinates here -- for now just use y-coordinates
		foreach(otherBndary; myblk.bc) {
		    foreach(j, vtx; otherBndary.vertices) {
			    vtx.pos[1].refx = vtx.pos[0].x;
			    vtx.pos[1].refy = vtx.pos[0].y;
			    vtx.pos[2].refx = vtx.pos[0].x;
			    vtx.pos[2].refy = vtx.pos[0].y;
		    }
		}
		size_t gtl;
		auto P0 = bndary.bezier.B[i].y;
		double dP = ETA;
		gtl = 1;
		bndary.bezier.B[i].refy = P0 + dP;
		foreach(j, vtx; bndary.vertices) {
		    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
		    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
		}
		inverse_distance_weighting(myblk, gtl);
		gtl = 2;
		bndary.bezier.B[i].refy = P0 - dP;
		foreach(j, vtx; bndary.vertices) {
		    vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
		    vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
		}
		inverse_distance_weighting(myblk, gtl);
		foreach(j, vtx; myblk.vertices) {
		    dXdD_T[i*ndim, j*ndim] = 0.0;
		    dXdD_T[i*ndim, j*ndim+1] = 0.0;
		    dXdD_T[i*ndim+1, j*ndim] = (vtx.pos[1].x - vtx.pos[2].x)/(2.0*dP);
		    dXdD_T[i*ndim+1, j*ndim+1] = (vtx.pos[1].y - vtx.pos[2].y)/(2.0*dP);
		}
		bndary.bezier.B[i].refy = P0;
	    }
	}
    }

    // ------------------------------
    // 3. Compute shape sensitivities
    // ------------------------------
    writeln("compute shape sensitivities");
    // temp matrix multiplication
    size_t ncells = myblk.cells.length;
    Matrix tempMatrix;
    tempMatrix = new Matrix(ndvars, ncells*np);
    for (int i = 0; i < ndvars; i++) {
        for (int j = 0; j < np*ncells; j++) {
            tempMatrix[i,j] = 0;
            for (int k = 0; k < nvertices*ndim; k++) {
                tempMatrix[i,j] += dXdD_T[i,k]*globaldRdXT[k,j];
            }
        }
    }
    double[] adjointGradients;
    adjointGradients.length = tempMatrix.nrows;
    dot(tempMatrix, psi, adjointGradients);    
    writeln("adjoint gradients = ", adjointGradients);

    // clear the sensitivity matrices from memory
    destroy(globaldRdXT);
    destroy(dXdD_T);
    destroy(tempMatrix);
    GC.minimize();
    
    // -----------------------------------------------------
    // Finite difference verification
    // -----------------------------------------------------

    if (withFiniteDiffVerification) {
        double[] finiteDiffGradients; double J1; double J0;
        
        // read original grid in
        ensure_directory_is_present(make_path_name!"grid-original"(0));
        string gridFileName = make_file_name!"grid-original"(jobName, myblk.id, 0, gridFileExt = "gz");
        myblk.read_new_underlying_grid(gridFileName);
        myblk.sync_vertices_from_underlying_grid(0);

        foreach (bndary; myblk.bc) {
	    if (designSurfaces.canFind(bndary.group)) {
		foreach (i; 0 .. bndary.bezier.B.length) {
		    writeln("computing finite difference gradient for variable ", i+1, " out of ", bndary.bezier.B.length, " variables");
		    // could loop over all x, y, z coordinates here -- for now just use y-coordinates
		    auto P0 = bndary.bezier.B[i].refy;
		    double dP = ETA;                
		    bndary.bezier.B[i].refy = P0 + dP;
		    J0 = finite_difference_grad(jobName, last_tindx, myblk, bndary, p_target);
		    bndary.bezier.B[i].refy = P0 - dP;
		    J1 = finite_difference_grad(jobName, last_tindx, myblk, bndary, p_target);
		    finiteDiffGradients ~= (J0 - J1)/(2.0*dP);
		    bndary.bezier.B[i].refy = P0;
		}
	    }
	}
        writeln("finite difference gradients = ", finiteDiffGradients);
	foreach (bndary; myblk.bc) {
	    if (designSurfaces.canFind(bndary.group)) {
		foreach (i; 0 .. bndary.bezier.B.length) {
		    double err = ((adjointGradients[i*2+1] - finiteDiffGradients[i])/finiteDiffGradients[i]) * 100.0;
		    writeln("% error for variable ", i+1, ": ", abs(err));
		}
	    }
	}
    }
    writeln("Simulation complete.");
}

double finite_difference_grad(string jobName, int last_tindx, FluidBlock blk, BoundaryCondition bndary, double[] p_target) {
    size_t gtl = 1;
    foreach(otherBndary; blk.bc) {
	foreach(j, vtx; otherBndary.vertices) {
	    vtx.pos[1].refx = vtx.pos[0].x;
	    vtx.pos[1].refy = vtx.pos[0].y;
	    vtx.pos[2].refx = vtx.pos[0].x;
	    vtx.pos[2].refy = vtx.pos[0].y;
	}
    }
    foreach(j, vtx; bndary.vertices) {
	vtx.pos[gtl].refx = bndary.bezier(bndary.ts[j]).x;
	vtx.pos[gtl].refy = bndary.bezier(bndary.ts[j]).y;
    }
    inverse_distance_weighting(blk, gtl);

    foreach(j, vtx; blk.vertices) {
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
    
    // read perturbed solution
    blk.read_solution(make_file_name!"flow"(jobName ~ "-perturb", blk.id, last_tindx, flowFileExt), false);
    
    // compute cost function
    double J = 0.0;
    foreach (i, cell; blk.cells) {
        J += 0.5*(cell.fs.gas.p - p_target[i])*(cell.fs.gas.p - p_target[i]);
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

void construct_flow_jacobian(FluidBlock blk, size_t ndim, size_t np, double EPSILON, double MU) {
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

void construct_inviscid_flow_jacobian_stencils(FluidBlock blk) {
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
    if (blk.myConfig.interpolation_order < 2) { 
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

            if (blk.grid_type == Grid_t.structured_grid) {

                // collect cells
                refs_unordered ~= pcell;
                pos_array[pcell.id] = refs_unordered.length-1;
                cell_ids ~= pcell.id;

                auto sblk = cast(SFluidBlock) blk;
                size_t[3] ijk = sblk.cell_id_to_ijk_indices(pcell.id);
                size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2]; 
                FVCell[] stencil_cells;
                stencil_cells ~= sblk.get_cell(i-2, j, k);
                stencil_cells ~= sblk.get_cell(i-1, j, k);
                stencil_cells ~= sblk.get_cell(i+1, j, k);
                stencil_cells ~= sblk.get_cell(i+2, j, k);
                stencil_cells ~= sblk.get_cell(i, j-2, k);
                stencil_cells ~= sblk.get_cell(i, j-1, k);
                stencil_cells ~= sblk.get_cell(i, j+1, k);
                stencil_cells ~= sblk.get_cell(i, j+2, k);
                foreach(cell; stencil_cells) {
                    // add cell to cell stencil
                    if (cell.id < ghost_cell_start_id) {
                        refs_unordered ~= cell;
                        pos_array[cell.id] = refs_unordered.length-1;
                        cell_ids ~= cell.id;
                    }
                }
                // collect faces
                stencil_cells = [];
                stencil_cells ~= sblk.get_cell(i, j, k);
                stencil_cells ~= sblk.get_cell(i-1, j, k);
                stencil_cells ~= sblk.get_cell(i+1, j, k);
                stencil_cells ~= sblk.get_cell(i, j+1, k);
                stencil_cells ~= sblk.get_cell(i, j-1, k);
                foreach(cell; stencil_cells) {
                    foreach(face; cell.iface) {
                        if (face_ids.canFind(face.id) == false && cell.id < ghost_cell_start_id)
                            pcell.jacobian_face_stencil ~= face; face_ids ~= face.id;
                    }
                }
            }
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

void construct_viscous_flow_jacobian_stencils(FluidBlock blk) {
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

        if (blk.grid_type == Grid_t.structured_grid) {
            // first loop around and store faces in stencil
            foreach(f; pcell.iface) {
                // store face
                pcell.jacobian_face_stencil ~= f; face_ids ~= f.id;
                
                // loop around neighbour
                if (f.left_cell.id != pcell.id && f.left_cell.id < ghost_cell_start_id) {
                    foreach(face; f.left_cell.iface) {
                        if (face_ids.canFind(face.id) == false )
                            { pcell.jacobian_face_stencil ~= face; face_ids ~= face.id; }
                    }
                }
                if (f.right_cell.id != pcell.id && f.right_cell.id < ghost_cell_start_id) {
                    foreach(face; f.right_cell.iface) {
                        if (face_ids.canFind(face.id) == false )
                            { pcell.jacobian_face_stencil ~= face; face_ids ~= face.id; }
                    }
                } else continue;
            }
            
            // now loop through face stencil and add left, and right cells
            foreach(f; pcell.jacobian_face_stencil) {
                // store (non-ghost) neighbour cells in cells stencil
                if (cell_ids.canFind(f.left_cell.id) == false && f.left_cell.id < ghost_cell_start_id) {
                    refs_unordered ~= f.left_cell;
                    pos_array[f.left_cell.id] = refs_unordered.length-1;
                    cell_ids ~= f.left_cell.id;
                }
                if (cell_ids.canFind(f.right_cell.id) == false && f.right_cell.id < ghost_cell_start_id) {
                    refs_unordered ~= f.right_cell;
                    pos_array[f.right_cell.id] = refs_unordered.length-1;
                    cell_ids ~= f.right_cell.id;
                } else continue;
            } // end foreach face
        }
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

void construct_mesh_jacobian_stencils(FluidBlock blk) {
    /++
     
     This is the stencil of surrounding cells (and faces) that are effected by a perturbation in a 
     vertex, for viscous flow.
     
     We must have the cells stored in id order, for efficient storage of the Jacobian.
     
     This stencil is a superposition of the flow Jacobian stencils for the cells that share the perturbed vertex.
     The general formula is to loop around the cells that share the perturbed vertex, and collect the cells, and faces
     from that cells flow Jacobian stencil.

     pvtx = perturbed vertex.
     
     ++/

    foreach(i, pvtx; blk.vertices) {
        FVCell[] refs_ordered;
        FVCell[] refs_unordered;
        size_t[size_t] pos_array; // used to identify where the cell is in the unordered list
        size_t[] cell_ids;
        size_t[] face_ids;

        foreach (cid; blk.cellIndexListPerVertex[pvtx.id]) {

            // collect cells from cells flow Jacobian stencil
            foreach(cell; blk.cells[cid].jacobian_cell_stencil) {
                if (cell_ids.canFind(cell.id) == false) {
                    refs_unordered ~= cell;
                    pos_array[cell.id] = refs_unordered.length-1;
                    cell_ids ~= cell.id;
                }
            }

            // collect faces from faces flow Jacobian stencil
            foreach(face; blk.cells[cid].jacobian_face_stencil) {
                if (face_ids.canFind(face.id) == false) {
                    pvtx.jacobian_face_stencil ~= face;
                    face_ids ~= face.id;
                }
            }
        }
        // finally sort ids, and store sorted cell references
        cell_ids.sort();
        foreach(id; cell_ids) {
            refs_ordered ~= refs_unordered[pos_array[id]];
        }
        pvtx.jacobian_cell_stencil ~= refs_ordered;
    }
}

void construct_mesh_jacobian(FluidBlock blk, size_t ndim, size_t np, double OMEGA, double MU) {
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

    // perform a flux computation
    blk.clear_fluxes_of_conserved_quantities();
    blk.applyPreReconAction(0.0, 0, 0);
    blk.convective_flux_phase0();
    blk.convective_flux_phase1();
    blk.applyPostConvFluxAction(0.0, 0, 0);
    if (GlobalConfig.viscous) {
        blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
        blk.flow_property_spatial_derivatives(0); 
        blk.estimate_turbulence_viscosity();
        blk.viscous_flux();
        blk.applyPostDiffFluxAction(0.0, 0, 0);
    }
    
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
        foreach(ci, c; vtx.jacobian_cell_stencil) {
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
        foreach (iface;  vtx.jacobian_face_stencil) {
            foreach (i; 0..iface.dFdU.length) {
                foreach (j; 0..iface.dFdU[i].length) {
                    iface.dFdU[i][j] = 0.0;
                }
            }
        }
    } // end foreach vtx
}

void compute_perturbed_flux(FluidBlock blk, FVCell[] cell_list, FVInterface[] iface_list, FVInterface[] ifaceP_list) {
    /++
     This method computes a perturbed flux at a given interface. 

     NB. The accuracy of the adjoint solver is largely dependent on how well this
     routine replicates the flow solver flux update procedures.
     ++/

    foreach(iface; iface_list) iface.F.clear_values();

    // currently apply bc's for every cell -- TODO: only update for perturbed boundary cells
    blk.applyPreReconAction(0.0, 0, 0);  // assume sim_time = 0.0, gtl = 0, ftl = 0
    
    // Convective flux update
    if (blk.grid_type == Grid_t.structured_grid) {
        foreach(iface; iface_list) { 
            // we need to cast an SFluidBlock here to reach some methods and data
            auto sblk = cast(SFluidBlock) blk;
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
    }
    else { // unstructured grid
	// compute new gradients for all cells in the stencil
        if (blk.myConfig.interpolation_order > 1) {
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

        // compute flux
        foreach(iface; iface_list) {
            auto ublk = cast(UFluidBlock) blk;
            ublk.lsq.interp_both(iface, 0, ublk.Lft, ublk.Rght); // gtl assumed 0
            iface.fs.copy_average_values_from(ublk.Lft, ublk.Rght);
            compute_interface_flux(ublk.Lft, ublk.Rght, iface, ublk.myConfig, ublk.omegaz);
        }
    }

    // currently apply bc's for every cell -- TODO: only update for perturbed boundary cells
    blk.applyPostConvFluxAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0

    if (GlobalConfig.viscous) {

        // we should have the least-squares weights up to date
        //foreach(iface; iface_list) {
        //    iface.grad.set_up_workspace_leastsq(iface.cloud_pos, iface.pos, false, iface.ws_grad);
        //}
        
        // currently apply bc's for every cell -- TODO: only update for perturbed boundary cells
        blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);

        
        // only for least-squares at faces
        foreach(iface; iface_list) {
            iface.grad.gradients_leastsq(iface.cloud_fs, iface.cloud_pos, iface.ws_grad); // blk.flow_property_spatial_derivatives(0); 
        }

        //writeln("BEFORE");
        //foreach(cell; cell_list) writeln(cell);
        // estimate turbulent viscosity
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
            iface.viscous_flux_calc(); // blk.viscous_flux();
        }
        
        // currently apply bc's for every cell -- TODO: only update for perturbed boundary cells
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
    codeStr ~= "compute_perturbed_flux(blk, cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePm);"; 
    codeStr ~= "cell.copy_values_from(cellOrig, CopyDataOption.all);";
    // ------------------ positive perturbation ------------------
    codeStr ~= "cell.fs."~varName~" += h;";
    if ( includeThermoUpdate ) {
        codeStr ~= "blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);";
    }
    codeStr ~= "compute_perturbed_flux(blk,  cell.jacobian_cell_stencil, cell.jacobian_face_stencil, ifacePp);"; 
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

string computeCellVolumeSensitivity(string varName)
{
    string codeStr;
    codeStr ~= "h = OMEGA;";
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
    codeStr ~= "h = OMEGA;";
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
    codeStr ~= "h = OMEGA;";
    codeStr ~= "foreach (i, iface; vtx.jacobian_face_stencil) {";
    codeStr ~= "ifaceOrig[i].copy_values_from(iface, CopyDataOption.all);";
    codeStr ~= "}";
    // ------------------ negative perturbation ------------------
    codeStr ~= "vtx."~varName~" -= h;";
    codeStr ~= "update_geometry(blk, vtx);";
    codeStr ~= "compute_perturbed_flux(blk, vtx.jacobian_cell_stencil, vtx.jacobian_face_stencil, ifacePm);";
    codeStr ~= "vtx."~varName~" += h;";
    // ------------------ positive perturbation ------------------
    codeStr ~= "vtx."~varName~" += h;";
    codeStr ~= "update_geometry(blk, vtx);";
    codeStr ~= "compute_perturbed_flux(blk, vtx.jacobian_cell_stencil, vtx.jacobian_face_stencil, ifacePp);"; 
    codeStr ~= "vtx."~varName~" -= h;";
    // ------------------ compute interface flux derivatives ------------------
    codeStr ~= "foreach (i, iface; vtx.jacobian_face_stencil) {";
    codeStr ~= "diff = ifacePp[i].F.mass - ifacePm[i].F.mass;";
    codeStr ~= "iface.dFdU[0][" ~ posInArray ~ "] = diff/(2.0*h);";         
    codeStr ~= "diff = ifacePp[i].F.momentum.x - ifacePm[i].F.momentum.x;";
    codeStr ~= "iface.dFdU[1][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp[i].F.momentum.y - ifacePm[i].F.momentum.y;";
    codeStr ~= "iface.dFdU[2][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "diff = ifacePp[i].F.total_energy - ifacePm[i].F.total_energy;";
    codeStr ~= "iface.dFdU[3][" ~ posInArray ~ "] = diff/(2.0*h);";
    codeStr ~= "}";
    // ------------------ restore original values ------------------
    codeStr ~= "update_geometry(blk, vtx);";
    codeStr ~= "foreach (i, iface; vtx.jacobian_face_stencil) {";
    codeStr ~= "iface.copy_values_from(ifaceOrig[i], CopyDataOption.all);";
    codeStr ~= "}";
    return codeStr;
}

void update_geometry(FluidBlock blk, FVVertex vtx) {

    // update cell geometry
    foreach (cid; blk.cellIndexListPerVertex[vtx.id])
        blk.cells[cid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);

    //foreach (cell; vtx.jacobian_cell_stencil)
    //  cell.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);
    
    // update face geometry
    foreach (fid; blk.faceIndexListPerVertex[vtx.id])
        blk.faces[fid].update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);

    //foreach (face; vtx.jacobian_face_stencil)
    //  face.update_2D_geometric_data(0, dedicatedConfig[blk.id].axisymmetric);
    
    // for USG grid, update least-squares weights
    if (blk.grid_type == Grid_t.unstructured_grid) {
        auto ublk = cast(UFluidBlock) blk;
        foreach (cell; vtx.jacobian_cell_stencil)
            ublk.compute_least_squares_setup_for_cell(cell, 0, false);
    }
    
    if (GlobalConfig.viscous) {
        // we should have the spatial derivatives least-squares weights up to date as well
        foreach(iface; vtx.jacobian_face_stencil) {
            iface.grad.set_up_workspace_leastsq(iface.cloud_pos, iface.pos, false, iface.ws_grad);
        }
    }
}

void construct_cost_function_sensitivity(FluidBlock blk, ref double[] dJdV, ref double[] p_target) {
    // for the moment we just have a hard-coded cost function
    // J = 0.5*integral[over domain] (p-p*)^2

    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;

    // target pressure distribution saved in file target.dat
    auto file = File("target.dat", "r");
    foreach(line; 0..ncells) {
        auto lineContent = file.readln().strip();
        auto tokens = lineContent.split();
        p_target ~= to!double(tokens[8]);
    }
    writeln("target pressure imported");

    // hand-differentiated code
    foreach(i, cell; blk.cells) {
        dJdV ~= 0.0; // rho
        dJdV ~= 0.0; // x-vel
        dJdV ~= 0.0; // y-vel
        dJdV ~= 0.5*(2.0*cell.fs.gas.p - 2.0*p_target[i]); // P
    }
}

double[] adjoint_solver(SMatrix globalJacobianT, double[] dJdV, FluidBlock blk, size_t np) {
    size_t ncells = blk.cells.length;
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

    // compute ILU[p] for preconditioning
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
    
    return psiN;
}

void write_adjoint_variables_to_file(FluidBlock blk, double[] psi, size_t np, string jobName) {
    size_t ncells = blk.cells.length;
    size_t nvertices = blk.vertices.length;
    // write out adjoint variables in VTK-format
    if (blk.grid_type == Grid_t.structured_grid) {
        auto sblk = cast(SFluidBlock) blk; 
        auto outFile = File("adjointVars.vtk", "w");
        outFile.writef("# vtk DataFile Version 3.0 \n");
        outFile.writef("%s \n", jobName);
        outFile.writef("ASCII \n");
        outFile.writef("DATASET STRUCTURED_GRID \n");
        outFile.writef("DIMENSIONS %d %d %d \n", sblk.imax, sblk.jmax, sblk.kmax+1);
        outFile.writef("POINTS %d double \n", nvertices);
        // write grid data
        foreach(i, vtx; blk.vertices) {
            outFile.writef("%.16f %.16f %.16f \n", vtx.pos[0].x, vtx.pos[0].y, vtx.pos[0].z); 
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
