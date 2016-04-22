/** simcore.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Rowan G. and Peter J.
 * First code: 2016-03-29
 */

module steadystate;

import std.math;
import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.datetime;
import std.parallelism;

import nm.smla;
import nm.bbla;

import globalconfig;
import globaldata;
import block;
import sblock;
import simcore;
import user_defined_source_terms;

//----------------------------------------------------------------------------

// We will use the flexible GMRES algorithm to solve for the Newton
// updated. We'll precondition the system with an inner GMRES solve
// on a system with low-order flux evaluation. This is an inner-outer
// method. The function here is the outer part of the method, and
// calls the inner GMRES iterations.

void flexGMRES_solve(double pseudo_sim_time, double dt, double residTol)
{
    Block blk = gasBlocks[0];
    int maxIters = 30;
    size_t iterCount;
    size_t m = to!size_t(maxIters);
    
    size_t n = blk.FU_o.length;
    // 0. Initialise some values in the pre-allocated storage
    blk.g0_o[] = 0.0;
    blk.g1_o[] = 0.0;
    blk.H0_o.zeros();
    blk.H1_o.zeros();
    blk.Gamma_o.eye();

    // 1. Evaluate r0, beta, v1
    // This actually involves doing all the steps necessary to populate dUdt[0].
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    blk.applyPreReconAction(pseudo_sim_time, 0, 0);
    blk.convective_flux();
    blk.applyPostConvFluxAction(pseudo_sim_time, 0, 0);

    double max_Frho = 0.0;
    double max_FrhoU = 0.0;
    double max_FrhoE = 0.0;
    int cellOffset = 0;    
    foreach (cell; blk.cells) {
	cell.add_inviscid_source_vector(0, 0.0);
	if (blk.myConfig.udf_source_terms) {
	    addUDFSourceTermsToCell(blk.myL, cell, 0, 
				    pseudo_sim_time, blk.myConfig.gmodel);
	}
	cell.time_derivatives(0, 0, false);
	blk.FU_o[cellOffset+0] = cell.dUdt[0].mass;
	max_Frho = max(max_Frho, fabs(cell.dUdt[0].mass));
	blk.FU_o[cellOffset+1] = cell.dUdt[0].momentum.x;
	max_FrhoU = max(max_FrhoU, fabs(cell.dUdt[0].momentum.x));
	blk.FU_o[cellOffset+2] = cell.dUdt[0].momentum.y;
	max_FrhoU = max(max_FrhoU, fabs(cell.dUdt[0].momentum.x));
	blk.FU_o[cellOffset+3] = cell.dUdt[0].total_energy;
	max_FrhoE = max(max_FrhoE, fabs(cell.dUdt[0].total_energy));
	cellOffset += 4;
    }
    // Probably need some guards against max values being 0.0
    //writefln("max_Frho= %20.16e", max_Frho);
    //writefln("max_FrhoU= %20.16e", max_FrhoU);
    //writefln("max_FrhoE= %20.16e", max_FrhoE);
    cellOffset = 0;
    foreach (cell; blk.cells) {
	blk.r0_o[cellOffset+0] = (1/max_Frho)*blk.FU_o[cellOffset+0];
	blk.r0_o[cellOffset+1] = (1/max_FrhoU)*blk.FU_o[cellOffset+1];
	blk.r0_o[cellOffset+2] = (1/max_FrhoU)*blk.FU_o[cellOffset+2];
	blk.r0_o[cellOffset+3] = (1/max_FrhoE)*blk.FU_o[cellOffset+3];
	cellOffset += 4;
    }

    auto beta = norm2(blk.r0_o);
    //writeln(format("beta= %20.12e", beta));
    blk.g0_o[0] = beta;
    foreach (k; 0 .. n) {
	blk.v_o[k] = blk.r0_o[k]/beta;
	blk.V_o[k,0] = blk.v_o[k];
    }

	    // 2. Begin iterations
    foreach (j; 0 .. maxIters) {
	writefln("==== ITERATION %d ====", j);
	iterCount = j+1;
	if ( blk.myConfig.interpolation_order > 1 ) {
	    GMRES_solve(pseudo_sim_time, dt, 10*residTol, blk.v_o, blk.z_o, max_Frho, max_FrhoU, max_FrhoE);
	}
	else {
	    // Skip the inner iterations.
	    blk.z_o[] = blk.v_o[];
	}
	// Save z vector in Z.
	foreach (k; 0 .. n) blk.Z_o[k,j] = blk.z_o[k];
	// Get w ready with (I/dt)p term.
	double dtInv = 1.0/dt;
	foreach (k; 0 .. n) {
	    blk.w_o[k] = dtInv*blk.z_o[k]; 
	}
	// Next do finite-difference perturbation to compute Jz
	blk.clear_fluxes_of_conserved_quantities();
	foreach (cell; blk.cells) cell.clear_source_vector();
	// Compute sigma, the perturbation constant using approach of
	// Cai et al (1995)
	double sigma = sqrt(n*double.epsilon);
	cellOffset = 0;
	foreach (cell; blk.cells) {
	    cell.U[1].copy_values_from(cell.U[0]);
	    cell.U[1].mass += sigma*max_Frho*blk.z_o[cellOffset+0];
	    cell.U[1].momentum.refx += sigma*max_FrhoU*blk.z_o[cellOffset+1];
	    cell.U[1].momentum.refy += sigma*max_FrhoU*blk.z_o[cellOffset+2];
	    cell.U[1].total_energy += sigma*max_FrhoE*blk.z_o[cellOffset+3];
	    cell.decode_conserved(0, 1, 0.0);
	    cellOffset += 4;
	}
	blk.applyPreReconAction(pseudo_sim_time, 0, 1);
	blk.convective_flux();
	blk.applyPostConvFluxAction(pseudo_sim_time, 0, 1);
	cellOffset = 0;
	foreach (cell; blk.cells) {
	    cell.add_inviscid_source_vector(1, 0.0);
	    if (blk.myConfig.udf_source_terms) {
		addUDFSourceTermsToCell(blk.myL, cell, 0, 
					pseudo_sim_time, blk.myConfig.gmodel);
	    }
	    cell.time_derivatives(0, 1, false);
	    blk.z_o[cellOffset+0] = (1.0/max_Frho)*((cell.dUdt[1].mass - blk.FU_o[cellOffset+0])/sigma);
	    blk.z_o[cellOffset+1] = (1.0/max_FrhoU)*((cell.dUdt[1].momentum.x - blk.FU_o[cellOffset+1])/sigma);
	    blk.z_o[cellOffset+2] = (1.0/max_FrhoU)*((cell.dUdt[1].momentum.y - blk.FU_o[cellOffset+2])/sigma);
	    blk.z_o[cellOffset+3] = (1.0/max_FrhoE)*((cell.dUdt[1].total_energy - blk.FU_o[cellOffset+3])/sigma);
	    cellOffset += 4;
	}
	// Finally update w.
	foreach (i; 0 .. n) {
	    blk.w_o[i] += blk.z_o[i];
	}
	// The remainder of the algorithm looks a lot like the
	// GMRES method in smla.d
	foreach (i; 0 .. j+1) {
	    foreach (k; 0 .. n ) blk.v_o[k] = blk.V_o[k,i]; // Extract column 'i'
	    blk.H0_o[i,j] = dot(blk.w_o, blk.v_o);
	    foreach (k; 0 .. n) blk.w_o[k] -= blk.H0_o[i,j]*blk.v_o[k]; 
	}
	blk.H0_o[j+1,j] = norm2(blk.w_o);
	
	foreach (k; 0 .. n) {
	    blk.v_o[k] = blk.w_o[k]/blk.H0_o[j+1,j];
	    blk.V_o[k,j+1] = blk.v_o[k];
	}

	// Build rotated Hessenberg progressively
	if ( j != 0 ) {
	    // Extract final column in H
	    foreach (i; 0 .. j+1) blk.h_o[i] = blk.H0_o[i,j];
	    // Rotate column by previous rotations (stored in Q0)
	    nm.bbla.dot(blk.Q0_o, j+1, j+1, blk.h_o, blk.hR_o);
	    // Place column back in H
	    foreach (i; 0 .. j+1) blk.H0_o[i,j] = blk.hR_o[i];
	}
	// Now form new Gamma
	blk.Gamma_o.eye();
	auto denom = sqrt(blk.H0_o[j,j]*blk.H0_o[j,j] + blk.H0_o[j+1,j]*blk.H0_o[j+1,j]);
	auto s_j = blk.H0_o[j+1,j]/denom; 
	auto c_j = blk.H0_o[j,j]/denom;
	blk.Gamma_o[j,j] = c_j; blk.Gamma_o[j,j+1] = s_j;
	blk.Gamma_o[j+1,j] = -s_j; blk.Gamma_o[j+1,j+1] = c_j;
	// Apply rotations
	nm.bbla.dot(blk.Gamma_o, j+2, j+2, blk.H0_o, j+1, blk.H1_o);
	nm.bbla.dot(blk.Gamma_o, j+2, j+2, blk.g0_o, blk.g1_o);
	// Accumulate Gamma rotations in Q.
	if ( j == 0 ) {
	    copy(blk.Gamma_o, blk.Q1_o);
	}
	else {
	    nm.bbla.dot(blk.Gamma_o, j+2, j+2, blk.Q0_o, j+2, blk.Q1_o);
	}
	// Get residual
	double resid = fabs(blk.g1_o[j+1]);
	writefln("residual= %20.12e", resid);
	if ( resid <= residTol ) {
	    m = j+1;
	    break;
	}
	
	// Prepare for next step
	copy(blk.H1_o, blk.H0_o);
	blk.g0_o[] = blk.g1_o[];
	copy(blk.Q1_o, blk.Q0_o);
    }
	
    if ( iterCount == maxIters )
	m = maxIters;
    
    // At end H := R up to row m
    //        g := gm up to row m
    upperSolve(blk.H1_o, to!int(m), blk.g1_o);
    nm.bbla.dot(blk.Z_o, n, m, blk.g1_o, blk.dU);
    // Rescale.
    cellOffset = 0;
    foreach (cell; blk.cells) {
	blk.dU[cellOffset+0] *= max_Frho;
	blk.dU[cellOffset+1] *= max_FrhoU;
	blk.dU[cellOffset+2] *= max_FrhoU;
	blk.dU[cellOffset+3] *= max_FrhoE;
	cellOffset += 4;
    }
    
} // end flexGMRES_solve

    // GMRES solve is the inner iterations, working with low-order Jacobian

void GMRES_solve(double pseudo_sim_time, double dt, double residTol, double[] v, double[] z,
		 double max_Frho, double max_FrhoU, double max_FrhoE)
{
    SBlock blk = cast(SBlock) gasBlocks[0];
    int interpOrderSave = blk.get_interpolation_order();
    int maxIters = 2;
    size_t iterCount;
    size_t m = to!size_t(maxIters);

    size_t n = v.length;
    // 0. Initialise some values in the pre-allocated storage
    blk.g0_i[] = 0.0;
    blk.g1_i[] = 0.0;
    blk.H0_i.zeros();
    blk.H1_i.zeros();
    blk.Gamma_i.eye();

    // 1. Evaluate r0, beta, v1
    // This actually involves doing all the steps necessary to populate dUdt[0]
    // with a low-order reconstruction.
    blk.set_interpolation_order(1);
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    blk.applyPreReconAction(pseudo_sim_time, 0, 0);
    blk.convective_flux();
    blk.applyPostConvFluxAction(pseudo_sim_time, 0, 0);

    int cellOffset = 0;    
    foreach (cell; blk.cells) {
	cell.add_inviscid_source_vector(0, 0.0);
	cell.time_derivatives(0, 0, false);
	blk.FU_i[cellOffset+0] = cell.dUdt[0].mass;
	blk.FU_i[cellOffset+1] = cell.dUdt[0].momentum.x;
	blk.FU_i[cellOffset+2] = cell.dUdt[0].momentum.y;
	blk.FU_i[cellOffset+3] = cell.dUdt[0].total_energy;
	cellOffset += 4;
    }

    cellOffset = 0;
    foreach (cell; blk.cells) {
	blk.r0_i[cellOffset+0] = (1/max_Frho)*blk.FU_i[cellOffset+0];
	blk.r0_i[cellOffset+1] = (1/max_FrhoU)*blk.FU_i[cellOffset+1];
	blk.r0_i[cellOffset+2] = (1/max_FrhoU)*blk.FU_i[cellOffset+2];
	blk.r0_i[cellOffset+3] = (1/max_FrhoE)*blk.FU_i[cellOffset+3];
	cellOffset += 4;
    }

    auto beta = norm2(blk.r0_i);
    blk.g0_i[0] = beta;
    foreach (k; 0 .. n) {
	blk.v_i[k] = blk.r0_i[k]/beta;
	blk.V_i[k,0] = blk.v_i[k];
    }

    // 2. Begin iterations
    foreach (j; 0 .. maxIters) {
	iterCount = j+1;
	// Get w ready with (I/dt)v term.
	double dtInv = 1.0/dt;
	foreach (i; 0 .. n) {
	    blk.w_i[i] = dtInv*blk.v_i[i]; 
	}
	// Next do finite-difference perturbation to compute JP^{-1}v
	blk.clear_fluxes_of_conserved_quantities();
	foreach (cell; blk.cells) cell.clear_source_vector();
	// Compute sigma, the perturbation constant using approach of
	// Cai et al (1995)
	double sigma = sqrt(n*double.epsilon);
	cellOffset = 0;
	foreach (cell; blk.cells) {
	    cell.U[1].copy_values_from(cell.U[0]);
	    cell.U[1].mass += sigma*max_Frho*blk.v_i[cellOffset+0];
	    cell.U[1].momentum.refx += sigma*max_FrhoU*blk.v_i[cellOffset+1];
	    cell.U[1].momentum.refy += sigma*max_FrhoU*blk.v_i[cellOffset+2];
	    cell.U[1].total_energy += sigma*max_FrhoE*blk.v_i[cellOffset+3];
	    cell.decode_conserved(0, 1, 0.0);
	    cellOffset += 4;
	}
	blk.applyPreReconAction(pseudo_sim_time, 0, 1);
	blk.convective_flux();
	blk.applyPostConvFluxAction(pseudo_sim_time, 0, 1);
	cellOffset = 0;
	foreach (cell; blk.cells) {
	    cell.add_inviscid_source_vector(1, 0.0);
	    cell.time_derivatives(0, 1, false);
	    blk.v_i[cellOffset+0] = (1.0/max_Frho)*((cell.dUdt[1].mass - blk.FU_i[cellOffset+0])/sigma);
	    blk.v_i[cellOffset+1] = (1.0/max_FrhoU)*((cell.dUdt[1].momentum.x - blk.FU_i[cellOffset+1])/sigma);
	    blk.v_i[cellOffset+2] = (1.0/max_FrhoU)*((cell.dUdt[1].momentum.y - blk.FU_i[cellOffset+2])/sigma);
	    blk.v_i[cellOffset+3] = (1.0/max_FrhoE)*((cell.dUdt[1].total_energy - blk.FU_i[cellOffset+3])/sigma);
	    cellOffset += 4;
	}

	// Finally update w.
	foreach (k; 0 .. n) {
	    blk.w_i[k] += blk.v_i[k];
	}

	// The remainder of the algorithm looks a lot like the
	// GMRES method in smla.d
	foreach (i; 0 .. j+1) {
	    foreach (k; 0 .. n ) blk.v_i[k] = blk.V_i[k,i]; // Extract column 'i'
	    blk.H0_i[i,j] = dot(blk.w_i, blk.v_i);
	    foreach (k; 0 .. n) blk.w_i[k] -= blk.H0_i[i,j]*blk.v_i[k]; 
	}
	blk.H0_i[j+1,j] = norm2(blk.w_i);
	
	foreach (k; 0 .. n) {
	    blk.v_i[k] = blk.w_i[k]/blk.H0_i[j+1,j];
	    blk.V_i[k,j+1] = blk.v_i[k];
	}

	// Build rotated Hessenberg progressively
	if ( j != 0 ) {
	    // Extract final column in H
	    foreach (i; 0 .. j+1) blk.h_i[i] = blk.H0_i[i,j];
	    // Rotate column by previous rotations (stored in Q0)
	    nm.bbla.dot(blk.Q0_i, j+1, j+1, blk.h_i, blk.hR_i);
	    // Place column back in H
	    foreach (i; 0 .. j+1) blk.H0_i[i,j] = blk.hR_i[i];
	}
	// Now form new Gamma
	blk.Gamma_i.eye();
	auto denom = sqrt(blk.H0_i[j,j]*blk.H0_i[j,j] + blk.H0_i[j+1,j]*blk.H0_i[j+1,j]);
	auto s_j = blk.H0_i[j+1,j]/denom; 
	auto c_j = blk.H0_i[j,j]/denom;
	blk.Gamma_i[j,j] = c_j; blk.Gamma_i[j,j+1] = s_j;
	blk.Gamma_i[j+1,j] = -s_j; blk.Gamma_i[j+1,j+1] = c_j;
	// Apply rotations
	nm.bbla.dot(blk.Gamma_i, j+2, j+2, blk.H0_i, j+1, blk.H1_i);
	nm.bbla.dot(blk.Gamma_i, j+2, j+2, blk.g0_i, blk.g1_i);
	// Accumulate Gamma rotations in Q.
	if ( j == 0 ) {
	    copy(blk.Gamma_i, blk.Q1_i);
	}
	else {
	    nm.bbla.dot(blk.Gamma_i, j+2, j+2, blk.Q0_i, j+2, blk.Q1_i);
	}
	// Get residual
	double resid = fabs(blk.g1_i[j+1]);
	if ( resid <= residTol ) {
	    writefln("Breaking.");
	    m = j+1;
	    break;
	}

	// Prepare for next step
	copy(blk.H1_i, blk.H0_i);
	blk.g0_i[] = blk.g1_i[];
	copy(blk.Q1_i, blk.Q0_i);
    }

    if ( iterCount == maxIters )
	m = maxIters;

    // At end H := R up to row m
    //        g := gm up to row m
    upperSolve(blk.H1_i, to!int(m), blk.g1_i);
    nm.bbla.dot(blk.V_i, n, m, blk.g1_i, z);

    // Restore the original interpolation order
    blk.set_interpolation_order(interpOrderSave);
} // end GMRES_solve

