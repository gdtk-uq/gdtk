/**
 * viscousflux.d
 * Viscous-Flux calculation, where the fluxes are driven by molecular-transport effects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-05-02: port essentials from Eilmer3 and refactor (a lot).
 *          2015-05-03: added gradient estimation for 2D flow
 */

module viscousflux;

import std.math;
import std.stdio;
import std.conv;
import geom;
import gas;
import flowstate;
import conservedquantities;
import fvcore;
import fvinterface;
import fvvertex;
import globalconfig;


// Stand-alone functions for the computation of gradients.

// [TODO] pass in cloud_pos and cloud_fs arrays directly 
// so that we're not limiting these functions to the vertices.

@nogc
void gradients_xy_div(ref FVVertex vtx, bool diffusion)
// Using the divergence theorem (I think), compute the average gradients
// for the flow conditions over a polygon in the xy-plane.
//     2-----1   1
//     |     |   |\
//     |     |   | \
//     3-----0   2--0
//  y
//  |
//  o--x
//
// Since we are embedding this function in a nominally 3D code,
// the z-coordinate derivatives are set to zero.
{
    // Number of corners in our polygon.
    size_t n = vtx.cloud_pos.length;
    // Compute our own estimate of *twice* the area in xy plane here.
    // Start with the contribution from the final segment of the bounding contour.
    double areaxy = (vtx.cloud_pos[0].x + vtx.cloud_pos[n-1].x) *
	(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    // Accumulate the contributions from the other segments.
    foreach (i; 0 .. n-1) {
	areaxy += (vtx.cloud_pos[i+1].x + vtx.cloud_pos[i].x) *
	    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
    }
    double area_inv = 1.0 / areaxy;
    //
    // Apply the divergence theorem to flow properties.
    //
    // Start with the contribution from the final segment of the bounding contour.
    double gradient_x = (vtx.cloud_fs[0].vel.x + vtx.cloud_fs[n-1].vel.x) *
	(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    double gradient_y = (vtx.cloud_fs[0].vel.x + vtx.cloud_fs[n-1].vel.x) *
	(vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
    // Accumulate the contributions from the other segments.
    foreach (i; 0 .. n-1) {
	gradient_x += (vtx.cloud_fs[i+1].vel.x + vtx.cloud_fs[i].vel.x) *
	    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
	gradient_y += (vtx.cloud_fs[i+1].vel.x + vtx.cloud_fs[i].vel.x) *
	    (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
    }
    // Compute average gradient and pack away.
    vtx.grad_vel[0][0] = gradient_x * area_inv;
    vtx.grad_vel[0][1] = -gradient_y * area_inv;
    vtx.grad_vel[0][2] = 0.0;
    //
    gradient_x = (vtx.cloud_fs[0].vel.y + vtx.cloud_fs[n-1].vel.y) *
	(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    gradient_y = (vtx.cloud_fs[0].vel.y + vtx.cloud_fs[n-1].vel.y) *
	(vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
    foreach (i; 0 .. n-1) {
	gradient_x += (vtx.cloud_fs[i+1].vel.y + vtx.cloud_fs[i].vel.y) *
	    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
	gradient_y += (vtx.cloud_fs[i+1].vel.y + vtx.cloud_fs[i].vel.y) *
	    (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
    }
    vtx.grad_vel[1][0] = gradient_x * area_inv;
    vtx.grad_vel[1][1] = -gradient_y * area_inv;
    vtx.grad_vel[1][2] = 0.0;
    //
    vtx.grad_vel[2][0] = 0.0;
    vtx.grad_vel[2][1] = 0.0;
    vtx.grad_vel[2][2] = 0.0;
    //
    gradient_x = (vtx.cloud_fs[0].gas.T[0] + vtx.cloud_fs[n-1].gas.T[0]) *
	(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    gradient_y = (vtx.cloud_fs[0].gas.T[0] + vtx.cloud_fs[n-1].gas.T[0]) *
	(vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
    foreach (i; 0 .. n-1) {
	gradient_x += (vtx.cloud_fs[i+1].gas.T[0] + vtx.cloud_fs[i].gas.T[0]) *
	    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
	gradient_y += (vtx.cloud_fs[i+1].gas.T[0] + vtx.cloud_fs[i].gas.T[0]) *
	    (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
    }
    vtx.grad_T.refx = gradient_x * area_inv;
    vtx.grad_T.refy = -gradient_y * area_inv;
    vtx.grad_T.refz = 0.0;
    //
    size_t nsp = vtx.cloud_fs[0].gas.massf.length;
    if (diffusion) {
	foreach(isp; 0 .. nsp) {
	    gradient_x = (vtx.cloud_fs[0].gas.massf[isp] + vtx.cloud_fs[n-1].gas.massf[isp]) *
		(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
	    gradient_y = (vtx.cloud_fs[0].gas.massf[isp] + vtx.cloud_fs[n-1].gas.massf[isp]) *
		(vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
	    foreach (i; 0 .. n-1) {
		gradient_x += (vtx.cloud_fs[i+1].gas.massf[isp] + vtx.cloud_fs[i].gas.massf[isp]) *
		    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
		gradient_y += (vtx.cloud_fs[i+1].gas.massf[isp] + vtx.cloud_fs[i].gas.massf[isp]) *
		    (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
	    }
	    vtx.grad_f[isp].refx = gradient_x * area_inv;
	    vtx.grad_f[isp].refy = -gradient_y * area_inv;
	    vtx.grad_f[isp].refz = 0.0;
	}
    } else {
	foreach(isp; 0 .. nsp) {
	    vtx.grad_f[isp].refx = 0.0;
	    vtx.grad_f[isp].refy = 0.0;
	    vtx.grad_f[isp].refz = 0.0;
	}
    }
    //
    gradient_x = (vtx.cloud_fs[0].tke + vtx.cloud_fs[n-1].tke) *
	(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    gradient_y = (vtx.cloud_fs[0].tke + vtx.cloud_fs[n-1].tke) *
	(vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
    foreach (i; 0 .. n-1) {
	gradient_x += (vtx.cloud_fs[i+1].tke + vtx.cloud_fs[i].tke) *
	    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
	gradient_y += (vtx.cloud_fs[i+1].tke + vtx.cloud_fs[i].tke) *
	    (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
    }
    vtx.grad_tke.refx = gradient_x * area_inv;
    vtx.grad_tke.refy = -gradient_y * area_inv;
    vtx.grad_tke.refz = 0.0;
    //
    gradient_x = (vtx.cloud_fs[0].omega + vtx.cloud_fs[n-1].omega) *
	(vtx.cloud_pos[0].y - vtx.cloud_pos[n-1].y);
    gradient_y = (vtx.cloud_fs[0].omega + vtx.cloud_fs[n-1].omega) *
	(vtx.cloud_pos[0].x - vtx.cloud_pos[n-1].x);
    foreach (i; 0 .. n-1) {
	gradient_x += (vtx.cloud_fs[i+1].omega + vtx.cloud_fs[i].omega) *
	    (vtx.cloud_pos[i+1].y - vtx.cloud_pos[i].y);
	gradient_y += (vtx.cloud_fs[i+1].omega + vtx.cloud_fs[i].omega) *
	    (vtx.cloud_pos[i+1].x - vtx.cloud_pos[i].x);
    }
    vtx.grad_omega.refx = gradient_x * area_inv;
    vtx.grad_omega.refy = -gradient_y * area_inv;
    vtx.grad_omega.refz = 0.0;
} // end gradients_xy_div()

@nogc
void computeInverse(int N)(ref double[2*N][N] c, double very_small_value=1.0e-16)
// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
{
    foreach(j; 0 .. N) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. N) {
	    if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
	}
	assert(abs(c[p][j]) > very_small_value, "matrix is essentially singular");
	if ( p != j ) { // Swap rows
	    foreach(col; 0 .. 2*N) {
		double tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
	    }
	}
	// Scale row j to get unity on the diagonal.
	double cjj = c[j][j];
	foreach(col; 0 .. 2*N) c[j][col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. N) {
	    if ( i == j ) continue;
	    double cij = c[i][j];
	    foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col]; 
	}
    } // end foreach j
} // end gaussJordanElimination()()

@nogc
void solveGradients(int N)(ref double[2*N][N] c, ref double[N] rhs, ref double[N] qgrad)
// Multiply right-hand-side by the inverse part of the augmented matrix.
{
    foreach(i; 0 .. N) {
	qgrad[i] = 0.0;
	foreach(j; 0 .. N) {
	    qgrad[i] += c[i][N+j] * rhs[j];
	}
    }
} // end solveGradients()()

@nogc
void gradients_xyz_leastsq(ref FVVertex vtx, bool diffusion)
// Fit a linear model to the cloud of flow-quantity points
// in order to extract approximations to the flow-field gradients.
// 3D
{
    double[6][3] xTx; // normal matrix, augmented to give 6 entries per row
    double[3] rhs, gradients;
    // Assemble and invert the normal matrix.
    // We'll reuse the resulting inverse for each flow-field quantity.
    size_t n = vtx.cloud_pos.length;
    double xx = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yy = 0.0;
    double yz = 0.0;
    double zz = 0.0;
    foreach (i; 1 .. n) {
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	xx += dx*dx; xy += dx*dy; xz += dx*dz;
	yy += dy*dy; yz += dy*dz; zz += dz*dz;
    }
    xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
    xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
    xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
    xTx[0][3] = 1.0; xTx[0][4] = 0.0; xTx[0][5] = 0.0;
    xTx[1][3] = 0.0; xTx[1][4] = 1.0; xTx[1][5] = 0.0;
    xTx[2][3] = 0.0; xTx[2][4] = 0.0; xTx[2][5] = 1.0;
    computeInverse!3(xTx);
    // x-velocity
    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
    foreach (i; 1 .. n) {
	double dvx = vtx.cloud_fs[i].vel.x - vtx.cloud_fs[0].vel.x;
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	rhs[0] += dx*dvx; rhs[1] += dy*dvx; rhs[2] += dz*dvx;
    }
    solveGradients!3(xTx, rhs, gradients);
    foreach (j; 0 .. 3) { vtx.grad_vel[0][j] = gradients[j]; }
    // y-velocity
    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
    foreach (i; 1 .. n) {
	double dvy = vtx.cloud_fs[i].vel.y - vtx.cloud_fs[0].vel.y;
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	rhs[0] += dx*dvy; rhs[1] += dy*dvy; rhs[2] += dz*dvy;
    }
    solveGradients!3(xTx, rhs, gradients);
    foreach (j; 0 .. 3) { vtx.grad_vel[1][j] = gradients[j]; }
    // z-velocity
    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
    foreach (i; 1 .. n) {
	double dvz = vtx.cloud_fs[i].vel.z - vtx.cloud_fs[0].vel.z;
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	rhs[0] += dx*dvz; rhs[1] += dy*dvz; rhs[2] += dz*dvz;
    }
    solveGradients!3(xTx, rhs, gradients);
    foreach (j; 0 .. 3) { vtx.grad_vel[2][j] = gradients[j]; }
    // T[0]
    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
    foreach (i; 1 .. n) {
	double dT = vtx.cloud_fs[i].gas.T[0] - vtx.cloud_fs[0].gas.T[0];
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	rhs[0] += dx*dT; rhs[1] += dy*dT; rhs[2] += dz*dT;
    }
    solveGradients!3(xTx, rhs, gradients);
    vtx.grad_T.refx = gradients[0];
    vtx.grad_T.refy = gradients[1];
    vtx.grad_T.refz = gradients[2];
    // massf
    size_t nsp = vtx.cloud_fs[0].gas.massf.length;
    if (diffusion) {
	foreach(isp; 0 .. nsp) {
	    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
	    foreach (i; 1 .. n) {
		double df = vtx.cloud_fs[i].gas.massf[isp] - vtx.cloud_fs[0].gas.massf[isp];
		double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
		double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
		double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
		rhs[0] += dx*df; rhs[1] += dy*df; rhs[2] += dz*df;
	    }
	    solveGradients!3(xTx, rhs, gradients);
	    vtx.grad_f[isp].refx = gradients[0];
	    vtx.grad_f[isp].refy = gradients[1];
	    vtx.grad_f[isp].refz = gradients[2];
	} // foreach isp
    } else {
	foreach(isp; 0 .. nsp) {
	    vtx.grad_f[isp].refx = 0.0;
	    vtx.grad_f[isp].refy = 0.0;
	    vtx.grad_f[isp].refz = 0.0;
	} // foreach isp
    }
    // tke
    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
    foreach (i; 1 .. n) {
	double dtke = vtx.cloud_fs[i].tke - vtx.cloud_fs[0].tke;
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	rhs[0] += dx*dtke; rhs[1] += dy*dtke; rhs[2] += dz*dtke;
    }
    solveGradients!3(xTx, rhs, gradients);
    vtx.grad_tke.refx = gradients[0];
    vtx.grad_tke.refy = gradients[1];
    vtx.grad_tke.refz = gradients[2];
    // omega
    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
    foreach (i; 1 .. n) {
	double domega = vtx.cloud_fs[i].omega - vtx.cloud_fs[0].omega;
	double dx = vtx.cloud_pos[i].x - vtx.cloud_pos[0].x;
	double dy = vtx.cloud_pos[i].y - vtx.cloud_pos[0].y;
	double dz = vtx.cloud_pos[i].z - vtx.cloud_pos[0].z;
	rhs[0] += dx*domega; rhs[1] += dy*domega; rhs[2] += dz*domega;
    }
    solveGradients!3(xTx, rhs, gradients);
    vtx.grad_omega.refx = gradients[0];
    vtx.grad_omega.refy = gradients[1];
    vtx.grad_omega.refz = gradients[2];
} // end gradients_xyz_leastsq()

@nogc
void gradients_xy_leastsq(ref FVVertex vtx, bool diffusion)
// Fit a linear model to the cloud of flow-quantity points
// in order to extract approximations to the flow-field gradients.
// 2D, built as a specialization of the 3D code.
// Experiment with taking differences about a middle point/value.
{
    double[4][2] xTx; // normal matrix, augmented to give 4 entries per row
    double[2] rhs, gradients;
    // Assemble and invert the normal matrix.
    // We'll reuse the resulting inverse for each flow-field quantity.
    size_t n = vtx.cloud_pos.length;
    double x_mid = 0.0;
    double y_mid = 0.0;
    foreach (i; 0 .. n) {
	x_mid += vtx.cloud_pos[i].x;
	y_mid += vtx.cloud_pos[i].y;
    }
    x_mid /= n;
    y_mid /= n;
    double xx = 0.0;
    double xy = 0.0;
    double yy = 0.0;
    foreach (i; 0 .. n) {
	double dx = vtx.cloud_pos[i].x - x_mid;
	double dy = vtx.cloud_pos[i].y - y_mid;
	xx += dx*dx; xy += dx*dy; yy += dy*dy;
    }
    xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = 1.0; xTx[0][3] = 0.0;
    xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = 0.0; xTx[1][3] = 1.0;
    computeInverse!2(xTx);
    // x-velocity
    double vx_mid = 0.0;
    foreach (i; 0 .. n) {
	vx_mid += vtx.cloud_fs[i].vel.x;
    }
    vx_mid /= n;
    rhs[0] = 0.0; rhs[1] = 0.0;
    foreach (i; 0 .. n) {
	double dvx = vtx.cloud_fs[i].vel.x - vx_mid;
	double dx = vtx.cloud_pos[i].x - x_mid;
	double dy = vtx.cloud_pos[i].y - y_mid;
	rhs[0] += dx*dvx; rhs[1] += dy*dvx;
    }
    solveGradients!2(xTx, rhs, gradients);
    vtx.grad_vel[0][0] = gradients[0];
    vtx.grad_vel[0][1] = gradients[1];
    vtx.grad_vel[0][2] = 0.0;
    // y-velocity
    double vy_mid = 0.0;
    foreach (i; 0 .. n) {
	vy_mid += vtx.cloud_fs[i].vel.y;
    }
    vy_mid /= n;
    rhs[0] = 0.0; rhs[1] = 0.0;
    foreach (i; 0 .. n) {
	double dvy = vtx.cloud_fs[i].vel.y - vy_mid;
	double dx = vtx.cloud_pos[i].x - x_mid;
	double dy = vtx.cloud_pos[i].y - y_mid;
	rhs[0] += dx*dvy; rhs[1] += dy*dvy;
    }
    solveGradients!2(xTx, rhs, gradients);
    vtx.grad_vel[1][0] = gradients[0];
    vtx.grad_vel[1][1] = gradients[1];
    vtx.grad_vel[1][2] = 0.0;
    // z-velocity
    vtx.grad_vel[2][0] = 0.0;
    vtx.grad_vel[2][1] = 0.0;
    vtx.grad_vel[2][2] = 0.0;
    // T[0]
    double T_mid = 0.0;
    foreach (i; 0 .. n) {
	T_mid += vtx.cloud_fs[i].gas.T[0];
    }
    T_mid /= n;
    rhs[0] = 0.0; rhs[1] = 0.0;
    foreach (i; 0 .. n) {
	double dT = vtx.cloud_fs[i].gas.T[0] - T_mid;
	double dx = vtx.cloud_pos[i].x - x_mid;
	double dy = vtx.cloud_pos[i].y - y_mid;
	rhs[0] += dx*dT; rhs[1] += dy*dT;
    }
    solveGradients!2(xTx, rhs, gradients);
    vtx.grad_T.refx = gradients[0];
    vtx.grad_T.refy = gradients[1];
    vtx.grad_T.refz = 0.0;
    // massf
    size_t nsp = vtx.cloud_fs[0].gas.massf.length;
    if (diffusion) {
	foreach(isp; 0 .. nsp) {
	    double massf_mid = 0.0;
	    foreach (i; 0 .. n) {
		massf_mid += vtx.cloud_fs[i].gas.massf[isp];
	    }
	    massf_mid /= n;
	    rhs[0] = 0.0; rhs[1] = 0.0;
	    foreach (i; 0 .. n) {
		double df = vtx.cloud_fs[i].gas.massf[isp] - massf_mid;
		double dx = vtx.cloud_pos[i].x - x_mid;
		double dy = vtx.cloud_pos[i].y - y_mid;
		rhs[0] += dx*df; rhs[1] += dy*df;
	    }
	    solveGradients!2(xTx, rhs, gradients);
	    vtx.grad_f[isp].refx = gradients[0];
	    vtx.grad_f[isp].refy = gradients[1];
	    vtx.grad_f[isp].refz = 0.0;
	} // foreach isp
    } else {
	foreach(isp; 0 .. nsp) {
	    vtx.grad_f[isp].refx = 0.0;
	    vtx.grad_f[isp].refy = 0.0;
	    vtx.grad_f[isp].refz = 0.0;
	} // foreach isp
    }
    // tke
    double tke_mid = 0.0;
    foreach (i; 0 .. n) {
	tke_mid += vtx.cloud_fs[i].tke;
    }
    tke_mid /= n;
    rhs[0] = 0.0; rhs[1] = 0.0;
    foreach (i; 0 .. n) {
	double dtke = vtx.cloud_fs[i].tke - tke_mid;
	double dx = vtx.cloud_pos[i].x - x_mid;
	double dy = vtx.cloud_pos[i].y - y_mid;
	rhs[0] += dx*dtke; rhs[1] += dy*dtke;
    }
    solveGradients!2(xTx, rhs, gradients);
    vtx.grad_tke.refx = gradients[0];
    vtx.grad_tke.refy = gradients[1];
    vtx.grad_tke.refz = 0.0;
    // omega
    double omega_mid = 0.0;
    foreach (i; 0 .. n) {
	omega_mid += vtx.cloud_fs[i].omega;
    }
    omega_mid /= n;
    rhs[0] = 0.0; rhs[1] = 0.0;
    foreach (i; 0 .. n) {
	double domega = vtx.cloud_fs[i].omega - omega_mid;
	double dx = vtx.cloud_pos[i].x - x_mid;
	double dy = vtx.cloud_pos[i].y - y_mid;
	rhs[0] += dx*domega; rhs[1] += dy*domega;
    }
    solveGradients!2(xTx, rhs, gradients);
    vtx.grad_omega.refx = gradients[0];
    vtx.grad_omega.refy = gradients[1];
    vtx.grad_omega.refz = 0.0;
} // end gradients_xy_leastsq()
