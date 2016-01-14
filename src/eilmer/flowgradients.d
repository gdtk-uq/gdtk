/**
 * flowgradients.d
 * Flow-gradient calculation, for use in the viscous fluxes, 
 * that are driven by molecular-transport effects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-05-02: port essentials from Eilmer3 and refactor (a lot).
 *          2015-05-03: added gradient estimation for 2D flow
 *          2015-12-27: viscous flux calculation sent to the FVInterface class
 */

module flowgradients;

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

class FlowGradients {
    // Spatial derivatives of the flow quantities
public:
    double[][] vel; 
    // velocity derivatives stored as a second-order tensor
    // [[du/dx du/dy du/dz]
    //  [dv/dx dv/dy dv/dz]
    //  [dw/dx dw/dy dw/dz]]
    Vector3[] massf; // mass fraction derivatives
    Vector3 T; // Temperature derivatives (static temperature only)
    Vector3 tke; // turbulence kinetic energy
    Vector3 omega; // pseudo vorticity for k-omega turbulence

    this(size_t nspecies)
    {
	vel.length = 3;
	foreach (ref elem; vel) elem.length = 3;
	massf.length = nspecies; 
    }

    this(const FlowGradients other)
    {
	vel.length = 3;
	foreach(i; 0 .. 3) vel[i] = other.vel[i].dup();
	foreach(item; other.massf) massf ~= Vector3(item);
	T = other.T;
	tke = other.tke;
	omega = other.omega;
    }

    @nogc
    void copy_values_from(const FlowGradients other)
    {
	foreach (i; 0 .. 3) {
	    foreach (j; 0 .. 3) vel[i][j] = other.vel[i][j];
	}
	foreach (isp; 0 .. massf.length) {
	    massf[isp].refx = other.massf[isp].x;
	    massf[isp].refy = other.massf[isp].y;
	    massf[isp].refz = other.massf[isp].z;
	}
	T.refx = other.T.x;
	T.refy = other.T.y;
	T.refz = other.T.z;
	tke.refx = other.tke.x;
	tke.refy = other.tke.y;
	tke.refz = other.tke.z;
	omega.refx = other.omega.x;
	omega.refy = other.omega.y;
	omega.refz = other.omega.z;
    }

    @nogc
    void accumulate_values_from(const FlowGradients other)
    {
	foreach (i; 0 .. 3) {
	    foreach (j; 0 .. 3) vel[i][j] += other.vel[i][j];
	}
	foreach (isp; 0 .. massf.length) {
	    massf[isp].refx += other.massf[isp].x;
	    massf[isp].refy += other.massf[isp].y;
	    massf[isp].refz += other.massf[isp].z;
	}
	T.refx += other.T.x;
	T.refy += other.T.y;
	T.refz += other.T.z;
	tke.refx += other.tke.x;
	tke.refy += other.tke.y;
	tke.refz += other.tke.z;
	omega.refx += other.omega.x;
	omega.refy += other.omega.y;
	omega.refz += other.omega.z;
    }

    @nogc
    void scale_values_by(double factor)
    {
	foreach (i; 0 .. 3) {
	    foreach (j; 0 .. 3) vel[i][j] *= factor;
	} 
	foreach (isp; 0 .. massf.length) massf[isp] *= factor; 
	T *= factor;
	tke *= factor;
	omega *= factor;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FlowGradients(";
	repr ~= ", vel=[";
	foreach (i; 0 .. vel.length) {
	    repr ~= "[" ~ to!string(vel[i][0]);
	    foreach (j; 1 .. vel[i].length) repr ~= ", " ~ to!string(vel[i][j]); 
	    repr ~= "],";
	}
	repr ~= "]";
	repr ~= ", massf=[" ~ to!string(massf[0]);
	foreach (i; 1 .. massf.length) repr ~= ", " ~ to!string(massf);
	repr ~= "]";
	repr ~= ", T=" ~ to!string(T);
	repr ~= ", tke=" ~ to!string(tke);
	repr ~= ", omega=" ~ to!string(omega);
	repr ~= ")";
	return to!string(repr);
    }

    @nogc
    void gradients_xy_div(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos, bool diffusion)
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
	size_t n = cloud_pos.length;
	// Compute our own estimate of *twice* the area in xy plane here.
	// Start with the contribution from the final segment of the bounding contour.
	double areaxy = (cloud_pos[0].x + cloud_pos[n-1].x) *
	    (cloud_pos[0].y - cloud_pos[n-1].y);
	// Accumulate the contributions from the other segments.
	foreach (i; 0 .. n-1) {
	    areaxy += (cloud_pos[i+1].x + cloud_pos[i].x) *
		(cloud_pos[i+1].y - cloud_pos[i].y);
	}
	double area_inv = 1.0 / areaxy;
	//
	// Apply the divergence theorem to flow variables, generating code for each.
	//
	double gradient_x, gradient_y;
	string codeForGradients(string qname)
	{
	    string code = "
            // Start with the contribution from the final segment of the bounding contour.
            gradient_x = (cloud_fs[0]."~qname~" + cloud_fs[n-1]."~qname~") *
                (cloud_pos[0].y - cloud_pos[n-1].y);
	    gradient_y = (cloud_fs[0]."~qname~" + cloud_fs[n-1]."~qname~") *
                (cloud_pos[0].x - cloud_pos[n-1].x);
            // Accumulate the contributions from the other segments.
	    foreach (i; 0 .. n-1) {
                gradient_x += (cloud_fs[i+1]."~qname~" + cloud_fs[i]."~qname~") *
      	            (cloud_pos[i+1].y - cloud_pos[i].y);
                gradient_y += (cloud_fs[i+1]."~qname~" + cloud_fs[i]."~qname~") *
                    (cloud_pos[i+1].x - cloud_pos[i].x);
	    }";
	    return code;
	}
	mixin(codeForGradients("vel.x"));
	vel[0][0] = gradient_x * area_inv;
	vel[0][1] = -gradient_y * area_inv;
	vel[0][2] = 0.0;
	//
	mixin(codeForGradients("vel.y"));
	vel[1][0] = gradient_x * area_inv;
	vel[1][1] = -gradient_y * area_inv;
	vel[1][2] = 0.0;
	//
	vel[2][0] = 0.0;
	vel[2][1] = 0.0;
	vel[2][2] = 0.0;
	//
	mixin(codeForGradients("gas.T[0]"));
	T.refx = gradient_x * area_inv;
	T.refy = -gradient_y * area_inv;
	T.refz = 0.0;
	//
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]"));
		massf[isp].refx = gradient_x * area_inv;
		massf[isp].refy = -gradient_y * area_inv;
		massf[isp].refz = 0.0;
	    }
	} else {
	    foreach(isp; 0 .. nsp) {
		massf[isp].refx = 0.0;
		massf[isp].refy = 0.0;
		massf[isp].refz = 0.0;
	    }
	}
	//
	mixin(codeForGradients("tke"));
	tke.refx = gradient_x * area_inv;
	tke.refy = -gradient_y * area_inv;
	tke.refz = 0.0;
	//
	mixin(codeForGradients("omega"));
	omega.refx = gradient_x * area_inv;
	omega.refy = -gradient_y * area_inv;
	omega.refz = 0.0;
    } // end gradients_xy_div()

    @nogc
    void gradients_xyz_leastsq(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos, bool diffusion)
    // Fit a linear model to the cloud of flow-quantity points
    // in order to extract approximations to the flow-field gradients.
    // 3D
    {
	double[6][3] xTx; // normal matrix, augmented to give 6 entries per row
	double[3] rhs, gradients;
	// Assemble and invert the normal matrix.
	// We'll reuse the resulting inverse for each flow-field quantity.
	size_t n = cloud_pos.length;
	double xx = 0.0;
	double xy = 0.0;
	double xz = 0.0;
	double yy = 0.0;
	double yz = 0.0;
	double zz = 0.0;
	foreach (i; 1 .. n) {
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
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
	    double dvx = cloud_fs[i].vel.x - cloud_fs[0].vel.x;
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
	    rhs[0] += dx*dvx; rhs[1] += dy*dvx; rhs[2] += dz*dvx;
	}
	solveGradients!3(xTx, rhs, gradients);
	foreach (j; 0 .. 3) { vel[0][j] = gradients[j]; }
	// y-velocity
	foreach (j; 0 .. 3) { rhs[j] = 0.0; }
	foreach (i; 1 .. n) {
	    double dvy = cloud_fs[i].vel.y - cloud_fs[0].vel.y;
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
	    rhs[0] += dx*dvy; rhs[1] += dy*dvy; rhs[2] += dz*dvy;
	}
	solveGradients!3(xTx, rhs, gradients);
	foreach (j; 0 .. 3) { vel[1][j] = gradients[j]; }
	// z-velocity
	foreach (j; 0 .. 3) { rhs[j] = 0.0; }
	foreach (i; 1 .. n) {
	    double dvz = cloud_fs[i].vel.z - cloud_fs[0].vel.z;
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
	    rhs[0] += dx*dvz; rhs[1] += dy*dvz; rhs[2] += dz*dvz;
	}
	solveGradients!3(xTx, rhs, gradients);
	foreach (j; 0 .. 3) { vel[2][j] = gradients[j]; }
	// T[0]
	foreach (j; 0 .. 3) { rhs[j] = 0.0; }
	foreach (i; 1 .. n) {
	    double dT = cloud_fs[i].gas.T[0] - cloud_fs[0].gas.T[0];
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
	    rhs[0] += dx*dT; rhs[1] += dy*dT; rhs[2] += dz*dT;
	}
	solveGradients!3(xTx, rhs, gradients);
	T.refx = gradients[0];
	T.refy = gradients[1];
	T.refz = gradients[2];
	// massf
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (diffusion) {
	    foreach(isp; 0 .. nsp) {
		foreach (j; 0 .. 3) { rhs[j] = 0.0; }
		foreach (i; 1 .. n) {
		    double df = cloud_fs[i].gas.massf[isp] - cloud_fs[0].gas.massf[isp];
		    double dx = cloud_pos[i].x - cloud_pos[0].x;
		    double dy = cloud_pos[i].y - cloud_pos[0].y;
		    double dz = cloud_pos[i].z - cloud_pos[0].z;
		    rhs[0] += dx*df; rhs[1] += dy*df; rhs[2] += dz*df;
		}
		solveGradients!3(xTx, rhs, gradients);
		massf[isp].refx = gradients[0];
		massf[isp].refy = gradients[1];
		massf[isp].refz = gradients[2];
	    } // foreach isp
	} else {
	    foreach(isp; 0 .. nsp) {
		massf[isp].refx = 0.0;
		massf[isp].refy = 0.0;
		massf[isp].refz = 0.0;
	    } // foreach isp
	}
	// tke
	foreach (j; 0 .. 3) { rhs[j] = 0.0; }
	foreach (i; 1 .. n) {
	    double dtke = cloud_fs[i].tke - cloud_fs[0].tke;
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
	    rhs[0] += dx*dtke; rhs[1] += dy*dtke; rhs[2] += dz*dtke;
	}
	solveGradients!3(xTx, rhs, gradients);
	tke.refx = gradients[0];
	tke.refy = gradients[1];
	tke.refz = gradients[2];
	// omega
	foreach (j; 0 .. 3) { rhs[j] = 0.0; }
	foreach (i; 1 .. n) {
	    double domega = cloud_fs[i].omega - cloud_fs[0].omega;
	    double dx = cloud_pos[i].x - cloud_pos[0].x;
	    double dy = cloud_pos[i].y - cloud_pos[0].y;
	    double dz = cloud_pos[i].z - cloud_pos[0].z;
	    rhs[0] += dx*domega; rhs[1] += dy*domega; rhs[2] += dz*domega;
	}
	solveGradients!3(xTx, rhs, gradients);
	omega.refx = gradients[0];
	omega.refy = gradients[1];
	omega.refz = gradients[2];
    } // end gradients_xyz_leastsq()

    @nogc
    void gradients_xy_leastsq(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos, bool diffusion)
    // Fit a linear model to the cloud of flow-quantity points
    // in order to extract approximations to the flow-field gradients.
    // 2D, built as a specialization of the 3D code.
    // Experiment with taking differences about a middle point/value.
    {
	double[4][2] xTx; // normal matrix, augmented to give 4 entries per row
	double[2] rhs, gradients;
	// Assemble and invert the normal matrix.
	// We'll reuse the resulting inverse for each flow-field quantity.
	size_t n = cloud_pos.length;
	double x_mid = 0.0;
	double y_mid = 0.0;
	foreach (i; 0 .. n) {
	    x_mid += cloud_pos[i].x;
	    y_mid += cloud_pos[i].y;
	}
	x_mid /= n;
	y_mid /= n;
	double xx = 0.0;
	double xy = 0.0;
	double yy = 0.0;
	foreach (i; 0 .. n) {
	    double dx = cloud_pos[i].x - x_mid;
	    double dy = cloud_pos[i].y - y_mid;
	    xx += dx*dx; xy += dx*dy; yy += dy*dy;
	}
	xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = 1.0; xTx[0][3] = 0.0;
	xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = 0.0; xTx[1][3] = 1.0;
	computeInverse!2(xTx);
	// x-velocity
	double q_mid;
	string codeForGradients(string qname)
	{
	    string code = "
            q_mid = 0.0;
	    foreach (i; 0 .. n) {
	        q_mid += cloud_fs[i]."~qname~";
	    }
	    q_mid /= n;
	    rhs[0] = 0.0; rhs[1] = 0.0;
	    foreach (i; 0 .. n) {
	        double dq = cloud_fs[i]."~qname~" - q_mid;
	        double dx = cloud_pos[i].x - x_mid;
	        double dy = cloud_pos[i].y - y_mid;
	        rhs[0] += dx*dq; rhs[1] += dy*dq;
	    }
	    solveGradients!2(xTx, rhs, gradients);";
	    return code;
	}
	mixin(codeForGradients("vel.x"));
	vel[0][0] = gradients[0];
	vel[0][1] = gradients[1];
	vel[0][2] = 0.0;
	// y-velocity
	mixin(codeForGradients("vel.y"));
	vel[1][0] = gradients[0];
	vel[1][1] = gradients[1];
	vel[1][2] = 0.0;
	// z-velocity
	vel[2][0] = 0.0;
	vel[2][1] = 0.0;
	vel[2][2] = 0.0;
	// T[0]
	mixin(codeForGradients("gas.T[0]"));
	T.refx = gradients[0];
	T.refy = gradients[1];
	T.refz = 0.0;
	// massf
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]"));
		massf[isp].refx = gradients[0];
		massf[isp].refy = gradients[1];
		massf[isp].refz = 0.0;
	    } // foreach isp
	} else {
	    foreach(isp; 0 .. nsp) {
		massf[isp].refx = 0.0;
		massf[isp].refy = 0.0;
		massf[isp].refz = 0.0;
	    } // foreach isp
	}
	// tke
	mixin(codeForGradients("tke"));
	tke.refx = gradients[0];
	tke.refy = gradients[1];
	tke.refz = 0.0;
	// omega
	mixin(codeForGradients("omega"));
	omega.refx = gradients[0];
	omega.refy = gradients[1];
	omega.refz = 0.0;
    } // end gradients_xy_leastsq()

} // end class FlowGradients


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
