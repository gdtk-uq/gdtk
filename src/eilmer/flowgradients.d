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
import nm.rsla;
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
	repr ~= "vel=[";
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
    void gradients_finitediff(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos, FVInterface iface, bool diffusion)
    {
	// Uses finite differences to calculate the gradient at the interface using the distance normal to the face.
	
	/+ This is the first attempt at a finite difference approach to calculating the interface gradients used in
	 + the viscous flux update. Currently the code correctly calculates gradients parallel to the face normal, 
	 + however does not capture the gradients perpendicular to the face normal. For example, looking at a 2D cell 
	 + face with it's normal aligned with the global x-direction, the code
	 + calculates du/dx and dv/dx correctly but does not capture dv/dy and du/dy.
	 + 
	 + TO_DO: FIX ME
         +/
	size_t count = cloud_pos.length;
	double nx = iface.n.x;
	double ny = iface.n.y;	
	double nz = iface.n.z;
	double gradn;

	string codeForGradients(string qname)
	{
	    string code = "
            gradn = 0.0;
            foreach (i; 0 .. count) {
                double dx = (iface.pos.x - cloud_pos[i].x);
	        double dy = (iface.pos.y - cloud_pos[i].y);
                double dz = (iface.pos.z - cloud_pos[i].z);
                double dn = dx*nx + dy*ny + dz*nz;
                double dqdn = (iface.fs."~qname~" - cloud_fs[i]."~qname~")/dn;
                gradn += dqdn;
            }
            gradn /= count;";
	    return code;
	}
	mixin(codeForGradients("vel.x"));
	vel[0][0] = gradn * nx;
	vel[0][1] = gradn * ny;
	vel[0][2] = gradn * nz;
	// y-velocity
	mixin(codeForGradients("vel.y"));
	vel[1][0] = gradn * nx;
	vel[1][1] = gradn * ny;
	vel[1][2] = gradn * nz;
	// z-velocity
	mixin(codeForGradients("vel.z"));
	vel[2][0] = gradn * nx;
	vel[2][1] = gradn * ny;
	vel[2][2] = gradn * nz;
	// T[0]
	mixin(codeForGradients("gas.T[0]"));
	T.refx = gradn * nx;
	T.refy = gradn * ny;
	T.refz = gradn * nz;
	// massf
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]"));
		massf[isp].refx = gradn * nx;
		massf[isp].refy = gradn * ny;
		massf[isp].refz = gradn * nz;
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
	tke.refx = gradn * nx;
	tke.refy = gradn * ny;
	tke.refz = gradn * nz;
	// omega
	mixin(codeForGradients("omega"));
	omega.refx = gradn * nx;
	omega.refy = gradn * ny;
	omega.refz = gradn * nz;	
    } // end gradients_finitediff
     
  @nogc
  void gradients_xyz_leastsq(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos, bool diffusion)
    // Fit a linear model to the cloud of flow-quantity points
    // in order to extract approximations to the flow-field gradients.
    // 3D
    // As for 2D, take differences about a middle point/value.
    {
	double[6][3] xTx; // normal matrix, augmented to give 6 entries per row
	double[3] rhs, gradients;
	// Assemble and invert the normal matrix.
	// We'll reuse the resulting inverse for each flow-field quantity.
	size_t n = cloud_pos.length;
	double x_mid = 0.0;
	double y_mid = 0.0;
	double z_mid = 0.0;
	foreach (i; 0 .. n) {
	    x_mid += cloud_pos[i].x;
	    y_mid += cloud_pos[i].y;
	    z_mid += cloud_pos[i].z;
	}
	x_mid /= n;
	y_mid /= n;
	z_mid /= n;
	double xx = 0.0;
	double xy = 0.0;
	double xz = 0.0;
	double yy = 0.0;
	double yz = 0.0;
	double zz = 0.0;
	foreach (i; 0 .. n) {
	    double dx = cloud_pos[i].x - x_mid;
	    double dy = cloud_pos[i].y - y_mid;
	    double dz = cloud_pos[i].z - z_mid;
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
	double q_mid;
	string codeForGradients(string qname)
	{
	    string code = "
            q_mid = 0.0;
	    foreach (i; 0 .. n) {
	        q_mid += cloud_fs[i]."~qname~";
	    }
	    q_mid /= n;
            foreach (j; 0 .. 3) { rhs[j] = 0.0; }
            foreach (i; 0 .. n) {
                double dq = cloud_fs[i]."~qname~" - q_mid;
                double dx = cloud_pos[i].x - x_mid;
                double dy = cloud_pos[i].y - y_mid;
                double dz = cloud_pos[i].z - z_mid;
                rhs[0] += dx*dq; rhs[1] += dy*dq; rhs[2] += dz*dq;
	    }
	    solveWithInverse!3(xTx, rhs, gradients);";
	    return code;
	}
	mixin(codeForGradients("vel.x"));
	foreach (j; 0 .. 3) { vel[0][j] = gradients[j]; }
	// y-velocity
	mixin(codeForGradients("vel.y"));
	foreach (j; 0 .. 3) { vel[1][j] = gradients[j]; }
	// z-velocity
	mixin(codeForGradients("vel.z"));	
	foreach (j; 0 .. 3) { vel[2][j] = gradients[j]; }
	// T[0]
	mixin(codeForGradients("gas.T[0]"));
	T.refx = gradients[0];
	T.refy = gradients[1];
	T.refz = gradients[2];
	// massf
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]"));
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
	mixin(codeForGradients("tke"));
	tke.refx = gradients[0];
	tke.refy = gradients[1];
	tke.refz = gradients[2];
	// omega
	mixin(codeForGradients("omega"));
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
	    double w = 1.0/(dx*dx+dy*dy);
	    xx += w*dx*dx; xy += w*dx*dy; yy += w*dy*dy;
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
                double w = 1.0/(dx*dx+dy*dy);
                rhs[0] += w*dx*dq; rhs[1] += w*dy*dq;
	    }
	    solveWithInverse!2(xTx, rhs, gradients);";
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
