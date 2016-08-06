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


class WLSQGradWorkspace {
    // A place to hold the intermediate results for assembling
    // the normal equations, to allow reuse of the values.
public:
    double[] weights2; // square of the weights for the linear constraint eqns.
    bool compute_about_mid;
    double[] dx, dy, dz;
    size_t loop_init; // starting index for loops: 0=compute_about_mid, 1=compute_about_[0]
    size_t n; // cloud_pos.length;
    double x0, y0, z0; // reference point may be a computed midpoint
    double[6][3] xTx; // normal matrix, augmented to give 6 entries per row

    this()
    {
	// don't need to do anything
    }

    this(const WLSQGradWorkspace other)
    {
	dx = other.dx.dup(); dy = other.dy.dup(); dz = other.dz.dup();
	weights2 = other.weights2.dup();
	compute_about_mid = other.compute_about_mid;
	loop_init = other.loop_init;
	n = other.n;
	x0 = other.x0; y0 = other.y0; z0 = other.z0;
	foreach (i; 0 .. 3) {
	    foreach (j; 0 .. 6) { xTx[i][j] = other.xTx[i][j]; }
	}
    }
} // end class WLSQGradWorkspace


class FlowGradients {
    // Spatial derivatives of the flow quantities
public:
    double[3][3] vel; 
    // velocity derivatives stored as a second-order tensor
    // [[du/dx du/dy du/dz]
    //  [dv/dx dv/dy dv/dz]
    //  [dw/dx dw/dy dw/dz]]
    double[3][] massf; // mass fraction derivatives
    double[3] T; // Temperature derivatives (static temperature only)
    double[3] tke; // turbulence kinetic energy
    double[3] omega; // pseudo vorticity for k-omega turbulence
private:
    LocalConfig myConfig;

public:
    this(ref LocalConfig myConfig)
    {
	this.myConfig = myConfig;
	massf.length = myConfig.gmodel.n_species;
    }

    this(const FlowGradients other)
    {
	foreach(i; 0 .. 3) vel[i][] = other.vel[i][];
	massf.length = other.massf.length;
	foreach(isp; 0 .. other.massf.length) { massf[isp][] = other.massf[isp][]; }
	T[] = other.T[];
	tke[] = other.tke[];
	omega[] = other.omega[];
    }

    @nogc
    void copy_values_from(const FlowGradients other)
    {
	foreach (i; 0 .. 3) { vel[i][] = other.vel[i][]; }
	foreach (isp; 0 .. other.massf.length) { massf[isp][] = other.massf[isp][]; }
	T[] = other.T[];
	tke[] = other.tke[];
	omega[] = other.omega[];
    }

    @nogc
    void accumulate_values_from(const FlowGradients other)
    {
	foreach (i; 0 .. 3) { vel[i][] += other.vel[i][]; }
	foreach (isp; 0 .. massf.length) { massf[isp][] += other.massf[isp][]; }
	T[] += other.T[];
	tke[] += other.tke[];
	omega[] += other.omega[];
    }

    @nogc
    void scale_values_by(double factor)
    {
	foreach (i; 0 .. 3) { vel[i][] *= factor; } 
	foreach (isp; 0 .. massf.length) { massf[isp][] *= factor; } 
	T[] *= factor;
	tke[] *= factor;
	omega[] *= factor;
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
    void gradients_xy_div(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos)
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
	T[0] = gradient_x * area_inv;
	T[1] = -gradient_y * area_inv;
	T[2] = 0.0;
	//
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (myConfig.diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]"));
		massf[isp][0] = gradient_x * area_inv;
		massf[isp][1] = -gradient_y * area_inv;
		massf[isp][2] = 0.0;
	    }
	} else {
	    foreach(isp; 0 .. nsp) {
		massf[isp][0] = 0.0;
		massf[isp][1] = 0.0;
		massf[isp][2] = 0.0;
	    }
	}
	//
	mixin(codeForGradients("tke"));
	tke[0] = gradient_x * area_inv;
	tke[1] = -gradient_y * area_inv;
	tke[2] = 0.0;
	//
	mixin(codeForGradients("omega"));
	omega[0] = gradient_x * area_inv;
	omega[1] = -gradient_y * area_inv;
	omega[2] = 0.0;
    } // end gradients_xy_div()

    void weights_leastsq(ref Vector3*[] cloud_pos, ref Vector3 pos,
			 bool compute_about_mid, ref WLSQGradWorkspace ws)
    // Calculate weights used in the least-squares gradient calculation.
    // The weights are calculated with the current interface/vertex as the reference point (pos).
    // NB. For the "faces" spatial location we are expecting the primary point (i.e. the face at which
    //     we are calculating the gradients) to be in the first cloud position. 
    {
	size_t n = cloud_pos.length;
	size_t loop_init;
	if (ws.weights2.length < n) ws.weights2.length = n;
	ws.compute_about_mid = compute_about_mid;
	if (compute_about_mid) {
	    // We will compute differences about a "mid" point that is assumed
	    // to be distinct from all points in the incoming cloud.
	    loop_init = 0; // All points count.
	} else {
	    loop_init = 1; // first point is reference for differences.
	    ws.weights2[0] = 0.0; // and doesn't enter into the sum itself.
	}
	ws.n = n;
	ws.loop_init = loop_init;
	double x0 = pos.x;
	double y0 = pos.y;
	if (myConfig.dimensions == 2) {
	    foreach (i; loop_init .. n) {
		double dx = cloud_pos[i].x - x0;
		double dy = cloud_pos[i].y - y0;
		ws.weights2[i] = 1.0/(dx*dx+dy*dy);
	    }
	} else { //3D
	    double z0 = pos.z;
	    foreach (i; loop_init .. n) {
		double dx = cloud_pos[i].x - x0;
		double dy = cloud_pos[i].y - y0;
		double dz = cloud_pos[i].z - z0;
		ws.weights2[i] = 1.0/(dx*dx+dy*dy+dz*dz);
	    }
	}
    } // end weights_leastsq()

    void set_up_workspace_for_gradients_xyz_leastsq(ref Vector3*[] cloud_pos,
						    ref WLSQGradWorkspace ws)
    {
	size_t loop_init = 0; // starting index for loops: 0=compute_about_mid, 1=compute_about_[0]
	size_t n = cloud_pos.length;
	// If computing about mid-point, calculate mid-point.
	double x0 = 0.0; double y0 = 0.0; double z0 = 0.0;
	if (ws.compute_about_mid) {
	    loop_init = 0;
	    foreach (i; loop_init .. n) {
		x0 += cloud_pos[i].x; y0 += cloud_pos[i].y; z0 += cloud_pos[i].z;
	    }
	    x0 /= n; y0 /= n; z0 /= n; // midpoint
	} else { // else use interface point (assumed to be in cloud position 0)
	    loop_init = 1;
	    x0 = cloud_pos[0].x; y0 = cloud_pos[0].y; z0 = cloud_pos[0].z;
	}
	assert(ws, "We are missing the workspace!");
	ws.n = n;
	ws.loop_init = loop_init;
	ws.x0 = x0; ws.y0 = y0; ws.z0 = z0;
	if (ws.dx.length < n) { ws.dx.length = n; }
	if (ws.dy.length < n) { ws.dy.length = n; }
	if (ws.dz.length < n) { ws.dz.length = n; }
	//
	// Assemble and invert the normal matrix.
	// We'll reuse the resulting inverse for each flow-field quantity.
	double xx = 0.0; double xy = 0.0; double xz = 0.0;
	double yy = 0.0; double yz = 0.0; double zz = 0.0;
	foreach (i; loop_init .. n) {
	    double dx = cloud_pos[i].x - x0;
	    double dy = cloud_pos[i].y - y0;
	    double dz = cloud_pos[i].z - z0;
	    xx += ws.weights2[i]*dx*dx; xy += ws.weights2[i]*dx*dy; xz += ws.weights2[i]*dx*dz;
	    yy += ws.weights2[i]*dy*dy; yz += ws.weights2[i]*dy*dz; zz += ws.weights2[i]*dz*dz;
	    ws.dx[i] = dx; ws.dy[i] = dy; ws.dz[i] = dz;
	}
	ws.xTx[0][0] = xx; ws.xTx[0][1] = xy; ws.xTx[0][2] = xz;
	ws.xTx[1][0] = xy; ws.xTx[1][1] = yy; ws.xTx[1][2] = yz;
	ws.xTx[2][0] = xz; ws.xTx[2][1] = yz; ws.xTx[2][2] = zz;
	ws.xTx[0][3] = 1.0; ws.xTx[0][4] = 0.0; ws.xTx[0][5] = 0.0;
	ws.xTx[1][3] = 0.0; ws.xTx[1][4] = 1.0; ws.xTx[1][5] = 0.0;
	ws.xTx[2][3] = 0.0; ws.xTx[2][4] = 0.0; ws.xTx[2][5] = 1.0;
	computeInverse!(3,3)(ws.xTx);
    } // end set_up_workspace_for_gradients_xyz_leastsq()

    void gradients_xyz_leastsq(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos,
			       ref WLSQGradWorkspace ws)
    // Fit a linear model to the cloud of flow-quantity points
    // in order to extract approximations to the flow-field gradients.
    // 3D
    // As for 2D, if using vertices spatial locations then take differences about a middle point/value.
    // Else take differences about the point at which the gradient is being calculated (i.e. the faces).
    {
	size_t loop_init = ws.loop_init;
	size_t n = ws.n;
	double[3] rhs;
	//
	double q0;
	string codeForGradients(string qname, string gname)
	{
	    string code = "
            if (ws.compute_about_mid) {
                q0 = 0.0;
                foreach (i; loop_init .. n) { q0 += cloud_fs[i]."~qname~"; }
	        q0 /= n;
            } else { 
                q0 = cloud_fs[0]."~qname~";
            }
	    foreach (j; 0 .. 3) { rhs[j] = 0.0; }
            foreach (i; loop_init .. n) {
                double dq = cloud_fs[i]."~qname~" - q0;
                rhs[0] += ws.weights2[i] * ws.dx[i] * dq;
                rhs[1] += ws.weights2[i] * ws.dy[i] * dq;
                rhs[2] += ws.weights2[i] * ws.dz[i] * dq;
	    }
	    solveWithInverse!(3,3)(ws.xTx, rhs, "~gname~");";
	    return code;
	}
	mixin(codeForGradients("vel.x", "vel[0]"));
	mixin(codeForGradients("vel.y", "vel[1]"));
	mixin(codeForGradients("vel.z", "vel[2]"));	
	mixin(codeForGradients("gas.T[0]", "T"));
	// massf
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (myConfig.diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]", "massf[isp]"));
	    }
	} else {
	    foreach(isp; 0 .. nsp) {
		massf[isp][0] = 0.0; massf[isp][1] = 0.0; massf[isp][2] = 0.0;
	    }
	}
	mixin(codeForGradients("tke", "tke"));
	mixin(codeForGradients("omega", "omega"));
    } // end gradients_xyz_leastsq()

    void set_up_workspace_for_gradients_xy_leastsq(ref Vector3*[] cloud_pos,
						   ref WLSQGradWorkspace ws)
    {
	size_t loop_init = 0; // starting index for loops: 0=compute_about_mid, 1=compute_about_[0]
	size_t n = cloud_pos.length;
	// If computing about mid-point, calculate mid-point.
	double x0 = 0.0; double y0 = 0.0; double z0 = 0.0;
	if (ws.compute_about_mid) {
	    loop_init = 0;
	    foreach (i; loop_init .. n) {
		x0 += cloud_pos[i].x; y0 += cloud_pos[i].y;
	    }
	    x0 /= n; y0 /= n; // midpoint
	} else { // else use interface point (assumed to be in cloud position 0)
	    loop_init = 1;
	    x0 = cloud_pos[0].x; y0 = cloud_pos[0].y;
	}
	assert(ws, "We are missing the workspace!");
	ws.n = n;
	ws.loop_init = loop_init;
	ws.x0 = x0; ws.y0 = y0; ws.z0 = z0;
	if (ws.dx.length < n) { ws.dx.length = n; }
	if (ws.dy.length < n) { ws.dy.length = n; }
	if (ws.dz.length < n) { ws.dz.length = n; }
	//
	// Assemble and invert the normal matrix.
	// We'll reuse the resulting inverse for each flow-field quantity.
	double xx = 0.0; double xy = 0.0; double yy = 0.0;
	foreach (i; loop_init .. n) {
	    double dx = cloud_pos[i].x - x0;
	    double dy = cloud_pos[i].y - y0;
	    xx += ws.weights2[i]*dx*dx; xy += ws.weights2[i]*dx*dy; yy += ws.weights2[i]*dy*dy;
	    ws.dx[i] = dx; ws.dy[i] = dy; ws.dz[i] = 0.0;
	}
	ws.xTx[0][0] = xx; ws.xTx[0][1] = xy; ws.xTx[0][2] = 1.0; ws.xTx[0][3] = 0.0;
	ws.xTx[1][0] = xy; ws.xTx[1][1] = yy; ws.xTx[1][2] = 0.0; ws.xTx[1][3] = 1.0;
	computeInverse!(2,3)(ws.xTx);
    } // end set_up_workspace_for_gradients_xy_leastsq()

    void gradients_xy_leastsq(ref FlowState[] cloud_fs, ref Vector3*[] cloud_pos,
			      ref WLSQGradWorkspace ws)
    // Fit a linear model to the cloud of flow-quantity points
    // in order to extract approximations to the flow-field gradients.
    // 2D, built as a specialization of the 3D code.
    // Experiment with taking differences about a middle point/value for the vertices spatial locations and
    // taking differences about the interface point for the faces spatial location.
    {
	size_t loop_init = ws.loop_init;
	size_t n = ws.n;
	double[3] rhs;
	double q0;
	string codeForGradients(string qname, string gname)
	{
	    string code = "
            if (ws.compute_about_mid) {
                q0 = 0.0;
                foreach (i; loop_init .. n) { q0 += cloud_fs[i]."~qname~"; }
	        q0 /= n;
            } else {
                q0 = cloud_fs[0]."~qname~";
            }
	    rhs[0] = 0.0; rhs[1] = 0.0;
	    foreach (i; loop_init .. n) {
	        double dq = cloud_fs[i]."~qname~" - q0;
                rhs[0] += ws.weights2[i] * ws.dx[i] * dq;
                rhs[1] += ws.weights2[i] * ws.dy[i] * dq;
	    }
	    solveWithInverse!(2,3)(ws.xTx, rhs, "~gname~");
            "~gname~"[2] = 0.0;";
	    return code;
	}
	mixin(codeForGradients("vel.x", "vel[0]"));
	mixin(codeForGradients("vel.y", "vel[1]"));
	vel[2][0] = 0.0; vel[2][1] = 0.0; vel[2][2] = 0.0; // z-velocity
	mixin(codeForGradients("gas.T[0]", "T"));
	//
	size_t nsp = cloud_fs[0].gas.massf.length;
	if (myConfig.diffusion) {
	    foreach(isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]", "massf[isp]"));
	    }
	} else {
	    foreach(isp; 0 .. nsp) {
		massf[isp][0] = 0.0; massf[isp][1] = 0.0; massf[isp][2] = 0.0;
	    }
	}
	mixin(codeForGradients("tke", "tke"));
	mixin(codeForGradients("omega", "omega"));
    } // end gradients_xy_leastsq()

} // end class FlowGradients
