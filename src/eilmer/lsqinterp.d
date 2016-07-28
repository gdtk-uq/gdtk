// lsqinterp.d
// Least-squares interpolation/reconstruction of flow field.
//

import std.math;
import std.stdio;
import std.algorithm;
import std.conv;
import nm.rsla;
import geom;
import gas;
import fvcore;
import globalconfig;
import flowstate;
import fvinterface;
import fvcell;


class LSQInterpWorkspace {
public:
    // A place to hold the intermediate results for assembling
    // the normal equations, to allow reuse of the values.
    double[] dx, dy, dz;
    double[6][3] xTx; // normal matrix, augmented to give 6 entries per row

    this()
    {
	// do nothing
    }

    this(ref LSQInterpWorkspace other)
    {
	// Cannot have const in constructor signature because of the duplication of references.
	dx = other.dx.dup(); dy = other.dy.dup(); dz = other.dz.dup();
	foreach (i; 0 .. 3) {
	    foreach (j; 0 .. 6) { xTx[i][j] = other.xTx[i][j]; }
	}
    }
    
    void assemble_and_invert_normal_matrix(FVCell[] cell_cloud, int dimensions, size_t gtl)
    {
	auto np = cell_cloud.length;
	dx.length = np; dy.length = np; dz.length = np;
	foreach (i; 1 .. np) {
	    Vector3 dr = cell_cloud[i].pos[gtl]; dr -= cell_cloud[0].pos[gtl];
	    dx[i] = dr.x; dy[i] = dr.y; dz[i] = dr.z;
	}	    
	// Prepare the normal matrix for the cloud and invert it.
	// Depending on how aligned the points are, there may be insufficient
	// variation in the local y- or z-positions to make a sensible 3D matrix.
	// We will form the sums and look at their relative sizes.
	double xx = 0.0; double xy = 0.0; double xz = 0.0;
	double yy = 0.0; double yz = 0.0;
	double zz = 0.0;
	foreach (i; 1 .. np) {
	    xx += dx[i]*dx[i]; xy += dx[i]*dy[i]; xz += dx[i]*dz[i];
	    yy += dy[i]*dy[i]; yz += dy[i]*dz[i]; zz += dz[i]*dz[i];
	}
	if (dimensions == 3) {
	    xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
	    xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
	    xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
	} else {
	    // dimensions == 2
	    xTx[0][0] =  xx; xTx[0][1] =  xy; xTx[0][2] = 0.0;
	    xTx[1][0] =  xy; xTx[1][1] =  yy; xTx[1][2] = 0.0;
	    xTx[2][0] = 0.0; xTx[2][1] = 0.0; xTx[2][2] = 1.0;
	}
	xTx[0][3] = 1.0; xTx[0][4] = 0.0; xTx[0][5] = 0.0;
	xTx[1][3] = 0.0; xTx[1][4] = 1.0; xTx[1][5] = 0.0;
	xTx[2][3] = 0.0; xTx[2][4] = 0.0; xTx[2][5] = 1.0;
	if (0 != computeInverse!(3,3)(xTx)) {
	    throw new FlowSolverException("Failed to invert LSQ normal matrix");
	    // Assume that the rows are linearly dependent 
	    // because the sample points are colinear.
	    // Maybe we could proceed by working as a single-dimensional interpolation.
	}
    } // end assemble_and_invert_normal_matrix()
} // end class LSQInterpWorkspace

class LSQInterpGradients {
    // These are the quantities that will be interpolated from cell centres 
    // to the Left and Right sides of interfaces.
    // We need to hold onto their gradients within cells.
public:
    double velx_x, velx_y, velx_z;
    double vely_x, vely_y, vely_z;
    double velz_x, velz_y, velz_z;
    double Bx_x, Bx_y, Bx_z;
    double By_x, By_y, By_z;
    double Bz_x, Bz_y, Bz_z;
    double psi_x, psi_y, psi_z;
    double tke_x, tke_y, tke_z;
    double omega_x, omega_y, omega_z;
    double[] massf_x, massf_y, massf_z;
    double rho_x, rho_y, rho_z;
    double p_x, p_y, p_z;
    double[] T_x, T_y, T_z;
    double[] e_x, e_y, e_z;

    this(size_t nsp, size_t nmodes)
    {
	massf_x.length = nsp; massf_y.length = nsp; massf_z.length = nsp;
	T_x.length = nmodes; T_y.length = nmodes; T_z.length = nmodes;
	e_x.length = nmodes; e_y.length = nmodes; e_z.length = nmodes;
    }

    this(ref const LSQInterpGradients other)
    {
	velx_x = other.velx_x; velx_y = other.velx_y; velx_z = other.velx_z;
	vely_x = other.vely_x; vely_y = other.vely_y; vely_z = other.vely_z;
	velz_x = other.velz_x; velz_y = other.velz_y; velz_z = other.velz_z;
	Bx_x = other.Bx_x; Bx_y = other.Bx_y; Bx_z = other.Bx_z;
	By_x = other.By_x; By_y = other.By_y; By_z = other.By_z;
	Bz_x = other.Bz_x; Bz_y = other.Bz_y; Bz_z = other.Bz_z;
	psi_x = other.psi_x; psi_y = other.psi_y; psi_z = other.psi_z;
	tke_x = other.tke_x; tke_y = other.tke_y; tke_z = other.tke_z;
	omega_x = other.omega_x; omega_y = other.omega_y; omega_z = other.omega_z;
	foreach(i; 0 .. other.massf_x.length) {
	    massf_x[i] = other.massf_x[i];
	    massf_y[i] = other.massf_y[i];
	    massf_z[i] = other.massf_z[i];
	}
	rho_x = other.rho_x; rho_y = other.rho_y; rho_z = other.rho_z;
	p_x = other.p_x; p_y = other.p_y; p_z = other.p_z;
	foreach(i; 0 .. other.T_x.length) {
	    T_x[i] = other.T_x[i];
	    T_y[i] = other.T_y[i];
	    T_z[i] = other.T_z[i];
	}
	foreach(i; 0 .. other.e_x.length) {
	    e_x[i] = other.e_x[i];
	    e_y[i] = other.e_y[i];
	    e_z[i] = other.e_z[i];
	}
    }

    void compute_lsq_values(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
	double[3] rhs, gradients;
	auto np = cell_cloud.length;
	// The following function to be used at compile time.
	string codeForGradients(string qname, string gradname)
	{
	    string code = "{
                double q0 = cell_cloud[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhs[j] = 0.0; }
                foreach (i; 1 .. np) {
                    double dq = cell_cloud[i].fs."~qname~" - q0;
                    rhs[0] += ws.dx[i]*dq; rhs[1] += ws.dy[i]*dq; rhs[2] += ws.dz[i]*dq;
                }
	        solveWithInverse!(3,3)(ws.xTx, rhs, gradients);
                "~gradname~"_x = gradients[0];
                "~gradname~"_y = gradients[1];
                "~gradname~"_z = gradients[2];
                }
                ";
	    return code;
	}
	string codeForGradientArrays(string qname, string gradname, string indexname)
	{
	    string code = "{
                double q0 = cell_cloud[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhs[j] = 0.0; }
                foreach (i; 1 .. np) {
                    double dq = cell_cloud[i].fs."~qname~" - q0;
                    rhs[0] += ws.dx[i]*dq; rhs[1] += ws.dy[i]*dq; rhs[2] += ws.dz[i]*dq;
                }
	        solveWithInverse!(3,3)(ws.xTx, rhs, gradients);
                "~gradname~"_x["~indexname~"] = gradients[0];
                "~gradname~"_y["~indexname~"] = gradients[1];
                "~gradname~"_z["~indexname~"] = gradients[2];
                }
                ";
	    return code;
	}
	// x-velocity
	mixin(codeForGradients("vel.x", "velx"));
	mixin(codeForGradients("vel.y", "vely"));
	mixin(codeForGradients("vel.z", "velz"));
	if (myConfig.MHD) {
	    mixin(codeForGradients("B.x", "Bx"));
	    mixin(codeForGradients("B.y", "By"));
	    mixin(codeForGradients("B.z", "Bz"));
	    if (myConfig.divergence_cleaning) {
		mixin(codeForGradients("psi", "psi"));
	    }
	}
	if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
	    mixin(codeForGradients("tke", "tke"));
	    mixin(codeForGradients("omega", "omega"));
	}
	auto nsp = myConfig.gmodel.n_species;
	if (nsp > 1) {
	    // Multiple species.
	    foreach (isp; 0 .. nsp) {
		mixin(codeForGradientArrays("gas.massf[isp]", "massf", "isp"));
	    }
	} else {
	    // Only one possible gradient value for a single species.
	    massf_x[0] = 0.0; massf_y[0] = 0.0; massf_z[0] = 0.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	auto nmodes = myConfig.gmodel.n_modes;
	final switch (myConfig.thermo_interpolator) {
	case InterpolateOption.pt: 
	    mixin(codeForGradients("gas.p", "p"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradientArrays("gas.T[imode]", "T", "imode"));
	    }
	    break;
	case InterpolateOption.rhoe:
	    mixin(codeForGradients("gas.rho", "rho"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradientArrays("gas.e[imode]", "e", "imode"));
	    }
	    break;
	case InterpolateOption.rhop:
	    mixin(codeForGradients("gas.rho", "rho"));
	    mixin(codeForGradients("gas.p", "p"));
	    break;
	case InterpolateOption.rhot: 
	    mixin(codeForGradients("gas.rho", "rho"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradientArrays("gas.T[imode]", "T", "imode"));
	    }
	    break;
	} // end switch thermo_interpolator
    } // end compute_lsq_gradients()
} // end class LSQInterpGradients


class LsqInterpolator {

private:
    LocalConfig myConfig;

public:
    this(LocalConfig myConfig) 
    {
	this.myConfig = myConfig;
    }

    @nogc double clip_to_limits(double q, double A, double B)
    // Returns q if q is between the values A and B, else
    // it returns the closer limit of the range [A,B].
    {
	double lower_limit = fmin(A, B);
	double upper_limit = fmax(A, B);
	return fmin(upper_limit, fmax(lower_limit, q));
    } // end clip_to_limits()

    @nogc void min_mod_limit(ref double a, ref double b)
    // A rough slope limiter.
    {
	if (a * b < 0.0) {
	    a = 0.0; b = 0.0;
	    return;
	}
	a = copysign(fmin(fabs(a), fabs(b)), a);
	b = a;
    }

    @nogc void van_albada_limit(ref double a, ref double b)
    // A smooth slope limiter.
    {
	immutable double eps = 1.0e-12;
	double s = (a*b + fabs(a*b))/(a*a + b*b + eps);
	a *= s;
	b *= s;
    }

    void venkatakrishan_limit(ref double[3] grad, FVCell cell, double U, double Umin, double Umax,
			      double dx, double dy, double dz)
    // A smooth slope limiter developed for unstructured grids (improvement on Barth limiter).
    // Reported to have 2nd order accuracy at the cost of exhibiting non-montonicity near shocks
    // -- small oscillations in regions of large gradients can occur.
    {
	immutable double w = 1.0e-12;
	immutable double K = 100.0;
	double a, b, numer, denom, s;
	double[] phi;
	double h;
	if (myConfig.dimensions == 3) h =  cbrt(dx*dy*dz);
	else h = sqrt(dx*dy);
	double eps = (K*h) * (K*h) * (K*h);
	if (cell.iface.length == 0) s = 1.0; // if ghost cell => do not limit
	else {
	    foreach (i, f; cell.iface) { // loop over gauss points
		Vector3 dr = f.pos; dr -= cell.pos[0];
		dr.transform_to_local_frame(f.n, f.t1, f.t2);
		double dxFace = dr.x; double dyFace = dr.y; double dzFace = dr.z;
		b = grad[0] * dxFace + grad[1] * dyFace;
		if (myConfig.dimensions == 3) b += grad[2] * dzFace;
		b = sgn(b) * (fabs(b) + w);
		if (b > 0.0) a = Umax - U;
		else if (b < 0.0) a = Umin - U;
		numer = (a*a + eps)*b + 2.0*b*b*a;
		denom = a*a + 2.0*b*b + a*b + eps;
		s = (1.0/b) * (numer/denom);
		if (b == 0.0) s = 1.0;
		phi ~= s;
	    }
	    s = phi[0];
	    foreach (i; 1..phi.length) { // take the mininum limiting factor 
		if (phi[i] < s) s = phi[i];
	    }
	}
	grad[0] *= s;
	grad[1] *= s;
	if (myConfig.dimensions == 3) grad[2] *= s;
    }

    void barth_limit(ref double[3] grad, FVCell cell, double U, double Umin, double Umax,
		     double dx, double dy, double dz)
    // A non-differentiable slope limiter developed for unstructured grids.
    // This limiter exhibits monotonicity, although it has poor convergence order.
    {
	double a, b, numer, denom, s;
	double[] phi;
	immutable double w = 1.0e-12;
	if (cell.iface.length == 0) s = 1.0; // if ghost cell => do not limit
	else {
	    foreach (i, f; cell.iface) { // loop over gauss points
		Vector3 dr = f.pos; dr -= cell.pos[0];
		dr.transform_to_local_frame(f.n, f.t1, f.t2);
		double dxFace = dr.x; double dyFace = dr.y; double dzFace = dr.z;
		b = grad[0] * dxFace + grad[1] * dyFace;
		if (myConfig.dimensions == 3) b += grad[2] * dzFace;
		if (b > 0.0) {
		    a = Umax - U;
		    phi ~= min(1.0, a/b);
		}
		else if (b < 0.0) {
		    a = Umin - U;
		    phi ~= min(1.0, a/b);
		}
		else {
		    phi ~= 1.0;
		}
	    }
	    s = phi[0];
	    foreach (i; 1..phi.length) { // take the mininum limiting factor
		if (phi[i] < s) s = phi[i];
	    }
	}
	grad[0] *= s;
	grad[1] *= s;
	if (myConfig.dimensions == 3) grad[2] *= s;
    }
    
    void interp_both(ref FVInterface IFace, size_t gtl, ref FlowState Lft, ref FlowState Rght)
    {
	auto gmodel = myConfig.gmodel;
	auto nsp = gmodel.n_species;
	auto nmodes = gmodel.n_modes;
	FVCell cL0 = IFace.left_cell;
	FVCell cR0 = IFace.right_cell;
	// Low-order reconstruction just copies data from adjacent FV_Cell.
	// Even for high-order reconstruction, we depend upon this copy for
	// the viscous-transport and diffusion coefficients.
	Lft.copy_values_from(cL0.fs);
	Rght.copy_values_from(cR0.fs);
	if (myConfig.interpolation_order > 1) {
	    // High-order reconstruction for some properties.
	    //
	    LSQInterpWorkspace wsL = cL0.ws;
	    LSQInterpWorkspace wsR = cR0.ws;
	    Vector3 dL = IFace.pos; dL -= cL0.pos[gtl]; // vector from left-cell-centre to face midpoint
	    Vector3 dR = IFace.pos; dR -= cR0.pos[gtl];
	    double mygradL_x, mygradL_y, mygradL_z;
	    double mygradR_x, mygradR_y, mygradR_z;
	    //
	    // Always reconstruct in the global frame of reference -- for now
	    //
	    // x-velocity
	    string codeForReconstruction(string qname, string gname, string tname)
	    {
		string code = "{
                double qL0 = cL0.fs."~qname~";
                // double qMinL = qL0;
                // double qMaxL = qL0;
                if (wsL) { // an active cell will have a workspace
                    // foreach (i; 1 .. cL0.cell_cloud.length) {
                    //     qMinL = min(qMinL, cL0.cell_cloud[i].fs."~qname~");
                    //     qMaxL = max(qMaxL, cL0.cell_cloud[i].fs."~qname~");
                    // }
                    mygradL_x = cL0.gradients."~gname~"_x;
                    mygradL_y = cL0.gradients."~gname~"_y;
                    mygradL_z = cL0.gradients."~gname~"_z;
                } else {
                    mygradL_x = 0.0;
                    mygradL_y = 0.0;
                    mygradL_z = 0.0;
                }
                double qR0 = cR0.fs."~qname~";
                // double qMinR = qR0;
                // double qMaxR = qR0;
                if (wsR) { // an active cell will have a workspace
                    // foreach (i; 1 .. cR0.cell_cloud.length) {
                    //     qMinR = min(qMinR, cR0.cell_cloud[i].fs."~qname~");
                    //     qMaxR = max(qMaxR, cR0.cell_cloud[i].fs."~qname~");
                    // }
                    mygradR_x = cR0.gradients."~gname~"_x;
                    mygradR_y = cR0.gradients."~gname~"_y;
                    mygradR_z = cR0.gradients."~gname~"_z;
                } else {
                    mygradR_x = 0.0;
                    mygradR_y = 0.0;
                    mygradR_z = 0.0;
                }
                if (myConfig.apply_limiter) {
                    // venkatakrishan_limit(cL0.gradients."~gname~", cL0, qL0, qMinL, qMaxL, cL0.iLength, cL0.jLength, cL0.kLength);
                    // venkatakrishan_limit(cR0.gradients."~gname~", cR0, qR0, qMinR, qMaxR, cR0.iLength, cR0.jLength, cR0.kLength);
                    van_albada_limit(mygradL_x, mygradR_x);
                    van_albada_limit(mygradL_y, mygradR_y);
                    van_albada_limit(mygradL_z, mygradR_z);
                }
                double qL = qL0 + dL.x*mygradL_x + dL.y*mygradL_y + dL.z*mygradL_z;
                double qR = qR0 + dR.x*mygradR_x + dR.y*mygradR_y + dR.z*mygradR_z;
                if (myConfig.extrema_clipping) {
                    Lft."~tname~" = clip_to_limits(qL, qL0, qR0);
                    Rght."~tname~" = clip_to_limits(qR, qL0, qR0);
                } else {
                    Lft."~tname~" = qL;
                    Rght."~tname~" = qR;
                }
                }
                ";
		return code;
	    }
	    string codeForReconstructionOfArray(string qname, string gname, string tname, string indexname)
	    {
		string code = "{
                double qL0 = cL0.fs."~qname~";
                // double qMinL = qL0;
                // double qMaxL = qL0;
                if (wsL) { // an active cell will have a workspace
                    // foreach (i; 1 .. cL0.cell_cloud.length) {
                    //     qMinL = min(qMinL, cL0.cell_cloud[i].fs."~qname~");
                    //     qMaxL = max(qMaxL, cL0.cell_cloud[i].fs."~qname~");
                    // }
                    mygradL_x = cL0.gradients."~gname~"_x["~indexname~"];
                    mygradL_y = cL0.gradients."~gname~"_y["~indexname~"];
                    mygradL_z = cL0.gradients."~gname~"_z["~indexname~"];
                } else {
                    mygradL_x = 0.0;
                    mygradL_y = 0.0;
                    mygradL_z = 0.0;
                }
                double qR0 = cR0.fs."~qname~";
                // double qMinR = qR0;
                // double qMaxR = qR0;
                if (wsR) { // an active cell will have a workspace
                    // foreach (i; 1 .. cR0.cell_cloud.length) {
                    //     qMinR = min(qMinR, cR0.cell_cloud[i].fs."~qname~");
                    //     qMaxR = max(qMaxR, cR0.cell_cloud[i].fs."~qname~");
                    // }
                    mygradR_x = cR0.gradients."~gname~"_x["~indexname~"];
                    mygradR_y = cR0.gradients."~gname~"_y["~indexname~"];
                    mygradR_z = cR0.gradients."~gname~"_z["~indexname~"];
                } else {
                    mygradR_x = 0.0;
                    mygradR_y = 0.0;
                    mygradR_z = 0.0;
                }
                if (myConfig.apply_limiter) {
                    // venkatakrishan_limit(cL0.gradients."~gname~", cL0, qL0, qMinL, qMaxL, cL0.iLength, cL0.jLength, cL0.kLength);
                    // venkatakrishan_limit(cR0.gradients."~gname~", cR0, qR0, qMinR, qMaxR, cR0.iLength, cR0.jLength, cR0.kLength);
                    van_albada_limit(mygradL_x, mygradR_x);
                    van_albada_limit(mygradL_y, mygradR_y);
                    van_albada_limit(mygradL_z, mygradR_z);
                }
                double qL = qL0 + dL.x*mygradL_x + dL.y*mygradL_y + dL.z*mygradL_z;
                double qR = qR0 + dR.x*mygradR_x + dR.y*mygradR_y + dR.z*mygradR_z;
                if (myConfig.extrema_clipping) {
                    Lft."~tname~" = clip_to_limits(qL, qL0, qR0);
                    Rght."~tname~" = clip_to_limits(qR, qL0, qR0);
                } else {
                    Lft."~tname~" = qL;
                    Rght."~tname~" = qR;
                }
                }
                ";
		return code;
	    }
	    mixin(codeForReconstruction("vel.x", "velx", "vel.refx"));
	    mixin(codeForReconstruction("vel.y", "vely", "vel.refy"));
	    mixin(codeForReconstruction("vel.z", "velz", "vel.refz"));
	    if (myConfig.MHD) {
		mixin(codeForReconstruction("B.x", "Bx", "B.refx"));
		mixin(codeForReconstruction("B.y", "By", "B.refy"));
		mixin(codeForReconstruction("B.z", "Bz", "B.refz"));
		if (myConfig.divergence_cleaning) {
		    mixin(codeForReconstruction("psi", "psi", "psi"));
		}
	    }
	    if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
		mixin(codeForReconstruction("tke", "tke", "tke"));
		mixin(codeForReconstruction("omega", "omega", "omega"));
	    }
	    if (nsp > 1) {
		// Multiple species.
		foreach (isp; 0 .. nsp) {
		    mixin(codeForReconstructionOfArray("gas.massf[isp]", "massf", "gas.massf[isp]", "isp"));
		}
		try {
		    scale_mass_fractions(Lft.gas.massf); 
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.gas.massf[] = IFace.left_cell.fs.gas.massf[];
		}
		try {
		    scale_mass_fractions(Rght.gas.massf);
		} catch (Exception e) {
		    writeln(e.msg);
		    Rght.gas.massf[] = IFace.right_cell.fs.gas.massf[];
		}
	    } else {
		// Only one possible mass-fraction value for a single species.
		Lft.gas.massf[0] = 1.0;
		Rght.gas.massf[0] = 1.0;
	    }
	    // Interpolate on two of the thermodynamic quantities, 
	    // and fill in the rest based on an EOS call. 
	    // If an EOS call fails, fall back to just copying cell-centre data.
	    // This does presume that the cell-centre data is valid. 
	    string codeForThermoUpdate(string funname)
	    {
		string code = "
		try {
		    gmodel.update_thermo_from_"~funname~"(Lft.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Lft.copy_values_from(IFace.left_cell.fs);
		}
		try {
		    gmodel.update_thermo_from_"~funname~"(Rght.gas);
		} catch (Exception e) {
		    writeln(e.msg);
		    Rght.copy_values_from(IFace.right_cell.fs);
		}
                ";
		return code;
	    }
	    final switch (myConfig.thermo_interpolator) {
	    case InterpolateOption.pt: 
		mixin(codeForReconstruction("gas.p", "p", "gas.p"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstructionOfArray("gas.T[imode]", "T", "gas.T[imode]", "imode"));
		}
		mixin(codeForThermoUpdate("pT"));
		break;
	    case InterpolateOption.rhoe:
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstructionOfArray("gas.e[imode]", "e", "gas.e[imode]", "imode"));
		}
		mixin(codeForThermoUpdate("rhoe"));
		break;
	    case InterpolateOption.rhop:
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho"));
		mixin(codeForReconstruction("gas.p", "p", "gas.p"));
		mixin(codeForThermoUpdate("rhop"));
		break;
	    case InterpolateOption.rhot: 
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstructionOfArray("gas.T[imode]", "T", "gas.T[imode]", "imode"));
		}
		mixin(codeForThermoUpdate("rhoT"));
		break;
	    } // end switch thermo_interpolator
	} // end of high-order reconstruction
    } // end interp_both()
    
} // end class LsqInterpolator

