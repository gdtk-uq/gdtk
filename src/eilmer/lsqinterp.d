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
    Vector3 velx, vely, velz;
    Vector3 Bx, By, Bz, psi;
    Vector3 tke, omega;
    Vector3[] massf;
    Vector3 rho, p;
    Vector3[] T, e;

    this(size_t nsp, size_t nmodes)
    {
	massf.length = nsp;
	T.length = nmodes;
	e.length = nmodes;
    }

    this(ref const LSQInterpGradients other)
    {
	velx.refx = other.velx.x; velx.refy = other.velx.y; velx.refz = other.velx.z;
	vely.refx = other.vely.x; vely.refy = other.vely.y; vely.refz = other.vely.z;
	velz.refx = other.velz.x; velz.refy = other.velz.y; velz.refz = other.velz.z;
	Bx.refx = other.Bx.x; Bx.refy = other.Bx.y; Bx.refz = other.Bx.z;
	By.refx = other.By.x; By.refy = other.By.y; By.refz = other.By.z;
	Bz.refx = other.Bz.x; Bz.refy = other.Bz.y; Bz.refz = other.Bz.z;
	psi.refx = other.psi.x; psi.refy = other.psi.y; psi.refz = other.psi.z;
	tke.refx = other.tke.x; tke.refy = other.tke.y; tke.refz = other.tke.z;
	omega.refx = other.omega.x; omega.refy = other.omega.y; omega.refz = other.omega.z;
	foreach(ref m; other.massf) { massf ~= Vector3(m); }
	rho.refx = other.rho.x; rho.refy = other.rho.y; rho.refz = other.rho.z;
	p.refx = other.p.x; p.refy = other.p.y; p.refz = other.p.z;
	foreach(ref Ti; other.T) { T ~= Vector3(Ti); }
	foreach(ref ei; other.e) { massf ~= Vector3(ei); }
    }

    void compute_lsq_values(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
	double[3] rhs, gradients;
	auto np = cell_cloud.length;
	// The following function to be used at compile time.
	string codeForGradients(string qname, string gname)
	{
	    string code = "{
                double q0 = cell_cloud[0].fs."~qname~";
                foreach (j; 0 .. 3) { rhs[j] = 0.0; }
                foreach (i; 1 .. np) {
                    double dq = cell_cloud[i].fs."~qname~" - q0;
                    rhs[0] += ws.dx[i]*dq; rhs[1] += ws.dy[i]*dq; rhs[2] += ws.dz[i]*dq;
                }
	        solveWithInverse!(3,3)(ws.xTx, rhs, gradients);
                "~gname~".refx = gradients[0];
                "~gname~".refy = gradients[1];
                "~gname~".refz = gradients[2];
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
		mixin(codeForGradients("gas.massf[isp]", "massf[isp]"));
	    }
	} else {
	    // Only one possible gradient value for a single species.
	    massf[0].refx = 0.0; massf[0].refy = 0.0; massf[0].refz = 0.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	auto nmodes = myConfig.gmodel.n_modes;
	final switch (myConfig.thermo_interpolator) {
	case InterpolateOption.pt: 
	    mixin(codeForGradients("gas.p", "p"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradients("gas.T[imode]", "T[imode]"));
	    }
	    break;
	case InterpolateOption.rhoe:
	    mixin(codeForGradients("gas.rho", "rho"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradients("gas.e[imode]", "e[imode]"));
	    }
	    break;
	case InterpolateOption.rhop:
	    mixin(codeForGradients("gas.rho", "rho"));
	    mixin(codeForGradients("gas.p", "p"));
	    break;
	case InterpolateOption.rhot: 
	    mixin(codeForGradients("gas.rho", "rho"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradients("gas.T[imode]", "T[imode]"));
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
	    Vector3 mygradL, mygradR;
	    //
	    // Always reconstruct in the global frame of reference -- for now
	    //
	    // x-velocity
	    string codeForReconstruction(string qname, string gname, string tname)
	    {
		string code = "{
                double qL0 = cL0.fs."~qname~";
                double qMinL = qL0;
                double qMaxL = qL0;
                if (wsL) { // an active cell will have a workspace
                    foreach (i; 1 .. cL0.cell_cloud.length) {
                        qMinL = min(qMinL, cL0.cell_cloud[i].fs."~qname~");
                        qMaxL = max(qMaxL, cL0.cell_cloud[i].fs."~qname~");
                    }
                    mygradL.refx = cL0.gradients."~gname~".x;
                    mygradL.refy = cL0.gradients."~gname~".y;
                    mygradL.refz = cL0.gradients."~gname~".z;
                } else {
                    mygradL.refx = 0.0;
                    mygradL.refy = 0.0;
                    mygradL.refz = 0.0;
                }
                double qR0 = cR0.fs."~qname~";
                double qMinR = qR0;
                double qMaxR = qR0;
                if (wsR) { // an active cell will have a workspace
                    foreach (i; 1 .. cR0.cell_cloud.length) {
                        qMinR = min(qMinR, cR0.cell_cloud[i].fs."~qname~");
                        qMaxR = max(qMaxR, cR0.cell_cloud[i].fs."~qname~");
                    }
                    mygradR.refx = cR0.gradients."~gname~".x;
                    mygradR.refy = cR0.gradients."~gname~".y;
                    mygradR.refz = cR0.gradients."~gname~".z;
                } else {
                    mygradR.refx = 0.0;
                    mygradR.refy = 0.0;
                    mygradR.refz = 0.0;
                }
                if (myConfig.apply_limiter) {
                    // venkatakrishan_limit(wsL.gradients."~gname~", cL0, qL0, qMinL, qMaxL, cL0.iLength, cL0.jLength, cL0.kLength);
                    // venkatakrishan_limit(wsR.gradients."~gname~", cR0, qR0, qMinR, qMaxR, cR0.iLength, cR0.jLength, cR0.kLength);
                    van_albada_limit(mygradL.refx, mygradR.refx);
                    van_albada_limit(mygradL.refy, mygradR.refy);
                    van_albada_limit(mygradL.refz, mygradR.refz);
                }
                double qL = qL0 + dL.dot(mygradL);
                double qR = qR0 + dR.dot(mygradR);
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
		    mixin(codeForReconstruction("gas.massf[isp]", "massf[isp]", "gas.massf[isp]"));
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
		    mixin(codeForReconstruction("gas.T[imode]", "T[imode]", "gas.T[imode]"));
		}
		mixin(codeForThermoUpdate("pT"));
		break;
	    case InterpolateOption.rhoe:
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.e[imode]", "e[imode]", "gas.e[imode]"));
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
		    mixin(codeForReconstruction("gas.T[imode]", "T[imode]", "gas.T[imode]"));
		}
		mixin(codeForThermoUpdate("rhoT"));
		break;
	    } // end switch thermo_interpolator
	} // end of high-order reconstruction
    } // end interp_both()
    
} // end class LsqInterpolator

