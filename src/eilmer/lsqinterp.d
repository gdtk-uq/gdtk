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

immutable size_t cloud_nmax = 12;

class LSQInterpWorkspace {
public:
    // A place to hold the intermediate results for computing
    // the least-squares model as a weighted sum of the flow data.
    double[cloud_nmax] wx, wy, wz; 

    this()
    {
	// do nothing
    }

    this(const LSQInterpWorkspace other)
    {
	wx[] = other.wx[]; wy[] = other.wy[]; wz[] = other.wz[];
    }
    
    void assemble_and_invert_normal_matrix(FVCell[] cell_cloud, int dimensions, size_t gtl)
    {
	auto np = cell_cloud.length;
	assert(np <= cloud_nmax, "Too many points in cloud.");
	double[cloud_nmax] dx, dy, dz;
	foreach (i; 1 .. np) {
	    dx[i] = cell_cloud[i].pos[gtl].x - cell_cloud[0].pos[gtl].x;
	    dy[i] = cell_cloud[i].pos[gtl].y - cell_cloud[0].pos[gtl].y;
	    if (dimensions == 3) {
		dz[i] = cell_cloud[i].pos[gtl].z - cell_cloud[0].pos[gtl].z;
	    } else {
		dz[i] = 0.0;
	    }
	}	    
	// Prepare the normal matrix for the cloud and invert it.
	if (dimensions == 3) {
	    double[6][3] xTx; // normal matrix, augmented to give 6 entries per row
	    double xx = 0.0; double xy = 0.0; double xz = 0.0;
	    double yy = 0.0; double yz = 0.0;
	    double zz = 0.0;
	    foreach (i; 1 .. np) {
		xx += dx[i]*dx[i]; xy += dx[i]*dy[i]; xz += dx[i]*dz[i];
		yy += dy[i]*dy[i]; yz += dy[i]*dz[i]; zz += dz[i]*dz[i];
	    }
	    xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
	    xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
	    xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
	    xTx[0][3] = 1.0; xTx[0][4] = 0.0; xTx[0][5] = 0.0;
	    xTx[1][3] = 0.0; xTx[1][4] = 1.0; xTx[1][5] = 0.0;
	    xTx[2][3] = 0.0; xTx[2][4] = 0.0; xTx[2][5] = 1.0;
	    if (0 != computeInverse!(3,3)(xTx)) {
		throw new FlowSolverException("Failed to invert LSQ normal matrix");
		// Assume that the rows are linearly dependent 
		// because the sample points are colinear.
		// Maybe we could proceed by working as a single-dimensional interpolation.
	    }
	    // Prepare final weights for later use in the reconstruction phase.
	    foreach (i; 1 .. np) {
		wx[i] = xTx[0][3]*dx[i] + xTx[0][4]*dy[i] + xTx[0][5]*dz[i];
		wy[i] = xTx[1][3]*dx[i] + xTx[1][4]*dy[i] + xTx[1][5]*dz[i];
		wz[i] = xTx[2][3]*dx[i] + xTx[2][4]*dy[i] + xTx[2][5]*dz[i];
	    }
	} else {
	    // dimensions == 2
	    double[4][2] xTx; // normal matrix, augmented to give 4 entries per row
	    double xx = 0.0; double xy = 0.0; double yy = 0.0;
	    foreach (i; 1 .. np) {
		xx += dx[i]*dx[i]; xy += dx[i]*dy[i]; yy += dy[i]*dy[i];
	    }
	    xTx[0][0] =  xx; xTx[0][1] =  xy;
	    xTx[1][0] =  xy; xTx[1][1] =  yy;
	    xTx[0][2] = 1.0; xTx[0][3] = 0.0;
	    xTx[1][2] = 0.0; xTx[1][3] = 1.0;
	    if (0 != computeInverse!(2,2)(xTx)) {
		throw new FlowSolverException("Failed to invert LSQ normal matrix");
		// Assume that the rows are linearly dependent 
		// because the sample points are colinear.
		// Maybe we could proceed by working as a single-dimensional interpolation.
	    }
	    // Prepare final weights for later use in the reconstruction phase.
	    foreach (i; 1 .. np) {
		wx[i] = xTx[0][2]*dx[i] + xTx[0][3]*dy[i];
		wy[i] = xTx[1][2]*dx[i] + xTx[1][3]*dy[i];
	    }
	}
    } // end assemble_and_invert_normal_matrix()
} // end class LSQInterpWorkspace

class LSQInterpGradients {
    // These are the quantities that will be interpolated from cell centres 
    // to the Left and Right sides of interfaces.
    // We need to hold onto their gradients within cells.
public:
    double[3] velx, vely, velz;
    double[3] Bx, By, Bz, psi;
    double[3] tke, omega;
    double[3][] massf;
    double[3] rho, p;
    double[3][] T, e;

    double velxPhi, velyPhi, velzPhi;
    double BxPhi, ByPhi, BzPhi, psiPhi;
    double tkePhi, omegaPhi;
    double[] massfPhi;
    double rhoPhi, pPhi;
    double[] TPhi, ePhi;

    double velxMax, velyMax, velzMax;
    double BxMax, ByMax, BzMax, psiMax;
    double tkeMax, omegaMax;
    double[] massfMax;
    double rhoMax, pMax;
    double[] TMax, eMax;
    double velxMin, velyMin, velzMin;
    double BxMin, ByMin, BzMin, psiMin;
    double tkeMin, omegaMin;
    double[] massfMin;
    double rhoMin, pMin;
    double[] TMin, eMin;

    this(size_t nsp, size_t nmodes)
    {
	massf.length = nsp;
	T.length = nmodes;
	e.length = nmodes;

	massfPhi.length = nsp;
	TPhi.length = nmodes;
	ePhi.length = nmodes;

	massfMax.length = nsp;
	TMax.length = nmodes;
	eMax.length = nmodes;

	massfMin.length = nsp;
	TMin.length = nmodes;
	eMin.length = nmodes;
    }

    this(ref const LSQInterpGradients other)
    {
	velx[] = other.velx[]; vely[] = other.vely[]; velz[] = other.velz[];
	Bx[] = other.Bx[]; By[] = other.By[]; Bz[] = other.Bz[]; psi[] = other.psi[];
	tke[] = other.tke[]; omega[] = other.omega[];
	massf.length = other.massf.length;
	foreach(i; 0 .. other.massf.length) { massf[i][] = other.massf[i][]; }
	rho[] = other.rho[]; p[] = other.p[];
	T.length = other.T.length; e.length = other.e.length;
	foreach(i; 0 .. other.T.length) { T[i][] = other.T[i][]; e[i][] = other.e[i][]; }

	velxPhi = other.velxPhi; velyPhi = other.velyPhi; velzPhi = other.velzPhi;
	BxPhi = other.BxPhi; ByPhi = other.ByPhi; BzPhi = other.BzPhi; psiPhi = other.psiPhi;
	tkePhi = other.tkePhi; omegaPhi = other.omegaPhi;
	massfPhi.length = other.massfPhi.length;
	foreach(i; 0 .. other.massf.length) { massfPhi[i] = other.massfPhi[i]; }
	rhoPhi = other.rhoPhi; pPhi = other.pPhi;
	TPhi.length = other.TPhi.length; ePhi.length = other.ePhi.length;
	foreach(i; 0 .. other.TPhi.length) { TPhi[i] = other.TPhi[i]; ePhi[i] = other.ePhi[i]; }

	velxMax = other.velxMax; velyMax = other.velyMax; velzMax = other.velzMax;
	BxMax = other.BxMax; ByMax = other.ByMax; BzMax = other.BzMax; psiMax = other.psiMax;
	tkeMax = other.tkeMax; omegaMax = other.omegaMax;
	massfMax.length = other.massfMax.length;
	foreach(i; 0 .. other.massf.length) { massfMax[i] = other.massfMax[i]; }
	rhoMax = other.rhoMax; pMax = other.pMax;
	TMax.length = other.TMax.length; eMax.length = other.eMax.length;
	foreach(i; 0 .. other.TMax.length) { TMax[i] = other.TMax[i]; eMax[i] = other.eMax[i]; }

	velxMin = other.velxMin; velyMin = other.velyMin; velzMin = other.velzMin;
	BxMin = other.BxMin; ByMin = other.ByMin; BzMin = other.BzMin; psiMin = other.psiMin;
	tkeMin = other.tkeMin; omegaMin = other.omegaMin;
	massfMin.length = other.massfMin.length;
	foreach(i; 0 .. other.massf.length) { massfMin[i] = other.massfMin[i]; }
	rhoMin = other.rhoMin; pMin = other.pMin;
	TMin.length = other.TMin.length; eMin.length = other.eMin.length;
	foreach(i; 0 .. other.TMin.length) { TMin[i] = other.TMin[i]; eMin[i] = other.eMin[i]; }
    }

    void barth_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
	size_t dimensions = myConfig.dimensions;
	double a, b, U, phi;
	// The following function to be used at compile time.
	string codeForLimits(string qname, string gname, string limFactorname, string qMaxname, string qMinname)
	{
	    string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (abs("~gname~"[0]) > 0.0 || abs("~gname~"[1]) > 0.0 || abs("~gname~"[2]) > 0.0) {
            foreach (i, f; cell_cloud[0].iface) {
                double dx = f.pos.x - cell_cloud[0].pos[0].x; 
                double dy = f.pos.y - cell_cloud[0].pos[0].y; 
                double dz = f.pos.z - cell_cloud[0].pos[0].z;
                double dxFace = dx*f.n.x + dy*f.n.y + dz*f.n.z;
                double dyFace = dx*f.t1.x + dy*f.t1.y + dz*f.t1.z;
                double dzFace = dx*f.t2.x + dy*f.t2.y + dz*f.t2.z;
                b = "~gname~"[0] * dxFace + "~gname~"[1] * dyFace;
		if (myConfig.dimensions == 3) b += "~gname~"[2] * dzFace;
		if (b > 0.0) {
		    a = "~qMaxname~" - U;
		    phi = min(phi, a/b);
		}
		else if (b < 0.0) {
		    a = "~qMinname~" - U;
		    phi = min(phi, a/b);
		}
	    }
            }
            "~limFactorname~" = phi;
            }
            ";
	    return code;
	}
	// x-velocity
	mixin(codeForLimits("vel.x", "velx", "velxPhi", "velxMax", "velxMin"));
	mixin(codeForLimits("vel.y", "vely", "velyPhi", "velyMax", "velyMin"));
	mixin(codeForLimits("vel.z", "velz", "velzPhi", "velzMax", "velzMin"));
	if (myConfig.MHD) {
	    mixin(codeForLimits("B.x", "Bx", "BxPhi", "BxMax", "BxMin"));
	    mixin(codeForLimits("B.y", "By", "ByPhi", "ByMax", "ByMin"));
	    mixin(codeForLimits("B.z", "Bz", "BzPhi", "BzMax", "BzMin"));
	    if (myConfig.divergence_cleaning) {
		mixin(codeForLimits("psi", "psi", "psiPhi", "psiMax", "psiMin"));
	    }
	}
	if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
	    mixin(codeForLimits("tke", "tke", "tkePhi", "tkeMax", "tkeMin"));
	    mixin(codeForLimits("omega", "omega", "omegaPhi", "omegaMax", "omegaMin"));
	}
	auto nsp = myConfig.gmodel.n_species;
	if (nsp > 1) {
	    // Multiple species.
	    foreach (isp; 0 .. nsp) {
		mixin(codeForLimits("gas.massf[isp]", "massf[isp]", "massfPhi[isp]", "massfMax[isp]", "massfMin[isp]"));
	    }
	} else {
	    // Only one possible gradient value for a single species.
	    massf[0][0] = 0.0; massf[0][1] = 0.0; massf[0][2] = 0.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	auto nmodes = myConfig.gmodel.n_modes;
	final switch (myConfig.thermo_interpolator) {
	case InterpolateOption.pt: 
	    mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForLimits("gas.T[imode]", "T[imode]", "TPhi[imode]", "TMax[imode]", "TMin[imode]"));
	    }
	    break;
	case InterpolateOption.rhoe:
	    mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForLimits("gas.e[imode]", "e[imode]", "ePhi[imode]", "eMax[imode]", "eMin[imode]"));
	    }
	    break;
	case InterpolateOption.rhop:
	    mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
	    mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
	    break;
	case InterpolateOption.rhot: 
	    mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForLimits("gas.T[imode]", "T[imode]", "TPhi[imode]", "TMax[imode]", "TMin[imode]"));
	    }
	    break;
	} // end switch thermo_interpolator
    } // end compute_lsq_gradients()

   void venkat_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
   {
	size_t dimensions = myConfig.dimensions;
	double a, b, U, phi, h, denom, numer, s;
	immutable double w = 1.0e-12;
        immutable double K = 100.0;
	if (myConfig.dimensions == 3) h =  cbrt(cell_cloud[0].iLength * cell_cloud[0].jLength * cell_cloud[0].kLength);  
        else h = sqrt(cell_cloud[0].iLength * cell_cloud[0].jLength);
	double eps = (K*h) * (K*h) * (K*h);
	// The following function to be used at compile time.
	string codeForLimits(string qname, string gname, string limFactorname, string qMaxname, string qMinname)
	{
	    string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (abs("~gname~"[0]) > 0.0 || abs("~gname~"[1]) > 0.0 || abs("~gname~"[2]) > 0.0) {
            foreach (i, f; cell_cloud[0].iface) {
                double dx = f.pos.x - cell_cloud[0].pos[0].x; 
                double dy = f.pos.y - cell_cloud[0].pos[0].y; 
                double dz = f.pos.z - cell_cloud[0].pos[0].z;
                double dxFace = dx*f.n.x + dy*f.n.y + dz*f.n.z;
                double dyFace = dx*f.t1.x + dy*f.t1.y + dz*f.t1.z;
                double dzFace = dx*f.t2.x + dy*f.t2.y + dz*f.t2.z;
                b = "~gname~"[0] * dxFace + "~gname~"[1] * dyFace;
		if (myConfig.dimensions == 3) b += "~gname~"[2] * dzFace;
		b = sgn(b) * (fabs(b) + w); 
                if (b > 0.0) a = "~qMaxname~" - U; 
                else if (b < 0.0) a = "~qMinname~" - U; 
                numer = (a*a + eps)*b + 2.0*b*b*a;
                denom = a*a + 2.0*b*b + a*b + eps;
                s = (1.0/b) * (numer/denom);                    
                phi = min(phi, s);
	    }
            }
            "~limFactorname~" = phi;
            }
            ";
	    return code;
	}
	// x-velocity
	mixin(codeForLimits("vel.x", "velx", "velxPhi", "velxMax", "velxMin"));
	mixin(codeForLimits("vel.y", "vely", "velyPhi", "velyMax", "velyMin"));
	mixin(codeForLimits("vel.z", "velz", "velzPhi", "velzMax", "velzMin"));
	if (myConfig.MHD) {
	    mixin(codeForLimits("B.x", "Bx", "BxPhi", "BxMax", "BxMin"));
	    mixin(codeForLimits("B.y", "By", "ByPhi", "ByMax", "ByMin"));
	    mixin(codeForLimits("B.z", "Bz", "BzPhi", "BzMax", "BzMin"));
	    if (myConfig.divergence_cleaning) {
		mixin(codeForLimits("psi", "psi", "psiPhi", "psiMax", "psiMin"));
	    }
	}
	if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
	    mixin(codeForLimits("tke", "tke", "tkePhi", "tkeMax", "tkeMin"));
	    mixin(codeForLimits("omega", "omega", "omegaPhi", "omegaMax", "omegaMin"));
	}
	auto nsp = myConfig.gmodel.n_species;
	if (nsp > 1) {
	    // Multiple species.
	    foreach (isp; 0 .. nsp) {
		mixin(codeForLimits("gas.massf[isp]", "massf[isp]", "massfPhi[isp]", "massfMax[isp]", "massfMin[isp]"));
	    }
	} else {
	    // Only one possible gradient value for a single species.
	    massf[0][0] = 0.0; massf[0][1] = 0.0; massf[0][2] = 0.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	auto nmodes = myConfig.gmodel.n_modes;
	final switch (myConfig.thermo_interpolator) {
	case InterpolateOption.pt: 
	    mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForLimits("gas.T[imode]", "T[imode]", "TPhi[imode]", "TMax[imode]", "TMin[imode]"));
	    }
	    break;
	case InterpolateOption.rhoe:
	    mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForLimits("gas.e[imode]", "e[imode]", "ePhi[imode]", "eMax[imode]", "eMin[imode]"));
	    }
	    break;
	case InterpolateOption.rhop:
	    mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
	    mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
	    break;
	case InterpolateOption.rhot: 
	    mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForLimits("gas.T[imode]", "T[imode]", "TPhi[imode]", "TMax[imode]", "TMin[imode]"));
	    }
	    break;
	} // end switch thermo_interpolator
    } // end compute_lsq_gradients()
    
    void compute_lsq_values(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
	size_t dimensions = myConfig.dimensions;
	auto np = cell_cloud.length;
	// The following function to be used at compile time.
	string codeForGradients(string qname, string gname, string qMaxname, string qMinname)
	{
	    string code = "{
                double q0 = cell_cloud[0].fs."~qname~";
                "~qMaxname~" = q0;
                "~qMinname~" = q0;
                "~gname~"[0] = 0.0; "~gname~"[1] = 0.0; "~gname~"[2] = 0.0;
                foreach (i; 1 .. np) {
                    double dq = cell_cloud[i].fs."~qname~" - q0;
                    "~gname~"[0] += ws.wx[i] * dq;
                    "~gname~"[1] += ws.wy[i] * dq;
                    if (dimensions == 3) { "~gname~"[2] += ws.wz[i] * dq; }
                    "~qMaxname~" = max("~qMaxname~", cell_cloud[i].fs."~qname~");
                    "~qMinname~" = min("~qMinname~", cell_cloud[i].fs."~qname~");
                }
                }
                ";
	    return code;
	}
	// x-velocity
	mixin(codeForGradients("vel.x", "velx", "velxMax", "velxMin"));
	mixin(codeForGradients("vel.y", "vely", "velyMax", "velyMin"));
	mixin(codeForGradients("vel.z", "velz", "velzMax", "velzMin"));
	if (myConfig.MHD) {
	    mixin(codeForGradients("B.x", "Bx", "BxMax", "BxMin"));
	    mixin(codeForGradients("B.y", "By", "ByMax", "ByMin"));
	    mixin(codeForGradients("B.z", "Bz", "BzMax", "BzMin"));
	    if (myConfig.divergence_cleaning) {
		mixin(codeForGradients("psi", "psi", "psiMax", "psiMin"));
	    }
	}
	if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
	    mixin(codeForGradients("tke", "tke", "tkeMax", "tkeMin"));
	    mixin(codeForGradients("omega", "omega", "omegaMax", "omegaMin"));
	}
	auto nsp = myConfig.gmodel.n_species;
	if (nsp > 1) {
	    // Multiple species.
	    foreach (isp; 0 .. nsp) {
		mixin(codeForGradients("gas.massf[isp]", "massf[isp]", "massfMax[isp]", "massfMin[isp]"));
	    }
	} else {
	    // Only one possible gradient value for a single species.
	    massf[0][0] = 0.0; massf[0][1] = 0.0; massf[0][2] = 0.0;
	}
	// Interpolate on two of the thermodynamic quantities, 
	// and fill in the rest based on an EOS call. 
	auto nmodes = myConfig.gmodel.n_modes;
	final switch (myConfig.thermo_interpolator) {
	case InterpolateOption.pt: 
	    mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradients("gas.T[imode]", "T[imode]", "TMax[imode]", "TMin[imode]"));
	    }
	    break;
	case InterpolateOption.rhoe:
	    mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradients("gas.e[imode]", "e[imode]", "eMax[imode]", "eMin[imode]"));
	    }
	    break;
	case InterpolateOption.rhop:
	    mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
	    mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
	    break;
	case InterpolateOption.rhot: 
	    mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
	    foreach (imode; 0 .. nmodes) {
		mixin(codeForGradients("gas.T[imode]", "T[imode]", "TMax[imode]", "TMin[imode]"));
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
	    double[3] mygradL, mygradR;
	    //
	    // Always reconstruct in the global frame of reference -- for now
	    //
	    // x-velocity
	    string codeForReconstruction(string qname, string gname, string tname, string lname)
	    {
		string code = "{
                double qL0 = cL0.fs."~qname~";
                double qMinL = qL0;
                double qMaxL = qL0;
                if (wsL) { // an active cell will have a workspace
                    mygradL[0] = cL0.gradients."~gname~"[0];
                    mygradL[1] = cL0.gradients."~gname~"[1];
                    mygradL[2] = cL0.gradients."~gname~"[2];
                } else {
                    // Guess that there is good data over on the other side.
                    mygradL[0] = cR0.gradients."~gname~"[0];
                    mygradL[1] = cR0.gradients."~gname~"[1];
                    mygradL[2] = cR0.gradients."~gname~"[2];
                }
                double qR0 = cR0.fs."~qname~";
                double qMinR = qR0;
                double qMaxR = qR0;
                if (wsR) { // an active cell will have a workspace
                    mygradR[0] = cR0.gradients."~gname~"[0];
                    mygradR[1] = cR0.gradients."~gname~"[1];
                    mygradR[2] = cR0.gradients."~gname~"[2];
                } else {
                    // Guess that there is good data over on the other side.
                    mygradR[0] = cL0.gradients."~gname~"[0];
                    mygradR[1] = cL0.gradients."~gname~"[1];
                    mygradR[2] = cL0.gradients."~gname~"[2];
                }
                if (myConfig.apply_limiter) {
                    final switch (myConfig.unstructured_limiter) {
                    case UnstructuredLimiter.van_albada:
                        van_albada_limit(mygradL[0], mygradR[0]);
                        van_albada_limit(mygradL[1], mygradR[1]);
                        van_albada_limit(mygradL[2], mygradR[2]);
                        break;
                    case UnstructuredLimiter.min_mod:
                        min_mod_limit(mygradL[0], mygradR[0]);
                        min_mod_limit(mygradL[1], mygradR[1]);
                        min_mod_limit(mygradL[2], mygradR[2]);
                        break;
                    case UnstructuredLimiter.barth:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat:
                        if (wsL) {
                            mygradL[0] *=  cL0.gradients."~lname~"; mygradL[1] *=  cL0.gradients."~lname~"; mygradL[2] *=  cL0.gradients."~lname~";
                        } else {                    
                            mygradL[0] *=  cR0.gradients."~lname~"; mygradL[1] *=  cR0.gradients."~lname~"; mygradL[2] *=  cR0.gradients."~lname~";
                        }
                        if (wsR) {
                            mygradR[0] *=  cR0.gradients."~lname~"; mygradR[1] *=  cR0.gradients."~lname~"; mygradR[2] *=  cR0.gradients."~lname~";
                        } else {                
                            mygradR[0] *=  cL0.gradients."~lname~"; mygradR[1] *=  cL0.gradients."~lname~"; mygradR[2] *=  cL0.gradients."~lname~";
                        }
                        break;
                    }
                }
                double qL = qL0 + dL.x*mygradL[0] + dL.y*mygradL[1];
                double qR = qR0 + dR.x*mygradR[0] + dR.y*mygradR[1];
                if (myConfig.dimensions == 3) {
                    qL += dL.z*mygradL[2];
                    qR += dR.z*mygradR[2];
                }
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
	    mixin(codeForReconstruction("vel.x", "velx", "vel.refx", "velxPhi"));
	    mixin(codeForReconstruction("vel.y", "vely", "vel.refy", "velyPhi"));
	    mixin(codeForReconstruction("vel.z", "velz", "vel.refz", "velzPhi"));
	    if (myConfig.MHD) {
		mixin(codeForReconstruction("B.x", "Bx", "B.refx", "BxPhi"));
		mixin(codeForReconstruction("B.y", "By", "B.refy", "ByPhi"));
		mixin(codeForReconstruction("B.z", "Bz", "B.refz", "BxPhi"));
		if (myConfig.divergence_cleaning) {
		    mixin(codeForReconstruction("psi", "psi", "psi", "psiPhi"));
		}
	    }
	    if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
		mixin(codeForReconstruction("tke", "tke", "tke", "tkePhi"));
		mixin(codeForReconstruction("omega", "omega", "omega", "omegaPhi"));
	    }
	    if (nsp > 1) {
		// Multiple species.
		foreach (isp; 0 .. nsp) {
		    mixin(codeForReconstruction("gas.massf[isp]", "massf[isp]", "gas.massf[isp]", "massfPhi[isp]"));
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
		mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "T[imode]", "gas.T[imode]", "TPhi[imode]"));
		}
		mixin(codeForThermoUpdate("pT"));
		break;
	    case InterpolateOption.rhoe:
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.e[imode]", "e[imode]", "gas.e[imode]", "ePhi[imode]"));
		}
		mixin(codeForThermoUpdate("rhoe"));
		break;
	    case InterpolateOption.rhop:
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
		mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
		mixin(codeForThermoUpdate("rhop"));
		break;
	    case InterpolateOption.rhot: 
		mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
		foreach (imode; 0 .. nmodes) {
		    mixin(codeForReconstruction("gas.T[imode]", "T[imode]", "gas.T[imode]", "TPhi[imode]"));
		}
		mixin(codeForThermoUpdate("rhoT"));
		break;
	    } // end switch thermo_interpolator
	} // end of high-order reconstruction
    } // end interp_both()
    
} // end class LsqInterpolator

