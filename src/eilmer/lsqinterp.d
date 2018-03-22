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

immutable size_t cloud_nmax = 33;
immutable double ESSENTIALLY_ZERO = 1.0e-50;

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
            double very_small_value = 1.0e-16*(normInf!(3,3)(xTx))^^3;
            if (0 != computeInverse!(3,3)(xTx, very_small_value)) {
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
            double very_small_value = 1.0e-16*(normInf!(2,2)(xTx))^^2;
            if (0 != computeInverse!(2,2)(xTx, very_small_value)) {
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
    double[3] T, u;
    double[3][] T_modes, u_modes;

    double velxPhi, velyPhi, velzPhi;
    double BxPhi, ByPhi, BzPhi, psiPhi;
    double tkePhi, omegaPhi;
    double[] massfPhi;
    double rhoPhi, pPhi;
    double TPhi, uPhi;
    double[] T_modesPhi, u_modesPhi;

    double velxMax, velyMax, velzMax;
    double BxMax, ByMax, BzMax, psiMax;
    double tkeMax, omegaMax;
    double[] massfMax;
    double rhoMax, pMax;
    double TMax, uMax;
    double[] T_modesMax, u_modesMax;
    double velxMin, velyMin, velzMin;
    double BxMin, ByMin, BzMin, psiMin;
    double tkeMin, omegaMin;
    double[] massfMin;
    double rhoMin, pMin;
    double TMin, uMin;
    double[] T_modesMin, u_modesMin;

    this(size_t nsp, size_t nmodes)
    {
        massf.length = nsp;
        T_modes.length = nmodes;
        u_modes.length = nmodes;

        massfPhi.length = nsp;
        T_modesPhi.length = nmodes;
        u_modesPhi.length = nmodes;

        massfMax.length = nsp;
        T_modesMax.length = nmodes;
        u_modesMax.length = nmodes;

        massfMin.length = nsp;
        T_modesMin.length = nmodes;
        u_modesMin.length = nmodes;
    }

    this(ref const LSQInterpGradients other)
    {
        this.copy_values_from(other);
    }

    void copy_values_from(ref const LSQInterpGradients other)
    {
        velx[] = other.velx[]; vely[] = other.vely[]; velz[] = other.velz[];
        Bx[] = other.Bx[]; By[] = other.By[]; Bz[] = other.Bz[]; psi[] = other.psi[];
        tke[] = other.tke[]; omega[] = other.omega[];
        massf.length = other.massf.length;
        foreach(i; 0 .. other.massf.length) { massf[i][] = other.massf[i][]; }
        rho[] = other.rho[]; p[] = other.p[];
        T[] = other.T[]; u[] = other.u[];
        T_modes.length = other.T_modes.length; u_modes.length = other.u_modes.length;
        foreach(i; 0 .. other.T_modes.length) {
            T_modes[i][] = other.T_modes[i][]; u_modes[i][] = other.u_modes[i][];
        }

        velxPhi = other.velxPhi; velyPhi = other.velyPhi; velzPhi = other.velzPhi;
        BxPhi = other.BxPhi; ByPhi = other.ByPhi; BzPhi = other.BzPhi; psiPhi = other.psiPhi;
        tkePhi = other.tkePhi; omegaPhi = other.omegaPhi;
        massfPhi.length = other.massfPhi.length;
        foreach(i; 0 .. other.massf.length) { massfPhi[i] = other.massfPhi[i]; }
        rhoPhi = other.rhoPhi; pPhi = other.pPhi;
        TPhi = other.TPhi; uPhi = other.uPhi;
        T_modesPhi.length = other.T_modesPhi.length; u_modesPhi.length = other.u_modesPhi.length;
        foreach(i; 0 .. other.T_modesPhi.length) {
            T_modesPhi[i] = other.T_modesPhi[i]; u_modesPhi[i] = other.u_modesPhi[i];
        }

        velxMax = other.velxMax; velyMax = other.velyMax; velzMax = other.velzMax;
        BxMax = other.BxMax; ByMax = other.ByMax; BzMax = other.BzMax; psiMax = other.psiMax;
        tkeMax = other.tkeMax; omegaMax = other.omegaMax;
        massfMax.length = other.massfMax.length;
        foreach(i; 0 .. other.massf.length) { massfMax[i] = other.massfMax[i]; }
        rhoMax = other.rhoMax; pMax = other.pMax;
        TMax = other.TMax; uMax = other.uMax;
        T_modesMax.length = other.T_modesMax.length; u_modesMax.length = other.u_modesMax.length;
        foreach(i; 0 .. other.T_modesMax.length) {
            T_modesMax[i] = other.T_modesMax[i]; u_modesMax[i] = other.u_modesMax[i];
        }

        velxMin = other.velxMin; velyMin = other.velyMin; velzMin = other.velzMin;
        BxMin = other.BxMin; ByMin = other.ByMin; BzMin = other.BzMin; psiMin = other.psiMin;
        tkeMin = other.tkeMin; omegaMin = other.omegaMin;
        massfMin.length = other.massfMin.length;
        foreach(i; 0 .. other.massf.length) { massfMin[i] = other.massfMin[i]; }
        rhoMin = other.rhoMin; pMin = other.pMin;
        TMin = other.TMin; uMin = other.uMin;
        T_modesMin.length = other.T_modesMin.length; u_modesMin.length = other.u_modesMin.length;
        foreach(i; 0 .. other.T_modesMin.length) {
            T_modesMin[i] = other.T_modesMin[i]; u_modesMin[i] = other.u_modesMin[i];
        }
    } // end copy_values_from()
    
    void barth_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
        size_t dimensions = myConfig.dimensions;
        double a, b, U, phi;
        immutable double w = 1.0e-12;
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname, string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (abs("~gname~"[0]) > ESSENTIALLY_ZERO || abs("~gname~"[1]) > ESSENTIALLY_ZERO || abs("~gname~"[2]) > ESSENTIALLY_ZERO) {
            foreach (i, f; cell_cloud[0].iface) {
                double dx = f.pos.x - cell_cloud[0].pos[0].x; 
                double dy = f.pos.y - cell_cloud[0].pos[0].y; 
                double dz = f.pos.z - cell_cloud[0].pos[0].z;
                b = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) b += "~gname~"[2] * dz; 
                b = sgn(b) * (fabs(b) + w);
                if (b >  ESSENTIALLY_ZERO) {
                    a = "~qMaxname~" - U;
                    phi = min(phi, a/b);
                }
                else if (b <  ESSENTIALLY_ZERO) {
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
                mixin(codeForLimits("gas.massf[isp]", "massf[isp]", "massfPhi[isp]",
                                    "massfMax[isp]", "massfMin[isp]"));
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
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                    "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                    "u_modesMax[imode]", "u_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                    "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    void mlp_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
        // The implementation of the MLP limiter follows the implementation in VULCAN
        // i.e. it uses the MLP approach, and the Van Leer limiting function
        // as outlined in the NASA publication on VULCAN's unstructured solver [White et al 2017]
        size_t dimensions = myConfig.dimensions;
        double eps, a, b, U, phi, h, denom, numer, s;
        immutable double w = 1.0e-12;
        if (myConfig.dimensions == 3) h =  cbrt(cell_cloud[0].volume[0]);  
        else h = sqrt(cell_cloud[0].volume[0]);
        eps = 1.0e-12;
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname, string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            foreach (i, vtx; cell_cloud[0].vtx) {
                double dx = vtx.pos[0].x - cell_cloud[0].pos[0].x; 
                double dy = vtx.pos[0].y - cell_cloud[0].pos[0].y; 
                double dz = vtx.pos[0].z - cell_cloud[0].pos[0].z;
                b = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) b += "~gname~"[2] * dz;
                b = sgn(b) * (fabs(b) + w); 
                if (b > 0.0) a = 0.5*(vtx.gradients."~qMaxname~" - U); 
                else if (b < 0.0) a = 0.5*(vtx.gradients."~qMinname~" - U); 
                numer = b*abs(a) + a*abs(b);
                denom = abs(a) + abs(b) + eps;
                s = (1.0/b) * (numer/denom);                    
                phi = min(phi, s);
                if (b == 0.0) phi = 1.0;
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
                mixin(codeForLimits("gas.massf[isp]", "massf[isp]", "massfPhi[isp]",
                                    "massfMax[isp]", "massfMin[isp]"));
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
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                    "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                    "u_modesMax[imode]", "u_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                    "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    void store_max_min_values_for_mlp_limiter(FVCell[] cell_cloud, ref LocalConfig myConfig)
    {
        // the MLP limiter stores the maximum, and minimum flowstate values that surround a node, at the node
        size_t dimensions = myConfig.dimensions;
        auto np = cell_cloud.length;
        // The following function to be used at compile time.
        string codeForGradients(string qname, string gname, string qMaxname, string qMinname)
        {
            string code = "{
                double q0 = cell_cloud[0].fs."~qname~";
                "~qMaxname~" = q0;
                "~qMinname~" = q0;
                foreach (i; 1 .. np) {
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
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                       "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.u", "u", "uMax", "uMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.u_modes[imode]", "u_modes[imode]",
                                       "u_modesMax[imode]", "u_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                       "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    
    void venkat_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig, size_t gtl=0)
   {
        size_t dimensions = myConfig.dimensions;
        double a, b, U, phi, h, denom, numer, s;
        immutable double w = 1.0e-12;
        immutable double K = 0.3;
        if (myConfig.dimensions == 3) h =  cbrt(cell_cloud[0].volume[gtl]);  
        else h = sqrt(cell_cloud[0].volume[gtl]);
        double eps = (K*h) * (K*h) * (K*h);
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname, string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (abs("~gname~"[0]) > ESSENTIALLY_ZERO || abs("~gname~"[1]) > ESSENTIALLY_ZERO || abs("~gname~"[2]) > ESSENTIALLY_ZERO) {
                foreach (i, f; cell_cloud[0].iface) {
                    double dx = f.pos.x - cell_cloud[0].pos[gtl].x; 
                    double dy = f.pos.y - cell_cloud[0].pos[gtl].y; 
                    double dz = f.pos.z - cell_cloud[0].pos[gtl].z;
                    b = "~gname~"[0] * dx + "~gname~"[1] * dy;
                    if (myConfig.dimensions == 3) b += "~gname~"[2] * dz;
                    b = sgn(b) * (fabs(b) + w); 
                    if ( fabs(b) > ESSENTIALLY_ZERO) {
                        if (b > ESSENTIALLY_ZERO) a = "~qMaxname~" - U; 
                        else if (b < ESSENTIALLY_ZERO) a = "~qMinname~" - U; 
                        numer = (a*a + eps)*b + 2.0*b*b*a;
                        denom = a*a + 2.0*b*b + a*b + eps;
                        s = (1.0/b) * (numer/denom);                    
                    }
                    else s = 1.0;
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
                mixin(codeForLimits("gas.massf[isp]", "massf[isp]", "massfPhi[isp]",
                                    "massfMax[isp]", "massfMin[isp]"));
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
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                    "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                    "u_modesMax[imode]", "u_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                    "T_modesMax[imode]", "T_modesMin[imode]"));
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
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                       "T_modesMax[imode]", "T_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.u", "u", "uMax", "uMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.u_modes[imode]", "u_modes[imode]",
                                       "u_modesMax[imode]", "u_modesMin[imode]"));
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                       "T_modesMax[imode]", "T_modesMin[imode]"));
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

    @nogc int get_interpolation_order()
    {
        return myConfig.interpolation_order;
    }

    @nogc void set_interpolation_order(int order)
    {
        myConfig.interpolation_order = order;
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
    // Interpolate the flow field quantities at the left- and right-side of the interface,
    // given information in both cells attached to this interface.
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
        // for some simulations we would like to have the boundaries to remain 1st order
        if (myConfig.suppress_reconstruction_at_boundaries && IFace.is_on_boundary) return;
        // else apply higher-order interpolation to all faces
        if (myConfig.interpolation_order > 1) {
            // High-order reconstruction for some properties.
            //
            LSQInterpWorkspace wsL = cL0.ws;
            LSQInterpWorkspace wsR = cR0.ws;
            // vector from left-cell-centre to face midpoint
            double dLx = IFace.pos.x - cL0.pos[gtl].x;
            double dLy = IFace.pos.y - cL0.pos[gtl].y;
            double dLz = IFace.pos.z - cL0.pos[gtl].z;
            double dRx = IFace.pos.x - cR0.pos[gtl].x;
            double dRy = IFace.pos.y - cR0.pos[gtl].y;
            double dRz = IFace.pos.z - cR0.pos[gtl].z;
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
                mygradL[0] = cL0.gradients."~gname~"[0];
                mygradL[1] = cL0.gradients."~gname~"[1];
                mygradL[2] = cL0.gradients."~gname~"[2];
                double qR0 = cR0.fs."~qname~";
                double qMinR = qR0;
                double qMaxR = qR0;
                mygradR[0] = cR0.gradients."~gname~"[0];
                mygradR[1] = cR0.gradients."~gname~"[1];
                mygradR[2] = cR0.gradients."~gname~"[2];
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
                    case UnstructuredLimiter.mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.barth:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat:
                        mygradL[0] *=  cL0.gradients."~lname~"; mygradL[1] *=  cL0.gradients."~lname~"; mygradL[2] *=  cL0.gradients."~lname~";
                        mygradR[0] *=  cR0.gradients."~lname~"; mygradR[1] *=  cR0.gradients."~lname~"; mygradR[2] *=  cR0.gradients."~lname~";
                        break;
                    }
                }
                double qL = qL0 + dLx*mygradL[0] + dLy*mygradL[1];
                double qR = qR0 + dRx*mygradR[0] + dRy*mygradR[1];
                if (myConfig.dimensions == 3) {
                    qL += dLz*mygradL[2];
                    qR += dRz*mygradR[2];
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
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                "gas.T_modes[imode]", "T_modesPhi[imode]"));
                }
                mixin(codeForThermoUpdate("pT"));
                break;
            case InterpolateOption.rhou:
                mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                mixin(codeForReconstruction("gas.u", "u", "gas.u", "uPhi"));
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                "gas.u_modes[imode]", "u_modesPhi[imode]"));
                }
                mixin(codeForThermoUpdate("rhou"));
                break;
            case InterpolateOption.rhop:
                mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
                mixin(codeForThermoUpdate("rhop"));
                break;
            case InterpolateOption.rhot: 
                mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                "gas.T_modes[imode]", "T_modesPhi[imode]"));
                }
                mixin(codeForThermoUpdate("rhoT"));
                break;
            } // end switch thermo_interpolator
        } // end of high-order reconstruction
    } // end interp_both()
    
} // end class LsqInterpolator

