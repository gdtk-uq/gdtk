// lsqinterp.d
// Least-squares interpolation/reconstruction of flow field.
//

import std.math;
import std.stdio;
import std.algorithm;
import std.conv;
import nm.complex;
import nm.number;
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
    number[cloud_nmax] wx, wy, wz; 

    this()
    {
        // do nothing
    }

    this(const LSQInterpWorkspace other)
    {
        wx[] = other.wx[]; wy[] = other.wy[]; wz[] = other.wz[];
    }

    @nogc
    void assemble_and_invert_normal_matrix(FVCell[] cell_cloud, int dimensions, size_t gtl)
    {
        auto np = cell_cloud.length;
        assert(np <= cloud_nmax, "Too many points in cloud.");
        number[cloud_nmax] dx, dy, dz;
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
            number[6][3] xTx; // normal matrix, augmented to give 6 entries per row
            number xx = 0.0; number xy = 0.0; number xz = 0.0;
            number yy = 0.0; number yz = 0.0;
            number zz = 0.0;
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
            double very_small_value = 1.0e-16*(normInf!(3,3,6,number)(xTx).re)^^3;
            if (0 != computeInverse!(3,3,6,number)(xTx, very_small_value)) {
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
            number[4][2] xTx; // normal matrix, augmented to give 4 entries per row
            number xx = 0.0; number xy = 0.0; number yy = 0.0;
            foreach (i; 1 .. np) {
                xx += dx[i]*dx[i]; xy += dx[i]*dy[i]; yy += dy[i]*dy[i];
            }
            xTx[0][0] =  xx; xTx[0][1] =  xy;
            xTx[1][0] =  xy; xTx[1][1] =  yy;
            xTx[0][2] = 1.0; xTx[0][3] = 0.0;
            xTx[1][2] = 0.0; xTx[1][3] = 1.0;
            double very_small_value = 1.0e-16*(normInf!(2,2,4,number)(xTx).re)^^2;
            if (0 != computeInverse!(2,2,4,number)(xTx, very_small_value)) {
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
    number[3] velx, vely, velz;
    number velxPhi, velyPhi, velzPhi;
    number velxMax, velyMax, velzMax;
    number velxMin, velyMin, velzMin;
    version(MHD) {
        number[3] Bx, By, Bz, psi;
        number BxPhi, ByPhi, BzPhi, psiPhi;
        number BxMax, ByMax, BzMax, psiMax;
        number BxMin, ByMin, BzMin, psiMin;
    }
    version(komega) {
        number[3] tke, omega;
        number tkePhi, omegaPhi;
        number tkeMax, omegaMax;
        number tkeMin, omegaMin;
    }
    version(multi_species_gas) {
        number[3][] massf;
        number[] massfPhi;
        number[] massfMax;
        number[] massfMin;
    }
    number[3] rho, p;
    number rhoPhi, pPhi;
    number rhoMax, pMax;
    number rhoMin, pMin;
    number[3] T, u;
    number TPhi, uPhi;
    number TMax, uMax;
    number TMin, uMin;
    version(multi_T_gas) {
        number[3][] T_modes, u_modes;
        number[] T_modesPhi, u_modesPhi;
        number[] T_modesMax, u_modesMax;
        number[] T_modesMin, u_modesMin;
    }


    this(size_t nsp, size_t nmodes)
    {
        version(multi_species_gas) {
            massf.length = nsp;
            massfPhi.length = nsp;
            massfMax.length = nsp;
            massfMin.length = nsp;
        }
        version(multi_T_gas) {
            T_modes.length = nmodes;
            u_modes.length = nmodes;
            T_modesPhi.length = nmodes;
            u_modesPhi.length = nmodes;
            T_modesMax.length = nmodes;
            u_modesMax.length = nmodes;
            T_modesMin.length = nmodes;
            u_modesMin.length = nmodes;
        }
    }

    this(ref const(LSQInterpGradients) other)
    {
        this.copy_values_from(other);
    }

    @nogc
    void copy_values_from(ref const(LSQInterpGradients) other)
    {
        velx[] = other.velx[]; vely[] = other.vely[]; velz[] = other.velz[];
        velxPhi = other.velxPhi; velyPhi = other.velyPhi; velzPhi = other.velzPhi;
        velxMax = other.velxMax; velyMax = other.velyMax; velzMax = other.velzMax;
        velxMin = other.velxMin; velyMin = other.velyMin; velzMin = other.velzMin;
        version(MHD) {
            Bx[] = other.Bx[]; By[] = other.By[]; Bz[] = other.Bz[]; psi[] = other.psi[];
            BxPhi = other.BxPhi; ByPhi = other.ByPhi; BzPhi = other.BzPhi; psiPhi = other.psiPhi;
            BxMax = other.BxMax; ByMax = other.ByMax; BzMax = other.BzMax; psiMax = other.psiMax;
            BxMin = other.BxMin; ByMin = other.ByMin; BzMin = other.BzMin; psiMin = other.psiMin;
        }
        version(komega) {
            tke[] = other.tke[]; omega[] = other.omega[];
            tkePhi = other.tkePhi; omegaPhi = other.omegaPhi;
            tkeMax = other.tkeMax; omegaMax = other.omegaMax;
            tkeMin = other.tkeMin; omegaMin = other.omegaMin;
        }
        version(multi_species_gas) {
            assert(massf.length == other.massf.length, "Mismatch in massf length");
            foreach(i; 0 .. other.massf.length) { massf[i][] = other.massf[i][]; }
            assert(massfPhi.length == other.massfPhi.length, "Mismatch in massfPhi length");
            foreach(i; 0 .. other.massf.length) { massfPhi[i] = other.massfPhi[i]; }
            assert(massfMax.length == other.massfMax.length, "Mismatch in massfMax length");
            foreach(i; 0 .. other.massf.length) { massfMax[i] = other.massfMax[i]; }
            assert(massfMin.length == other.massfMin.length, "Mismatch in massfMin length");
            foreach(i; 0 .. other.massf.length) { massfMin[i] = other.massfMin[i]; }
        }
        rho[] = other.rho[]; p[] = other.p[];
        rhoPhi = other.rhoPhi; pPhi = other.pPhi;
        rhoMax = other.rhoMax; pMax = other.pMax;
        rhoMin = other.rhoMin; pMin = other.pMin;
        T[] = other.T[]; u[] = other.u[];
        TPhi = other.TPhi; uPhi = other.uPhi;
        TMax = other.TMax; uMax = other.uMax;
        TMin = other.TMin; uMin = other.uMin;
        version(multi_T_gas) {
            assert(T_modes.length == other.T_modes.length, "Mismatch in T_modes length");
            assert(u_modes.length == other.u_modes.length, "Mismatch in u_modes length");
            foreach(i; 0 .. other.T_modes.length) {
                T_modes[i][] = other.T_modes[i][]; u_modes[i][] = other.u_modes[i][];
            }
            assert(T_modesPhi.length == other.T_modesPhi.length, "Mismatch in T_modesPhi length");
            assert(u_modesPhi.length == other.u_modesPhi.length, "Mismatch in u_modesPhi length");
            foreach(i; 0 .. other.T_modesPhi.length) {
                T_modesPhi[i] = other.T_modesPhi[i]; u_modesPhi[i] = other.u_modesPhi[i];
            }
            assert(T_modesMax.length == other.T_modesMax.length, "Mismatch in T_modesMax length");
            assert(u_modesMax.length == other.u_modesMax.length, "Mismatch in u_modesMax length");
            foreach(i; 0 .. other.T_modesMax.length) {
                T_modesMax[i] = other.T_modesMax[i]; u_modesMax[i] = other.u_modesMax[i];
            }
            assert(T_modesMin.length == other.T_modesMin.length, "Mismatch in T_modesMin length");
            assert(u_modesMin.length == other.u_modesMin.length, "Mismatch in u_modesMin length");
            foreach(i; 0 .. other.T_modesMin.length) {
                T_modesMin[i] = other.T_modesMin[i]; u_modesMin[i] = other.u_modesMin[i];
            }
        }
    } // end copy_values_from()

    @nogc
    void barth_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
        size_t dimensions = myConfig.dimensions;
        number a, b, U, phi;
        immutable double w = 1.0e-12;
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (fabs("~gname~"[0]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[1]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[2]) > ESSENTIALLY_ZERO) {
            foreach (i, f; cell_cloud[0].iface) {
                number dx = f.pos.x - cell_cloud[0].pos[0].x; 
                number dy = f.pos.y - cell_cloud[0].pos[0].y; 
                number dz = f.pos.z - cell_cloud[0].pos[0].z;
                b = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) b += "~gname~"[2] * dz; 
                b = copysign(((fabs(b) + w)), b);
                if (b >  ESSENTIALLY_ZERO) {
                    a = "~qMaxname~" - U;
                    phi = fmin(phi, a/b);
                }
                else if (b <  ESSENTIALLY_ZERO) {
                    a = "~qMinname~" - U;
                    phi = fmin(phi, a/b);
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
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForLimits("B.x", "Bx", "BxPhi", "BxMax", "BxMin"));
                mixin(codeForLimits("B.y", "By", "ByPhi", "ByMax", "ByMin"));
                mixin(codeForLimits("B.z", "Bz", "BzPhi", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForLimits("psi", "psi", "psiPhi", "psiMax", "psiMin"));
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                mixin(codeForLimits("tke", "tke", "tkePhi", "tkeMax", "tkeMin"));
                mixin(codeForLimits("omega", "omega", "omegaPhi", "omegaMax", "omegaMin"));
            }
        }
        version(multi_species_gas) {
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
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        auto nmodes = myConfig.gmodel.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    @nogc
    void mlp_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                   ref LocalConfig myConfig, size_t gtl=0)
    {
        // The implementation of the MLP limiter follows the implementation in VULCAN
        // i.e. it uses the MLP approach, and the Van Leer limiting function
        // as outlined in the NASA publication on VULCAN's unstructured solver [White et al 2017]
        size_t dimensions = myConfig.dimensions;
        number a, b, U, phi, denom, numer, s;
        immutable double eps = 1.0e-12;
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (fabs("~gname~"[0]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[1]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[2]) > ESSENTIALLY_ZERO) {
                foreach (i, vtx; cell_cloud[0].vtx) {
                    number dx = vtx.pos[gtl].x - cell_cloud[0].pos[gtl].x; 
                    number dy = vtx.pos[gtl].y - cell_cloud[0].pos[gtl].y; 
                    number dz = vtx.pos[gtl].z - cell_cloud[0].pos[gtl].z;
                    a = "~gname~"[0] * dx + "~gname~"[1] * dy;
                    if (myConfig.dimensions == 3) a += "~gname~"[2] * dz;
                    a = copysign(((fabs(a) + eps)), a); 
                    if (fabs(a) > ESSENTIALLY_ZERO) {
                        b = (a > 0.0) ? vtx.gradients."~qMaxname~" - U: vtx.gradients."~qMinname~" - U;
                        numer = b*fabs(a) + a*fabs(b);
                        denom = fabs(a) + fabs(b) + eps;
                        s = (1.0/a) * (numer/denom);          
                    } else {
                        s = 1.0;
                    }
                    phi = fmin(phi, s);
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
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForLimits("B.x", "Bx", "BxPhi", "BxMax", "BxMin"));
                mixin(codeForLimits("B.y", "By", "ByPhi", "ByMax", "ByMin"));
                mixin(codeForLimits("B.z", "Bz", "BzPhi", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForLimits("psi", "psi", "psiPhi", "psiMax", "psiMin"));
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                mixin(codeForLimits("tke", "tke", "tkePhi", "tkeMax", "tkeMin"));
                mixin(codeForLimits("omega", "omega", "omegaPhi", "omegaMax", "omegaMin"));
            }
        }
        version(multi_species_gas) {
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
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        auto nmodes = myConfig.gmodel.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    @nogc
    void store_max_min_values_for_mlp_limiter(FVCell[] cell_cloud, ref LocalConfig myConfig)
    {
        // the MLP limiter stores the maximum and minimum flowstate values
        // that surround a node, at the node
        size_t dimensions = myConfig.dimensions;
        auto np = cell_cloud.length;
        // The following function to be used at compile time.
        string codeForGradients(string qname, string gname,
                                string qMaxname, string qMinname)
        {
            string code = "{
                number q0 = cell_cloud[0].fs."~qname~";
                "~qMaxname~" = q0;
                "~qMinname~" = q0;
                foreach (i; 1 .. np) {
                    "~qMaxname~" = fmax("~qMaxname~", cell_cloud[i].fs."~qname~");
                    "~qMinname~" = fmin("~qMinname~", cell_cloud[i].fs."~qname~");
                }
                }
                ";
            return code;
        }
        // x-velocity
        mixin(codeForGradients("vel.x", "velx", "velxMax", "velxMin"));
        mixin(codeForGradients("vel.y", "vely", "velyMax", "velyMin"));
        mixin(codeForGradients("vel.z", "velz", "velzMax", "velzMin"));
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForGradients("B.x", "Bx", "BxMax", "BxMin"));
                mixin(codeForGradients("B.y", "By", "ByMax", "ByMin"));
                mixin(codeForGradients("B.z", "Bz", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForGradients("psi", "psi", "psiMax", "psiMin"));
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                mixin(codeForGradients("tke", "tke", "tkeMax", "tkeMin"));
                mixin(codeForGradients("omega", "omega", "omegaMax", "omegaMin"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.gmodel.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForGradients("gas.massf[isp]", "massf[isp]",
                                           "massfMax[isp]", "massfMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                massf[0][0] = 0.0; massf[0][1] = 0.0; massf[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        auto nmodes = myConfig.gmodel.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                           "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.u", "u", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.u_modes[imode]", "u_modes[imode]",
                                           "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                           "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    @nogc
    void venkat_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                      ref LocalConfig myConfig, size_t gtl=0)
    {
        size_t dimensions = myConfig.dimensions;
        number a, b, U, phi, h, denom, numer, s;
        immutable double w = 1.0e-12;
        double K = myConfig.venkat_K_value;
        if (myConfig.dimensions == 3) {
            h = cell_cloud[0].volume[gtl]^^(1.0/3.0); // cbrt  
        } else {
            h = sqrt(cell_cloud[0].volume[gtl]);
        }
        number eps = (K*h) * (K*h) * (K*h);
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            if (fabs("~gname~"[0]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[1]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[2]) > ESSENTIALLY_ZERO) {
                foreach (i, f; cell_cloud[0].iface) {
                    number dx = f.pos.x - cell_cloud[0].pos[gtl].x; 
                    number dy = f.pos.y - cell_cloud[0].pos[gtl].y; 
                    number dz = f.pos.z - cell_cloud[0].pos[gtl].z;
                    b = "~gname~"[0] * dx + "~gname~"[1] * dy;
                    if (myConfig.dimensions == 3) { b += "~gname~"[2] * dz; }
                    b = copysign(((fabs(b) + w)), b); 
                    if (fabs(b) > ESSENTIALLY_ZERO) {
                        a = (b > 0.0) ? "~qMaxname~" - U: "~qMinname~" - U;
                        numer = (a*a + eps)*b + 2.0*b*b*a;
                        denom = a*a + 2.0*b*b + a*b + eps;
                        s = (1.0/b) * (numer/denom);                    
                    } else {
                        s = 1.0;
                    }
                    phi = fmin(phi, s);
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
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForLimits("B.x", "Bx", "BxPhi", "BxMax", "BxMin"));
                mixin(codeForLimits("B.y", "By", "ByPhi", "ByMax", "ByMin"));
                mixin(codeForLimits("B.z", "Bz", "BzPhi", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForLimits("psi", "psi", "psiPhi", "psiMax", "psiMin"));
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                mixin(codeForLimits("tke", "tke", "tkePhi", "tkeMax", "tkeMin"));
                mixin(codeForLimits("omega", "omega", "omegaPhi", "omegaMax", "omegaMin"));
            }
        }
        version(multi_species_gas) {
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
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        auto nmodes = myConfig.gmodel.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    @nogc
    void heuristic_van_albada_limit(FVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                                    ref LocalConfig myConfig, size_t gtl=0)
    {
        size_t dimensions = myConfig.dimensions;
        number phi, s, a, b;
        number h, dP, dPx, dPy, dPz, dx1, dy1, dz1, dx2, dy2, dz2;
        immutable double eps = 1.0e-12;
        // The following function to be used at compile time.
        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{
            phi = 1.0;
            if (fabs("~gname~"[0]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[1]) > ESSENTIALLY_ZERO ||
                fabs("~gname~"[2]) > ESSENTIALLY_ZERO) {
                foreach (i, f; cell_cloud[0].iface) {
                    // heuristic pressure limiter value
                    if (f.left_cell && f.right_cell && f.left_cell.is_interior && f.right_cell.is_interior) {
                        dx1 = f.pos.x - f.left_cell.pos[gtl].x; 
                        dy1 = f.pos.y - f.left_cell.pos[gtl].y; 
                        dz1 = f.pos.z - f.left_cell.pos[gtl].z;
                        dx2 = f.pos.x - f.right_cell.pos[gtl].x; 
                        dy2 = f.pos.y - f.right_cell.pos[gtl].y; 
                        dz2 = f.pos.z - f.right_cell.pos[gtl].z;
                        dPx = dx1*f.left_cell.gradients.p[0] - dx2*f.right_cell.gradients.p[0];
                        dPy = dy1*f.left_cell.gradients.p[1] - dy2*f.right_cell.gradients.p[1];
                        dPz = dz1*f.left_cell.gradients.p[2] - dz2*f.right_cell.gradients.p[2];
                        dP = sqrt(dPx*dPx + dPy*dPy + dPz*dPz);
                        h = 1.0 - tanh(dP/fmin(f.left_cell.fs.gas.p, f.right_cell.fs.gas.p));
                    } else if (f.left_cell && f.left_cell.is_interior) {
                        dx1 = f.pos.x - f.left_cell.pos[gtl].x; 
                        dy1 = f.pos.y - f.left_cell.pos[gtl].y; 
                        dz1 = f.pos.z - f.left_cell.pos[gtl].z;
                        dPx = dx1*f.left_cell.gradients.p[0];
                        dPy = dy1*f.left_cell.gradients.p[1];
                        dPz = dz1*f.left_cell.gradients.p[2];
                        dP = sqrt(dPx*dPx + dPy*dPy + dPz*dPz);
                        h = 1.0 - tanh(dP/f.left_cell.fs.gas.p);
                    } else if (f.right_cell && f.right_cell.is_interior) {
                        dx1 = f.pos.x - f.right_cell.pos[gtl].x; 
                        dy1 = f.pos.y - f.right_cell.pos[gtl].y; 
                        dz1 = f.pos.z - f.right_cell.pos[gtl].z;
                        dPx = dx1*f.right_cell.gradients.p[0];
                        dPy = dy1*f.right_cell.gradients.p[1];
                        dPz = dz1*f.right_cell.gradients.p[2];
                        dP = sqrt(dPx*dPx + dPy*dPy + dPz*dPz);
                        h = 1.0 - tanh(dP/f.right_cell.fs.gas.p);
                    }
                    // Van Albada limiter value
                    if (f.left_cell && f.right_cell && f.left_cell.is_interior && f.right_cell.is_interior) {
                        a = sqrt(f.left_cell.gradients."~gname~"[0]^^2 + 
                                 f.left_cell.gradients."~gname~"[1]^^2 + 
                                 f.left_cell.gradients."~gname~"[2]^^2);
                        b = sqrt(f.right_cell.gradients."~gname~"[0]^^2 + 
                                 f.right_cell.gradients."~gname~"[1]^^2 + 
                                 f.right_cell.gradients."~gname~"[2]^^2);
                        s = h * (a*b + fabs(a*b))/(a*a + b*b + eps);
                    } else if (f.left_cell && f.left_cell.is_interior) {
                        s = h * 1.0;
                    } else if (f.right_cell && f.right_cell.is_interior) {
                        s = h * 1.0;  
                    }
                    phi = fmin(phi, s);
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
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForLimits("B.x", "Bx", "BxPhi", "BxMax", "BxMin"));
                mixin(codeForLimits("B.y", "By", "ByPhi", "ByMax", "ByMin"));
                mixin(codeForLimits("B.z", "Bz", "BzPhi", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForLimits("psi", "psi", "psiPhi", "psiMax", "psiMin"));
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                mixin(codeForLimits("tke", "tke", "tkePhi", "tkeMax", "tkeMin"));
                mixin(codeForLimits("omega", "omega", "omegaPhi", "omegaMax", "omegaMin"));
            }
        }
        version(multi_species_gas) {
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
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        auto nmodes = myConfig.gmodel.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.u", "u", "uPhi", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.p", "p", "pPhi", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForLimits("gas.rho", "rho", "rhoPhi", "rhoMax", "rhoMin"));
            mixin(codeForLimits("gas.T", "T", "TPhi", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.T_modes[imode]", "T_modes[imode]", "T_modesPhi[imode]",
                                        "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
    } // end compute_lsq_gradients()

    @nogc
    void compute_lsq_values(FVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                            ref LocalConfig myConfig)
    {
        size_t dimensions = myConfig.dimensions;
        auto np = cell_cloud.length;
        // The following function to be used at compile time.
        string codeForGradients(string qname, string gname,
                                string qMaxname, string qMinname)
        {
            string code = "{
                number q0 = cell_cloud[0].fs."~qname~";
                "~qMaxname~" = q0;
                "~qMinname~" = q0;
                "~gname~"[0] = 0.0; "~gname~"[1] = 0.0; "~gname~"[2] = 0.0;
                foreach (i; 1 .. np) {
                    number dq = cell_cloud[i].fs."~qname~" - q0;
                    "~gname~"[0] += ws.wx[i] * dq;
                    "~gname~"[1] += ws.wy[i] * dq;
                    if (dimensions == 3) { "~gname~"[2] += ws.wz[i] * dq; }
                    "~qMaxname~" = fmax("~qMaxname~", cell_cloud[i].fs."~qname~");
                    "~qMinname~" = fmin("~qMinname~", cell_cloud[i].fs."~qname~");
                }
                }
                ";
            return code;
        }
        // x-velocity
        mixin(codeForGradients("vel.x", "velx", "velxMax", "velxMin"));
        mixin(codeForGradients("vel.y", "vely", "velyMax", "velyMin"));
        mixin(codeForGradients("vel.z", "velz", "velzMax", "velzMin"));
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForGradients("B.x", "Bx", "BxMax", "BxMin"));
                mixin(codeForGradients("B.y", "By", "ByMax", "ByMin"));
                mixin(codeForGradients("B.z", "Bz", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForGradients("psi", "psi", "psiMax", "psiMin"));
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                mixin(codeForGradients("tke", "tke", "tkeMax", "tkeMin"));
                mixin(codeForGradients("omega", "omega", "omegaMax", "omegaMin"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.gmodel.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForGradients("gas.massf[isp]", "massf[isp]",
                                           "massfMax[isp]", "massfMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                massf[0][0] = 0.0; massf[0][1] = 0.0; massf[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        auto nmodes = myConfig.gmodel.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                           "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.u", "u", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.u_modes[imode]", "u_modes[imode]",
                                           "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            break;
        case InterpolateOption.rhot: 
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.T", "T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]",
                                           "T_modesMax[imode]", "T_modesMin[imode]"));
                }
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

    @nogc number clip_to_limits(number q, number A, number B)
    // Returns q if q is between the values A and B, else
    // it returns the closer limit of the range [A,B].
    {
        number lower_limit = fmin(A, B);
        number upper_limit = fmax(A, B);
        return fmin(upper_limit, fmax(lower_limit, q));
    } // end clip_to_limits()

    @nogc void min_mod_limit(ref number a, ref number b)
    // A rough slope limiter.
    {
        if (a * b < 0.0) {
            a = 0.0; b = 0.0;
            return;
        }
        a = copysign(fmin(fabs(a), fabs(b)), a);
        b = a;
    }

    @nogc void van_albada_limit(ref number a, ref number b)
    // A smooth slope limiter.
    {
        immutable double eps = 1.0e-12;
        number s = (a*b + fabs(a*b))/(a*a + b*b + eps);
        a *= s;
        b *= s;
    }

    @nogc
    void interp_both(ref FVInterface IFace, size_t gtl,
                     ref FlowState Lft, ref FlowState Rght,
                     bool allow_high_order_interpolation)
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
        // enforce first order reconstruction for cells that capure the shocks
        if (GlobalConfig.suppress_reconstruction_at_captured_shocks)
            if ( Lft.S || Rght.S ) return;
        // else apply higher-order interpolation to all faces
        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            // High-order reconstruction for some properties.
            //
            LSQInterpWorkspace wsL = cL0.ws;
            LSQInterpWorkspace wsR = cR0.ws;
            // vector from left-cell-centre to face midpoint
            number dLx = IFace.pos.x - cL0.pos[gtl].x;
            number dLy = IFace.pos.y - cL0.pos[gtl].y;
            number dLz = IFace.pos.z - cL0.pos[gtl].z;
            number dRx = IFace.pos.x - cR0.pos[gtl].x;
            number dRy = IFace.pos.y - cR0.pos[gtl].y;
            number dRz = IFace.pos.z - cR0.pos[gtl].z;
            number[3] mygradL, mygradR;
            //
            // Always reconstruct in the global frame of reference -- for now
            //
            // x-velocity
            string codeForReconstruction(string qname, string gname,
                                         string tname, string lname)
            {
                string code = "{
                number qL0 = cL0.fs."~qname~";
                number qMinL = qL0;
                number qMaxL = qL0;
                mygradL[0] = cL0.gradients."~gname~"[0];
                mygradL[1] = cL0.gradients."~gname~"[1];
                mygradL[2] = cL0.gradients."~gname~"[2];
                number qR0 = cR0.fs."~qname~";
                number qMinR = qR0;
                number qMaxR = qR0;
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
                    case UnstructuredLimiter.heuristic_van_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat:
                        mygradL[0] *= cL0.gradients."~lname~";
                        mygradL[1] *= cL0.gradients."~lname~";
                        mygradL[2] *= cL0.gradients."~lname~";
                        mygradR[0] *= cR0.gradients."~lname~";
                        mygradR[1] *= cR0.gradients."~lname~";
                        mygradR[2] *= cR0.gradients."~lname~";
                        break;
                    }
                }
                number qL = qL0 + dLx*mygradL[0] + dLy*mygradL[1];
                number qR = qR0 + dRx*mygradR[0] + dRy*mygradR[1];
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
            version(MHD) {
                if (myConfig.MHD) {
                    mixin(codeForReconstruction("B.x", "Bx", "B.refx", "BxPhi"));
                    mixin(codeForReconstruction("B.y", "By", "B.refy", "ByPhi"));
                    mixin(codeForReconstruction("B.z", "Bz", "B.refz", "BxPhi"));
                    if (myConfig.divergence_cleaning) {
                        mixin(codeForReconstruction("psi", "psi", "psi", "psiPhi"));
                    }
                }
            }
            version(komega) {
                if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                    mixin(codeForReconstruction("tke", "tke", "tke", "tkePhi"));
                    mixin(codeForReconstruction("omega", "omega", "omega", "omegaPhi"));
                }
            }
            version(multi_species_gas) {
                if (nsp > 1) {
                    // Multiple species.
                    foreach (isp; 0 .. nsp) {
                        mixin(codeForReconstruction("gas.massf[isp]", "massf[isp]",
                                                    "gas.massf[isp]", "massfPhi[isp]"));
                    }
                    try {
                        scale_mass_fractions(Lft.gas.massf); 
                    } catch (Exception e) {
                        debug { writeln(e.msg); }
                        Lft.gas.massf[] = IFace.left_cell.fs.gas.massf[];
                    }
                    try {
                        scale_mass_fractions(Rght.gas.massf);
                    } catch (Exception e) {
                        debug { writeln(e.msg); }
                        Rght.gas.massf[] = IFace.right_cell.fs.gas.massf[];
                    }
                } else {
                    // Only one possible mass-fraction value for a single species.
                    Lft.gas.massf[0] = 1.0;
                    Rght.gas.massf[0] = 1.0;
                }
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
                    debug { writeln(e.msg); }
                    Lft.copy_values_from(IFace.left_cell.fs);
                }
                try {
                    gmodel.update_thermo_from_"~funname~"(Rght.gas);
                } catch (Exception e) {
                    debug { writeln(e.msg); }
                    Rght.copy_values_from(IFace.right_cell.fs);
                }
                ";
                return code;
            }
            final switch (myConfig.thermo_interpolator) {
            case InterpolateOption.pt: 
                mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                version(multi_T_gas) {
                    foreach (imode; 0 .. nmodes) {
                        mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                    "gas.T_modes[imode]", "T_modesPhi[imode]"));
                    }
                }
                mixin(codeForThermoUpdate("pT"));
                break;
            case InterpolateOption.rhou:
                mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                mixin(codeForReconstruction("gas.u", "u", "gas.u", "uPhi"));
                version(multi_T_gas) {
                    foreach (imode; 0 .. nmodes) {
                        mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                    "gas.u_modes[imode]", "u_modesPhi[imode]"));
                    }
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
                version(multi_T_gas) {
                    foreach (imode; 0 .. nmodes) {
                        mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                    "gas.T_modes[imode]", "T_modesPhi[imode]"));
                    }
                }
                mixin(codeForThermoUpdate("rhoT"));
                break;
            } // end switch thermo_interpolator
        } // end of high-order reconstruction
    } // end interp_both()

    @nogc
    void interp_right(ref FVInterface IFace, size_t gtl, ref FlowState Rght,
                      bool allow_high_order_interpolation)
    // Interpolate the flow field quantities at the right-side of the interface,
    // given information in right-cell attached to this interface.
    {
        auto gmodel = myConfig.gmodel;
        auto nsp = gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        FVCell cR0 = IFace.right_cell;
        // Low-order reconstruction just copies data from adjacent FV_Cell.
        // Even for high-order reconstruction, we depend upon this copy for
        // the viscous-transport and diffusion coefficients.
        Rght.copy_values_from(cR0.fs);
        // [TODO] 2018-06-10 higher-order interpolation.
        return;
    }

    @nogc
    void interp_left(ref FVInterface IFace, size_t gtl, ref FlowState Lft,
                     bool allow_high_order_interpolation)
    // Interpolate the flow field quantities at the left-side of the interface,
    // given information in the left-cell attached to this interface.
    {
        auto gmodel = myConfig.gmodel;
        auto nsp = gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        FVCell cL0 = IFace.left_cell;
        // Low-order reconstruction just copies data from adjacent FV_Cell.
        // Even for high-order reconstruction, we depend upon this copy for
        // the viscous-transport and diffusion coefficients.
        Lft.copy_values_from(cL0.fs);
        // [TODO] 2018-06-10 higher-order interpolation.
        return;
    }

} // end class LsqInterpolator

