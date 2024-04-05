// lsqinterp.d
// Least-squares interpolation/reconstruction of flow field.
//

module lsqinterp;

import std.math;
import std.stdio;
import std.algorithm;
import std.conv;
import ntypes.complex;
import nm.number;
import nm.rsla;
import nm.limiters;

import geom;
import gas;
import globalconfig;
import flowstate;
import fvinterface;
import lmr.fluidfvcell;

immutable size_t cloud_nmax = 112; // ??? WHY? [TODO] Add some comment for this selection please
immutable double ESSENTIALLY_ZERO = 1.0e-50;


// TODO: These objects violate RAII. Is there a good reason for it? (NNG 30/05/22)
struct LSQInterpWorkspace {
public:
    // A place to hold the intermediate results for computing
    // the least-squares model as a weighted sum of the flow data.
    number[cloud_nmax] wx, wy, wz;

    this(in LSQInterpWorkspace other)
    {
        wx[] = other.wx[]; wy[] = other.wy[]; wz[] = other.wz[];
    }

    @nogc
    void assemble_and_invert_normal_matrix(FluidFVCell[] cell_cloud, int dimensions, size_t gtl)
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
} // end struct LSQInterpWorkspace

struct LSQInterpGradients {
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
    version(turbulence) {
        number[3][2] turb;
        number[2] turbPhi;
        number[2] turbMax;
        number[2] turbMin;
    }
    version(multi_species_gas) {
        number[3][] rho_s;
        number[] rho_sPhi;
        number[] rho_sMax;
        number[] rho_sMin;
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


    this(size_t nsp, size_t nmodes, size_t nturb)
    {
        version(multi_species_gas) {
            rho_s.length    = nsp;
            rho_sPhi.length = nsp;
            rho_sMax.length = nsp;
            rho_sMin.length = nsp;
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

        // The user has the option to postprocess the limiter values for visualisation,
        // so we initialize these values to -1.0.
        // We need to do this because:
        // (1) in a simulation, we do not calculate a limiter value for every thermodynamic variable,
        //     (note that this means if a variable has limiter values of -1.0 in a solution, then we did
        //     not calculate or use that particular limiter value in that simulation)
        // and
        // (2) old versions of Paraview throw a segfault error when plotting NaNs
        //
        // KAD 22-10-2022
        //
        velxPhi = -1.0; velyPhi = -1.0; velzPhi = -1.0;
        version(MHD) {
            BxPhi = -1.0; ByPhi = -1.0; BzPhi = -1.0; psiPhi = -1.0;
        }
        version(turbulence) {
            foreach (ref val; turbPhi) { val = -1.0; }
        }
        version(multi_species_gas) {
            foreach (ref val; rho_sPhi) { val = -1.0; }
        }
        rhoPhi = -1.0; pPhi = -1.0;
        TPhi = -1.0; uPhi = -1.0;
        version(multi_T_gas) {
            foreach (ref val; u_modesPhi) { val = -1.0; }
            foreach (ref val; T_modesPhi) { val = -1.0; }
        }
    }

    this(ref const(LSQInterpGradients) other)
    {
        size_t nsp=0;
        size_t nmodes=0;
        size_t nturb=0;
        version(multi_species_gas) { nsp = other.rho_s.length; }
        version(multi_T_gas) { nmodes = other.T_modes.length;}
        version(turbulence) { nturb = other.turb.length;}
        this(nsp, nmodes, nturb);
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
        version(turbulence) {
            assert(turb.length == other.turb.length, "Mismatch in turb length");
            foreach(i; 0 .. turb.length){
                turb[i][] = other.turb[i][];
                turbPhi[i] = other.turbPhi[i];
                turbMax[i] = other.turbMax[i];
                turbMin[i] = other.turbMin[i];
            }
        }
        version(multi_species_gas) {
            assert(rho_s.length == other.rho_s.length, "Mismatch in rho_s length");
            foreach(i; 0 .. other.rho_s.length) { rho_s[i][] = other.rho_s[i][]; }
            assert(rho_sPhi.length == other.rho_sPhi.length, "Mismatch in rho_sPhi length");
            foreach(i; 0 .. other.rho_s.length) { rho_sPhi[i] = other.rho_sPhi[i]; }
            assert(rho_sMax.length == other.rho_sMax.length, "Mismatch in rho_sMax length");
            foreach(i; 0 .. other.rho_s.length) { rho_sMax[i] = other.rho_sMax[i]; }
            assert(rho_sMin.length == other.rho_sMin.length, "Mismatch in rho_sMin length");
            foreach(i; 0 .. other.rho_s.length) { rho_sMin[i] = other.rho_sMin[i]; }
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
    void barth_limit(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws, ref LocalConfig myConfig)
    {
        // This is the classic Barth and Jespersen limiter from ref. [1], refer to ref. [2] for implementation details.
        //
        // references:
        // [1] Barth TJ, Jespersen DC.
        //     The design and application of upwind schemes on unstructured meshes.
        //     AIAA Paper 89-0366; 1989
        // [2] Blazek J.
        //     CFD principles and applications
        //     pg. 156, Third edition, 2007
        //
        number delp, delm, U, phi;
        immutable double eps = 1.0e-25;
        string codeForLimits(string qname, string gname,
                             string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            foreach (i, f; cell_cloud[0].iface) {
                number dx = f.pos.x - cell_cloud[0].pos[0].x;
                number dy = f.pos.y - cell_cloud[0].pos[0].y;
                number dz = f.pos.z - cell_cloud[0].pos[0].z;
                delm = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) { delm += "~gname~"[2] * dz; }
                delm += eps;
                if (delm >= 0.0) {
                    delp = "~qMaxname~" - U;
                } else {
                    delp = "~qMinname~" - U;
                }
                phi = fmin(phi, delp/delm);
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
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForLimits("turb[it]","turb[it]","turbPhi[it]","turbMax[it]","turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForLimits("gas.rho_s[isp]", "rho_s[isp]", "rho_sPhi[isp]",
                                        "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
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
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
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
    } // end barth_limit()

    @nogc
    void venkat_mlp_limit(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                          bool apply_heuristic_pressure_limiter, ref LocalConfig myConfig, size_t gtl=0)
    {
        // This is an implementation of the multi-dimensional limiting process from ref. [1],
        // which  uses the Venkatakrishnan limiter function from ref. [2].
        //
        // Note that we use non-dimensional quantities in the limiter function to improve the
        // effect of the smoothing parameter.
        //
        // references:
        // [1] J. S. Park and C. Kim
        //     Multi-dimensional limiting process for ﬁnite volume methods on unstructured grids
        //     Computers & Fluids, vol. 65, pg. 8-24 (2012)
        // [2] V. Venkatakrishnan
        //     Convergence to steady state solutions of the Euler equations on unstructured grids with limiters.
        //     Journal of Computational Physics vol.118 pp120-130 (1995)
        //

        number delp, delm, delu, theta, U, Umax, Umin, phi, h, phi_f, nondim;
        immutable double K = myConfig.smooth_limiter_coeff;
        if (myConfig.dimensions == 3) {
            h = pow(cell_cloud[0].volume[gtl], 1.0/3.0);
        } else {
            h = sqrt(cell_cloud[0].volume[gtl]);
        }
        number eps2;

        // Park heuristic pressure limiter
        number phi_hp = 1.0;
        if (apply_heuristic_pressure_limiter) {
            park_limit(cell_cloud, ws, myConfig);
            phi_hp = velxPhi; // we choose velxPhi here since it will always be set regardless of the thermo_interpolator
        }

        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{ ";
            if ( qname == "vel.x" || qname == "vel.y" || qname == "vel.z" ) {
                code ~= "number velMax = sqrt(velxMax*velxMax+velyMax*velyMax+velzMax*velzMax);
                         number velMin = sqrt(velxMin*velxMin+velyMin*velyMin+velzMin*velzMin);
                         nondim = 0.5*fabs(velMax+velMin) + 1.0e-25;";
            } else if ( qname == "gas.rho_s[isp]" ) {
                code ~= "nondim = 1.0;";
            } else {
                code ~= "nondim = 0.5*fabs("~qMaxname~" + "~qMinname~");";
            }
            code ~= "
            U = cell_cloud[0].fs."~qname~";
            Umin = U;
            Umax = U;
            foreach (ci; cell_cloud) {
                foreach (cj; ci.cell_cloud) {
                    Umin = fmin(Umin, cj.fs."~qname~");
                    Umax = fmax(Umax, cj.fs."~qname~");
                }
            }
            delu = (Umax-Umin)/nondim;
            theta = delu/(K*pow(h, 1.5));
            eps2 = (K*delu*delu)/(1.0+theta);
            eps2 += 1.0e-25; // prevent division by zero
            phi = 1.0;
            foreach (i, v; cell_cloud[0].vtx) {
                number dx = v.pos[gtl].x - cell_cloud[0].pos[gtl].x;
                number dy = v.pos[gtl].y - cell_cloud[0].pos[gtl].y;
                number dz = v.pos[gtl].z - cell_cloud[0].pos[gtl].z;
                delm = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) { delm += "~gname~"[2] * dz; }
                delp = (delm >= 0.0) ? Umax - U: Umin - U;
                delp /= nondim;
                delm /= nondim;
                if (delm == 0.0) {
                    phi_f = 1.0;
                } else {
                    phi_f = (delp*delp + 2.0*delp*delm + eps2)/(delp*delp + 2.0*delm*delm + delp*delm + eps2);
                }
                phi = fmin(phi, phi_f);
            }
            "~limFactorname~" = phi*phi_hp;
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
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForLimits("turb[it]","turb[it]","turbPhi[it]","turbMax[it]","turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForLimits("gas.rho_s[isp]", "rho_s[isp]", "rho_sPhi[isp]",
                                        "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
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
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
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
    } // end venkat_mlp_limit()

    @nogc
    void van_albada_limit(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                          bool apply_heuristic_pressure_limiter, ref LocalConfig myConfig, size_t gtl=0)
    {
        // This is the classic Van Albada limiter from ref. [1] implemented in the unstructured grid
        // format detailed in ref. [2]. The original limiter function returned the limited slope directly,
        // ref. [3] presents the limiter form as eq. B.11, which returns a value between 0 and 1.
        //
        // Note that we use non-dimensional quantities in the limiter function for proper application of the
        // effect of the smoothing parameter.
        //
        // references:
        // [1] G. D. van Albada, B. van Leer and W. W. Roberts
        //     A Comparative Study of Computational Methods in Cosmic Gas Dynamics
        //     Astron. Astrophys., Vol. 108, No. 1, 1982, pp. 76–84
        // [2] J. A. White, R. A. Baurle, B. J. Pase, S. C. Spiegel and H. Nishikawa
        //     Geometrically flexible and efficient flow analysis of high speed vehicles via domain decomposition
        //     2017 JANNAF paper
        // [3] H. Nishikawa
        //     New Unstructured-Grid Limiter Functions
        //     2022 AIAA SciTech forum
        //
        number delp, delm, U, phi, phi_f, h, nondim;
        immutable double eps = myConfig.smooth_limiter_coeff;

        // Park heuristic pressure limiter
        number phi_hp = 1.0;
        if (apply_heuristic_pressure_limiter) {
            park_limit(cell_cloud, ws, myConfig);
            phi_hp = velxPhi; // we choose velxPhi here since it will always be set regardless of the thermo_interpolator
        }

        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{ ";
            if ( qname == "vel.x" || qname == "vel.y" || qname == "vel.z" ) {
                code ~= "number velMax = sqrt(velxMax*velxMax+velyMax*velyMax+velzMax*velzMax);
                         number velMin = sqrt(velxMin*velxMin+velyMin*velyMin+velzMin*velzMin);
                         nondim = 0.5*fabs(velMax+velMin) + 1.0e-25;";
            } else if ( qname == "gas.rho_s[isp]" ) {
                code ~= "nondim = 1.0;";
            } else {
                code ~= "nondim = 0.5*fabs("~qMaxname~" + "~qMinname~");";
            }
            code ~= "
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            foreach (i, f; cell_cloud[0].iface) {
                number dx = f.pos.x - cell_cloud[0].pos[gtl].x;
                number dy = f.pos.y - cell_cloud[0].pos[gtl].y;
                number dz = f.pos.z - cell_cloud[0].pos[gtl].z;
                delm = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) { delm += "~gname~"[2] * dz; }
                delp = (delm >= 0.0) ? 0.5*("~qMaxname~" - U): 0.5*("~qMinname~" - U);
                delp /= nondim;
                delm /= nondim;
                if (delm == 0.0) {
                    phi_f = 1.0;
                } else {
                    phi_f = (2.0*(delm*delp+eps*eps))/(delp*delp+delm*delm+2.0*eps*eps);
                }
                phi = fmin(phi, phi_f);
            }
            "~limFactorname~" = phi*phi_hp;
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
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForLimits("turb[it]","turb[it]","turbPhi[it]","turbMax[it]","turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForLimits("gas.rho_s[isp]", "rho_s[isp]", "rho_sPhi[isp]",
                                        "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
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
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
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
    } // end van_albada_limit()

    @nogc
    void venkat_limit(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                      bool apply_heuristic_pressure_limiter, ref LocalConfig myConfig, size_t gtl=0)
    {
        // This is the classic Venkatakrishnan limiter from ref. [1], refer to ref. [2] for implementation details.
        // We use the new definition for eps2 as given in ref. [3].
        // We have found that it provides a less noisy triggering of the limiter compared to the original
        // definition provided in ref. [1] of eps2 = (K*h) * (K*h) * (K*h).
        //
        // Note that we use non-dimensional quantities in the limiter function to improve the
        // effect of the smoothing parameter.
        //
        // references:
        // [1] V. Venkatakrishnan
        //     Convergence to steady state solutions of the Euler equations on unstructured grids with limiters.
        //     Journal of Computational Physics vol.118 pp120-130 (1995)
        // [2] Blazek J.
        //     CFD principles and applications
        //     pg. 156, Third edition, 2007
        // [3] J. S. Park and C. Kim
        //     Multi-dimensional limiting process for ﬁnite volume methods on unstructured grids
        //     Computers & Fluids, vol. 65, pg. 8-24 (2012)
        //
        number delp, delm, delu, U, theta, phi, h, phi_f, nondim;
        immutable double K = myConfig.smooth_limiter_coeff;
        if (myConfig.dimensions == 3) {
            h = pow(cell_cloud[0].volume[gtl], 1.0/3.0);
        } else {
            h = sqrt(cell_cloud[0].volume[gtl]);
        }
        number eps2;

        // Park heuristic pressure limiter
        number phi_hp = 1.0;
        if (apply_heuristic_pressure_limiter) {
            park_limit(cell_cloud, ws, myConfig);
            phi_hp = velxPhi; // we choose velxPhi here since it will always be set regardless of the thermo_interpolator
        }

        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{ ";
            if ( qname == "vel.x" || qname == "vel.y" || qname == "vel.z" ) {
                code ~= "number velMax = sqrt(velxMax*velxMax+velyMax*velyMax+velzMax*velzMax);
                         number velMin = sqrt(velxMin*velxMin+velyMin*velyMin+velzMin*velzMin);
                         nondim = 0.5*fabs(velMax+velMin) + 1.0e-25;";
            } else if ( qname == "gas.rho_s[isp]" ) {
                code ~= "nondim = 1.0;";
            } else {
                code ~= "nondim = 0.5*fabs("~qMaxname~" + "~qMinname~");";
            }
            code ~= "
            U = cell_cloud[0].fs."~qname~";
            delu = ("~qMaxname~"-"~qMinname~")/nondim;
            theta = delu/(K*pow(h, 1.5));
            eps2 = (K*delu*delu)/(1.0+theta);
            eps2 += 1.0e-25; // prevent division by zero
            phi = 1.0;
            foreach (i, f; cell_cloud[0].iface) {
                number dx = f.pos.x - cell_cloud[0].pos[gtl].x;
                number dy = f.pos.y - cell_cloud[0].pos[gtl].y;
                number dz = f.pos.z - cell_cloud[0].pos[gtl].z;
                delm = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) { delm += "~gname~"[2] * dz; }
                delp = (delm >= 0.0) ? "~qMaxname~" - U: "~qMinname~" - U;
                delp /= nondim;
                delm /= nondim;
                if (delm == 0.0) {
                    phi_f = 1.0;
                } else {
                    phi_f = (delp*delp + 2.0*delp*delm + eps2)/(delp*delp + 2.0*delm*delm + delp*delm + eps2);
                }
                phi = fmin(phi, phi_f);
            }
            "~limFactorname~" = phi*phi_hp;
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
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForLimits("turb[it]","turb[it]","turbPhi[it]","turbMax[it]","turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForLimits("gas.rho_s[isp]", "rho_s[isp]", "rho_sPhi[isp]",
                                        "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
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
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
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
    } // end venkat_limit()

    @nogc
    void nishikawa_limit(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                         bool apply_heuristic_pressure_limiter, ref LocalConfig myConfig, size_t gtl=0)
    {
        // This is the R3 limiter from ref. [1].
        //
        // Note that we use non-dimensional quantities in the limiter function to improve the
        // effect of the smoothing parameter.
        //
        // references:
        // [1] H. Nishikawa
        //     New Unstructured-Grid Limiter Functions
        //     2022 AIAA SciTech Forum
        //
        number a, b, S, delp, delm, U, phi, h, phi_f, nondim;
        immutable double K = myConfig.smooth_limiter_coeff;
        if (myConfig.dimensions == 3) {
            h = pow(cell_cloud[0].volume[gtl], 1.0/3.0);
        } else {
            h = sqrt(cell_cloud[0].volume[gtl]);
        }
        number eps = (K*h)*(K*h)*(K*h)*(K*h);

        // Park heuristic pressure limiter
        number phi_hp = 1.0;
        if (apply_heuristic_pressure_limiter) {
            park_limit(cell_cloud, ws, myConfig);
            phi_hp = velxPhi; // we choose velxPhi here since it will always be set regardless of the thermo_interpolator
        }

        string codeForLimits(string qname, string gname, string limFactorname,
                             string qMaxname, string qMinname)
        {
            string code = "{ ";
            if ( qname == "vel.x" || qname == "vel.y" || qname == "vel.z" ) {
                code ~= "number velMax = sqrt(velxMax*velxMax+velyMax*velyMax+velzMax*velzMax);
                         number velMin = sqrt(velxMin*velxMin+velyMin*velyMin+velzMin*velzMin);
                         nondim = 0.5*fabs(velMax+velMin) + 1.0e-25;";
            } else if ( qname == "gas.rho_s[isp]" ) {
                code ~= "nondim = 1.0;";
            } else {
                code ~= "nondim = 0.5*fabs("~qMaxname~" + "~qMinname~");";
            }
            code ~= "
            U = cell_cloud[0].fs."~qname~";
            phi = 1.0;
            foreach (i, f; cell_cloud[0].iface) {
                number dx = f.pos.x - cell_cloud[0].pos[gtl].x;
                number dy = f.pos.y - cell_cloud[0].pos[gtl].y;
                number dz = f.pos.z - cell_cloud[0].pos[gtl].z;
                delm = "~gname~"[0] * dx + "~gname~"[1] * dy;
                if (myConfig.dimensions == 3) { delm += "~gname~"[2] * dz; }
                delp = (delm >= 0.0) ? "~qMaxname~" - U: "~qMinname~" - U;
                delp /= nondim;
                a = fabs(delp);
                delm /= nondim;
                b = fabs(delm);
                if (a > 2*b) {
                    phi_f = 1.0;
                } else {
                    S = 4*b*b;
                    phi_f = (a*a*a + eps + a*S)/(a*a*a + eps + b*(delp*delp + S));
                }
                phi = fmin(phi, phi_f);
            }
            "~limFactorname~" = phi*phi_hp;
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
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForLimits("turb[it]","turb[it]","turbPhi[it]","turbMax[it]","turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForLimits("gas.rho_s[isp]", "rho_s[isp]", "rho_sPhi[isp]",
                                        "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
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
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForLimits("gas.u_modes[imode]", "u_modes[imode]", "u_modesPhi[imode]",
                                        "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
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
    } // end nishikawa_limit()

    @nogc
    void park_limit(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws,
                    ref LocalConfig myConfig, size_t gtl=0)
    {
        // Pressure-based heuristic limiter
        // Implementation details from
        //     M. A. Park
        //     Anisotropic Output-Based Adaptation with Tetrahedral Cut Cells for Compressible Flows
        //     Thesis @ Massachusetts Institute of Technology, 2008
        // Modified by NNG to use the harmonic mean of the s value computed at each face,
        // instead of just the minimum s value. This smooths out the limiter spatially and
        // removes a source of noise in the reconstruction process, while still prioritising
        // small values of s in the averaging process. (09/08/22)

        FluidFVCell ncell;
        number pmin;
        number phi = 0.0;
        number n = 0.0;
        foreach (i, f; cell_cloud[0].iface) {
            if (f.left_cell.id == cell_cloud[0].id) { ncell = f.right_cell; }
            else { ncell = f.left_cell; }
            number dx1 = f.pos.x - cell_cloud[0].pos[gtl].x;
            number dy1 = f.pos.y - cell_cloud[0].pos[gtl].y;
            number dz1 = f.pos.z - cell_cloud[0].pos[gtl].z;
            //number dx2 = f.pos.x - ncell.pos[gtl].x;
            //number dy2 = f.pos.y - ncell.pos[gtl].y;
            //number dz2 = f.pos.z - ncell.pos[gtl].z;
            //number dpx = dx1*cell_cloud[0].gradients.p[0] - dx2*ncell.gradients.p[0];
            //number dpy = dy1*cell_cloud[0].gradients.p[1] - dy2*ncell.gradients.p[1];
            //number dpz = dz1*cell_cloud[0].gradients.p[2] - dz2*ncell.gradients.p[2];
            // this step is a modification on the original algorithm since we don't have cell gradients from neighbouring blocks at this point
            number dpx = dx1*cell_cloud[0].gradients.p[0] - (0.5*(cell_cloud[0].fs.gas.p + ncell.fs.gas.p) - ncell.fs.gas.p);
            number dpy = dy1*cell_cloud[0].gradients.p[1] - (0.5*(cell_cloud[0].fs.gas.p + ncell.fs.gas.p) - ncell.fs.gas.p);
            number dpz = dz1*cell_cloud[0].gradients.p[2] - (0.5*(cell_cloud[0].fs.gas.p + ncell.fs.gas.p) - ncell.fs.gas.p);
            number dp = dpx*dpx + dpy*dpy;
            if (myConfig.dimensions == 3) { dp += dpz*dpz; }
            dp = sqrt(dp);
            pmin = fmin(cell_cloud[0].fs.gas.p, ncell.fs.gas.p);
            number s = 1.0-tanh(dp/pmin);
            phi += 1.0/fmax(s, 1e-16);
            n += 1.0;
        }
        phi = n/phi;

        // the limiter value for each variable is set to the pressure-based value
        velxPhi = phi; velyPhi = phi; velzPhi = phi;
        version(MHD) {
            if (myConfig.MHD) {
                BxPhi = phi; ByPhi = phi; BzPhi = phi;
                if (myConfig.divergence_cleaning) { psiPhi = phi; }
            }
        }
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) { turbPhi[it] = phi; }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) { rho_sPhi[isp] = phi; }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }

        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            pPhi = phi; TPhi = phi;
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) { T_modesPhi[imode] = phi; }
            }
            break;
        case InterpolateOption.rhou:
            rhoPhi = phi; uPhi = phi;
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) { u_modesPhi[imode] = phi; }
            }
            break;
        case InterpolateOption.rhop:
            rhoPhi = phi; pPhi = phi;
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) { u_modesPhi[imode] = phi; }
            }
            break;
        case InterpolateOption.rhot:
            rhoPhi = phi; TPhi = phi;
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) { T_modesPhi[imode] = phi;
                }
            }
            break;
        } // end switch thermo_interpolator

        return;
    } // end park_limit()

    @nogc
    void store_max_min_values_for_compact_stencil(FluidFVCell[] cell_cloud, ref LocalConfig myConfig)
    {
        // Some of the limiters require the maximum and minimum flowstate variables in the reconstruction stencil

        size_t dimensions = myConfig.dimensions;
        auto np = cell_cloud.length;

        // The following function to be used at compile time.
        string codeForMaxMin(string qname, string qMaxname, string qMinname)
        {
            string code = "{
                "~qMaxname~" = cell_cloud[0].fs."~qname~";
                "~qMinname~" = cell_cloud[0].fs."~qname~";
                foreach (i; 1 .. np) {
                    "~qMaxname~" = fmax("~qMaxname~", cell_cloud[i].fs."~qname~");
                    "~qMinname~" = fmin("~qMinname~", cell_cloud[i].fs."~qname~");
                }
            }
            ";
            return code;
        }

        mixin(codeForMaxMin("vel.x", "velxMax", "velxMin"));
        mixin(codeForMaxMin("vel.y", "velyMax", "velyMin"));
        mixin(codeForMaxMin("vel.z", "velzMax", "velzMin"));
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForMaxMin("B.x", "BxMax", "BxMin"));
                mixin(codeForMaxMin("B.y", "ByMax", "ByMin"));
                mixin(codeForMaxMin("B.z", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForMaxMin("psi", "psiMax", "psiMin"));
                }
            }
        }
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForMaxMin("turb[it]", "turbMax[it]", "turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForMaxMin("gas.rho_s[isp]", "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            mixin(codeForMaxMin("gas.p", "pMax", "pMin"));
            mixin(codeForMaxMin("gas.T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.T_modes[imode]", "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForMaxMin("gas.rho", "rhoMax", "rhoMin"));
            mixin(codeForMaxMin("gas.u", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.u_modes[imode]", "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForMaxMin("gas.rho", "rhoMax", "rhoMin"));
            mixin(codeForMaxMin("gas.p", "pMax", "pMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.u_modes[imode]", "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhot:
            mixin(codeForMaxMin("gas.rho", "rhoMax", "rhoMin"));
            mixin(codeForMaxMin("gas.T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.T_modes[imode]", "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
        return;
    } // end store_max_min_values_for_compact_stencil()

    @nogc
    void store_max_min_values_for_extended_stencil(FluidFVCell[] cell_cloud, ref LocalConfig myConfig)
    {
        // the MLP limiter is unique in that it requires the max/min values from a stencil that includes all cells attached to the vertices of a cell

        size_t dimensions = myConfig.dimensions;
        auto np = cell_cloud.length;

        // The following function to be used at compile time.
        string codeForMaxMin(string qname, string qMaxname, string qMinname)
        {
            string code = "{
                "~qMaxname~" = cell_cloud[0].fs."~qname~";
                "~qMinname~" = cell_cloud[0].fs."~qname~";
                foreach (i, v; cell_cloud[0].vtx) {
                    foreach (c; v.cell_cloud) {
                        "~qMaxname~" = fmax("~qMaxname~", c.fs."~qname~");
                        "~qMinname~" = fmin("~qMinname~", c.fs."~qname~");
                    }
                }
            }
            ";
            return code;
        }

        mixin(codeForMaxMin("vel.x", "velxMax", "velxMin"));
        mixin(codeForMaxMin("vel.y", "velyMax", "velyMin"));
        mixin(codeForMaxMin("vel.z", "velzMax", "velzMin"));
        version(MHD) {
            if (myConfig.MHD) {
                mixin(codeForMaxMin("B.x", "BxMax", "BxMin"));
                mixin(codeForMaxMin("B.y", "ByMax", "ByMin"));
                mixin(codeForMaxMin("B.z", "BzMax", "BzMin"));
                if (myConfig.divergence_cleaning) {
                    mixin(codeForMaxMin("psi", "psiMax", "psiMin"));
                }
            }
        }
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForMaxMin("turb[it]", "turbMax[it]", "turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForMaxMin("gas.rho_s[isp]", "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            mixin(codeForMaxMin("gas.p", "pMax", "pMin"));
            mixin(codeForMaxMin("gas.T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.T_modes[imode]", "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhou:
            mixin(codeForMaxMin("gas.rho", "rhoMax", "rhoMin"));
            mixin(codeForMaxMin("gas.u", "uMax", "uMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.u_modes[imode]", "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForMaxMin("gas.rho", "rhoMax", "rhoMin"));
            mixin(codeForMaxMin("gas.p", "pMax", "pMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.u_modes[imode]", "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
            break;
        case InterpolateOption.rhot:
            mixin(codeForMaxMin("gas.rho", "rhoMax", "rhoMin"));
            mixin(codeForMaxMin("gas.T", "TMax", "TMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForMaxMin("gas.T_modes[imode]", "T_modesMax[imode]", "T_modesMin[imode]"));
                }
            }
            break;
        } // end switch thermo_interpolator
        return;
    } // end store_max_min_values_for_extended_stencil()

    @nogc
    void compute_lsq_values(FluidFVCell[] cell_cloud, ref LSQInterpWorkspace ws,
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
                "~gname~"[0] = 0.0; "~gname~"[1] = 0.0; "~gname~"[2] = 0.0;
                foreach (i; 1 .. np) {
                    number dq = cell_cloud[i].fs."~qname~" - q0;
                    "~gname~"[0] += ws.wx[i] * dq;
                    "~gname~"[1] += ws.wy[i] * dq;
                    if (dimensions == 3) { "~gname~"[2] += ws.wz[i] * dq; }
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
        version(turbulence) {
            foreach (it; 0 .. myConfig.turb_model.nturb) {
                mixin(codeForGradients("turb[it]","turb[it]","turbMax[it]","turbMin[it]"));
            }
        }
        version(multi_species_gas) {
            auto nsp = myConfig.n_species;
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    mixin(codeForGradients("gas.rho_s[isp]", "rho_s[isp]",
                                           "rho_sMax[isp]", "rho_sMin[isp]"));
                }
            } else {
                // Only one possible gradient value for a single species.
                rho_s[0][0] = 0.0; rho_s[0][1] = 0.0; rho_s[0][2] = 0.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        auto nmodes = myConfig.n_modes;
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
            // the Park limiter needs the pressure gradient
            if (myConfig.unstructured_limiter == UnstructuredLimiter.park ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvan_albada ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvenkat ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hnishikawa) {
                mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            }
            break;
        case InterpolateOption.rhop:
            mixin(codeForGradients("gas.rho", "rho", "rhoMax", "rhoMin"));
            mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            version(multi_T_gas) {
                foreach (imode; 0 .. nmodes) {
                    mixin(codeForGradients("gas.u_modes[imode]", "u_modes[imode]",
                                           "u_modesMax[imode]", "u_modesMin[imode]"));
                }
            }
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
            // the Park limiter needs the pressure gradient
            if (myConfig.unstructured_limiter == UnstructuredLimiter.park ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvan_albada ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvenkat ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hnishikawa) {
                mixin(codeForGradients("gas.p", "p", "pMax", "pMin"));
            }
            break;
        } // end switch thermo_interpolator
        return;
    } // end compute_lsq_gradients()

} // end class LSQInterpGradients



    @nogc
    void interp_both(LocalConfig myConfig, ref FVInterface IFace, size_t gtl,
                     ref FlowState Lft, ref FlowState Rght,
                     bool allow_high_order_interpolation)
    // Interpolate the flow field quantities at the left- and right-side of the interface,
    // given information in both cells attached to this interface.
    {
        auto gmodel = myConfig.gmodel;
        auto nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        auto nturb = myConfig.turb_model.nturb;
        FluidFVCell cL0 = IFace.left_cell;
        FluidFVCell cR0 = IFace.right_cell;
        // Low-order reconstruction just copies data from adjacent FV_Cell.
        // Even for high-order reconstruction, we depend upon this copy for
        // the viscous-transport and diffusion coefficients.
        Lft.copy_values_from(cL0.fs);
        Rght.copy_values_from(cR0.fs);
        IFace.fs.copy_average_values_from(Lft, Rght);
        // Enforce first order reconstruction for cells that capure the shocks,
        if (myConfig.suppress_reconstruction_at_shocks && ((Lft.S == 1.0) || (Rght.S == 1.0))) { return; }
        // else apply higher-order interpolation to all faces.
        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            // High-order reconstruction for some properties.
            //
            LSQInterpWorkspace* wsL = cL0.ws;
            LSQInterpWorkspace* wsR = cR0.ws;
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
                    case UnstructuredLimiter.svan_albada:
                        double eps = myConfig.smooth_limiter_coeff/100.0;
                        number qqL = fmax(1e-12, fabs(qL0));
                        number qqR = fmax(1e-12, fabs(qR0));
                        number qq = fmax(qqL, qqR);

                        number dqL = -4.0*dLx*mygradL[0] + -4.0*dLy*mygradL[1];
                        number dqR = -4.0*dRx*mygradR[0] + -4.0*dRy*mygradR[1];
                        if (myConfig.dimensions == 3) {
                            dqL += -4.0*dLz*mygradL[2];
                            dqR += -4.0*dRz*mygradR[2];
                        }

                        number qL1 = qR0 + dqL;
                        number qR1 = qL0 + dqR;

                        number delLminus = (qL0 - qL1);
                        number del = (qR0 - qL0);
                        number delRplus = (qR1 - qR0);

                        // val Albada limiter, modified with the non-dimensional
                        // smoothing parameter qq*eps. See notes 15/11/22 (NNG)
                        number sL = (delLminus*del + fabs(delLminus*del) + qq*eps) /
                                    (delLminus*delLminus + del*del + qq*eps);
                        number sR = (del*delRplus + fabs(del*delRplus) + qq*eps) /
                                    (del*del + delRplus*delRplus + qq*eps);

                        mygradL[0] *= sL;
                        mygradL[1] *= sL;
                        mygradL[2] *= sL;

                        mygradR[0] *= sR;
                        mygradR[1] *= sR;
                        mygradR[2] *= sR;
                        break;
                    case UnstructuredLimiter.min_mod:
                        min_mod_limit(mygradL[0], mygradR[0]);
                        min_mod_limit(mygradL[1], mygradR[1]);
                        min_mod_limit(mygradL[2], mygradR[2]);
                        break;
                    case UnstructuredLimiter.barth:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.park:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvan_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.van_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hnishikawa:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.nishikawa:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvenkat_mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat_mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvenkat:
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
            mixin(codeForReconstruction("vel.x", "velx", "vel.x", "velxPhi"));
            mixin(codeForReconstruction("vel.y", "vely", "vel.y", "velyPhi"));
            mixin(codeForReconstruction("vel.z", "velz", "vel.z", "velzPhi"));
            version(MHD) {
                if (myConfig.MHD) {
                    mixin(codeForReconstruction("B.x", "Bx", "B.x", "BxPhi"));
                    mixin(codeForReconstruction("B.y", "By", "B.y", "ByPhi"));
                    mixin(codeForReconstruction("B.z", "Bz", "B.z", "BxPhi"));
                    if (myConfig.divergence_cleaning) {
                        mixin(codeForReconstruction("psi", "psi", "psi", "psiPhi"));
                    }
                }
            }
            version(turbulence) {
                if (myConfig.allow_reconstruction_for_turbulent_variables) {
                    foreach (it; 0 .. myConfig.turb_model.nturb) {
                        mixin(codeForReconstruction("turb[it]","turb[it]","turb[it]","turbPhi[it]"));
                    }
                } else {
                    Lft.turb[] = IFace.left_cell.fs.turb[];
                    Rght.turb[] = IFace.right_cell.fs.turb[];
                }
            }
            version(multi_species_gas) {
                if (nsp > 1) {
                    // Reconstruct species densities
                    if (myConfig.allow_reconstruction_for_species) {
                        foreach (isp; 0 .. nsp) {
                            mixin(codeForReconstruction("gas.rho_s[isp]", "rho_s[isp]",
                                                        "gas.rho_s[isp]", "rho_sPhi[isp]"));
                        }
                    } else {
                        Lft.gas.rho_s[] = IFace.left_cell.fs.gas.rho_s[];
                        Rght.gas.rho_s[] = IFace.right_cell.fs.gas.rho_s[];
                    }
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
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                        "gas.T_modes[imode]", "T_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.T_modes[] = IFace.left_cell.fs.gas.T_modes[];
                        Rght.gas.T_modes[] = IFace.right_cell.fs.gas.T_modes[];
                    }
                }
                mixin(codeForThermoUpdate("pT"));
                version(multi_species_gas) {
                    if (nsp > 1) {
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp]  = Lft.gas.rho_s[isp]/Lft.gas.rho;
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        Lft.gas.massf[0]  = 1.0;
                        Rght.gas.massf[0] = 1.0;
                    }
                } else {
                    Lft.gas.massf[0]  = 1.0;
                    Rght.gas.massf[0] = 1.0;
                }
                break;
            case InterpolateOption.rhou:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_L = 0.0;
                        number rho_R = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_L += Lft.gas.rho_s[isp];
                            rho_R += Rght.gas.rho_s[isp];
                        }
                        Lft.gas.rho  = rho_L;
                        Rght.gas.rho = rho_R;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp] = Lft.gas.rho_s[isp]/Lft.gas.rho;
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.u", "u", "gas.u", "uPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                        "gas.u_modes[imode]", "u_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.u_modes[] = IFace.left_cell.fs.gas.u_modes[];
                        Rght.gas.u_modes[] = IFace.right_cell.fs.gas.u_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhou"));
                break;
            case InterpolateOption.rhop:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_L = 0.0;
                        number rho_R = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_L += Lft.gas.rho_s[isp];
                            rho_R += Rght.gas.rho_s[isp];
                        }
                        Lft.gas.rho  = rho_L;
                        Rght.gas.rho = rho_R;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp] = Lft.gas.rho_s[isp]/Lft.gas.rho;
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                        "gas.u_modes[imode]", "u_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.u_modes[] = IFace.left_cell.fs.gas.u_modes[];
                        Rght.gas.u_modes[] = IFace.right_cell.fs.gas.u_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhop"));
                break;
            case InterpolateOption.rhot:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_L = 0.0;
                        number rho_R = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_L += Lft.gas.rho_s[isp];
                            rho_R += Rght.gas.rho_s[isp];
                        }
                        Lft.gas.rho  = rho_L;
                        Rght.gas.rho = rho_R;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp] = Lft.gas.rho_s[isp]/Lft.gas.rho;
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                        "gas.T_modes[imode]", "T_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.T_modes[] = IFace.left_cell.fs.gas.T_modes[];
                        Rght.gas.T_modes[] = IFace.right_cell.fs.gas.T_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhoT"));
                break;
            } // end switch thermo_interpolator
        } // end of high-order reconstruction
        return;
    } // end interp_both()

    @nogc
    void interp_right(LocalConfig myConfig, ref FVInterface IFace, size_t gtl, ref FlowState Rght,
                      bool allow_high_order_interpolation)
    // Interpolate the flow field quantities at the right-side of the interface,
    // given information in right-cell attached to this interface.
    {
        auto gmodel = myConfig.gmodel;
        auto nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        auto nturb = myConfig.turb_model.nturb;
        FluidFVCell cR0 = IFace.right_cell;
        // Low-order reconstruction just copies data from adjacent FV_Cell.
        // Even for high-order reconstruction, we depend upon this copy for
        // the viscous-transport and diffusion coefficients.
        Rght.copy_values_from(cR0.fs);
        IFace.fs.copy_values_from(Rght);
        // Enforce first order reconstruction for cells that capure the shocks,
        if (myConfig.suppress_reconstruction_at_shocks && (Rght.S == 1.0)) { return; }
        // else apply higher-order interpolation to all faces.
        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            // High-order reconstruction for some properties.
            //
            LSQInterpWorkspace* wsR = cR0.ws;
            // vector from left-cell-centre to face midpoint
            number dRx = IFace.pos.x - cR0.pos[gtl].x;
            number dRy = IFace.pos.y - cR0.pos[gtl].y;
            number dRz = IFace.pos.z - cR0.pos[gtl].z;
            number[3] mygradR;
            //
            // Always reconstruct in the global frame of reference -- for now
            // At a boundary, don't apply limiter unless extrema_clipping also.
            //
            string codeForReconstruction(string qname, string gname,
                                         string tname, string lname)
            {
                string code = "{
                number qR0 = cR0.fs."~qname~";
                number qMinR = qR0;
                number qMaxR = qR0;
                mygradR[0] = cR0.gradients."~gname~"[0];
                mygradR[1] = cR0.gradients."~gname~"[1];
                mygradR[2] = cR0.gradients."~gname~"[2];
                if (myConfig.apply_limiter && myConfig.extrema_clipping) {
                    final switch (myConfig.unstructured_limiter) {
                    case UnstructuredLimiter.svan_albada:
                        break;
                    case UnstructuredLimiter.min_mod:
                        break;
                    case UnstructuredLimiter.barth:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.park:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvan_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.van_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hnishikawa:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.nishikawa:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvenkat_mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat_mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvenkat:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat:
                        mygradR[0] *= cR0.gradients."~lname~";
                        mygradR[1] *= cR0.gradients."~lname~";
                        mygradR[2] *= cR0.gradients."~lname~";
                        break;
                    }
                }
                number qR = qR0 + dRx*mygradR[0] + dRy*mygradR[1];
                if (myConfig.dimensions == 3) { qR += dRz*mygradR[2]; }
                Rght."~tname~" = qR;
                }
                ";
                return code;
            }
            mixin(codeForReconstruction("vel.x", "velx", "vel.x", "velxPhi"));
            mixin(codeForReconstruction("vel.y", "vely", "vel.y", "velyPhi"));
            mixin(codeForReconstruction("vel.z", "velz", "vel.z", "velzPhi"));
            version(MHD) {
                if (myConfig.MHD) {
                    mixin(codeForReconstruction("B.x", "Bx", "B.x", "BxPhi"));
                    mixin(codeForReconstruction("B.y", "By", "B.y", "ByPhi"));
                    mixin(codeForReconstruction("B.z", "Bz", "B.z", "BxPhi"));
                    if (myConfig.divergence_cleaning) {
                        mixin(codeForReconstruction("psi", "psi", "psi", "psiPhi"));
                    }
                }
            }
            version(turbulence) {
                if (myConfig.allow_reconstruction_for_turbulent_variables) {
                    foreach (it; 0 .. myConfig.turb_model.nturb) {
                        mixin(codeForReconstruction("turb[it]","turb[it]","turb[it]","turbPhi[it]"));
                    }
                } else {
                    Rght.turb[] = IFace.right_cell.fs.turb[];
                }
            }
            version(multi_species_gas) {
                if (nsp > 1) {
                    // Reconstruct species densities
                    if (myConfig.allow_reconstruction_for_species) {
                        foreach (isp; 0 .. nsp) {
                            mixin(codeForReconstruction("gas.rho_s[isp]", "rho_s[isp]",
                                                        "gas.rho_s[isp]", "rho_sPhi[isp]"));
                        }
                    } else {
                        Rght.gas.rho_s[] = IFace.right_cell.fs.gas.rho_s[];
                    }
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
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                        "gas.T_modes[imode]", "T_modesPhi[imode]"));
                        }
                    } else {
                        Rght.gas.T_modes[] = IFace.right_cell.fs.gas.T_modes[];
                    }
                }
                mixin(codeForThermoUpdate("pT"));
                version(multi_species_gas) {
                    if (nsp > 1) {
                        foreach (isp; 0 .. nsp) {
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        Rght.gas.massf[0] = 1.0;
                    }
                } else {
                    Rght.gas.massf[0] = 1.0;
                }
                break;
            case InterpolateOption.rhou:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_R = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_R += Rght.gas.rho_s[isp];
                        }
                        Rght.gas.rho = rho_R;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.u", "u", "gas.u", "uPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                        "gas.u_modes[imode]", "u_modesPhi[imode]"));
                        }
                    } else {
                        Rght.gas.u_modes[] = IFace.right_cell.fs.gas.u_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhou"));
                break;
            case InterpolateOption.rhop:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_R = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_R += Rght.gas.rho_s[isp];
                        }
                        Rght.gas.rho = rho_R;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                        "gas.u_modes[imode]", "u_modesPhi[imode]"));
                        }
                    } else {
                        Rght.gas.u_modes[] = IFace.right_cell.fs.gas.u_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhop"));
                break;
            case InterpolateOption.rhot:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_R = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_R += Rght.gas.rho_s[isp];
                        }
                        Rght.gas.rho = rho_R;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Rght.gas.massf[isp] = Rght.gas.rho_s[isp]/Rght.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Rght.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                        "gas.T_modes[imode]", "T_modesPhi[imode]"));
                        }
                    } else {
                        Rght.gas.T_modes[] = IFace.right_cell.fs.gas.T_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhoT"));
                break;
            } // end switch thermo_interpolator
        } // end of high-order reconstruction
        return;
    } // end interp_right()

    @nogc
    void interp_left(LocalConfig myConfig, ref FVInterface IFace, size_t gtl, ref FlowState Lft,
                     bool allow_high_order_interpolation)
    // Interpolate the flow field quantities at the left-side of the interface,
    // given information in the left-cell attached to this interface.
    {
        auto gmodel = myConfig.gmodel;
        auto nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        auto nturb = myConfig.turb_model.nturb;
        FluidFVCell cL0 = IFace.left_cell;
        // Low-order reconstruction just copies data from adjacent FV_Cell.
        // Even for high-order reconstruction, we depend upon this copy for
        // the viscous-transport and diffusion coefficients.
        Lft.copy_values_from(cL0.fs);
        IFace.fs.copy_values_from(Lft);
        // Enforce first order reconstruction for cells that capure the shocks,
        if (myConfig.suppress_reconstruction_at_shocks && (Lft.S == 1.0)) { return; }
        // else apply higher-order interpolation to all faces.
        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            // High-order reconstruction for some properties.
            //
            LSQInterpWorkspace* wsL = cL0.ws;
            // vector from left-cell-centre to face midpoint
            number dLx = IFace.pos.x - cL0.pos[gtl].x;
            number dLy = IFace.pos.y - cL0.pos[gtl].y;
            number dLz = IFace.pos.z - cL0.pos[gtl].z;
            number[3] mygradL;
            //
            // Always reconstruct in the global frame of reference -- for now
            // At a boundary, don't apply limiter unless extrema_clipping also.
            //
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
                if (myConfig.apply_limiter && myConfig.extrema_clipping) {
                    final switch (myConfig.unstructured_limiter) {
                    case UnstructuredLimiter.svan_albada:
                        break;
                    case UnstructuredLimiter.min_mod:
                        break;
                    case UnstructuredLimiter.barth:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.park:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvan_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.van_albada:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hnishikawa:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.nishikawa:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvenkat_mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat_mlp:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.hvenkat:
                        goto case UnstructuredLimiter.venkat;
                    case UnstructuredLimiter.venkat:
                        mygradL[0] *= cL0.gradients."~lname~";
                        mygradL[1] *= cL0.gradients."~lname~";
                        mygradL[2] *= cL0.gradients."~lname~";
                        break;
                    }
                }
                number qL = qL0 + dLx*mygradL[0] + dLy*mygradL[1];
                if (myConfig.dimensions == 3) { qL += dLz*mygradL[2]; }
                Lft."~tname~" = qL;
                }
                ";
                return code;
            }
            mixin(codeForReconstruction("vel.x", "velx", "vel.x", "velxPhi"));
            mixin(codeForReconstruction("vel.y", "vely", "vel.y", "velyPhi"));
            mixin(codeForReconstruction("vel.z", "velz", "vel.z", "velzPhi"));
            version(MHD) {
                if (myConfig.MHD) {
                    mixin(codeForReconstruction("B.x", "Bx", "B.x", "BxPhi"));
                    mixin(codeForReconstruction("B.y", "By", "B.y", "ByPhi"));
                    mixin(codeForReconstruction("B.z", "Bz", "B.z", "BxPhi"));
                    if (myConfig.divergence_cleaning) {
                        mixin(codeForReconstruction("psi", "psi", "psi", "psiPhi"));
                    }
                }
            }
            version(turbulence) {
                if (myConfig.allow_reconstruction_for_turbulent_variables) {
                    foreach (it; 0 .. myConfig.turb_model.nturb) {
                        mixin(codeForReconstruction("turb[it]","turb[it]","turb[it]","turbPhi[it]"));
                    }
                } else {
                    Lft.turb[] = IFace.left_cell.fs.turb[];
                }
            }
            version(multi_species_gas) {
                if (nsp > 1) {
                    // Reconstruct species densities
                    if (myConfig.allow_reconstruction_for_species) {
                        foreach (isp; 0 .. nsp) {
                            mixin(codeForReconstruction("gas.rho_s[isp]", "rho_s[isp]",
                                                        "gas.rho_s[isp]", "rho_sPhi[isp]"));
                        }
                    } else {
                        Lft.gas.rho_s[] = IFace.left_cell.fs.gas.rho_s[];
                    }
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
                ";
                return code;
            }
            final switch (myConfig.thermo_interpolator) {
            case InterpolateOption.pt:
                mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                        "gas.T_modes[imode]", "T_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.T_modes[] = IFace.left_cell.fs.gas.T_modes[];
                    }
                }
                mixin(codeForThermoUpdate("pT"));
                version(multi_species_gas) {
                    if (nsp > 1) {
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp]  = Lft.gas.rho_s[isp]/Lft.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                        }
                    } else {
                        Lft.gas.massf[0]  = 1.0;
                    }
                } else {
                    Lft.gas.massf[0]  = 1.0;
                }
                break;
            case InterpolateOption.rhou:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_L = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_L += Lft.gas.rho_s[isp];
                        }
                        Lft.gas.rho  = rho_L;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp] = Lft.gas.rho_s[isp]/Lft.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.u", "u", "gas.u", "uPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                        "gas.u_modes[imode]", "u_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.u_modes[] = IFace.left_cell.fs.gas.u_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhou"));
                break;
            case InterpolateOption.rhop:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_L = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_L += Lft.gas.rho_s[isp];
                        }
                        Lft.gas.rho  = rho_L;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp] = Lft.gas.rho_s[isp]/Lft.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.p", "p", "gas.p", "pPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.u_modes[imode]", "u_modes[imode]",
                                                        "gas.u_modes[imode]", "u_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.u_modes[] = IFace.left_cell.fs.gas.u_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhop"));
                break;
            case InterpolateOption.rhot:
                version(multi_species_gas) {
                    if (nsp > 1) {
                        // compute total density as a sum of species densities
                        number rho_L = 0.0;
                        foreach (isp; 0 .. nsp) {
                            rho_L += Lft.gas.rho_s[isp];
                        }
                        Lft.gas.rho  = rho_L;
                        // compute mass fractions from total density and species densities
                        foreach (isp; 0 .. nsp) {
                            Lft.gas.massf[isp] = Lft.gas.rho_s[isp]/Lft.gas.rho;
                        }
                        if (myConfig.scale_species_after_reconstruction) {
                            scale_mass_fractions(Lft.gas.massf);
                        }
                    } else {
                        mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                    }
                } else {
                    mixin(codeForReconstruction("gas.rho", "rho", "gas.rho", "rhoPhi"));
                }
                mixin(codeForReconstruction("gas.T", "T", "gas.T", "TPhi"));
                version(multi_T_gas) {
                    if (myConfig.allow_reconstruction_for_energy_modes) {
                        foreach (imode; 0 .. nmodes) {
                            mixin(codeForReconstruction("gas.T_modes[imode]", "T_modes[imode]",
                                                        "gas.T_modes[imode]", "T_modesPhi[imode]"));
                        }
                    } else {
                        Lft.gas.T_modes[] = IFace.left_cell.fs.gas.T_modes[];
                    }
                }
                mixin(codeForThermoUpdate("rhoT"));
                break;
            } // end switch thermo_interpolator
        } // end of high-order reconstruction
        return;
    } // end interp_left()

