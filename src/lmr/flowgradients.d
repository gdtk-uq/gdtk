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

module lmr.flowgradients;

import std.conv;
import std.math;
import std.stdio;

import gas;
import geom;
import nm.number;
import nm.rsla;
import ntypes.complex;

import lmr.conservedquantities;
import lmr.flowstate;
import lmr.fvinterface;
import lmr.fvvertex;
import lmr.globalconfig;
import lmr.mass_diffusion;

immutable size_t cloud_nmax = 12;

struct WLSQGradWorkspace {
public:
    // A place to hold the intermediate results for computing
    // the least-squares model as a weighted sum of the flow data.
    number[cloud_nmax] wx, wy, wz;
    bool compute_about_mid;
    size_t loop_init; // starting index for loops:
    // 0=compute_about_mid, 1=compute_about_[0]
    size_t n; // cloud_pos.length;

    this(in WLSQGradWorkspace other)
    {
        wx[] = other.wx[]; wy[] = other.wy[]; wz[] = other.wz[];
        compute_about_mid = other.compute_about_mid;
        loop_init = other.loop_init;
        n = other.n;
    }
} // end struct WLSQGradWorkspace


struct FlowGradients {
    // Spatial derivatives of the flow quantities
public:
    number[3][3] vel;
    // velocity derivatives stored as a second-order tensor
    // [[du/dx du/dy du/dz]
    //  [dv/dx dv/dy dv/dz]
    //  [dw/dx dw/dy dw/dz]]
    version(multi_species_gas) {
        number[3][] massf; // mass fraction derivatives
    }
    number[3] T; // Temperature derivatives
    version(multi_T_gas) {
        number[3][] T_modes;
    }
    version(turbulence) {
        number[3][2] turb; // turbulence primitive variables
    }

    number[3] p_e; // electron pressure

    this(in LocalConfig myConfig)
    {
        // Note that the gradient values are initialized to zero.
        // This is because some of the updates (e.g. the point-implicit update)
        // will pick up some gradient values before thay have been computed and
        // add them to a tally for taking an average.
        // On the first iteration, having initial Nan values would poison
        // that average but zero values are tolerable.
        // On later iterations, the gradient value from a previous iteration is fine.
        foreach (i; 0 .. 3) {
            foreach (j; 0 .. 3) { vel[i][j] = to!number(0.0); }
            T[i] = to!number(0.0);

            if (myConfig.gmodel.is_plasma){
                p_e[i] = to!number(0.0);
            }
        }
        version(multi_species_gas) {
            massf.length = myConfig.n_species;
            foreach (i; 0 .. myConfig.n_species) {
                foreach (j; 0 .. 3) { massf[i][j] = 0.0; }
            }
        }
        version(multi_T_gas) {
            T_modes.length = myConfig.n_modes;
            foreach (i; 0 .. myConfig.n_modes) {
                foreach (j; 0 .. 3) { T_modes[i][j] = to!number(0.0); }
            }
        }
        version(turbulence) {
            foreach (i; 0 .. myConfig.turb_model.nturb) {
                foreach (j; 0 .. 3) { turb[i][j] = to!number(0.0); }
            }
        }
    }

    this(in FlowGradients other)
    {
        foreach(i; 0 .. 3) vel[i][] = other.vel[i][];
        version(multi_species_gas) {
            massf.length = other.massf.length;
            foreach(isp; 0 .. other.massf.length) { massf[isp][] = other.massf[isp][]; }
        }
        T[] = other.T[];
        version(multi_T_gas) {
            T_modes.length = other.T_modes.length;
            foreach (imode; 0 .. other.T_modes.length) { T_modes[imode][] = other.T_modes[imode][]; }
        }
        version(turbulence) {
            foreach(i; 0 .. turb.length) turb[i][] = other.turb[i][];
        }
    }

    @nogc
    void copy_values_from(in FlowGradients other)
    {
        foreach (i; 0 .. 3) { vel[i][] = other.vel[i][]; }
        version(multi_species_gas) {
            foreach (isp; 0 .. other.massf.length) { massf[isp][] = other.massf[isp][]; }
        }
        T[] = other.T[];
        version(multi_T_gas) {
            foreach (imode; 0 .. other.T_modes.length) { T_modes[imode][] = other.T_modes[imode][]; }
        }
        version(turbulence) {
            foreach(i; 0 .. turb.length) turb[i][] = other.turb[i][];
        }
    }

    @nogc
    void accumulate_values_from(in FlowGradients other)
    {
        foreach (i; 0 .. 3) {
            foreach (j; 0 .. 3) vel[i][j] += other.vel[i][j];
        }
        version(multi_species_gas) {
            foreach (isp; 0 .. massf.length) { 
                foreach (j; 0 .. 3) massf[isp][j] += other.massf[isp][j];
            }
        }
        foreach (j; 0 .. 3) T[j] += other.T[j];
        version(multi_T_gas) {
            foreach (imode; 0 .. T_modes.length) { 
                foreach (j; 0 .. 3) T_modes[imode][j] += other.T_modes[imode][j];
            }
        }
        version(turbulence) {
            foreach(i; 0 .. turb.length) {
                foreach (j; 0 .. 3) turb[i][j] += other.turb[i][j];
            }
        }
    }

    @nogc
    void accumulate_values_from(in FlowGradients other, number factor)
    {
        foreach (i; 0 .. 3) {
            foreach (j; 0 .. 3) vel[i][j] += other.vel[i][j] * factor;
        }
        version(multi_species_gas) {
            foreach (isp; 0 .. massf.length) {
                foreach (j; 0 .. 3) massf[isp][j] += other.massf[isp][j] * factor;
            }
        }
        foreach (j; 0 .. 3) T[j] += other.T[j] * factor;
        version(multi_T_gas) {
            foreach (imode; 0 .. T_modes.length) {
                foreach (j; 0 .. 3) T_modes[imode][j] += other.T_modes[imode][j] * factor;
            }
        }
        version(turbulence) {
            foreach(i; 0 .. turb.length) {
                foreach (j; 0 .. 3) turb[i][j] += other.turb[i][j] * factor;
            }
        }
    }

    @nogc
    void scale_values_by(number factor)
    {
        foreach (i; 0 .. 3) { vel[i][] *= factor; }
        version(multi_species_gas) {
            foreach (isp; 0 .. massf.length) { massf[isp][] *= factor; }
        }
        T[] *= factor;
        version(multi_T_gas) {
            foreach (imode; 0 .. T_modes.length) { T_modes[imode][] *= factor; }
        }
        version(turbulence) {
            foreach(i; 0 .. turb.length) turb[i][] *= factor;
        }
    }

    string toString() const
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
        version(multi_species_gas) {
            repr ~= ", massf=[" ~ to!string(massf[0]);
            foreach (i; 1 .. massf.length) repr ~= ", " ~ to!string(massf);
            repr ~= "]";
        }
        repr ~= ", T=" ~ to!string(T);
        version(multi_T_gas) {
            repr ~= ", T_modes=[";
            if (T_modes.length > 0) repr ~= to!string(T_modes[0]);
            foreach (i; 1 .. T_modes.length) repr ~= ", " ~ to!string(T_modes);
            repr ~= "]";
        }
        version(turbulence) {
            repr ~= ", turb=[";
            if (turb.length>0) repr ~= to!string(turb[0]);
            foreach (i; 1 .. turb.length) repr ~= ", " ~ to!string(turb[i]);
            repr ~= "]";
        }
        repr ~= ")";
        return to!string(repr);
    }

    @nogc
    void gradients_xy_div(in LocalConfig myConfig, ref FlowState*[] cloud_fs, ref Vector3*[] cloud_pos)
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
        number areaxy = (cloud_pos[0].x + cloud_pos[n-1].x) *
            (cloud_pos[0].y - cloud_pos[n-1].y);
        // Accumulate the contributions from the other segments.
        foreach (i; 0 .. n-1) {
            areaxy += (cloud_pos[i+1].x + cloud_pos[i].x) *
                (cloud_pos[i+1].y - cloud_pos[i].y);
        }
        number area_inv = 1.0 / areaxy;
        //
        // Apply the divergence theorem to flow variables, generating code for each.
        //
        number gradient_x, gradient_y;
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
        mixin(codeForGradients("gas.T"));
        T[0] = gradient_x * area_inv;
        T[1] = -gradient_y * area_inv;
        T[2] = 0.0;
        //
        version(multi_T_gas) {
            size_t nmodes = cloud_fs[0].gas.T_modes.length;
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]"));
                T_modes[imode][0] = gradient_x * area_inv;
                T_modes[imode][1] = -gradient_y * area_inv;
                T_modes[imode][2] = 0.0;
            }
        }
        //
        version(multi_species_gas) {
            size_t nsp = cloud_fs[0].gas.massf.length;
            if (myConfig.turb_model.isTurbulent ||
                myConfig.mass_diffusion_model != MassDiffusionModel.none) {
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
        }
        //
        version(turbulence) {
            foreach(tidx; 0 .. myConfig.turb_model.nturb){
                mixin(codeForGradients("turb[tidx]"));
                turb[tidx][0] = gradient_x * area_inv;
                turb[tidx][1] = -gradient_y * area_inv;
                turb[tidx][2] = 0.0;
            }
            //
        }
    } // end gradients_xy_div()

    @nogc
    void set_up_workspace_leastsq_via_normal(in LocalConfig myConfig, ref Vector3*[] cloud_pos, ref Vector3 pos,
                                             bool compute_about_mid, ref WLSQGradWorkspace ws)
    {
        assert(&ws !is null, "We are missing the workspace!");
        size_t n = cloud_pos.length;
        assert(n <= cloud_nmax, "Too many points in cloud.");
        number[cloud_nmax] weights2;
        ws.compute_about_mid = compute_about_mid;
        size_t loop_init;
        if (compute_about_mid) {
            // We will compute differences about a "mid" point that is assumed
            // to be distinct from all points in the incoming cloud.
            loop_init = 0; // All points count.
        } else {
            loop_init = 1; // first point is reference for differences.
            weights2[0] = 0.0; // and doesn't enter into the sum itself.
        }
        if (n < 3) {
            string msg = "Not enough points in cloud.";
            debug {
                import std.format;
                msg ~= format("\n  pos=%s", pos);
                msg ~= format("\n  compute_about_mid=%s, loop_init=%d, n=%d", compute_about_mid, loop_init, n);
                foreach (i; 0 .. n) { msg ~= format("\n  cloud_pos[%d]=%s", i, cloud_pos[i]); }
            }
            throw new FlowSolverException(msg);
        }
        ws.n = n;
        ws.loop_init = loop_init;
        //
        // Calculate weights used in the least-squares gradient calculation.
        // These are the square of the weights on the original linear constraint eqns.
        // These weights are calculated with the current interface/vertex
        // as the reference point (pos).
        // For the "faces" spatial location we are expecting the primary point
        // (i.e. the face at which we are calculating the gradients) to be in
        // the first cloud position.
        number x0 = pos.x; number y0 = pos.y; number z0 = pos.z;
        bool apply_weighting = false;
        if (myConfig.viscous_least_squares_type == ViscousLeastSquaresType.weighted_normal) {
            apply_weighting = true;
        }
        if (myConfig.dimensions == 2) {
            foreach (i; loop_init .. n) {
                number dx = cloud_pos[i].x - x0;
                number dy = cloud_pos[i].y - y0;
                if (apply_weighting) {
                    weights2[i] = 1.0/(dx*dx+dy*dy);
                } else {
                    weights2[i] = 1.0;
                }
            }
        } else { //3D
            foreach (i; loop_init .. n) {
                number dx = cloud_pos[i].x - x0;
                number dy = cloud_pos[i].y - y0;
                number dz = cloud_pos[i].z - z0;
                if (apply_weighting) {
                    weights2[i] = 1.0/(dx*dx+dy*dy+dz*dz);
                } else {
                    weights2[i] = 1.0;
                }
            }
        }
        // If computing about a mid-point, calculate that mid-point.
        if (compute_about_mid) {
            x0 = 0.0; y0 = 0.0; z0 = 0.0;
            foreach (i; 0 .. n) {
                x0 += cloud_pos[i].x; y0 += cloud_pos[i].y; z0 += cloud_pos[i].z;
            }
            x0 /= n; y0 /= n; z0 /= n; // midpoint
        } else { // else use the primary point (assumed to be in cloud position 0)
            x0 = cloud_pos[0].x; y0 = cloud_pos[0].y; z0 = cloud_pos[0].z;
        }
        number[cloud_nmax] dx, dy, dz;
        //
        // Assemble and invert the normal matrix.
        // We'll reuse the resulting inverse for each flow-field quantity.
        if (myConfig.dimensions == 3) {
            number[6][3] xTx; // normal matrix, augmented to give 6 entries per row
            number xx = 0.0; number xy = 0.0; number xz = 0.0;
            number yy = 0.0; number yz = 0.0; number zz = 0.0;
            foreach (i; loop_init .. n) {
                dx[i] = cloud_pos[i].x - x0;
                dy[i] = cloud_pos[i].y - y0;
                dz[i] = cloud_pos[i].z - z0;
                xx += weights2[i]*dx[i]*dx[i];
                xy += weights2[i]*dx[i]*dy[i];
                xz += weights2[i]*dx[i]*dz[i];
                yy += weights2[i]*dy[i]*dy[i];
                yz += weights2[i]*dy[i]*dz[i];
                zz += weights2[i]*dz[i]*dz[i];
            }
            xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
            xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
            xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
            xTx[0][3] = 1.0; xTx[0][4] = 0.0; xTx[0][5] = 0.0;
            xTx[1][3] = 0.0; xTx[1][4] = 1.0; xTx[1][5] = 0.0;
            xTx[2][3] = 0.0; xTx[2][4] = 0.0; xTx[2][5] = 1.0;
            double very_small_value = 1.0e-32*(normInf!(3,3,6,number)(xTx).re)^^3;
            if (0 != computeInverse!(3,3,6,number)(xTx, very_small_value)) {
                string msg = "Failed to invert LSQ normal matrix";
                debug {
                    import std.format;
                    msg ~= format("\n  pos=%s", pos);
                    msg ~= format("\n  compute_about_mid=%s, loop_init=%d, n=%d", compute_about_mid, loop_init, n);
                    foreach (i; 0 .. n) { msg ~= format("\n  cloud_pos[%d]=%s", i, cloud_pos[i]); }
                }
                throw new FlowSolverException(msg);
            }
            // Prepare final weights for later use in the reconstruction phase.
            foreach (i; loop_init .. n) {
                ws.wx[i] = xTx[0][3]*dx[i] + xTx[0][4]*dy[i] + xTx[0][5]*dz[i];
                ws.wx[i] *= weights2[i];
                ws.wy[i] = xTx[1][3]*dx[i] + xTx[1][4]*dy[i] + xTx[1][5]*dz[i];
                ws.wy[i] *= weights2[i];
                ws.wz[i] = xTx[2][3]*dx[i] + xTx[2][4]*dy[i] + xTx[2][5]*dz[i];
                ws.wz[i] *= weights2[i];
            }
        } else {
            // dimensions == 2
            number[4][2] xTx; // normal matrix, augmented to give 4 entries per row
            number xx = 0.0; number xy = 0.0; number yy = 0.0;
            foreach (i; loop_init .. n) {
                dx[i] = cloud_pos[i].x - x0;
                dy[i] = cloud_pos[i].y - y0;
                xx += weights2[i]*dx[i]*dx[i];
                xy += weights2[i]*dx[i]*dy[i];
                yy += weights2[i]*dy[i]*dy[i];
            }
            xTx[0][0] = xx; xTx[0][1] = xy;
            xTx[1][0] = xy; xTx[1][1] = yy;
            xTx[0][2] = 1.0; xTx[0][3] = 0.0;
            xTx[1][2] = 0.0; xTx[1][3] = 1.0;
            double very_small_value = 1.0e-32*(normInf!(2,2,4,number)(xTx).re)^^2;
            if (0 != computeInverse!(2,2,4,number)(xTx, very_small_value)) {
                string msg = "Failed to invert LSQ normal matrix";
                debug {
                    import std.format;
                    msg ~= format("\n  pos=%s", pos);
                    msg ~= format("\n  compute_about_mid=%s, loop_init=%d, n=%d", compute_about_mid, loop_init, n);
                    foreach (i; 0 .. n) { msg ~= format("\n  cloud_pos[%d]=%s", i, cloud_pos[i]); }
                }
                throw new FlowSolverException(msg);
            }
            // Prepare final weights for later use in the reconstruction phase.
            foreach (i; loop_init .. n) {
                ws.wx[i] = xTx[0][2]*dx[i] + xTx[0][3]*dy[i];
                ws.wx[i] *= weights2[i];
                ws.wy[i] = xTx[1][2]*dx[i] + xTx[1][3]*dy[i];
                ws.wy[i] *= weights2[i];
                ws.wz[i] = 0.0;
            }
        }
    } // end set_up_workspace_leastsq_via_normal()

    @nogc
    void set_up_workspace_leastsq_via_qr_factorization(in LocalConfig myConfig, ref Vector3*[] cloud_pos, ref Vector3 pos,
                                                       bool compute_about_mid, ref WLSQGradWorkspace ws)
    {
        /**
         *
         * The following algorithm is taken from ref. [1].
         * We use the inverse distance weighting as suggested in ref. [2] with a weighting coefficient of 0.75.
         *
         * references:
         * [1] A. Haselbacher and J. Blazek,
         *     Accurate and efficient discretization of Navier-Stokes equations on mixed grids,
         *     AIAA Journal, vol. 38, no. 11, 2000
         * [2] H. Nishikawa and J. White,
         *     An efficient cell-centered finite-volume method with face-averaged nodal-gradients for triangular grids,
         *     Journal of Computational Physics, vol. 411, 2020
         *
         * @author: Kyle A. Damm (2024-07-25)
         *
         **/

        assert(&ws !is null, "We are missing the workspace!");
        size_t n = cloud_pos.length;
        assert(n <= cloud_nmax, "Too many points in cloud.");
        number[cloud_nmax] theta;
        ws.compute_about_mid = compute_about_mid;
        size_t loop_init;
        if (compute_about_mid) {
            // We will compute differences about a "mid" point that is assumed
            // to be distinct from all points in the incoming cloud.
            loop_init = 0; // All points count.
        } else {
            loop_init = 1; // first point is reference for differences.
            theta[0] = 0.0; // and doesn't enter into the sum itself.
        }
        if (n < 3) {
            string msg = "Not enough points in cloud.";
            debug {
                import std.format;
                msg ~= format("\n  pos=%s", pos);
                msg ~= format("\n  compute_about_mid=%s, loop_init=%d, n=%d", compute_about_mid, loop_init, n);
                foreach (i; 0 .. n) { msg ~= format("\n  cloud_pos[%d]=%s", i, cloud_pos[i]); }
            }
            throw new FlowSolverException(msg);
        }
        ws.n = n;
        ws.loop_init = loop_init;

        // set the weighting coefficient (this is w from eq. 21 from ref. [2] with p = 0.75)
        number x0 = pos.x; number y0 = pos.y; number z0 = pos.z;
        bool apply_weighting = false;
        if (myConfig.viscous_least_squares_type == ViscousLeastSquaresType.weighted_qr) {
            apply_weighting = true;
        }
        if (apply_weighting) {
            if (myConfig.dimensions == 2) {
                foreach (i; loop_init .. n) {
                    number dx = cloud_pos[i].x - x0;
                    number dy = cloud_pos[i].y - y0;
                    theta[i] = 1.0/pow(sqrt(dx*dx+dy*dy), 0.75);
                }
            } else { //3D
                foreach (i; loop_init .. n) {
                    number dx = cloud_pos[i].x - x0;
                    number dy = cloud_pos[i].y - y0;
                    number dz = cloud_pos[i].z - z0;
                    theta[i] = 1.0/pow(sqrt(dx*dx+dy*dy+dz*dz), 0.75);
                }
            }
        } else {
            foreach (i; loop_init .. n) {
                theta[i] = 1.0;
            }
        }

        // If computing about a mid-point, calculate that mid-point.
        if (compute_about_mid) {
            x0 = 0.0; y0 = 0.0; z0 = 0.0;
            foreach (i; 0 .. n) {
                x0 += cloud_pos[i].x; y0 += cloud_pos[i].y; z0 += cloud_pos[i].z;
            }
            x0 /= n; y0 /= n; z0 /= n; // midpoint
        } else { // else use the primary point (assumed to be in cloud position 0)
            x0 = cloud_pos[0].x; y0 = cloud_pos[0].y; z0 = cloud_pos[0].z;
        }

        // compute dx, dy, dz from the A matrix in eq. 12/13 from ref. [1] (incorporating the weighting coefficients)
        number[cloud_nmax] dx, dy, dz;
        foreach (i; loop_init .. n) {
            dx[i] = theta[i]*(cloud_pos[i].x - x0);
            dy[i] = theta[i]*(cloud_pos[i].y - y0);
            if (myConfig.dimensions == 3) {
                dz[i] = theta[i]*(cloud_pos[i].z - z0);
            } else {
                dz[i] = 0.0;
            }
        }

        if (myConfig.dimensions == 2) {
            // compute the vector of weights from eq. 12/13 from ref. [1]
            number r11 = 0.0;
            foreach(i; loop_init .. n) { r11 += pow(dx[i], 2); }
            r11 = sqrt(r11);

            number r12 = 0.0;
            foreach(i; loop_init .. n) { r12 += dx[i]*dy[i]; }
            r12 /= r11;

            number r22 = 0.0;
            foreach(i; loop_init .. n) { r22 += pow(dy[i], 2); }
            r22 = sqrt(r22 - pow(r12,2));

            foreach (i; loop_init .. n) {
                number alpha1 = dx[i]*pow(r11,-2);
                number alpha2 = (dy[i]-(r12/r11)*dx[i])*pow(r22,-2);
                ws.wx[i] = theta[i]*(alpha1 - (r12/r11)*alpha2);
                ws.wy[i] = theta[i]*(alpha2);
            }
        } else {
            // compute the vector of weights from eq. 12/13 from ref. [1]
            number r11 = 0.0;
            foreach(i; loop_init .. n) { r11 += pow(dx[i], 2); }
            r11 = sqrt(r11);

            number r12 = 0.0;
            foreach(i; loop_init .. n) { r12 += dx[i]*dy[i]; }
            r12 /= r11;

            number r22 = 0.0;
            foreach(i; loop_init .. n) { r22 += pow(dy[i], 2); }
            r22 = sqrt(r22 - pow(r12,2));

            number r13 = 0.0;
            foreach(i; loop_init .. n) { r13 += dx[i]*dz[i]; }
            r13 /= r11;

            number r23 = 0.0;
            foreach(i; loop_init .. n) { r23 += dy[i]*dz[i]; }
            r23 = (1.0/r22) * (r23 - r12*r13);

            number r33 = 0.0;
            foreach(i; loop_init .. n) { r33 += pow(dz[i], 2); }
            r33 = sqrt(r33-(pow(r13,2)+pow(r23,2)));

            number beta = (r12*r23-r13*r22)/(r11*r22);

            foreach (i; loop_init .. n) {
                number alpha1 = dx[i]*pow(r11,-2);
                number alpha2 = (dy[i]-(r12/r11)*dx[i])*pow(r22,-2);
                number alpha3 = (dz[i]-(r23/r22)*dy[i]+beta*dx[i])*pow(r33,-2);
                ws.wx[i] = theta[i]*(alpha1 - (r12/r11)*alpha2 + beta*alpha3);
                ws.wy[i] = theta[i]*(alpha2 - (r23/r22)*alpha3);
                ws.wz[i] = theta[i]*(alpha3);
            }
        }
    } // end set_up_workspace_leastsq_via_qr_factorization()

    @nogc
    void gradients_leastsq(in LocalConfig myConfig, ref FlowState*[] cloud_fs, ref WLSQGradWorkspace ws)
    // Evaluate the gradients using the precomputed weights.
    {
        debug {
            if (&ws is null) {
                throw new Error("Uninitialized workspace for least-squares gradient calculation.");
            }
        }
        size_t loop_init = ws.loop_init;
        size_t n = ws.n;
        size_t dimensions = myConfig.dimensions;
        //
        number q0;
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
            "~gname~"[0] = 0.0; "~gname~"[1] = 0.0; "~gname~"[2] = 0.0;
            foreach (i; loop_init .. n) {
                number dq = cloud_fs[i]."~qname~" - q0;
                "~gname~"[0] += ws.wx[i] * dq;
                "~gname~"[1] += ws.wy[i] * dq;
                if (dimensions == 3) { "~gname~"[2] += ws.wz[i] * dq; }
            }";
            return code;
        }
        mixin(codeForGradients("vel.x", "vel[0]"));
        mixin(codeForGradients("vel.y", "vel[1]"));
        if (dimensions == 3) {
            mixin(codeForGradients("vel.z", "vel[2]"));
        } else {
            // 2D z-velocity
            vel[2][0] = 0.0; vel[2][1] = 0.0; vel[2][2] = 0.0;
        }
        mixin(codeForGradients("gas.T", "T"));

        if (myConfig.gmodel.is_plasma) {
            mixin(codeForGradients("gas.p_e", "p_e"));
        }
        version(multi_T_gas) {
            // T_modes
            size_t nmodes = cloud_fs[0].gas.T_modes.length;
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]"));
            }
        }
        version(multi_species_gas) {
            // massf
            size_t nsp = cloud_fs[0].gas.massf.length;
            if (myConfig.turb_model.isTurbulent ||
                myConfig.mass_diffusion_model != MassDiffusionModel.none) {
                foreach(isp; 0 .. nsp) {
                    mixin(codeForGradients("gas.massf[isp]", "massf[isp]"));
                }
            } else {
                foreach(isp; 0 .. nsp) {
                    massf[isp][0] = 0.0; massf[isp][1] = 0.0; massf[isp][2] = 0.0;
                }
            }
        }
        version(turbulence) {
            foreach(tidx; 0 .. myConfig.turb_model.nturb){
                mixin(codeForGradients("turb[tidx]", "turb[tidx]"));
            }
        }
    } // end gradients_leastsq()

    @nogc
    void gradients_at_cells_leastsq(in FlowState fs, FlowState[] fss, size_t[] cloud_idxs, in WLSQGradWorkspace ws,
                                    size_t n, bool is3d, size_t nsp, size_t nmodes, size_t nturb, bool doSpecies)
    /*
        Faster variant of the least-squares gradient code, for working with the dense data structures.
        @author: NNG (Feb 25)
    */
    {
        number q0;
        string codeForGradients(string qname, string gname)
        {
            string code = "
            q0 = fs."~qname~";
            "~gname~"[0] = 0.0; "~gname~"[1] = 0.0; "~gname~"[2] = 0.0;
            foreach (i; 0 .. n) {
                size_t io = cloud_idxs[i];
                number dq = fss[io]."~qname~" - q0;
                "~gname~"[0] += ws.wx[i+1] * dq;
                "~gname~"[1] += ws.wy[i+1] * dq;
                if (is3d) { "~gname~"[2] += ws.wz[i+1] * dq; }
            }";
            return code;
        }
        mixin(codeForGradients("vel.x", "vel[0]"));
        mixin(codeForGradients("vel.y", "vel[1]"));
        if (is3d) {
            mixin(codeForGradients("vel.z", "vel[2]"));
        } else {
            // 2D z-velocity
            vel[2][0] = 0.0; vel[2][1] = 0.0; vel[2][2] = 0.0;
        }
        mixin(codeForGradients("gas.T", "T"));
        version(multi_T_gas) {
            // T_modes
            foreach (imode; 0 .. nmodes) {
                mixin(codeForGradients("gas.T_modes[imode]", "T_modes[imode]"));
            }
        }
        version(multi_species_gas) {
            if (doSpecies) {
                foreach(isp; 0 .. nsp) {
                    mixin(codeForGradients("gas.massf[isp]", "massf[isp]"));
                }
            } else {
                foreach(isp; 0 .. nsp) {
                    massf[isp][0] = 0.0; massf[isp][1] = 0.0; massf[isp][2] = 0.0;
                }
            }
        }
        version(turbulence) {
            foreach(tidx; 0 .. nturb){
                mixin(codeForGradients("turb[tidx]", "turb[tidx]"));
            }
        }
    } // end gradients_leastsq()

} // end struct FlowGradients
