// onedinterp.d
// One-dimensional interpolation/reconstruction of flow field.
//
// See MBCNS workbook 2000/2 page 36 (26-Jan-2001) for formulation.
// and MBCNS workbook 2005/Apr page 36 for new index labels

module onedinterp;

import std.math;
import std.stdio;
import nm.complex;
import nm.number;

import gas;
import fvcore;
import globalconfig;
import flowstate;
import fvinterface;
import fvcell;
import limiters;

immutable double epsilon_van_albada = 1.0e-12;

class OneDInterpolator {

private:
    // The following variables must be set to appropriate values
    // by xxxx_prepare() before their use in interp_xxxx_scalar().
    number aL0;
    number aR0;
    number lenL0_;
    number lenR0_;
    number two_over_lenL0_plus_lenL1;
    number two_over_lenR0_plus_lenL0;
    number two_over_lenR1_plus_lenR0;
    number two_lenL0_plus_lenL1;
    number two_lenR0_plus_lenR1;
    number w0, w1;
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

    //------------------------------------------------------------------------------

    @nogc void l2r2_prepare(number lenL1, number lenL0, number lenR0, number lenR1)
    // Set up intermediate data that depends only on the cell geometry.
    // It will remain constant when reconstructing the different scalar fields
    // over the same set of cells.
    {
        lenL0_ = lenL0;
        lenR0_ = lenR0;
        aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
        two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
        two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
        two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
    } // end l2r2_prepare()

    @nogc void interp_l2r2_scalar(number qL1, number qL0, number qR0, number qR1,
                                  ref number qL, ref number qR)
    {
        number delLminus, del, delRplus, sL, sR;
        // Set up differences and limiter values.
        delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            sL = (delLminus*del + fabs(delLminus*del)) / 
                (delLminus*delLminus + del*del + epsilon_van_albada);
            sR = (del*delRplus + fabs(del*delRplus)) / 
                (del*del + delRplus*delRplus + epsilon_van_albada);
        } else {
            // Use unlimited high-order reconstruction.
            sL = 1.0;
            sR = 1.0;
        }
        // The actual high-order reconstruction, possibly limited.
        qL = qL0 + sL * aL0 * (del * two_lenL0_plus_lenL1 + delLminus * lenR0_);
        qR = qR0 - sR * aR0 * (delRplus * lenL0_ + del * two_lenR0_plus_lenR1);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            qL = clip_to_limits(qL, qL0, qR0);
            qR = clip_to_limits(qR, qL0, qR0);
        }
    } // end of interp_l2r2_scalar()

    //------------------------------------------------------------------------------

    @nogc void l2r1_prepare(number lenL1, number lenL0, number lenR0)
    {
        lenL0_ = lenL0;
        lenR0_ = lenR0;
        aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
        two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
        linear_interp_prepare(lenL0, lenR0);
    } // end l2r1_prepare()

    @nogc void interp_l2r1_scalar(number qL1, number qL0, number qR0, ref number qL, ref number qR)
    {
        number delLminus, del, sL;
        delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        if (myConfig.apply_limiter) {
            sL = (delLminus*del + fabs(delLminus*del)) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
        } else {
            sL = 1.0;
        }
        qL = qL0 + sL * aL0 * (del * two_lenL0_plus_lenL1 + delLminus * lenR0_);
        if (myConfig.apply_limiter && (delLminus*del < 0.0)) {
            qR = qR0;
        } else {
            qR = weight_scalar(qL0, qR0);
        }
        if (myConfig.extrema_clipping) {
            qL = clip_to_limits(qL, qL0, qR0);
            // qR is already clipped inside weight_scalar().
        }
    } // end of interp_l2r1_scalar()

    @nogc void l1r2_prepare(number lenL0, number lenR0, number lenR1)
    {
        lenL0_ = lenL0;
        lenR0_ = lenR0;
        aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
        two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
        linear_interp_prepare(lenL0, lenR0);
    } // end l1r2_prepare()

    @nogc void interp_l1r2_scalar(number qL0, number qR0, number qR1, ref number qL, ref number qR)
    {
        number del, delRplus, sR;
        del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        if (myConfig.apply_limiter) {
            sR = (del*delRplus + fabs(del*delRplus)) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        } else {
            sR = 1.0;
        }
        qR = qR0 - sR * aR0 * (delRplus * lenL0_ + del * two_lenR0_plus_lenR1);
        if (myConfig.apply_limiter && (delRplus*del < 0.0)) {
            qL = qL0;
        } else {
            qL = weight_scalar(qL0, qR0);
        }
        if (myConfig.extrema_clipping) {
            // qL is already clipped inside weight_scalar().
            qR = clip_to_limits(qR, qL0, qR0);
        }
    } // end of interp_l1r2_scalar()

    //------------------------------------------------------------------------------

    @nogc void linear_extrap_prepare(number len0, number len1)
    {
        // Set up weights for a linear combination if q0 and q1
        // assuming that we are extrapolating past q0 to the boundary len0/2 away.
        w0 = (2.0*len0 + len1)/(len0+len1);
        w1 = -len0/(len0+len1);
    }

    @nogc void linear_interp_prepare(number len0, number len1)
    {
        // Set up weights for a linear combination if q0 and q1
        // assuming that we are interpolating to the cell-face between
        // the cell-centre points.
        w0 = len1/(len0+len1);
        w1 = len0/(len0+len1);
    }

    @nogc number weight_scalar(number q0, number q1)
    {
        // The weights for interpolation or extrapolation may be used.
        number q = q0*w0 + q1*w1;
        if (myConfig.extrema_clipping) { q = clip_to_limits(q, q0, q1); }
        return q;
    }

    //------------------------------------------------------------------------------
    
    @nogc
    void interp_l2r2(ref FVInterface IFace,
                     ref FVCell cL1, ref FVCell cL0, ref FVCell cR0, ref FVCell cR1,
                     number cL1Length, number cL0Length, 
                     number cR0Length, number cR1Length, 
                     ref FlowState Lft, ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = (myConfig.sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // Paul Petrie-Repar and Jason Qin have noted that the velocity needs
            // to be reconstructed in the interface-local frame of reference so that
            // the normal velocities are not messed up for mirror-image at walls.
            // PJ 21-feb-2012
            cL1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        l2r2_prepare(cL1Length, cL0Length, cR0Length, cR1Length);
        interp_l2r2_scalar(cL1.fs.vel.x, cL0.fs.vel.x, cR0.fs.vel.x, cR1.fs.vel.x,
                           Lft.vel.refx, Rght.vel.refx);
        interp_l2r2_scalar(cL1.fs.vel.y, cL0.fs.vel.y, cR0.fs.vel.y, cR1.fs.vel.y,
                           Lft.vel.refy, Rght.vel.refy);
        interp_l2r2_scalar(cL1.fs.vel.z, cL0.fs.vel.z, cR0.fs.vel.z, cR1.fs.vel.z,
                           Lft.vel.refz, Rght.vel.refz);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l2r2_scalar(cL1.fs.B.x, cL0.fs.B.x, cR0.fs.B.x, cR1.fs.B.x,
                                   Lft.B.refx, Rght.B.refx);
                interp_l2r2_scalar(cL1.fs.B.y, cL0.fs.B.y, cR0.fs.B.y, cR1.fs.B.y,
                                   Lft.B.refy, Rght.B.refy);
                interp_l2r2_scalar(cL1.fs.B.z, cL0.fs.B.z, cR0.fs.B.z, cR1.fs.B.z,
                                   Lft.B.refz, Rght.B.refz);
                if (myConfig.divergence_cleaning) {
                    interp_l2r2_scalar(cL1.fs.psi, cL0.fs.psi, cR0.fs.psi, cR1.fs.psi,
                                       Lft.psi, Rght.psi);
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                interp_l2r2_scalar(cL1.fs.tke, cL0.fs.tke, cR0.fs.tke, cR1.fs.tke,
                                   Lft.tke, Rght.tke);
                interp_l2r2_scalar(cL1.fs.omega, cL0.fs.omega, cR0.fs.omega, cR1.fs.omega,
                                   Lft.omega, Rght.omega);
            }
        }
        auto gL1 = &(cL1.fs.gas); // Avoid construction of another object.
        auto gL0 = &(cL0.fs.gas);
        auto gR0 = &(cR0.fs.gas);
        auto gR1 = &(cR1.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    interp_l2r2_scalar(gL1.massf[isp], gL0.massf[isp], gR0.massf[isp], gR1.massf[isp],
                                       Lft.gas.massf[isp], Rght.gas.massf[isp]);
                }
                try {
                    scale_mass_fractions(Lft.gas.massf); 
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Lft.gas.massf[] = gL0.massf[];
                }
                try {
                    scale_mass_fractions(Rght.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Rght.gas.massf[] = gR0.massf[];
                }
            } else {
                // Only one possible mass-fraction value for a single species.
                Lft.gas.massf[0] = 1.0;
                Rght.gas.massf[0] = 1.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            interp_l2r2_scalar(gL1.p, gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p);
            interp_l2r2_scalar(gL1.T, gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           gR1.T_modes[i], Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.T_modes[i] = gL0.T_modes[i];
                        Rght.gas.T_modes[i] = gR0.T_modes[i];
                    }
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
            break;
        case InterpolateOption.rhou:
            interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r2_scalar(gL1.u, gL0.u, gR0.u, gR1.u, Lft.gas.u, Rght.gas.u);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                           gR1.u_modes[i], Lft.gas.u_modes[i], Rght.gas.u_modes[i]);
                    }
                } else {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.u_modes[i] = gL0.u_modes[i];
                        Rght.gas.u_modes[i] = gR0.u_modes[i];
                    }
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
            break;
        case InterpolateOption.rhop:
            interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r2_scalar(gL1.p, gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p);
            mixin(codeForThermoUpdateBoth("rhop"));
            break;
        case InterpolateOption.rhot: 
            interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r2_scalar(gL1.T, gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           gR1.T_modes[i], Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.T_modes[i] = gL0.T_modes[i];
                        Rght.gas.T_modes[i] = gR0.T_modes[i];
                    }
                }
            }
            mixin(codeForThermoUpdateBoth("rhoT"));
            break;
        } // end switch thermo_interpolator
        if (myConfig.interpolate_in_local_frame) {
            // Undo the transformation made earlier. PJ 21-feb-2012
            Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } // end interp_l2r2()

    @nogc
    void interp_l2r1(ref FVInterface IFace,
                     ref FVCell cL1, ref FVCell cL0, ref FVCell cR0,
                     number cL1Length, number cL0Length, number cR0Length,
                     ref FlowState Lft, ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = (myConfig.sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cL1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        l2r1_prepare(cL1Length, cL0Length, cR0Length);
        interp_l2r1_scalar(cL1.fs.vel.x, cL0.fs.vel.x, cR0.fs.vel.x, Lft.vel.refx, Rght.vel.refx);
        interp_l2r1_scalar(cL1.fs.vel.y, cL0.fs.vel.y, cR0.fs.vel.y, Lft.vel.refy, Rght.vel.refy);
        interp_l2r1_scalar(cL1.fs.vel.z, cL0.fs.vel.z, cR0.fs.vel.z, Lft.vel.refz, Rght.vel.refz);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l2r1_scalar(cL1.fs.B.x, cL0.fs.B.x, cR0.fs.B.x, Lft.B.refx, Rght.B.refx);
                interp_l2r1_scalar(cL1.fs.B.y, cL0.fs.B.y, cR0.fs.B.y, Lft.B.refy, Rght.B.refy);
                interp_l2r1_scalar(cL1.fs.B.z, cL0.fs.B.z, cR0.fs.B.z, Lft.B.refz, Rght.B.refz);
                if (myConfig.divergence_cleaning) {
                    interp_l2r1_scalar(cL1.fs.psi, cL0.fs.psi, cR0.fs.psi, Lft.psi, Rght.psi);
                }
            }
        }
        version(komega) {
            if ( myConfig.turbulence_model == TurbulenceModel.k_omega ) {
                interp_l2r1_scalar(cL1.fs.tke, cL0.fs.tke, cR0.fs.tke, Lft.tke, Rght.tke);
                interp_l2r1_scalar(cL1.fs.omega, cL0.fs.omega, cR0.fs.omega, Lft.omega, Rght.omega);
            }
        }
        auto gL1 = &(cL1.fs.gas); auto gL0 = &(cL0.fs.gas); auto gR0 = &(cR0.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    interp_l2r1_scalar(gL1.massf[isp], gL0.massf[isp], gR0.massf[isp],
                                       Lft.gas.massf[isp], Rght.gas.massf[isp]);
                }
                try {
                    scale_mass_fractions(Lft.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Lft.gas.massf[] = gL0.massf[];
                }
                try {
                    scale_mass_fractions(Rght.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Rght.gas.massf[] = gR0.massf[];
                }
            } else {
                // Only one possible mass-fraction value for a single species.
                Lft.gas.massf[0] = 1.0;
                Rght.gas.massf[0] = 1.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            interp_l2r1_scalar(gL1.p, gL0.p, gR0.p, Lft.gas.p, Rght.gas.p);
            interp_l2r1_scalar(gL1.T, gL0.T, gR0.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    interp_l2r1_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                       Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
            break;
        case InterpolateOption.rhou:
            interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r1_scalar(gL1.u, gL0.u, gR0.u, Lft.gas.u, Rght.gas.u);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    interp_l2r1_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                       Lft.gas.u_modes[i], Rght.gas.u_modes[i]);
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
            break;
        case InterpolateOption.rhop:
            interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r1_scalar(gL1.p, gL0.p, gR0.p, Lft.gas.p, Rght.gas.p);
            mixin(codeForThermoUpdateBoth("rhop"));
            break;
        case InterpolateOption.rhot: 
            interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r1_scalar(gL1.T, gL0.T, gR0.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    interp_l2r1_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                       Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateBoth("rhoT"));
            break;
        } // end switch thermo_interpolator
        if (myConfig.interpolate_in_local_frame) {
            // Undo the transformation made earlier.
            Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } // end interp_l2r1()

    @nogc
    void interp_l1r2(ref FVInterface IFace,
                     ref FVCell cL0, ref FVCell cR0, ref FVCell cR1,
                     number cL0Length, number cR0Length, number cR1Length,
                     ref FlowState Lft, ref FlowState Rght)
    // Reconstruct flow properties at an interface from cells L0,R0,R1.
    //
    // This is essentially a one-dimensional interpolation process.  It needs only
    // the cell-average data and the lengths of the cells in the interpolation direction.
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = (myConfig.sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        l1r2_prepare(cL0Length, cR0Length, cR1Length);
        interp_l1r2_scalar(cL0.fs.vel.x, cR0.fs.vel.x, cR1.fs.vel.x, Lft.vel.refx, Rght.vel.refx);
        interp_l1r2_scalar(cL0.fs.vel.y, cR0.fs.vel.y, cR1.fs.vel.y, Lft.vel.refy, Rght.vel.refy);
        interp_l1r2_scalar(cL0.fs.vel.z, cR0.fs.vel.z, cR1.fs.vel.z, Lft.vel.refz, Rght.vel.refz);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l1r2_scalar(cL0.fs.B.x, cR0.fs.B.x, cR1.fs.B.x, Lft.B.refx, Rght.B.refx);
                interp_l1r2_scalar(cL0.fs.B.y, cR0.fs.B.y, cR1.fs.B.y, Lft.B.refy, Rght.B.refy);
                interp_l1r2_scalar(cL0.fs.B.z, cR0.fs.B.z, cR1.fs.B.z, Lft.B.refz, Rght.B.refz);
                if (myConfig.divergence_cleaning) {
                    interp_l1r2_scalar(cL0.fs.psi, cR0.fs.psi, cR1.fs.psi, Lft.psi, Rght.psi);
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                interp_l1r2_scalar(cL0.fs.tke, cR0.fs.tke, cR1.fs.tke, Lft.tke, Rght.tke);
                interp_l1r2_scalar(cL0.fs.omega, cR0.fs.omega, cR1.fs.omega, Lft.omega, Rght.omega);
            }
        }
        auto gL0 = &(cL0.fs.gas); auto gR0 = &(cR0.fs.gas); auto gR1 = &(cR1.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    interp_l1r2_scalar(gL0.massf[isp], gR0.massf[isp], gR1.massf[isp],
                                       Lft.gas.massf[isp], Rght.gas.massf[isp]);
                }
                try {
                    scale_mass_fractions(Lft.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Lft.gas.massf[] = gL0.massf[];
                }
                try {
                    scale_mass_fractions(Rght.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Rght.gas.massf[] = gR0.massf[];
                }
            } else {
                // Only one possible mass-fraction value for a single species.
                Lft.gas.massf[0] = 1.0;
                Rght.gas.massf[0] = 1.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            interp_l1r2_scalar(gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p);
            interp_l1r2_scalar(gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    interp_l1r2_scalar(gL0.T_modes[i], gR0.T_modes[i], gR1.T_modes[i],
                                       Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
            break;
        case InterpolateOption.rhou:
            interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l1r2_scalar(gL0.u, gR0.u, gR1.u, Lft.gas.u, Rght.gas.u);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    interp_l1r2_scalar(gL0.u_modes[i], gR0.u_modes[i], gR1.u_modes[i],
                                       Lft.gas.u_modes[i], Rght.gas.u_modes[i]);
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
            break;
        case InterpolateOption.rhop:
            interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l1r2_scalar(gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p);
            mixin(codeForThermoUpdateBoth("rhop"));
            break;
        case InterpolateOption.rhot: 
            interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l1r2_scalar(gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    interp_l1r2_scalar(gL0.T_modes[i], gR0.T_modes[i], gR1.T_modes[i],
                                       Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateBoth("rhoT"));
            break;
        } // end switch thermo_interpolator
        if (myConfig.interpolate_in_local_frame) {
            // Undo the transformation made earlier.
            Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } // end interp_l1r2()

    @nogc
    void interp_l2r0(ref FVInterface IFace,
                     ref FVCell cL1, ref FVCell cL0,
                     number cL1Length, number cL0Length,
                     ref FlowState Lft, ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = (myConfig.sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cL1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        if (myConfig.extrema_clipping) {
            // Not much that we can do with linear extrapolation
            // that doesn't produce an new extreme value.
            w0 = 1.0; w1 = 0.0;
        } else {
            linear_extrap_prepare(cL0Length, cL1Length);
        }
        Lft.vel.refx = weight_scalar(cL0.fs.vel.x, cL1.fs.vel.x);
        Lft.vel.refy = weight_scalar(cL0.fs.vel.y, cL1.fs.vel.y);
        Lft.vel.refz = weight_scalar(cL0.fs.vel.z, cL1.fs.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                Lft.B.refx = weight_scalar(cL0.fs.B.x, cL1.fs.B.x);
                Lft.B.refy = weight_scalar(cL0.fs.B.y, cL1.fs.B.y);
                Lft.B.refz = weight_scalar(cL0.fs.B.z, cL1.fs.B.z);
                if (myConfig.divergence_cleaning) {
                    Lft.psi = weight_scalar(cL0.fs.psi, cL1.fs.psi);
                }
            }
        }
        version(komega) {
            if ( myConfig.turbulence_model == TurbulenceModel.k_omega ) {
                Lft.tke = weight_scalar(cL0.fs.tke, cL1.fs.tke);
                Lft.omega = weight_scalar(cL0.fs.omega, cL1.fs.omega);
            }
        }
        auto gL1 = &(cL1.fs.gas); auto gL0 = &(cL0.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    Lft.gas.massf[isp] = weight_scalar(gL0.massf[isp], gL1.massf[isp]);
                }
                try {
                    scale_mass_fractions(Lft.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Lft.gas.massf[] = gL0.massf[];
                }
            } else {
                // Only one possible mass-fraction value for a single species.
                Lft.gas.massf[0] = 1.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            Lft.gas.p = weight_scalar(gL0.p, gL1.p);
            Lft.gas.T = weight_scalar(gL0.T, gL1.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    Lft.gas.T_modes[i] = weight_scalar(gL0.T_modes[i], gL1.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateLft("pT"));
            break;
        case InterpolateOption.rhou:
            Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            Lft.gas.u = weight_scalar(gL0.u, gL1.u);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    Lft.gas.u_modes[i] = weight_scalar(gL0.u_modes[i], gL1.u_modes[i]);
                }
            }
            mixin(codeForThermoUpdateLft("rhou"));
            break;
        case InterpolateOption.rhop:
            Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            Lft.gas.p = weight_scalar(gL0.p, gL1.p);
            mixin(codeForThermoUpdateLft("rhop"));
            break;
        case InterpolateOption.rhot: 
            Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            Lft.gas.T = weight_scalar(gL0.T, gL1.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    Lft.gas.T_modes[i] = weight_scalar(gL0.T_modes[i], gL1.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateLft("rhoT"));
            break;
        } // end switch thermo_interpolator
        if (myConfig.interpolate_in_local_frame) {
            // Undo the transformation made earlier.
            Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } // end interp_l2r0()

    @nogc
    void interp_l0r2(ref FVInterface IFace,
                     ref FVCell cR0, ref FVCell cR1,
                     number cR0Length, number cR1Length,
                     ref FlowState Lft, ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = (myConfig.sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
        auto nmodes = gmodel.n_modes;
        //
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        if (myConfig.extrema_clipping) {
            // Not much that we can do with linear extrapolation
            // that doesn't produce an new extreme value.
            w0 = 1.0; w1 = 0.0;
        } else {
            linear_extrap_prepare(cR0Length, cR1Length);
        }
        Rght.vel.refx = weight_scalar(cR0.fs.vel.x, cR1.fs.vel.x);
        Rght.vel.refy = weight_scalar(cR0.fs.vel.y, cR1.fs.vel.y);
        Rght.vel.refz = weight_scalar(cR0.fs.vel.z, cR1.fs.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                Rght.B.refx = weight_scalar(cR0.fs.B.x, cR1.fs.B.x);
                Rght.B.refy = weight_scalar(cR0.fs.B.y, cR1.fs.B.y);
                Rght.B.refz = weight_scalar(cR0.fs.B.z, cR1.fs.B.z);
                if (myConfig.divergence_cleaning) {
                    Rght.psi = weight_scalar(cR0.fs.psi, cR1.fs.psi);
                }
            }
        }
        version(komega) {
            if (myConfig.turbulence_model == TurbulenceModel.k_omega) {
                Rght.tke = weight_scalar(cR0.fs.tke, cR1.fs.tke);
                Rght.omega = weight_scalar(cR0.fs.omega, cR1.fs.omega);
            }
        }
        auto gR0 = &(cR0.fs.gas); auto gR1 = &(cR1.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                foreach (isp; 0 .. nsp) {
                    Rght.gas.massf[isp] = weight_scalar(gR0.massf[isp], gR1.massf[isp]);
                }
                try {
                    scale_mass_fractions(Rght.gas.massf);
                } catch(Exception e) {
                    debug { writeln(e.msg); }
                    Rght.gas.massf[] = gR0.massf[];
                }
            } else {
                // Only one possible mass-fraction value for a single species.
                Rght.gas.massf[0] = 1.0;
            }
        }
        // Interpolate on two of the thermodynamic quantities, 
        // and fill in the rest based on an EOS call. 
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt: 
            Rght.gas.p = weight_scalar(gR0.p, gR1.p);
            Rght.gas.T = weight_scalar(gR0.T, gR1.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    Rght.gas.T_modes[i] = weight_scalar(gR0.T_modes[i], gR1.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateRght("pT"));
            break;
        case InterpolateOption.rhou:
            Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            Rght.gas.u = weight_scalar(gR0.u, gR1.u);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    Rght.gas.u_modes[i] = weight_scalar(gR0.u_modes[i], gR1.u_modes[i]);
                }
            }
            mixin(codeForThermoUpdateRght("rhou"));
            break;
        case InterpolateOption.rhop:
            Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            Rght.gas.p = weight_scalar(gR0.p, gR1.p);
            mixin(codeForThermoUpdateRght("rhop"));
            break;
        case InterpolateOption.rhot: 
            Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            Rght.gas.T = weight_scalar(gR0.T, gR1.T);
            version(multi_T_gas) {
                foreach (i; 0 .. nmodes) {
                    Rght.gas.T_modes[i] = weight_scalar(gR0.T_modes[i], gR1.T_modes[i]);
                }
            }
            mixin(codeForThermoUpdateRght("rhoT"));
            break;
        } // end switch thermo_interpolator
        if (myConfig.interpolate_in_local_frame) {
            // Undo the transformation made earlier.
            Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } // end interp_l0r2()

} // end class OneDInterpolator

//------------------------------------------------------------------------------

// Helper functions for code generation.
// If an EOS call fails, fall back to just copying cell-centre data.
// This does presume that the cell-centre data is valid. 
string codeForThermoUpdateLft(string funname)
{
    string code = "
        try {
            gmodel.update_thermo_from_"~funname~"(Lft.gas);
        } catch (Exception e) {
            debug { writeln(e.msg); }
            Lft.copy_values_from(cL0.fs);
        }
        ";
    return code;
}

string codeForThermoUpdateRght(string funname)
{
    string code = "
        try {
            gmodel.update_thermo_from_"~funname~"(Rght.gas);
        } catch (Exception e) {
            debug { writeln(e.msg); }
            Rght.copy_values_from(cR0.fs);
        }
        ";
    return code;
}

string codeForThermoUpdateBoth(string funname)
{
    return codeForThermoUpdateLft(funname) ~ codeForThermoUpdateRght(funname);
}
