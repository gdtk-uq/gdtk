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
import std.conv;

import gas;
import globalconfig;
import flowstate;
import fvinterface;
import fvcell;
import limiters;

immutable double epsilon_van_albada = 1.0e-12;


template GetCellVariable(string name, string location)
{
    const char[] GetCellVariable = 
    "pragma(inline) @nogc void "~name~"(size_t n)(FVCell[] cells, number[] scalars)
    {
        static foreach(i; 0 .. n) {
            scalars[i] = cells[i]."~location~";
        }
    }";
}

template GetCellArrayVariable(string name, string location)
{
    const char[] GetCellArrayVariable = 
    "pragma(inline) @nogc void "~name~"(size_t n)(FVCell[] cells, size_t index, number[] scalars)
    {
        static foreach(i; 0 .. n) {
            scalars[i] = cells[i]."~location~"[index];
        }
    }";
}

mixin(GetCellVariable!("GetVelX", "fs.vel.x"));
mixin(GetCellVariable!("GetVelY", "fs.vel.y"));
mixin(GetCellVariable!("GetVelZ", "fs.vel.z"));

version(MHD) {
mixin(GetCellVariable!("GetBx", "fs.B.x"));
mixin(GetCellVariable!("GetBy", "fs.B.y"));
mixin(GetCellVariable!("GetBz", "fs.B.z"));
mixin(GetCellVariable!("GetPsi", "fs.psi"));
}

version(turbulence) {
mixin(GetCellArrayVariable!("GetTurb", "fs.turb"));
}

version(multi_species_gas) {
mixin(GetCellArrayVariable!("GetMassF", "fs.gas.massf"));
}

mixin(GetCellVariable!("GetP", "fs.gas.p"));
mixin(GetCellVariable!("GetT", "fs.gas.T"));

version(multi_T_gas) {
mixin(GetCellArrayVariable!("GetTMode", "fs.gas.T_modes"));
mixin(GetCellArrayVariable!("GetUMode", "fs.gas.u_modes"));
}

mixin(GetCellVariable!("GetRho", "fs.gas.rho"));
mixin(GetCellVariable!("GetU", "fs.gas.u"));

//------------------------------------------------------------------------------

// Helper functions for code generation.
// If an EOS call fails, fall back to just copying cell-centre data.
// This does presume that the cell-centre data is valid.
string codeForThermoUpdateLft2(string funname)
{
    string code = "
        static if (numL > 0) {
            try {
                gmodel.update_thermo_from_"~funname~"(Lft.gas);
            } catch (Exception e) {
                debug { writeln(e.msg); }
                Lft.copy_values_from(cL[0].fs);
            }
        }
        ";
    return code;
}

string codeForThermoUpdateRght2(string funname)
{
    string code = "
        static if (numR > 0) {
            try {
                gmodel.update_thermo_from_"~funname~"(Rght.gas);
            } catch (Exception e) {
                debug { writeln(e.msg); }
                Rght.copy_values_from(cR[0].fs);
            }
        }
        ";
    return code;
}

string codeForThermoUpdateBoth2(string funname)
{
    return codeForThermoUpdateLft2(funname) ~ codeForThermoUpdateRght2(funname);
}


@nogc void copy_flowstate(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
{
    FVCell cL0 = (f.left_cells.length > 0) ? f.left_cells[0] : f.right_cells[0];
    FVCell cR0 = (f.right_cells.length > 0) ? f.right_cells[0]: f.left_cells[0];
    Lft.copy_values_from(cL0.fs);
    Rght.copy_values_from(cR0.fs);
}

@nogc number weight_scalar(number q0, number q1, number w0, number w1, bool extrema_clipping)
{
    // The weights for interpolation or extrapolation may be used.
    number q = q0*w0 + q1*w1;
    if (extrema_clipping) { q = clip_to_limits(q, q0, q1); }
    return q;
}

//=============================================================================

string codeForInterpolation()
{

    string code = "

    FVCell[numL] cL;
    FVCell[numR] cR;

    static foreach (i; 0 .. numL) {
        cL[i] = f.left_cells[i];
    }

    static foreach (i; 0 .. numR) {
        cR[i] = f.right_cells[i];
    }

    number[numW] weight;
    prep_calculate_scalar(f, cL, cR, weight);

    auto gmodel = myConfig.gmodel;
    uint nsp = myConfig.n_species;
    auto nmodes = myConfig.n_modes;
    // High-order reconstruction for some properties.
    if (myConfig.interpolate_in_local_frame) {
        // Paul Petrie-Repar and Jason Qin have noted that the velocity needs
        // to be reconstructed in the interface-local frame of reference so that
        // the normal velocities are not messed up for mirror-image at walls.
        // PJ 21-feb-2012
        static foreach (i; 0 .. numL) {
            cL[i].fs.vel.transform_to_local_frame(f.n, f.t1, f.t2);
        }
        static foreach (i; 0 .. numR) {
            cR[i].fs.vel.transform_to_local_frame(f.n, f.t1, f.t2);
        }
    }

    number[numL] sL;
    number[numR] sR;

    GetVelX!(numL)(cL, sL); GetVelX!(numR)(cR, sR); calculate_scalar(Lft.vel.x, Rght.vel.x, sL, sR, weight);
    GetVelY!(numL)(cL, sL); GetVelY!(numR)(cR, sR); calculate_scalar(Lft.vel.y, Rght.vel.y, sL, sR, weight);
    GetVelZ!(numL)(cL, sL); GetVelZ!(numR)(cR, sR); calculate_scalar(Lft.vel.z, Rght.vel.z, sL, sR, weight);

    version(MHD) {
        if (myConfig.MHD) {
            GetBx!(numL)(cL,sL); GetBx!(numR)(cR, sR); calculate_scalar(Lft.B.x, Rght.B.x, sL, sR, weight);
            GetBy!(numL)(cL,sL); GetBy!(numR)(cR, sR); calculate_scalar(Lft.B.y, Rght.B.y, sL, sR, weight);
            GetBz!(numL)(cL,sL); GetBz!(numR)(cR, sR); calculate_scalar(Lft.B.z, Rght.B.z, sL, sR, weight);
            if (myConfig.divergence_cleaning) {
                GetPsi!(numL)(cL,sL); GetPsi!(numR)(cR, sR); calculate_scalar(Lft.psi, Rght.psi, sL, sR, weight);
            }
        }
    }
    version(turbulence) {
        foreach (it; 0 .. myConfig.turb_model.nturb){
            GetTurb!(numL)(cL, it, sL); GetTurb!(numR)(cR, it, sR); 
            calculate_scalar(Lft.turb[it], Rght.turb[it], sL, sR, weight);
        }
    }

    version(multi_species_gas) {
        if (nsp > 1) {
            // Multiple species.
            if (myConfig.allow_reconstruction_for_species) {
                foreach (isp; 0 .. nsp) {
                    GetMassF!(numL)(cL, isp, sL); GetMassF!(numR)(cR, isp, sR); 
                    calculate_scalar(Lft.gas.massf[isp], Rght.gas.massf[isp], sL, sR, weight);

                }
                static if (numL > 0) {
                    try {
                        scale_mass_fractions(Lft.gas.massf);
                    } catch(Exception e) {
                        debug { writeln(e.msg); }
                        Lft.gas.massf[] = cL[0].fs.gas.massf[];
                    }
                }
                static if (numR > 0) {
                    try {
                        scale_mass_fractions(Rght.gas.massf);
                    } catch(Exception e) {
                        debug { writeln(e.msg); }
                        Rght.gas.massf[] = cR[0].fs.gas.massf[];
                    }
                }
            } else {
                static if (numL > 0) {
                    Lft.gas.massf[] = cL[0].fs.gas.massf[];
                }
                static if (numR > 0) {
                    Rght.gas.massf[] = cR[0].fs.gas.massf[];
                }
            }
        } else {
            // Only one possible mass-fraction value for a single species.
            static if (numL > 0) {
                Lft.gas.massf[0] = 1.0;
            }
            static if (numR > 0) {
                Rght.gas.massf[0] = 1.0;
            }
        }
    }
    // Interpolate on two of the thermodynamic quantities,
    // and fill in the rest based on an EOS call.
    final switch (myConfig.thermo_interpolator) {
    case InterpolateOption.pt:
        GetP!(numL)(cL,sL); GetP!(numR)(cR, sR); calculate_scalar(Lft.gas.p, Rght.gas.p, sL, sR, weight);
        GetT!(numL)(cL,sL); GetT!(numR)(cR, sR); calculate_scalar(Lft.gas.T, Rght.gas.T, sL, sR, weight);

        version(multi_T_gas) {
            if (myConfig.allow_reconstruction_for_energy_modes) {
                foreach (i; 0 .. nmodes) {
                    GetTMode!(numL)(cL, i, sL); GetTMode!(numR)(cR, i, sR); 
                    calculate_scalar(Lft.gas.T_modes[i], Rght.gas.T_modes[i], sL, sR, weight);
                }
            } else {
                static if (numL > 0) {
                    Lft.gas.T_modes[] = cL[0].fs.gas.T_modes[];
                }
                static if (numR > 0) {
                    Rght.gas.T_modes[] = cR[0].fs.gas.T_modes[];
                }
            }
        }
        mixin(codeForThermoUpdateBoth2(\"pT\"));
        break;
    case InterpolateOption.rhou:
        GetRho!(numL)(cL,sL); GetRho!(numR)(cR, sR); calculate_scalar(Lft.gas.rho, Rght.gas.rho, sL, sR, weight);
        GetU!(numL)(cL,sL); GetU!(numR)(cR, sR); calculate_scalar(Lft.gas.u, Rght.gas.u, sL, sR, weight);
        version(multi_T_gas) {
            if (myConfig.allow_reconstruction_for_energy_modes) {
                foreach (i; 0 .. nmodes) {
                    GetUMode!(numL)(cL, i, sL); GetUMode!(numR)(cR, i, sR); 
                    calculate_scalar(Lft.gas.u_modes[i], Rght.gas.u_modes[i], sL, sR, weight);
                }
            } else {
                static if (numL > 0) {
                    Lft.gas.u_modes[] = cL[0].fs.gas.u_modes[];
                }
                static if (numR > 0) {
                    Rght.gas.u_modes[] = cR[0].fs.gas.u_modes[];
                }
            }
        }
        mixin(codeForThermoUpdateBoth2(\"rhou\"));
        break;
    case InterpolateOption.rhop:
        GetRho!(numL)(cL,sL); GetRho!(numR)(cR, sR); calculate_scalar(Lft.gas.rho, Rght.gas.rho, sL, sR, weight);
        GetP!(numL)(cL,sL); GetP!(numR)(cR, sR); calculate_scalar(Lft.gas.p, Rght.gas.p, sL, sR, weight);
        mixin(codeForThermoUpdateBoth2(\"rhop\"));
        break;
    case InterpolateOption.rhot:
        GetRho!(numL)(cL,sL); GetRho!(numR)(cR, sR); calculate_scalar(Lft.gas.rho, Rght.gas.rho, sL, sR, weight);
        GetT!(numL)(cL,sL); GetT!(numR)(cR, sR); calculate_scalar(Lft.gas.T, Rght.gas.T, sL, sR, weight);
        version(multi_T_gas) {
            if (myConfig.allow_reconstruction_for_energy_modes) {
                foreach (i; 0 .. nmodes) {
                    GetTMode!(numL)(cL, i, sL); GetTMode!(numR)(cR, i, sR); 
                    calculate_scalar(Lft.gas.T_modes[i], Rght.gas.T_modes[i], sL, sR, weight);
                }
            } else {
                static if (numL > 0) {
                    Lft.gas.T_modes[] = cL[0].fs.gas.T_modes[];
                }
                static if (numR > 0) {
                    Rght.gas.T_modes[] = cR[0].fs.gas.T_modes[];
                }
            }
        }
        mixin(codeForThermoUpdateBoth2(\"rhoT\"));
        break;
    } // end switch thermo_interpolator
    if (myConfig.interpolate_in_local_frame) {
        // Undo the transformation made earlier. PJ 21-feb-2012
        static foreach (i; 0 .. numL) {
            cL[i].fs.vel.transform_to_global_frame(f.n, f.t1, f.t2);
        }
        static foreach (i; 0 .. numR) {
            cR[i].fs.vel.transform_to_global_frame(f.n, f.t1, f.t2);
        }

        static if (numL > 0) {
            Lft.vel.transform_to_global_frame(f.n, f.t1, f.t2);
        }
        static if (numR > 0) {
            Rght.vel.transform_to_global_frame(f.n, f.t1, f.t2);
        }

    }";

    return code;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L1R0_O1(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
{
    Lft.copy_values_from(f.left_cells[0].fs);
    Rght.copy_values_from(f.left_cells[0].fs);

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L0R1_O1(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
{
    Rght.copy_values_from(f.right_cells[0].fs);
    Lft.copy_values_from(f.right_cells[0].fs);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L1R1_O1(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
{
    Rght.copy_values_from(f.right_cells[0].fs);
    Lft.copy_values_from(f.left_cells[0].fs);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


@nogc void reconstruct_L0R2_O2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 0;
    enum numR = 2;
    enum numW = 2;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight)
    {

        const number len0 = cR[0].lengths[f.idir];
        const number len1 = cR[1].lengths[f.idir];

        // Set up weights for a linear combination if q0 and q1
        // assuming that we are extrapolating past q0 to the boundary len0/2 away.
        weight[0] = (2.0*len0 + len1)/(len0+len1);
        weight[1] = -len0/(len0+len1);
    }

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] w)
    {
        number q = sR[0]*w[0] + sR[1]*w[1];
        if (myConfig.extrema_clipping) { q = clip_to_limits(q, sR[0], sR[1]); }
        Rght = q;
    }


    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L0R2_O2

@nogc void reconstruct_L0R2_O2_smooth(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 0;
    enum numR = 2;
    enum numW = 0;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight){}

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] w)
    {
        number q = 1.5*(sR[0] - sR[1]);
        if (myConfig.extrema_clipping) { q = clip_to_limits(q, sR[0], sR[1]); }
        Rght = q;
    }


    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L0R2_O2_smooth


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


@nogc void reconstruct_L2R0_O2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 2;
    enum numR = 0;
    enum numW = 2;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight)
    {

        const number len0 = cL[0].lengths[f.idir];
        const number len1 = cL[1].lengths[f.idir];

        // Set up weights for a linear combination if q0 and q1
        // assuming that we are extrapolating past q0 to the boundary len0/2 away.
        weight[0] = (2.0*len0 + len1)/(len0+len1);
        weight[1] = -len0/(len0+len1);
    }

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] w)
    {
        number q = sL[0]*w[0] + sL[1]*w[1];
        if (myConfig.extrema_clipping) { q = clip_to_limits(q, sL[0], sL[1]); }
        Lft = q;
    }


    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L2R0_O2

@nogc void reconstruct_L2R0_O2_smooth(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 2;
    enum numR = 0;
    enum numW = 0;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight){}

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] w)
    {
        number q = 1.5*(sL[0] - sL[1]);
        if (myConfig.extrema_clipping) { q = clip_to_limits(q, sL[0], sL[1]); }
        Lft = q;
    }


    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L2R0_O2_smooth



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L1R2_O2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 1;
    enum numR = 2;
    enum numW = 7;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight)
    {
        const number lenL0 = cL[0].lengths[f.idir];
        const number lenR0 = cR[0].lengths[f.idir];
        const number lenR1 = cR[1].lengths[f.idir];

        weight[0] = lenR0/(lenL0+lenR0);
        weight[1] = lenL0/(lenL0+lenR0);
        weight[2] = 2.0 / (lenR0 + lenL0);
        weight[3] = 2.0 / (lenR1 + lenR0);
        weight[4] = (2.0*lenR0 + lenR1);
        weight[5] = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        weight[6] = lenL0;


    }

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {

        number del = (sR[0] - sL[0]) * weight[2];
        number delRplus = (sR[1] - sR[0]) * weight[3];
        number slopeR = 1.0;
        if (myConfig.apply_limiter) {
            slopeR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        Rght = sR[0] - slopeR * weight[5] * (delRplus * weight[6] + del * weight[4]);

        if (myConfig.apply_limiter && (delRplus*del < 0.0)) {
            Lft = sL[0];
        } else {
            Lft = weight_scalar(sL[0], sR[0], weight[0], weight[1], myConfig.extrema_clipping);
        }
        if (myConfig.extrema_clipping) {
            // Lft is already clipped inside weight_scalar().
            Rght = clip_to_limits(Rght, sL[0], sR[0]);
        }
    }

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L1R2_O2

@nogc void reconstruct_L1R2_O2_smooth(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 1;
    enum numR = 2;
    enum numW = 7;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight){}

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {

        number del = sR[0] - sL[0];
        number delRplus = sR[1] - sR[0];
        number slopeR = 1.0;
        if (myConfig.apply_limiter) {
            slopeR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        Rght = sR[0] - slopeR * 0.125 * (delRplus + del * 3);

        if (myConfig.apply_limiter && (delRplus*del < 0.0)) {
            Lft = sL[0];
        } else {
            Lft = weight_scalar(sL[0], sR[0], to!number(0.5), to!number(0.5), myConfig.extrema_clipping);
        }
        if (myConfig.extrema_clipping) {
            // Lft is already clipped inside weight_scalar().
            Rght = clip_to_limits(Rght, sL[0], sR[0]);
        }
    }

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L1R2_O2_smooth

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L2R1_O2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 2;
    enum numR = 1;
    enum numW = 7;


    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight)
    {
        const number lenL0 = cL[0].lengths[f.idir];
        const number lenL1 = cL[1].lengths[f.idir];
        const number lenR0 = cR[0].lengths[f.idir];

        weight[0] = lenR0/(lenL0+lenR0);
        weight[1] = lenL0/(lenL0+lenR0);
        weight[2] = 2.0 / (lenL0 + lenL1);
        weight[3] = 2.0 / (lenR0 + lenL0);
        weight[4] = (2.0*lenL0 + lenL1);
        weight[5] = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        weight[6] = lenR0;


    } // end prep_calculate_scalar()

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {

        number delLminus = (sL[0] - sL[1]) * weight[2];
        number del = (sR[0] - sL[0]) * weight[3];
        number slopeL = 1.0;
        if (myConfig.apply_limiter) {
            slopeL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
        }
        Lft = sL[0] + slopeL * weight[5] * (del * weight[4] + delLminus * weight[6]);
        if (myConfig.apply_limiter && (delLminus*del < 0.0)) {
            Rght = sR[0];
        } else {
            Rght = weight_scalar(sL[0], sR[0], weight[0], weight[1], myConfig.extrema_clipping);
        }
        if (myConfig.extrema_clipping) {
            Lft = clip_to_limits(Lft, sL[0], sR[0]);
            // Rght is already clipped inside weight_scalar().
        }
    } 

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());
} // end reconstruct_L2R1_O2

@nogc void reconstruct_L2R1_O2_smooth(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 2;
    enum numR = 1;
    enum numW = 7;


    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight){}

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {
        number delLminus = sL[0] - sL[1];
        number del = sR[0] - sL[0];
        number slopeL = 1.0;
        if (myConfig.apply_limiter) {
            slopeL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
        }
        Lft = sL[0] + slopeL * 0.125 * (del * 3 + delLminus);
        if (myConfig.apply_limiter && (delLminus*del < 0.0)) {
            Rght = sR[0];
        } else {
            Rght = weight_scalar(sL[0], sR[0], to!number(0.5), to!number(0.5), myConfig.extrema_clipping);
        }
        if (myConfig.extrema_clipping) {
            Lft = clip_to_limits(Lft, sL[0], sR[0]);
            // Rght is already clipped inside weight_scalar().
        }
    } 

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());
} // end reconstruct_L2R1_O2_smooth

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L2R2_O2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 2;
    enum numR = 2;
    enum numW = 9;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight)
    {
        const number lenL0 = cL[0].lengths[f.idir];
        const number lenL1 = cL[1].lengths[f.idir];
        const number lenR0 = cR[0].lengths[f.idir];
        const number lenR1 = cR[1].lengths[f.idir];

        weight[0] = lenL0;
        weight[1] = lenR0;
        weight[2] = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        weight[3] = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        weight[4] = 2.0 / (lenL0 + lenL1);
        weight[5] = 2.0 / (lenR0 + lenL0);
        weight[6] = 2.0 / (lenR1 + lenR0);
        weight[7] = (2.0*lenL0 + lenL1);
        weight[8] = (2.0*lenR0 + lenR1);
    }


    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {

        // Set up differences and limiter values.
        number delLminus = (sL[0] - sL[1]) * weight[4];
        number del = (sR[0] - sL[0]) * weight[5];
        number delRplus = (sR[1] - sR[0]) * weight[6];
        // Presume unlimited high-order reconstruction.
        number slopeL = 1.0;
        number slopeR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            slopeL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
            slopeR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        Lft = sL[0] + slopeL * weight[2] * (del * weight[7] + delLminus * weight[1]);
        Rght = sR[0] - slopeR * weight[3] * (delRplus * weight[0] + del * weight[8]);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            Lft = clip_to_limits(Lft, sL[0], sR[0]);
            Rght = clip_to_limits(Rght, sL[0], sR[0]);
        }
    }

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());
} // end reconstruct_L2R2_O2

@nogc void reconstruct_L2R2_O2_smooth(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 2;
    enum numR = 2;
    enum numW = 9;

    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight){}


    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {
        // Set up differences and limiter values.
        number delLminus = sL[0] - sL[1];
        number del = sR[0] - sL[0];
        number delRplus = sR[1] - sR[0];
        // Presume unlimited high-order reconstruction.
        number slopeL = 1.0;
        number slopeR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            slopeL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
            slopeR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        Lft = sL[0] + slopeL * 0.125 * (del * 3 + delLminus);
        Rght = sR[0] - slopeR * 0.125 * (delRplus + del * 3);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            Lft = clip_to_limits(Lft, sL[0], sR[0]);
            Rght = clip_to_limits(Rght, sL[0], sR[0]);
        }
    }

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());
} // end reconstruct_L2R2_O2_smooth

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@nogc void reconstruct_L3R3_O3(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 3;
    enum numR = 3;
    enum numW = 10;


    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight)
    {
        const number lenL0 = cL[0].lengths[f.idir];
        const number lenL1 = cL[1].lengths[f.idir];
        const number lenL2 = cL[2].lengths[f.idir];

        const number lenR0 = cR[0].lengths[f.idir];
        const number lenR1 = cR[1].lengths[f.idir];
        const number lenR2 = cR[2].lengths[f.idir];

        // Positions of the cell centres relative to interpolation point.
        const number xL0 = -(0.5*lenL0);
        const number xL1 = -(lenL0 + 0.5*lenL1);
        const number xL2 = -(lenL0 + lenL1 + 0.5*lenL2);
        const number xR0 = 0.5*lenR0;
        const number xR1 = lenR0 + 0.5*lenR1;
        const number xR2 = lenR0 + lenR1 + 0.5*lenR2;
        // Weights for Lagrangian interpolation at x=0.
        weight[0] = xL1*xL0*xR0*xR1/((xL2-xL1)*(xL2-xL0)*(xL2-xR0)*(xL2-xR1));
        weight[1] = xL2*xL0*xR0*xR1/((xL1-xL2)*(xL1-xL0)*(xL1-xR0)*(xL1-xR1));
        weight[2] = xL2*xL1*xR0*xR1/((xL0-xL2)*(xL0-xL1)*(xL0-xR0)*(xL0-xR1));
        weight[3] = xL2*xL1*xL0*xR1/((xR0-xL2)*(xR0-xL1)*(xR0-xL0)*(xR0-xR1));
        weight[4] = xL2*xL1*xL0*xR0/((xR1-xL2)*(xR1-xL1)*(xR1-xL0)*(xR1-xR0));
        weight[5] = xL0*xR0*xR1*xR2/((xL1-xL0)*(xL1-xR0)*(xL1-xR1)*(xL1-xR2));
        weight[6] = xL1*xR0*xR1*xR2/((xL0-xL1)*(xL0-xR0)*(xL0-xR1)*(xL0-xR2));
        weight[7] = xL1*xL0*xR1*xR2/((xR0-xL1)*(xR0-xL0)*(xR0-xR1)*(xR0-xR2));
        weight[8] = xL1*xL0*xR0*xR2/((xR1-xL1)*(xR1-xL0)*(xR1-xR0)*(xR1-xR2));
        weight[9] = xL1*xL0*xR0*xR1/((xR2-xL1)*(xR2-xL0)*(xR2-xR0)*(xR2-xR1));

    } // end l3r3_prepare()

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {

        number qL0 = sL[0]; number qR0 = sR[0];
        const number qL1 = sL[1]; const number qR1 = sL[1];
        const number qL2 = sL[2]; const number qR2 = sL[2];

        // Set up differences and limiter values.
        number delLminus = (qL0 - qL1);
        number del = (qR0 - qL0);
        number delRplus = (qR1 - qR0);
        // Presume unlimited high-order reconstruction.
        number slopeL = 1.0;
        number slopeR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            slopeL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
            slopeR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        Lft = qL0 + slopeL * (weight[0]*qL2 + weight[1]*qL1 + (weight[2]-1.0)*qL0 + weight[3]*qR0 + weight[4]*qR1);
        Rght = qR0 + slopeR * (weight[5]*qL1 + weight[6]*qL0 + (weight[7]-1.0)*qR0 + weight[8]*qR1 + weight[9]*qR2);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            Lft = clip_to_limits(Lft, qL0, qR0);
            Rght = clip_to_limits(Rght, qL0, qR0);
        }
    } // end of calculate_scalar()

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

}

@nogc void reconstruct_L3R3_O3_smooth(ref FVInterface f, ref FlowState Lft, ref FlowState Rght, LocalConfig myConfig)
{

    enum numL = 3;
    enum numR = 3;
    enum numW = 10;


    pragma(inline) @nogc void prep_calculate_scalar(ref FVInterface f, ref FVCell[numL] cL, ref FVCell[numR] cR, ref number[numW] weight){}

    pragma(inline) @nogc void calculate_scalar(ref number Lft, ref number Rght, ref number[numL] sL, ref number[numR] sR, ref number[numW] weight)
    {
        const double wL_L2 =  0.0234375;
        const double wL_L1 =  -0.15625;
        const double wL_L0 =  0.703125;
        const double wL_R0 =  0.46875;
        const double wL_R1 =  -0.0390625;
        const double wR_L1 =  -0.0390625;
        const double wR_L0 =  0.46875;
        const double wR_R0 =  0.703125;
        const double wR_R1 =  -0.15625;
        const double wR_R2 =  0.0234375;


        number qL0 = sL[0]; number qR0 = sR[0];
        const number qL1 = sL[1]; const number qR1 = sL[1];
        const number qL2 = sL[2]; const number qR2 = sL[2];

        // Set up differences and limiter values.
        number delLminus = (qL0 - qL1);
        number del = (qR0 - qL0);
        number delRplus = (qR1 - qR0);
        // Presume unlimited high-order reconstruction.
        number slopeL = 1.0;
        number slopeR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            slopeL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
            slopeR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        Lft = qL0 + slopeL * (wL_L2*qL2 + wL_L1*qL1 + (wL_L0-1.0)*qL0 + wL_R0*qR0 + wL_R1*qR1);
        Rght = qR0 + slopeR * (wR_L1*qL1 + wR_L0*qL0 + (wR_R0-1.0)*qR0 + wR_R1*qR1 + wR_R2*qR2);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            Lft = clip_to_limits(Lft, qL0, qR0);
            Rght = clip_to_limits(Rght, qL0, qR0);
        }
    } // end of calculate_scalar()

    copy_flowstate(f, Lft, Rght);

    mixin(codeForInterpolation());

} // end reconstruct_L3R3_O3_smooth

//=============================================================================

/**
 * base class for one dimensional reconstruction with adaptive cell counts
 */
 class OneDInterpolator {
    public:
    LocalConfig myConfig;
    
    this(){}
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


    @nogc abstract void interp(ref FVInterface f, ref FlowState Lft, ref FlowState Rght);

    @nogc abstract void interp_l0r2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght);


} // end class OneDInterpolator

class IrregularOneDInterpolator : OneDInterpolator {

    public:

    this(LocalConfig myConfig)
    {
        this.myConfig = myConfig;
    }


    @nogc final override void interp(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    {
        const size_t nL = f.left_cells.length;
        const size_t nR = f.right_cells.length;

        if (nL >= 3) {
            if (nR >= 3) {
                reconstruct_L3R3_O3(f, Lft, Rght, myConfig);
            } else if (nR == 2) {
                reconstruct_L2R2_O2(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L2R1_O2(f, Lft, Rght, myConfig);
            } else {
                reconstruct_L2R0_O2(f, Lft, Rght, myConfig);
            }
        } else if (nL == 2) {
            if (nR >= 2) {
                reconstruct_L2R2_O2(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L2R1_O2(f, Lft, Rght, myConfig);
            } else {
                reconstruct_L2R0_O2(f, Lft, Rght, myConfig);
            }
        } else if (nL == 1) {
            if (nR >= 2) {
                reconstruct_L1R2_O2(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L1R1_O1(f, Lft, Rght);
            } else {
                reconstruct_L1R0_O1(f, Lft, Rght);
            }
        } else if (nL == 0) {
            if (nR >= 2) {
                reconstruct_L0R2_O2(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L0R1_O1(f, Lft, Rght);
            } else {
                throw new Error("Stencils not suitable for reconstruction.");
            }
        }
    }

    @nogc final override void interp_l0r2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    {
        reconstruct_L0R2_O2(f, Lft, Rght, myConfig);
    }

} // end class IrregularOneDInterpolator


class RegularOneDInterpolator : OneDInterpolator {

    public:

    this(LocalConfig myConfig)
    {
        this.myConfig = myConfig;
    }


    @nogc final override void interp(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    {
        const size_t nL = f.left_cells.length;
        const size_t nR = f.right_cells.length;

        if (nL >= 3) {
            if (nR >= 3) {
                reconstruct_L3R3_O3(f, Lft, Rght, myConfig);
            } else if (nR == 2) {
                reconstruct_L2R2_O2_smooth(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L2R1_O2_smooth(f, Lft, Rght, myConfig);
            } else {
                reconstruct_L2R0_O2_smooth(f, Lft, Rght, myConfig);
            }
        } else if (nL == 2) {
            if (nR >= 2) {
                reconstruct_L2R2_O2_smooth(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L2R1_O2_smooth(f, Lft, Rght, myConfig);
            } else {
                reconstruct_L2R0_O2_smooth(f, Lft, Rght, myConfig);
            }
        } else if (nL == 1) {
            if (nR >= 2) {
                reconstruct_L1R2_O2_smooth(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L1R1_O1(f, Lft, Rght);
            } else {
                reconstruct_L1R0_O1(f, Lft, Rght);
            }
        } else if (nL == 0) {
            if (nR >= 2) {
                reconstruct_L0R2_O2_smooth(f, Lft, Rght, myConfig);
            } else if (nR == 1) {
                reconstruct_L0R1_O1(f, Lft, Rght);
            } else {
                throw new Error("Stencils not suitable for reconstruction.");
            }
        }
    }

    @nogc final override void interp_l0r2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    {
        reconstruct_L0R2_O2_smooth(f, Lft, Rght, myConfig);
    }

} // end class RegularOneDInterpolator


//=============================================================================

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


class LegacyOneDInterpolator  : OneDInterpolator {

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
    number wL_L2, wL_L1, wL_L0, wL_R0, wL_R1;
    number wR_L1, wR_L0, wR_R0, wR_R1, wR_R2;
    LocalConfig myConfig;

public:
    this(LocalConfig myConfig)
    {
        this.myConfig = myConfig;
    }

    @nogc override
    void  interp(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    // Top-level interpolation function delegates to the specific functions, below.
    {
        FVCell cL0 = (f.left_cells.length > 0) ? f.left_cells[0] : f.right_cells[0];
        FVCell cR0 = (f.right_cells.length > 0) ? f.right_cells[0]: f.left_cells[0];
        Lft.copy_values_from(cL0.fs);
        Rght.copy_values_from(cR0.fs);

        if (myConfig.interpolation_order == 3) {
            // A form of high-order flux calculation built on
            // reconstruction via Lagrangian interpolation across a 6-cell stencil.
            // It is a bit of an experiment...
            if (f.left_cells.length < 3 || f.right_cells.length < 3) {
                throw new Error("Stencils not suitable for third-order interpolation.");
            }
            cL0 = f.left_cells[0]; auto cL1 = f.left_cells[1]; auto cL2 = f.left_cells[2];
            cR0 = f.right_cells[0]; auto cR1 = f.right_cells[1]; auto cR2 = f.right_cells[2];
            number lL0, lL1, lL2, lR0, lR1, lR2;
            final switch (f.idir) {
            case IndexDirection.i:
                lL0 = cL0.iLength; lL1 = cL1.iLength; lL2 = cL2.iLength;
                lR0 = cR0.iLength; lR1 = cR1.iLength; lR2 = cR2.iLength;
                break;
            case IndexDirection.j:
                lL0 = cL0.jLength; lL1 = cL1.jLength; lL2 = cL2.jLength;
                lR0 = cR0.jLength; lR1 = cR1.jLength; lR2 = cR2.jLength;
                break;
            case IndexDirection.k:
                lL0 = cL0.kLength; lL1 = cL1.kLength; lL2 = cL2.kLength;
                lR0 = cR0.kLength; lR1 = cR1.kLength; lR2 = cR2.kLength;
                break;
            case IndexDirection.none:
                throw new Error("Invalid index direction.");
            }
            interp_l3r3(f, cL2, cL1, cL0, cR0, cR1, cR2, lL2, lL1, lL0, lR0, lR1, lR2, Lft, Rght);
        } else {
            // Eilmer's classic not-so-high-order reconstruction
            // that attempts to do piecewise-parabolic reconstruction
            // if there are sufficient cells in the stencils.
            if (f.left_cells.length == 0 && f.right_cells.length >= 2) {
                cR0 = f.right_cells[0]; auto cR1 = f.right_cells[1];
                number lR0, lR1;
                final switch (f.idir) {
                case IndexDirection.i:
                    lR0 = cR0.iLength; lR1 = cR1.iLength;
                    break;
                case IndexDirection.j:
                    lR0 = cR0.jLength; lR1 = cR1.jLength;
                    break;
                case IndexDirection.k:
                    lR0 = cR0.kLength; lR1 = cR1.kLength;
                    break;
                case IndexDirection.none:
                    throw new Error("Invalid index direction.");
                }
                interp_l0r2(f, cR0, cR1, lR0, lR1, Rght);
            } else if (f.left_cells.length == 1 && f.right_cells.length >= 2) {
                cL0 = f.left_cells[0];
                cR0 = f.right_cells[0]; auto cR1 = f.right_cells[1];
                number lL0, lR0, lR1;
                final switch (f.idir) {
                case IndexDirection.i:
                    lL0 = cL0.iLength; lR0 = cR0.iLength; lR1 = cR1.iLength;
                    break;
                case IndexDirection.j:
                    lL0 = cL0.jLength; lR0 = cR0.jLength; lR1 = cR1.jLength;
                    break;
                case IndexDirection.k:
                    lL0 = cL0.kLength; lR0 = cR0.kLength; lR1 = cR1.kLength;
                    break;
                case IndexDirection.none:
                    throw new Error("Invalid index direction.");
                }
                interp_l1r2(f, cL0, cR0, cR1, lL0, lR0, lR1, Lft, Rght);
            } else if (f.left_cells.length >= 2 && f.right_cells.length == 1) {
                cL0 = f.left_cells[0]; auto cL1 = f.left_cells[1];
                cR0 = f.right_cells[0];
                number lL0, lL1, lR0;
                final switch (f.idir) {
                case IndexDirection.i:
                    lL0 = cL0.iLength; lL1 = cL1.iLength; lR0 = cR0.iLength;
                    break;
                case IndexDirection.j:
                    lL0 = cL0.jLength; lL1 = cL1.jLength; lR0 = cR0.jLength;
                    break;
                case IndexDirection.k:
                    lL0 = cL0.kLength; lL1 = cL1.kLength; lR0 = cR0.kLength;
                    break;
                case IndexDirection.none:
                    throw new Error("Invalid index direction.");
                }
                interp_l2r1(f, cL1, cL0, cR0, lL1, lL0, lR0, Lft, Rght);
            } else if (f.left_cells.length >= 2 && f.right_cells.length == 0) {
                cL0 = f.left_cells[0]; auto cL1 = f.left_cells[1];
                number lL0, lL1;
                final switch (f.idir) {
                case IndexDirection.i:
                    lL0 = cL0.iLength; lL1 = cL1.iLength;
                    break;
                case IndexDirection.j:
                    lL0 = cL0.jLength; lL1 = cL1.jLength;
                    break;
                case IndexDirection.k:
                    lL0 = cL0.kLength; lL1 = cL1.kLength;
                    break;
                case IndexDirection.none:
                    throw new Error("Invalid index direction.");
                }
                interp_l2r0(f, cL1, cL0, cL1.iLength, cL0.iLength, Lft);
            } else if (f.left_cells.length >= 2 && f.right_cells.length >= 2) {
                // General symmetric reconstruction.
                cL0 = f.left_cells[0]; auto cL1 = f.left_cells[1];
                cR0 = f.right_cells[0]; auto cR1 = f.right_cells[1];
                number lL0, lL1, lR0, lR1;
                final switch (f.idir) {
                case IndexDirection.i:
                    lL0 = cL0.iLength; lL1 = cL1.iLength;
                    lR0 = cR0.iLength; lR1 = cR1.iLength;
                    break;
                case IndexDirection.j:
                    lL0 = cL0.jLength; lL1 = cL1.jLength;
                    lR0 = cR0.jLength; lR1 = cR1.jLength;
                    break;
                case IndexDirection.k:
                    lL0 = cL0.kLength; lL1 = cL1.kLength;
                    lR0 = cR0.kLength; lR1 = cR1.kLength;
                    break;
                case IndexDirection.none:
                    throw new Error("Invalid index direction.");
                }
                interp_l2r2(f, cL1, cL0, cR0, cR1, lL1, lL0, lR0, lR1, Lft, Rght);
            } else {
                throw new Error("Stencils not suitable for standard interpolation.");
            }
        }
    } // end interp()

    //------------------------------------------------------------------------------

    @nogc void l3r3_prepare(number lenL2, number lenL1, number lenL0,
                            number lenR0, number lenR1, number lenR2)
    // Set up intermediate data (Lagrange interpolation weights)
    // that depend only on the cell geometry.
    // They will remain constant when reconstructing the different scalar fields
    // over the same set of cells.
    {
        // Positions of the cell centres relative to interpolation point.
        number xL0 = -(0.5*lenL0);
        number xL1 = -(lenL0 + 0.5*lenL1);
        number xL2 = -(lenL0 + lenL1 + 0.5*lenL2);
        number xR0 = 0.5*lenR0;
        number xR1 = lenR0 + 0.5*lenR1;
        number xR2 = lenR0 + lenR1 + 0.5*lenR2;
        // Weights for Lagrangian interpolation at x=0.
        wL_L2 = xL1*xL0*xR0*xR1/((xL2-xL1)*(xL2-xL0)*(xL2-xR0)*(xL2-xR1));
        wL_L1 = xL2*xL0*xR0*xR1/((xL1-xL2)*(xL1-xL0)*(xL1-xR0)*(xL1-xR1));
        wL_L0 = xL2*xL1*xR0*xR1/((xL0-xL2)*(xL0-xL1)*(xL0-xR0)*(xL0-xR1));
        wL_R0 = xL2*xL1*xL0*xR1/((xR0-xL2)*(xR0-xL1)*(xR0-xL0)*(xR0-xR1));
        wL_R1 = xL2*xL1*xL0*xR0/((xR1-xL2)*(xR1-xL1)*(xR1-xL0)*(xR1-xR0));
        wR_L1 = xL0*xR0*xR1*xR2/((xL1-xL0)*(xL1-xR0)*(xL1-xR1)*(xL1-xR2));
        wR_L0 = xL1*xR0*xR1*xR2/((xL0-xL1)*(xL0-xR0)*(xL0-xR1)*(xL0-xR2));
        wR_R0 = xL1*xL0*xR1*xR2/((xR0-xL1)*(xR0-xL0)*(xR0-xR1)*(xR0-xR2));
        wR_R1 = xL1*xL0*xR0*xR2/((xR1-xL1)*(xR1-xL0)*(xR1-xR0)*(xR1-xR2));
        wR_R2 = xL1*xL0*xR0*xR1/((xR2-xL1)*(xR2-xL0)*(xR2-xR0)*(xR2-xR1));
    } // end l3r3_prepare()

    @nogc void interp_l3r3_scalar(number qL2, number qL1, number qL0,
                                  number qR0, number qR1, number qR2,
                                  ref number qL, ref number qR)
    {
        // Set up differences and limiter values.
        number delLminus = (qL0 - qL1);
        number del = (qR0 - qL0);
        number delRplus = (qR1 - qR0);
        // Presume unlimited high-order reconstruction.
        number sL = 1.0;
        number sR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            sL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
            sR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        qL = qL0 + sL * (wL_L2*qL2 + wL_L1*qL1 + (wL_L0-1.0)*qL0 + wL_R0*qR0 + wL_R1*qR1);
        qR = qR0 + sR * (wR_L1*qL1 + wR_L0*qL0 + (wR_R0-1.0)*qR0 + wR_R1*qR1 + wR_R2*qR2);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            qL = clip_to_limits(qL, qL0, qR0);
            qR = clip_to_limits(qR, qL0, qR0);
        }
    } // end of interp_l3r3_scalar()

    pragma(inline, true)
    @nogc void l2r2_prepare(number lenL1, number lenL0, number lenR0, number lenR1)
    // Set up intermediate data that depend only on the cell geometry.
    // They will remain constant when reconstructing the different scalar fields
    // over the same set of cells.
    // For piecewise parabolic reconstruction, see PJ workbook notes Jan 2001.
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

    pragma(inline, true)
    @nogc void interp_l2r2_scalar(number qL1, number qL0, number qR0, number qR1,
                                  ref number qL, ref number qR)
    {
        // Set up differences and limiter values.
        number delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        number del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        number delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        // Presume unlimited high-order reconstruction.
        number sL = 1.0;
        number sR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            sL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
            sR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
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
        number delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        number del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        number sL = 1.0;
        if (myConfig.apply_limiter) {
            sL = (delLminus*del + fabs(delLminus*del) + epsilon_van_albada) /
                (delLminus*delLminus + del*del + epsilon_van_albada);
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
        number del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        number delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        number sR = 1.0;
        if (myConfig.apply_limiter) {
            sR = (del*delRplus + fabs(del*delRplus) + epsilon_van_albada) /
                (del*del + delRplus*delRplus + epsilon_van_albada);
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
    void interp_l3r3(ref FVInterface IFace,
                     ref FVCell cL2, ref FVCell cL1, ref FVCell cL0,
                     ref FVCell cR0, ref FVCell cR1, ref FVCell cR2,
                     number cL2Length, number cL1Length, number cL0Length,
                     number cR0Length, number cR1Length, number cR2Length,
                     ref FlowState Lft, ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // Paul Petrie-Repar and Jason Qin have noted that the velocity needs
            // to be reconstructed in the interface-local frame of reference so that
            // the normal velocities are not messed up for mirror-image at walls.
            // PJ 21-feb-2012
            cL2.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR2.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        l3r3_prepare(cL2Length, cL1Length, cL0Length, cR0Length, cR1Length, cR2Length);
        interp_l3r3_scalar(cL2.fs.vel.x, cL1.fs.vel.x, cL0.fs.vel.x,
                           cR0.fs.vel.x, cR1.fs.vel.x, cR2.fs.vel.x,
                           Lft.vel.x, Rght.vel.x);
        interp_l3r3_scalar(cL2.fs.vel.y, cL1.fs.vel.y, cL0.fs.vel.y,
                           cR0.fs.vel.y, cR1.fs.vel.y, cR2.fs.vel.y,
                           Lft.vel.y, Rght.vel.y);
        interp_l3r3_scalar(cL2.fs.vel.z, cL1.fs.vel.z, cL0.fs.vel.z,
                           cR0.fs.vel.z, cR1.fs.vel.z, cR2.fs.vel.z,
                           Lft.vel.z, Rght.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l3r3_scalar(cL2.fs.B.x, cL1.fs.B.x, cL0.fs.B.x,
                                   cR0.fs.B.x, cR1.fs.B.x, cR2.fs.B.x,
                                   Lft.B.x, Rght.B.x);
                interp_l3r3_scalar(cL2.fs.B.y, cL1.fs.B.y, cL0.fs.B.y,
                                   cR0.fs.B.y, cR1.fs.B.y, cR2.fs.B.y,
                                   Lft.B.y, Rght.B.y);
                interp_l3r3_scalar(cL2.fs.B.z, cL1.fs.B.z, cL0.fs.B.z,
                                   cR0.fs.B.z, cR1.fs.B.z, cR2.fs.B.z,
                                   Lft.B.z, Rght.B.z);
                if (myConfig.divergence_cleaning) {
                    interp_l3r3_scalar(cL2.fs.psi, cL1.fs.psi, cL0.fs.psi,
                                       cR0.fs.psi, cR1.fs.psi, cR2.fs.psi,
                                       Lft.psi, Rght.psi);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                interp_l3r3_scalar(cL2.fs.turb[it], cL1.fs.turb[it], cL0.fs.turb[it],
                                   cR0.fs.turb[it], cR1.fs.turb[it], cR2.fs.turb[it],
                                   Lft.turb[it], Rght.turb[it]);
                }
        }
        auto gL2 = &(cL2.fs.gas); // Avoid construction of another object.
        auto gL1 = &(cL1.fs.gas);
        auto gL0 = &(cL0.fs.gas);
        auto gR0 = &(cR0.fs.gas);
        auto gR1 = &(cR1.fs.gas);
        auto gR2 = &(cR2.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        interp_l3r3_scalar(gL2.massf[isp], gL1.massf[isp], gL0.massf[isp],
                                           gR0.massf[isp], gR1.massf[isp], gR2.massf[isp],
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
                    Lft.gas.massf[] = gL0.massf[];
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
            interp_l3r3_scalar(gL2.p, gL1.p, gL0.p, gR0.p, gR1.p, gR2.p, Lft.gas.p, Rght.gas.p);
            interp_l3r3_scalar(gL2.T, gL1.T, gL0.T, gR0.T, gR1.T, gR2.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.T_modes[i], gL1.T_modes[i], gL0.T_modes[i],
                                           gR0.T_modes[i], gR1.T_modes[i], gR2.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
            break;
        case InterpolateOption.rhou:
            interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l3r3_scalar(gL2.u, gL1.u, gL0.u, gR0.u, gR1.u, gR2.u, Lft.gas.u, Rght.gas.u);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.u_modes[i], gL1.u_modes[i], gL0.u_modes[i],
                                           gR0.u_modes[i], gR1.u_modes[i], gR2.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i]);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
            break;
        case InterpolateOption.rhop:
            interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l3r3_scalar(gL2.p, gL1.p, gL0.p, gR0.p, gR1.p, gR2.p, Lft.gas.p, Rght.gas.p);
            mixin(codeForThermoUpdateBoth("rhop"));
            break;
        case InterpolateOption.rhot:
            interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l3r3_scalar(gL2.T, gL1.T, gL0.T, gR0.T, gR1.T, gR2.T, Lft.gas.T, Rght.gas.T);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.T_modes[i], gL1.T_modes[i], gL0.T_modes[i],
                                           gR0.T_modes[i], gR1.T_modes[i], gR2.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhoT"));
            break;
        } // end switch thermo_interpolator
        if (myConfig.interpolate_in_local_frame) {
            // Undo the transformation made earlier. PJ 21-feb-2012
            Lft.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            Rght.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL2.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
            cR2.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } // end interp_l3r3()

    //------------------------------------------------------------------------------

    @nogc
    void interp_l2r2(ref FVInterface IFace,
                     ref FVCell cL1, ref FVCell cL0, ref FVCell cR0, ref FVCell cR1,
                     number cL1Length, number cL0Length,
                     number cR0Length, number cR1Length,
                     ref FlowState Lft, ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
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
                           Lft.vel.x, Rght.vel.x);
        interp_l2r2_scalar(cL1.fs.vel.y, cL0.fs.vel.y, cR0.fs.vel.y, cR1.fs.vel.y,
                           Lft.vel.y, Rght.vel.y);
        interp_l2r2_scalar(cL1.fs.vel.z, cL0.fs.vel.z, cR0.fs.vel.z, cR1.fs.vel.z,
                           Lft.vel.z, Rght.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l2r2_scalar(cL1.fs.B.x, cL0.fs.B.x, cR0.fs.B.x, cR1.fs.B.x,
                                   Lft.B.x, Rght.B.x);
                interp_l2r2_scalar(cL1.fs.B.y, cL0.fs.B.y, cR0.fs.B.y, cR1.fs.B.y,
                                   Lft.B.y, Rght.B.y);
                interp_l2r2_scalar(cL1.fs.B.z, cL0.fs.B.z, cR0.fs.B.z, cR1.fs.B.z,
                                   Lft.B.z, Rght.B.z);
                if (myConfig.divergence_cleaning) {
                    interp_l2r2_scalar(cL1.fs.psi, cL0.fs.psi, cR0.fs.psi, cR1.fs.psi,
                                       Lft.psi, Rght.psi);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    interp_l2r2_scalar(cL1.fs.turb[it], cL0.fs.turb[it], cR0.fs.turb[it],
                                       cR1.fs.turb[it], Lft.turb[it], Rght.turb[it]);
                }
        }
        auto gL1 = &(cL1.fs.gas); // Avoid construction of another object.
        auto gL0 = &(cL0.fs.gas);
        auto gR0 = &(cR0.fs.gas);
        auto gR1 = &(cR1.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                if (myConfig.allow_reconstruction_for_species) {
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
                    Lft.gas.massf[] = gL0.massf[];
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
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
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
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
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
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
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
        uint nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cL1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        l2r1_prepare(cL1Length, cL0Length, cR0Length);
        interp_l2r1_scalar(cL1.fs.vel.x, cL0.fs.vel.x, cR0.fs.vel.x, Lft.vel.x, Rght.vel.x);
        interp_l2r1_scalar(cL1.fs.vel.y, cL0.fs.vel.y, cR0.fs.vel.y, Lft.vel.y, Rght.vel.y);
        interp_l2r1_scalar(cL1.fs.vel.z, cL0.fs.vel.z, cR0.fs.vel.z, Lft.vel.z, Rght.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l2r1_scalar(cL1.fs.B.x, cL0.fs.B.x, cR0.fs.B.x, Lft.B.x, Rght.B.x);
                interp_l2r1_scalar(cL1.fs.B.y, cL0.fs.B.y, cR0.fs.B.y, Lft.B.y, Rght.B.y);
                interp_l2r1_scalar(cL1.fs.B.z, cL0.fs.B.z, cR0.fs.B.z, Lft.B.z, Rght.B.z);
                if (myConfig.divergence_cleaning) {
                    interp_l2r1_scalar(cL1.fs.psi, cL0.fs.psi, cR0.fs.psi, Lft.psi, Rght.psi);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    interp_l2r1_scalar(cL1.fs.turb[it], cL0.fs.turb[it], cR0.fs.turb[it],
                                       Lft.turb[it], Rght.turb[it]);
                }
        }
        auto gL1 = &(cL1.fs.gas); auto gL0 = &(cL0.fs.gas); auto gR0 = &(cR0.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                if (myConfig.allow_reconstruction_for_species) {
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
                    Lft.gas.massf[] = gL0.massf[];
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
            break;
        case InterpolateOption.rhou:
            interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l2r1_scalar(gL1.u, gL0.u, gR0.u, Lft.gas.u, Rght.gas.u);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i]);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
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
        uint nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        l1r2_prepare(cL0Length, cR0Length, cR1Length);
        interp_l1r2_scalar(cL0.fs.vel.x, cR0.fs.vel.x, cR1.fs.vel.x, Lft.vel.x, Rght.vel.x);
        interp_l1r2_scalar(cL0.fs.vel.y, cR0.fs.vel.y, cR1.fs.vel.y, Lft.vel.y, Rght.vel.y);
        interp_l1r2_scalar(cL0.fs.vel.z, cR0.fs.vel.z, cR1.fs.vel.z, Lft.vel.z, Rght.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l1r2_scalar(cL0.fs.B.x, cR0.fs.B.x, cR1.fs.B.x, Lft.B.x, Rght.B.x);
                interp_l1r2_scalar(cL0.fs.B.y, cR0.fs.B.y, cR1.fs.B.y, Lft.B.y, Rght.B.y);
                interp_l1r2_scalar(cL0.fs.B.z, cR0.fs.B.z, cR1.fs.B.z, Lft.B.z, Rght.B.z);
                if (myConfig.divergence_cleaning) {
                    interp_l1r2_scalar(cL0.fs.psi, cR0.fs.psi, cR1.fs.psi, Lft.psi, Rght.psi);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    interp_l1r2_scalar(cL0.fs.turb[it], cR0.fs.turb[it], cR1.fs.turb[it],
                                       Lft.turb[it], Rght.turb[it]);
                }
        }
        auto gL0 = &(cL0.fs.gas); auto gR0 = &(cR0.fs.gas); auto gR1 = &(cR1.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                if (myConfig.allow_reconstruction_for_species) {
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
                    Lft.gas.massf[] = gL0.massf[];
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.T_modes[i], gR0.T_modes[i], gR1.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
            break;
        case InterpolateOption.rhou:
            interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho);
            interp_l1r2_scalar(gL0.u, gR0.u, gR1.u, Lft.gas.u, Rght.gas.u);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.u_modes[i], gR0.u_modes[i], gR1.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i]);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.T_modes[i], gR0.T_modes[i], gR1.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
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
                     ref FlowState Lft)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        if (myConfig.extrema_clipping) {
            // Not much that we can do with linear extrapolation
            // that doesn't produce an new extreme value.
            // Let the copy, made by the caller, stand.
            return;
        }
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cL1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cL0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        linear_extrap_prepare(cL0Length, cL1Length);
        Lft.vel.x = weight_scalar(cL0.fs.vel.x, cL1.fs.vel.x);
        Lft.vel.y = weight_scalar(cL0.fs.vel.y, cL1.fs.vel.y);
        Lft.vel.z = weight_scalar(cL0.fs.vel.z, cL1.fs.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                Lft.B.x = weight_scalar(cL0.fs.B.x, cL1.fs.B.x);
                Lft.B.y = weight_scalar(cL0.fs.B.y, cL1.fs.B.y);
                Lft.B.z = weight_scalar(cL0.fs.B.z, cL1.fs.B.z);
                if (myConfig.divergence_cleaning) {
                    Lft.psi = weight_scalar(cL0.fs.psi, cL1.fs.psi);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    Lft.turb[it] = weight_scalar(cL0.fs.turb[it], cL1.fs.turb[it]);
                }
        }
        auto gL1 = &(cL1.fs.gas); auto gL0 = &(cL0.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                if (myConfig.allow_reconstruction_for_species) {
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.T_modes[i] = weight_scalar(gL0.T_modes[i], gL1.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateLft("pT"));
            break;
        case InterpolateOption.rhou:
            Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            Lft.gas.u = weight_scalar(gL0.u, gL1.u);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.u_modes[i] = weight_scalar(gL0.u_modes[i], gL1.u_modes[i]);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.T_modes[i] = weight_scalar(gL0.T_modes[i], gL1.T_modes[i]);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
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
                     ref FlowState Rght)
    {
        auto gmodel = myConfig.gmodel;
        uint nsp = myConfig.n_species;
        auto nmodes = myConfig.n_modes;
        //
        if (myConfig.extrema_clipping) {
            // Not much that we can do with linear extrapolation
            // that doesn't produce an new extreme value.
            // Let the copy, made by the caller, stand.
            return;
        }
        // High-order reconstruction for some properties.
        if (myConfig.interpolate_in_local_frame) {
            // In the interface-local frame.
            cR0.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            cR1.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        }
        linear_extrap_prepare(cR0Length, cR1Length);
        Rght.vel.x = weight_scalar(cR0.fs.vel.x, cR1.fs.vel.x);
        Rght.vel.y = weight_scalar(cR0.fs.vel.y, cR1.fs.vel.y);
        Rght.vel.z = weight_scalar(cR0.fs.vel.z, cR1.fs.vel.z);
        version(MHD) {
            if (myConfig.MHD) {
                Rght.B.x = weight_scalar(cR0.fs.B.x, cR1.fs.B.x);
                Rght.B.y = weight_scalar(cR0.fs.B.y, cR1.fs.B.y);
                Rght.B.z = weight_scalar(cR0.fs.B.z, cR1.fs.B.z);
                if (myConfig.divergence_cleaning) {
                    Rght.psi = weight_scalar(cR0.fs.psi, cR1.fs.psi);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    Rght.turb[it] = weight_scalar(cR0.fs.turb[it], cR1.fs.turb[it]);
                }
        }
        auto gR0 = &(cR0.fs.gas); auto gR1 = &(cR1.fs.gas);
        version(multi_species_gas) {
            if (nsp > 1) {
                // Multiple species.
                if (myConfig.allow_reconstruction_for_species) {
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Rght.gas.T_modes[i] = weight_scalar(gR0.T_modes[i], gR1.T_modes[i]);
                    }
                } else {
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateRght("pT"));
            break;
        case InterpolateOption.rhou:
            Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            Rght.gas.u = weight_scalar(gR0.u, gR1.u);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Rght.gas.u_modes[i] = weight_scalar(gR0.u_modes[i], gR1.u_modes[i]);
                    }
                } else {
                    Rght.gas.u_modes[] = gR0.u_modes[];
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
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Rght.gas.T_modes[i] = weight_scalar(gR0.T_modes[i], gR1.T_modes[i]);
                    }
                } else {
                    Rght.gas.T_modes[] = gR0.T_modes[];
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

    @nogc final override void interp_l0r2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    {
        auto cR0 = f.right_cells[0]; auto cR1 = f.right_cells[1];
        number lR0, lR1;
        final switch (f.idir) {
        case IndexDirection.i:
            lR0 = cR0.iLength; lR1 = cR1.iLength;
            break;
        case IndexDirection.j:
            lR0 = cR0.jLength; lR1 = cR1.jLength;
            break;
        case IndexDirection.k:
            lR0 = cR0.kLength; lR1 = cR1.kLength;
            break;
        case IndexDirection.none:
            throw new Error("Invalid index direction.");
        }
        interp_l0r2(f, cR0, cR1, lR0, lR1, Rght);
    }

} // end class LegacyOneDInterpolator
