// onedinterp.d
// One-dimensional interpolation/reconstruction of flow field.
//
// See MBCNS workbook 2000/2 page 36 (26-Jan-2001) for formulation.
// and MBCNS workbook 2005/Apr page 36 for new index labels

module onedinterp;

import std.math;
import std.stdio;
import std.conv;
import ntypes.complex;
import nm.number;
import nm.limiters;

import gas;
import globalconfig;
import flowstate;
import fvinterface;
import lmr.fluidfvcell;

//------------------------------------------------------------------------------

/+
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
+/

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

struct L2R2InterpData {
    number lenL0_;
    number lenR0_;
    number aL0;
    number aR0;
    number two_over_lenL0_plus_lenL1;
    number two_over_lenR0_plus_lenL0;
    number two_over_lenR1_plus_lenR0;
    number two_lenL0_plus_lenL1;
    number two_lenR0_plus_lenR1;

    @nogc
    void set(number lenL1, number lenL0, number lenR0, number lenR1){
        lenL0_ = lenL0;
        lenR0_ = lenR0;
        aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
        two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
        two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
        two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
    }
}

struct L3R3InterpData {
    number wL_L2, wL_L1, wL_L0, wL_R0, wL_R1;
    number wR_L1, wR_L0, wR_R0, wR_R1, wR_R2;

    @nogc
    void set(number lenL2, number lenL1, number lenL0,
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
    }
}

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
    number wL_L2, wL_L1, wL_L0, wL_R0, wL_R1;
    number wR_L1, wR_L0, wR_R0, wR_R1, wR_R2;
    LocalConfig myConfig;

public:
    @disable this();

    this(LocalConfig myConfig)
    {
        this.myConfig = myConfig;
    }

    // The steady-state solver wants to override the order of interpolation
    // in various parts of its calculation.
    @nogc int get_interpolation_order()
    {
        return myConfig.interpolation_order;
    }

    @nogc void set_interpolation_order(int order)
    {
        myConfig.interpolation_order = order;
    }

    @nogc
    void  interp(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
    // Top-level interpolation function delegates to the specific functions, below.
    {
        FluidFVCell cL0 = (f.left_cells.length > 0) ? f.left_cells[0] : f.right_cells[0];
        FluidFVCell cR0 = (f.right_cells.length > 0) ? f.right_cells[0]: f.left_cells[0];
        Lft.copy_values_from(cL0.fs);
        Rght.copy_values_from(cR0.fs);

        number beta = 1.0;
        bool one_sided_interp = (f.left_cells.length == 0 || f.right_cells.length == 0);
        if (myConfig.apply_heuristic_pressure_based_limiting && !one_sided_interp) {
            // Optional pressure-based augmentation factor on the base limiter value,
            // intended for use on problems with strong bow shocks.
            // The equations are taken from pg. 33-34 of
            //     Simulation and Dynamics of Hypersonic Turbulent Combustion
            //     N. Gibbons, University of Queensland, 2019
            //
            number pmin = cL0.fs.gas.p;
            number pmax = cL0.fs.gas.p;
            foreach (arr; [f.left_cells, f.right_cells]) {
                foreach (cell; arr) {
                    pmin = fmin(pmin, cell.fs.gas.p);
                    pmax = fmax(pmax, cell.fs.gas.p);
                }
            }
            number alpha = fabs(pmax-pmin)/pmin;
            beta = 1.0/(1.0+alpha*alpha*alpha*alpha);
        }

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
            interp_l3r3(f, cL2, cL1, cL0, cR0, cR1, cR2, lL2, lL1, lL0, lR0, lR1, lR2, Lft, Rght, beta);
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
                interp_l1r2(f, cL0, cR0, cR1, lL0, lR0, lR1, Lft, Rght, beta);
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
                interp_l2r1(f, cL1, cL0, cR0, lL1, lL0, lR0, Lft, Rght, beta);
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
                interp_l2r2(f, cL1, cL0, cR0, cR1, lL1, lL0, lR0, lR1, Lft, Rght, beta);
            } else {
                throw new Error("Stencils not suitable for standard interpolation.");
            }
        }
    } // end interp()

    //------------------------------------------------------------------------------

    @nogc
    void l3r3_prepare(number lenL2, number lenL1, number lenL0,
                      number lenR0, number lenR1, number lenR2)
    // Set up intermediate data (Lagrange interpolation weights)
    // that depend only on the cell geometry.
    // They will remain constant when reconstructing the different scalar fields
    // over the same set of cells.
    {
        pragma(inline, true);
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

    @nogc
    void interp_l3r3_scalar(number qL2, number qL1, number qL0,
                            number qR0, number qR1, number qR2,
                            ref number qL, ref number qR, number beta)
    {
        pragma(inline, true);
        // Set up differences and limiter values.
        number delLminus = (qL0 - qL1);
        number del = (qR0 - qL0);
        number delRplus = (qR1 - qR0);
        // Presume unlimited high-order reconstruction.
        number sL = 1.0;
        number sR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            sL = (delLminus*del + fabs(delLminus*del) + myConfig.epsilon_van_albada) /
                (delLminus*delLminus + del*del + myConfig.epsilon_van_albada);
            sR = (del*delRplus + fabs(del*delRplus) + myConfig.epsilon_van_albada) /
                (del*del + delRplus*delRplus + myConfig.epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        qL = qL0 + beta * sL * (wL_L2*qL2 + wL_L1*qL1 + (wL_L0-1.0)*qL0 + wL_R0*qR0 + wL_R1*qR1);
        qR = qR0 + beta * sR * (wR_L1*qL1 + wR_L0*qL0 + (wR_R0-1.0)*qR0 + wR_R1*qR1 + wR_R2*qR2);
        if (myConfig.extrema_clipping) {
            // An extra limiting filter to ensure that we do not compute new extreme values.
            // This was introduced to deal with very sharp transitions in species.
            qL = clip_to_limits(qL, qL0, qR0);
            qR = clip_to_limits(qR, qL0, qR0);
        }
    } // end of interp_l3r3_scalar()

    @nogc
    void l2r2_prepare(number lenL1, number lenL0, number lenR0, number lenR1)
    // Set up intermediate data that depend only on the cell geometry.
    // They will remain constant when reconstructing the different scalar fields
    // over the same set of cells.
    // For piecewise parabolic reconstruction, see PJ workbook notes Jan 2001.
    {
        pragma(inline, true);
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

    @nogc
    void interp_l2r2_scalar(number qL1, number qL0, number qR0, number qR1,
                            ref number qL, ref number qR, number beta)
    {
        pragma(inline, true);
        // Set up differences and limiter values.
        number delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        number del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        number delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        // Presume unlimited high-order reconstruction.
        number sL = 1.0;
        number sR = 1.0;
        if (myConfig.apply_limiter) {
            // val Albada limiter as per Ian Johnston's thesis.
            sL = (delLminus*del + fabs(delLminus*del) + myConfig.epsilon_van_albada) /
                (delLminus*delLminus + del*del + myConfig.epsilon_van_albada);
            sR = (del*delRplus + fabs(del*delRplus) + myConfig.epsilon_van_albada) /
                (del*del + delRplus*delRplus + myConfig.epsilon_van_albada);
        }
        // The actual high-order reconstruction, possibly limited.
        qL = qL0 + beta * sL * aL0 * (del * two_lenL0_plus_lenL1 + delLminus * lenR0_);
        qR = qR0 - beta * sR * aR0 * (delRplus * lenL0_ + del * two_lenR0_plus_lenR1);
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
        pragma(inline, true);
        lenL0_ = lenL0;
        lenR0_ = lenR0;
        aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
        two_over_lenL0_plus_lenL1 = 2.0 / (lenL0 + lenL1);
        two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        two_lenL0_plus_lenL1 = (2.0*lenL0 + lenL1);
        linear_interp_prepare(lenL0, lenR0);
    } // end l2r1_prepare()

    @nogc void interp_l2r1_scalar(number qL1, number qL0, number qR0, ref number qL, ref number qR, number beta)
    {
        pragma(inline, true);
        number delLminus = (qL0 - qL1) * two_over_lenL0_plus_lenL1;
        number del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        number sL = 1.0;
        if (myConfig.apply_limiter) {
            sL = (delLminus*del + fabs(delLminus*del) + myConfig.epsilon_van_albada) /
                (delLminus*delLminus + del*del + myConfig.epsilon_van_albada);
        }
        qL = qL0 + beta * sL * aL0 * (del * two_lenL0_plus_lenL1 + delLminus * lenR0_);
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
        pragma(inline, true);
        lenL0_ = lenL0;
        lenR0_ = lenR0;
        aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
        two_over_lenR0_plus_lenL0 = 2.0 / (lenR0 + lenL0);
        two_over_lenR1_plus_lenR0 = 2.0 / (lenR1 + lenR0);
        two_lenR0_plus_lenR1 = (2.0*lenR0 + lenR1);
        linear_interp_prepare(lenL0, lenR0);
    } // end l1r2_prepare()

    @nogc void interp_l1r2_scalar(number qL0, number qR0, number qR1, ref number qL, ref number qR, number beta)
    {
        pragma(inline, true);
        number del = (qR0 - qL0) * two_over_lenR0_plus_lenL0;
        number delRplus = (qR1 - qR0) * two_over_lenR1_plus_lenR0;
        number sR = 1.0;
        if (myConfig.apply_limiter) {
            sR = (del*delRplus + fabs(del*delRplus) + myConfig.epsilon_van_albada) /
                (del*del + delRplus*delRplus + myConfig.epsilon_van_albada);
        }
        qR = qR0 - beta * sR * aR0 * (delRplus * lenL0_ + del * two_lenR0_plus_lenR1);
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
        pragma(inline, true);
        // Set up weights for a linear combination if q0 and q1
        // assuming that we are extrapolating past q0 to the boundary len0/2 away.
        w0 = (2.0*len0 + len1)/(len0+len1);
        w1 = -len0/(len0+len1);
    }

    @nogc void linear_interp_prepare(number len0, number len1)
    {
        pragma(inline, true);
        // Set up weights for a linear combination if q0 and q1
        // assuming that we are interpolating to the cell-face between
        // the cell-centre points.
        w0 = len1/(len0+len1);
        w1 = len0/(len0+len1);
    }

    @nogc number weight_scalar(number q0, number q1)
    {
        pragma(inline, true);
        // The weights for interpolation or extrapolation may be used.
        number q = q0*w0 + q1*w1;
        if (myConfig.extrema_clipping) { q = clip_to_limits(q, q0, q1); }
        return q;
    }

    //------------------------------------------------------------------------------

    @nogc
    void interp_l3r3(ref FVInterface IFace,
                     ref FluidFVCell cL2, ref FluidFVCell cL1, ref FluidFVCell cL0,
                     ref FluidFVCell cR0, ref FluidFVCell cR1, ref FluidFVCell cR2,
                     number cL2Length, number cL1Length, number cL0Length,
                     number cR0Length, number cR1Length, number cR2Length,
                     ref FlowState Lft, ref FlowState Rght, number beta)
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
                           Lft.vel.x, Rght.vel.x, beta);
        interp_l3r3_scalar(cL2.fs.vel.y, cL1.fs.vel.y, cL0.fs.vel.y,
                           cR0.fs.vel.y, cR1.fs.vel.y, cR2.fs.vel.y,
                           Lft.vel.y, Rght.vel.y, beta);
        interp_l3r3_scalar(cL2.fs.vel.z, cL1.fs.vel.z, cL0.fs.vel.z,
                           cR0.fs.vel.z, cR1.fs.vel.z, cR2.fs.vel.z,
                           Lft.vel.z, Rght.vel.z, beta);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l3r3_scalar(cL2.fs.B.x, cL1.fs.B.x, cL0.fs.B.x,
                                   cR0.fs.B.x, cR1.fs.B.x, cR2.fs.B.x,
                                   Lft.B.x, Rght.B.x, beta);
                interp_l3r3_scalar(cL2.fs.B.y, cL1.fs.B.y, cL0.fs.B.y,
                                   cR0.fs.B.y, cR1.fs.B.y, cR2.fs.B.y,
                                   Lft.B.y, Rght.B.y, beta);
                interp_l3r3_scalar(cL2.fs.B.z, cL1.fs.B.z, cL0.fs.B.z,
                                   cR0.fs.B.z, cR1.fs.B.z, cR2.fs.B.z,
                                   Lft.B.z, Rght.B.z, beta);
                if (myConfig.divergence_cleaning) {
                    interp_l3r3_scalar(cL2.fs.psi, cL1.fs.psi, cL0.fs.psi,
                                       cR0.fs.psi, cR1.fs.psi, cR2.fs.psi,
                                       Lft.psi, Rght.psi, beta);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                interp_l3r3_scalar(cL2.fs.turb[it], cL1.fs.turb[it], cL0.fs.turb[it],
                                   cR0.fs.turb[it], cR1.fs.turb[it], cR2.fs.turb[it],
                                   Lft.turb[it], Rght.turb[it], beta);
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
                // Reconstruct species densities
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        interp_l3r3_scalar(gL2.rho_s[isp], gL1.rho_s[isp], gL0.rho_s[isp],
                                           gR0.rho_s[isp], gR1.rho_s[isp], gR2.rho_s[isp],
                                           Lft.gas.rho_s[isp], Rght.gas.rho_s[isp], beta);
                    }
                } else {
                    Lft.gas.rho_s[]  = gL0.rho_s[];
                    Rght.gas.rho_s[] = gR0.rho_s[];
                }
            }
        }

        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            interp_l3r3_scalar(gL2.p, gL1.p, gL0.p, gR0.p, gR1.p, gR2.p, Lft.gas.p, Rght.gas.p, beta);
            interp_l3r3_scalar(gL2.T, gL1.T, gL0.T, gR0.T, gR1.T, gR2.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.T_modes[i], gL1.T_modes[i], gL0.T_modes[i],
                                           gR0.T_modes[i], gR1.T_modes[i], gR2.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
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
                    interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l3r3_scalar(gL2.u, gL1.u, gL0.u, gR0.u, gR1.u, gR2.u, Lft.gas.u, Rght.gas.u, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.u_modes[i], gL1.u_modes[i], gL0.u_modes[i],
                                           gR0.u_modes[i], gR1.u_modes[i], gR2.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
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
                    interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l3r3_scalar(gL2.p, gL1.p, gL0.p, gR0.p, gR1.p, gR2.p, Lft.gas.p, Rght.gas.p, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.u_modes[i], gL1.u_modes[i], gL0.u_modes[i],
                                           gR0.u_modes[i], gR1.u_modes[i], gR2.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhop"));
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
                    interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l3r3_scalar(gL2.rho, gL1.rho, gL0.rho, gR0.rho, gR1.rho, gR2.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l3r3_scalar(gL2.T, gL1.T, gL0.T, gR0.T, gR1.T, gR2.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l3r3_scalar(gL2.T_modes[i], gL1.T_modes[i], gL0.T_modes[i],
                                           gR0.T_modes[i], gR1.T_modes[i], gR2.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
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
                     ref FluidFVCell cL1, ref FluidFVCell cL0, ref FluidFVCell cR0, ref FluidFVCell cR1,
                     number cL1Length, number cL0Length,
                     number cR0Length, number cR1Length,
                     ref FlowState Lft, ref FlowState Rght, number beta)
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
                           Lft.vel.x, Rght.vel.x, beta);
        interp_l2r2_scalar(cL1.fs.vel.y, cL0.fs.vel.y, cR0.fs.vel.y, cR1.fs.vel.y,
                           Lft.vel.y, Rght.vel.y, beta);
        interp_l2r2_scalar(cL1.fs.vel.z, cL0.fs.vel.z, cR0.fs.vel.z, cR1.fs.vel.z,
                           Lft.vel.z, Rght.vel.z, beta);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l2r2_scalar(cL1.fs.B.x, cL0.fs.B.x, cR0.fs.B.x, cR1.fs.B.x,
                                   Lft.B.x, Rght.B.x, beta);
                interp_l2r2_scalar(cL1.fs.B.y, cL0.fs.B.y, cR0.fs.B.y, cR1.fs.B.y,
                                   Lft.B.y, Rght.B.y, beta);
                interp_l2r2_scalar(cL1.fs.B.z, cL0.fs.B.z, cR0.fs.B.z, cR1.fs.B.z,
                                   Lft.B.z, Rght.B.z, beta);
                if (myConfig.divergence_cleaning) {
                    interp_l2r2_scalar(cL1.fs.psi, cL0.fs.psi, cR0.fs.psi, cR1.fs.psi,
                                       Lft.psi, Rght.psi, beta);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    interp_l2r2_scalar(cL1.fs.turb[it], cL0.fs.turb[it], cR0.fs.turb[it],
                                       cR1.fs.turb[it], Lft.turb[it], Rght.turb[it], beta);
                }
        }
        auto gL1 = &(cL1.fs.gas); // Avoid construction of another object.
        auto gL0 = &(cL0.fs.gas);
        auto gR0 = &(cR0.fs.gas);
        auto gR1 = &(cR1.fs.gas);

        version(multi_species_gas) {
            if (nsp > 1) {
                // Reconstruct species densities
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        interp_l2r2_scalar(gL1.rho_s[isp], gL0.rho_s[isp], gR0.rho_s[isp], gR1.rho_s[isp],
                                           Lft.gas.rho_s[isp], Rght.gas.rho_s[isp], beta);
                    }
                } else {
                    Lft.gas.rho_s[]  = gL0.rho_s[];
                    Rght.gas.rho_s[] = gR0.rho_s[];
                }
            }
        }

        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            interp_l2r2_scalar(gL1.p, gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p, beta);
            interp_l2r2_scalar(gL1.T, gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           gR1.T_modes[i], Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
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
                    interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l2r2_scalar(gL1.u, gL0.u, gR0.u, gR1.u, Lft.gas.u, Rght.gas.u, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                           gR1.u_modes[i], Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
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
                    interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l2r2_scalar(gL1.p, gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                           gR1.u_modes[i], Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhop"));
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
                    interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l2r2_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l2r2_scalar(gL1.T, gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r2_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           gR1.T_modes[i], Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
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
                     ref FluidFVCell cL1, ref FluidFVCell cL0, ref FluidFVCell cR0,
                     number cL1Length, number cL0Length, number cR0Length,
                     ref FlowState Lft, ref FlowState Rght, number beta)
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
        interp_l2r1_scalar(cL1.fs.vel.x, cL0.fs.vel.x, cR0.fs.vel.x, Lft.vel.x, Rght.vel.x, beta);
        interp_l2r1_scalar(cL1.fs.vel.y, cL0.fs.vel.y, cR0.fs.vel.y, Lft.vel.y, Rght.vel.y, beta);
        interp_l2r1_scalar(cL1.fs.vel.z, cL0.fs.vel.z, cR0.fs.vel.z, Lft.vel.z, Rght.vel.z, beta);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l2r1_scalar(cL1.fs.B.x, cL0.fs.B.x, cR0.fs.B.x, Lft.B.x, Rght.B.x, beta);
                interp_l2r1_scalar(cL1.fs.B.y, cL0.fs.B.y, cR0.fs.B.y, Lft.B.y, Rght.B.y, beta);
                interp_l2r1_scalar(cL1.fs.B.z, cL0.fs.B.z, cR0.fs.B.z, Lft.B.z, Rght.B.z, beta);
                if (myConfig.divergence_cleaning) {
                    interp_l2r1_scalar(cL1.fs.psi, cL0.fs.psi, cR0.fs.psi, Lft.psi, Rght.psi, beta);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    interp_l2r1_scalar(cL1.fs.turb[it], cL0.fs.turb[it], cR0.fs.turb[it],
                                       Lft.turb[it], Rght.turb[it], beta);
                }
        }
        auto gL1 = &(cL1.fs.gas); auto gL0 = &(cL0.fs.gas); auto gR0 = &(cR0.fs.gas);

        version(multi_species_gas) {
            if (nsp > 1) {
                // Reconstruct species densities
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        interp_l2r1_scalar(gL1.rho_s[isp], gL0.rho_s[isp], gR0.rho_s[isp],
                                           Lft.gas.rho_s[isp], Rght.gas.rho_s[isp], beta);
                    }
                } else {
                    Lft.gas.rho_s[]  = gL0.rho_s[];
                    Rght.gas.rho_s[] = gR0.rho_s[];
                }
            }
        }

        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            interp_l2r1_scalar(gL1.p, gL0.p, gR0.p, Lft.gas.p, Rght.gas.p, beta);
            interp_l2r1_scalar(gL1.T, gL0.T, gR0.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
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
                    interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l2r1_scalar(gL1.u, gL0.u, gR0.u, Lft.gas.u, Rght.gas.u, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
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
                    interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l2r1_scalar(gL1.p, gL0.p, gR0.p, Lft.gas.p, Rght.gas.p, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.u_modes[i], gL0.u_modes[i], gR0.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhop"));
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
                    interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l2r1_scalar(gL1.rho, gL0.rho, gR0.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l2r1_scalar(gL1.T, gL0.T, gR0.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l2r1_scalar(gL1.T_modes[i], gL0.T_modes[i], gR0.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
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
                     ref FluidFVCell cL0, ref FluidFVCell cR0, ref FluidFVCell cR1,
                     number cL0Length, number cR0Length, number cR1Length,
                     ref FlowState Lft, ref FlowState Rght, number beta)
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
        interp_l1r2_scalar(cL0.fs.vel.x, cR0.fs.vel.x, cR1.fs.vel.x, Lft.vel.x, Rght.vel.x, beta);
        interp_l1r2_scalar(cL0.fs.vel.y, cR0.fs.vel.y, cR1.fs.vel.y, Lft.vel.y, Rght.vel.y, beta);
        interp_l1r2_scalar(cL0.fs.vel.z, cR0.fs.vel.z, cR1.fs.vel.z, Lft.vel.z, Rght.vel.z, beta);
        version(MHD) {
            if (myConfig.MHD) {
                interp_l1r2_scalar(cL0.fs.B.x, cR0.fs.B.x, cR1.fs.B.x, Lft.B.x, Rght.B.x, beta);
                interp_l1r2_scalar(cL0.fs.B.y, cR0.fs.B.y, cR1.fs.B.y, Lft.B.y, Rght.B.y, beta);
                interp_l1r2_scalar(cL0.fs.B.z, cR0.fs.B.z, cR1.fs.B.z, Lft.B.z, Rght.B.z, beta);
                if (myConfig.divergence_cleaning) {
                    interp_l1r2_scalar(cL0.fs.psi, cR0.fs.psi, cR1.fs.psi, Lft.psi, Rght.psi, beta);
                }
            }
        }
        version(turbulence) {
                foreach (it; 0 .. myConfig.turb_model.nturb){
                    interp_l1r2_scalar(cL0.fs.turb[it], cR0.fs.turb[it], cR1.fs.turb[it],
                                       Lft.turb[it], Rght.turb[it], beta);
                }
        }
        auto gL0 = &(cL0.fs.gas); auto gR0 = &(cR0.fs.gas); auto gR1 = &(cR1.fs.gas);

        version(multi_species_gas) {
            if (nsp > 1) {
                // Reconstruct species densities
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        interp_l1r2_scalar(gL0.rho_s[isp], gR0.rho_s[isp], gR1.rho_s[isp],
                                           Lft.gas.rho_s[isp], Rght.gas.rho_s[isp], beta);
                    }
                } else {
                    Lft.gas.rho_s[]  = gL0.rho_s[];
                    Rght.gas.rho_s[] = gR0.rho_s[];
                }
            }
        }

        // Interpolate on two of the thermodynamic quantities,
        // and fill in the rest based on an EOS call.
        final switch (myConfig.thermo_interpolator) {
        case InterpolateOption.pt:
            interp_l1r2_scalar(gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p, beta);
            interp_l1r2_scalar(gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.T_modes[i], gR0.T_modes[i], gR1.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
                    }
                } else {
                    Lft.gas.T_modes[] = gL0.T_modes[];
                    Rght.gas.T_modes[] = gR0.T_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("pT"));
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
                    interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l1r2_scalar(gL0.u, gR0.u, gR1.u, Lft.gas.u, Rght.gas.u, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.u_modes[i], gR0.u_modes[i], gR1.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhou"));
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
                    interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l1r2_scalar(gL0.p, gR0.p, gR1.p, Lft.gas.p, Rght.gas.p, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.u_modes[i], gR0.u_modes[i], gR1.u_modes[i],
                                           Lft.gas.u_modes[i], Rght.gas.u_modes[i], beta);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateBoth("rhop"));
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
                    interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
                }
            } else {
                interp_l1r2_scalar(gL0.rho, gR0.rho, gR1.rho, Lft.gas.rho, Rght.gas.rho, beta);
            }
            interp_l1r2_scalar(gL0.T, gR0.T, gR1.T, Lft.gas.T, Rght.gas.T, beta);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        interp_l1r2_scalar(gL0.T_modes[i], gR0.T_modes[i], gR1.T_modes[i],
                                           Lft.gas.T_modes[i], Rght.gas.T_modes[i], beta);
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
                     ref FluidFVCell cL1, ref FluidFVCell cL0,
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
                // Reconstruct species densities
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        Lft.gas.rho_s[isp] = weight_scalar(gL0.rho_s[isp], gL1.rho_s[isp]);
                    }
                } else {
                    Lft.gas.rho_s[]  = gL0.rho_s[];
                }
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
                    Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
                }
            } else {
                Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            }
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
                    Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
                }
            } else {
                Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            }
            Lft.gas.p = weight_scalar(gL0.p, gL1.p);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Lft.gas.u_modes[i] = weight_scalar(gL0.u_modes[i], gL1.u_modes[i]);
                    }
                } else {
                    Lft.gas.u_modes[] = gL0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateLft("rhop"));
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
                    Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
                }
            } else {
                Lft.gas.rho = weight_scalar(gL0.rho, gL1.rho);
            }
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
                     ref FluidFVCell cR0, ref FluidFVCell cR1,
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
                // Reconstruct species densities
                if (myConfig.allow_reconstruction_for_species) {
                    foreach (isp; 0 .. nsp) {
                        Rght.gas.rho_s[isp] = weight_scalar(gR0.rho_s[isp], gR1.rho_s[isp]);
                    }
                } else {
                    Rght.gas.rho_s[] = gR0.rho_s[];
                }
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
                    Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
                }
            } else {
                Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            }
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
                    Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
                }
            } else {
                Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            }
            Rght.gas.p = weight_scalar(gR0.p, gR1.p);
            version(multi_T_gas) {
                if (myConfig.allow_reconstruction_for_energy_modes) {
                    foreach (i; 0 .. nmodes) {
                        Rght.gas.u_modes[i] = weight_scalar(gR0.u_modes[i], gR1.u_modes[i]);
                    }
                } else {
                    Rght.gas.u_modes[] = gR0.u_modes[];
                }
            }
            mixin(codeForThermoUpdateRght("rhop"));
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
                    Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
                }
            } else {
                Rght.gas.rho = weight_scalar(gR0.rho, gR1.rho);
            }
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

    // The following variant of interpolation with right-cells-only is used in
    // the shock-fitting calculation, over in module grid_motion_shock_fitting.d.
    @nogc
    void interp_l0r2(ref FVInterface f, ref FlowState Lft, ref FlowState Rght)
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

} // end class OneDInterpolator
