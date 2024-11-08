// interp.d -- Part of the Puffin and Lorikeet flow calculators.
//
// PA Jacobs
// 2022-02-02
//
module interp;

import std.format;
import std.math;

import geom;
import gas;
import gasdyn.gasflow;
import config;
import flow;
import face;
import cell;


@nogc
void interp_l2r2(Face2D f, ref FlowState2D fsL, ref FlowState2D fsR,
                 GasModel gmodel, bool clipFlag=false)
// Reconstruct flow states fsL,fsR, at the middle interface for a stencil of 4 cell-centred values.
{
    auto fsL1 = &(f.left_cells[1].fs); auto fsL0 = &(f.left_cells[0].fs);
    auto fsR0 = &(f.right_cells[0].fs); auto fsR1 = &(f.right_cells[1].fs);
    auto gasL1 = &(fsL1.gas); auto gasL0 = &(fsL0.gas);
    auto gasR0 = &(fsR0.gas); auto gasR1 = &(fsR1.gas);
    auto gasL = &(fsL.gas); auto gasR = &(fsR.gas);
    // First-order reconstruction is just a copy from the nearest cell centre.
    gasL.copy_values_from(*gasL0);
    gasR.copy_values_from(*gasR0);
    // We will interpolate only some properties.
    interp_l2r2_scalar(gasL1.rho, gasL0.rho, gasR0.rho, gasR1.rho, gasL.rho, gasR.rho, clipFlag);
    interp_l2r2_scalar(gasL1.u, gasL0.u, gasR0.u, gasR1.u, gasL.u, gasR.u, clipFlag);
    gmodel.update_thermo_from_rhou(*gasL);
    gmodel.update_sound_speed(*gasL);
    gmodel.update_thermo_from_rhou(*gasR);
    gmodel.update_sound_speed(*gasR);
    //
    interp_l2r2_scalar(fsL1.vel.x, fsL0.vel.x, fsR0.vel.x, fsR1.vel.x, fsL.vel.x, fsR.vel.x, clipFlag);
    interp_l2r2_scalar(fsL1.vel.y, fsL0.vel.y, fsR0.vel.y, fsR1.vel.y, fsL.vel.y, fsR.vel.y, clipFlag);
    return;
} // end interp_l2r2()

@nogc
void interp_l2r2_scalar(double qL1, double qL0, double qR0, double qR1,
                        ref double qL, ref double qR, bool clipFlag=false)
// Reconstruct values, qL,qR, at the middle interface for a stencil of 4 cell-centred values.
// Assume equal cell widths.
{
    import nm.limiters;
    // Set up differences and limiter values.
    double delLminus = (qL0 - qL1);
    double del = (qR0 - qL0);
    double delRplus = (qR1 - qR0);
    double sL = van_albada_limit1(delLminus, del);
    double sR = van_albada_limit1(del, delRplus);
    // The actual high-order reconstruction, possibly limited.
    qL = qL0 + sL * 0.125 * (3.0*del + delLminus);
    qR = qR0 - sR * 0.125 * (delRplus + 3.0*del);
    if (clipFlag) {
        // An extra limiting filter to ensure that we do not compute new extreme values.
        // This was introduced to deal with very sharp transitions in species.
        qL = clip_to_limits(qL, qL0, qR0);
        qR = clip_to_limits(qR, qL0, qR0);
    }
} // end of interp_l2r2_scalar()
