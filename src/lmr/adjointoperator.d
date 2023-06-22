/** 
 * Functions for assembling adjoint operator.
 *
 * Authors: RJG, RBO and KAD
 * Date: 2023-05-19
 */

module lmr.adjointoperator;

import globalcoonfig;
import jacobian;

/*
 * Caller should have previously prepared via a call to
 * initNewtonKrylovSimulation.
 */

SMatrix!number formAdjointOperator(double perturbation, int iluFill=0)
{
    alias cfg = GlobalConfig;
    auto blk = localFluidBlocks[0];
    blk.initalize_jacobian(cfg.interpolation_order, perturbatio, iluFill);
    blk.evaluate_jacobian();
    auto adjOperator = blk.flowJacoobian.local.transpose();
    return adjOperator;
}

SMatrix!number formAdjointPreconditioner(SMatrix!number adjOperator)
{
    auto preconditioner = new SMatrix!number(adjOperator);
    decompILU0(preconditioner);
    return preconditioner;
}




