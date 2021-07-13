/* jacobian.d
 *
 * Numerical Jacobian class
 *
 * Author: Kyle Damm
 * First code: 2020-11-03
 */

module jacobian;

import std.conv;

import nm.complex;
import nm.number;
import nm.smla;
import nm.bbla;

import globalconfig;
import fvcell;
import fvinterface;


class FlowJacobian {

public:

    bool wrtConserved = true;
    size_t spatial_order;
    SMatrix!number local;       // residual sensitivities of local cells
                                // to perturbations of local cells
    SMatrix!number external;    // residual sensitivities of cells in neighbour blocks
                                // to perturbations of local cells

    number[][] dudU;            // used in the Jacobian boundary conditions to store the
                                // the sensitivity of the ghost cell flow states to the internal
                                // cell flow states
    number eps;                 // numerical perturbation parameter

    // indices used for the sparse matrix representation of the Jacobian
    size_t aa_idx = 0;
    size_t ja_idx = 0;
    size_t ia_idx = 0;

    number[] diagonal;
    Matrix!number D;
    Matrix!number Dinv;
    Matrix!number LDinv;
    number[] x;
    this (size_t dimensions, size_t nConserved, int spatial_order, size_t nentry, size_t ncells)
    {
        this.spatial_order = (spatial_order > 1) ? 2 : 1;
        import std.stdio;
        version(complex_numbers) { eps = complex(0.0, 1.0e-30); }
        else { eps = 1.0e-10; }

        size_local_matrix(nConserved, nentry, ncells);

        dudU.length = nConserved;
        foreach (ref val; dudU) {val.length = nConserved; }
    }

    void prepare_crs_indexes()
    {
        /*
          Prepares the indices used to index into the pre-sized sparse matrix
          representation of the Jacobian (Compressed Row Storage (CRS) format).
          This method needs to be called each time the Jacobian is either first
          evaluated or re-evaluated.
         */

        // zero the inidices
        aa_idx = 0; ja_idx = 0; ia_idx = 0;

        // the first entry will always be filled, let's prepare for this entry
        local.ia[ia_idx] = 0;
        ia_idx += 1;
    } // end prepare_crs_indexes()

    void size_local_matrix(size_t nConserved, size_t nentry, size_t ncells)
    {
        // reserve memory for the local entries
        size_t size = nentry * nConserved * nConserved;
        local = new SMatrix!number();
	local.aa.length = size;
	local.ja.length = size;
	local.ia.length = ncells * nConserved + 1;
        D = new Matrix!number(nConserved,nConserved);
        Dinv = new Matrix!number(nConserved,nConserved);
        LDinv = new Matrix!number(nConserved,nConserved);
        diagonal.length = ncells*nConserved*nConserved;
        x.length = ncells*nConserved;
    } // end size_local_matrix()

    void prepare_jacobi_preconditioner(FVCell[] cells, double dt, size_t ncells, size_t nConserved, size_t fill_in_level = 0)
    {
        /*
          This method prepares the flow Jacobian for use as a Jacobi precondition matrix.
          It does this in 3 steps:

          1. multiply Jacobian by -1 to match our mathematical formulation of the implicit problem
          2. add 1/dt to the diagonal
          3. invert the block diagonals

          TODO: move this operation to a more appropriate location.
         */

        foreach ( ref entry; local.aa) { entry *= -1; }
        number dtInv;
        foreach (i; 0..ncells) {
            foreach (j; 0..nConserved) {
                if (GlobalConfig.with_local_time_stepping) {
                    FVCell cell = cells[i];
                    dtInv = 1.0/cell.dt_local;
                } else {
                    dtInv = 1.0/dt;
                }
                ulong idx = i*nConserved + j;
                local[idx,idx] = local[idx,idx] + dtInv;
            }
        }

        nm.smla.invert_diagonal(local, to!int(nConserved), D, Dinv);
    } // end prepare_jacobi_preconditioner()

    void prepare_ilu_preconditioner(FVCell[] cells, double dt, size_t ncells, size_t nConserved, size_t fill_in_level = 0)
    {
        /*
          This method prepares the flow Jacobian for use as a precondition matrix.
          It does this in 3 steps:

          1. multiply Jacobian by -1 to match our mathematical formulation of the implicit problem
          2. add 1/dt to the diagonal
          3. perform an ILU decomposition

          TODO: move this operation to a more appropriate location.
         */

        foreach ( ref entry; local.aa) { entry *= -1; }
        number dtInv;
        foreach (i; 0..ncells) {
            foreach (j; 0..nConserved) {
                if (GlobalConfig.with_local_time_stepping) {
                    FVCell cell = cells[i];
                    dtInv = 1.0/cell.dt_local;
                } else {
                    dtInv = 1.0/dt;
                }        // TODO: local-time-stepping
                ulong idx = i*nConserved + j;
                local[idx,idx] = local[idx,idx] + dtInv;
            }
        }

        if (fill_in_level > 0) { decompILUp(local, to!int(fill_in_level)); }
        else { decompILU0(local); }

    } // end prepare_preconditioner()

    void prepare_sgs_preconditioner(FVCell[] cells, double dt, size_t ncells, size_t nConserved, size_t fill_in_level = 0)
    {
        /*
          This method prepares the flow Jacobian for use as a precondition matrix.
          It does this in 3 steps:

          1. multiply Jacobian by -1 to match our mathematical formulation of the implicit problem
          2. add 1/dt to the diagonal

          TODO: move this operation to a more appropriate location.
         */

        foreach ( ref entry; local.aa) { entry *= -1; }
        number dtInv;
        foreach (i; 0..ncells) {
            foreach (j; 0..nConserved) {
                if (GlobalConfig.with_local_time_stepping) {
                    FVCell cell = cells[i];
                    dtInv = 1.0/cell.dt_local;
                } else {
                    dtInv = 1.0/dt;
                }
                ulong idx = i*nConserved + j;
                local[idx,idx] = local[idx,idx] + dtInv;
            }
        }

        foreach (k; 0..ncells) {
            foreach (i; 0..nConserved) {
                int idx = to!int(k*nConserved + i);
                foreach (j; 0..nConserved) {
                    int jdx = to!int(k*nConserved + j);
                    diagonal[k*nConserved*nConserved + i*nConserved + j] = local[idx,jdx];
                }
            }
        }

        invert_diagonal(local, to!int(nConserved), D, Dinv);

    } // end prepare_sgs_preconditioner()

    void prepare_sgsr_preconditioner(FVCell[] cells, double dt, size_t ncells, size_t nConserved, size_t fill_in_level = 0)
    {
        /*
          This method prepares the flow Jacobian for use as a precondition matrix.
          It does this in 3 steps:

          1. multiply Jacobian by -1 to match our mathematical formulation of the implicit problem
          2. add 1/dt to the diagonal
          3. perform an ILU decomposition

          we also scale the LU decomposition needed for the transpose solve method
          TODO: move this operation to a more appropriate location.
         */

        foreach ( ref entry; local.aa) { entry *= -1; }
        number dtInv;
        foreach (i; 0..ncells) {
            foreach (j; 0..nConserved) {
                if (GlobalConfig.with_local_time_stepping) {
                    FVCell cell = cells[i];
                    dtInv = 1.0/cell.dt_local;
                } else {
                    dtInv = 1.0/dt;
                }        // TODO: local-time-stepping
                ulong idx = i*nConserved + j;
                local[idx,idx] = local[idx,idx] + dtInv;
            }
        }

        invert_diagonal(local, to!int(nConserved), D, Dinv);

    } // end prepare_sgsr_preconditioner()

} // end class FlowJacobian
