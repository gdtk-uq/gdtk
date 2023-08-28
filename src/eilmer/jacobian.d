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
    SMatrix!double local;       // residual sensitivities of local cells
                                // to perturbations of local cells
    SMatrix!double external;    // residual sensitivities of cells in neighbour blocks
                                // to perturbations of local cells

    double[][] dudU;            // used in the Jacobian boundary conditions to store the
                                // the sensitivity of the ghost cell flow states to the internal
                                // cell flow states
    number eps;                 // numerical perturbation parameter

    // indices used for the sparse matrix representation of the Jacobian
    size_t aa_idx = 0;
    size_t ja_idx = 0;
    size_t ia_idx = 0;

    Matrix!double D;
    Matrix!double Dinv;
    this (double sigma, size_t dimensions, size_t nConserved, int spatial_order, size_t nentry, size_t ncells)
    {
        this.spatial_order = (spatial_order > 1) ? 2 : 1;
        import std.stdio;
        eps = sigma;

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
        local = new SMatrix!double();
	local.aa.length = size;
	local.ja.length = size;
	local.ia.length = ncells * nConserved + 1;
        D = new Matrix!double(nConserved,nConserved);
        Dinv = new Matrix!double(nConserved,nConserved);
    } // end size_local_matrix()

    void augment_with_dt(FVCell[] cells, double dt, size_t ncells, size_t nConserved, bool dual_time_stepping, int temporal_order, double dt_physical, double dt_physical_old)
    {
        /*
          This method augments the Jacobian by adding the inverse pseudo-time term in the form A = 1/dt - dR/dU
          NB. we multiply the Jacobian by -1 to match our mathematical formulation of the implicit problem
         */

        foreach ( ref entry; local.aa) { entry *= -1; }
        double dtInv;
        foreach (i; 0..ncells) {
            foreach (j; 0..nConserved) {
                if (GlobalConfig.with_local_time_stepping) {
                    FVCell cell = cells[i];
                    dtInv = 1.0/cell.dt_local;
                } else {
                    dtInv = 1.0/dt;
                }
                if (dual_time_stepping) {
                    if (temporal_order == 1) {
                        dtInv = dtInv + 1.0/dt_physical;
                    } else {
                        double dtInv_physical  = 1.0/dt_physical + 1.0/dt_physical_old - dt_physical/(dt_physical_old*(dt_physical+dt_physical_old));
                        dtInv = dtInv + dtInv_physical;
                    }
                }
                ulong idx = i*nConserved + j;
                local[idx,idx] = local[idx,idx] + dtInv;
            }
        }

    } // end augment_with_dt()

} // end class FlowJacobian
