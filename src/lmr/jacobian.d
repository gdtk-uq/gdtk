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

    Matrix!number D;
    Matrix!number Dinv;
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
        local = new SMatrix!number();
	local.aa.length = size;
	local.ja.length = size;
	local.ia.length = ncells * nConserved + 1;
        D = new Matrix!number(nConserved,nConserved);
        Dinv = new Matrix!number(nConserved,nConserved);
    } // end size_local_matrix()

    void augment_with_dt(FVCell[] cells, size_t nConserved)
    {
        /*
          This method augments the Jacobian by adding the inverse pseudo-time term in the form A = 1/dt - dR/dU
          NB. we multiply the Jacobian by -1 to match our mathematical formulation of the implicit problem
         */

        foreach ( ref entry; local.aa) { entry *= -1; }
        number dtInv;
        foreach (i, cell; cells) {
            foreach (j; 0..nConserved) {
                dtInv = 1.0/cell.dt_local;
                size_t idx = i*nConserved + j;
                local[idx,idx] = local[idx,idx] + dtInv;
            }
        }

    } // end augment_with_dt()

} // end class FlowJacobian
