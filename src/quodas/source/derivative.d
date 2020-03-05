// derivative.d
module derivative;
import std.stdio;
import std.conv;
import std.math;
import config;
import number;
import block;
import finite_volume;
import linalg;

class Derivative {
public:
    // derivatives
    double[][] jac_V;           // Jacobian w.r.t. primitive (V) variables
    double[][] jac_U;           // Jacobian w.r.t. conservative (U) variables
    double[][] transform_V2U;   // transform matrix to convert between V and U Jacobians
    double[][] jac_VT;          // jac_V transpose
    double[][] jac_UT;          // jac_U transpose
    double[][] jac_inv;         // inverse of Jacobian matrix
    double[][] resid_sens;      // sensitivity of residual vector to the design variables
    double[] obj_sens;          // sensitivty of cost function (J) w.r.t. flow variables
    double[][] dRdq;
    double[][] dqdQ;
    double[][] Lext;

    FVInterface[] ifaces_cpy;
    FVCell[] cells_cpy;
    
    this(Params params, size_t ncells, size_t nifaces) {
        size_t nc = params.ncells;
        size_t np = nprimitive;
        size_t nd = ndesign;

        prep_matrix(this.jac_V, nc*np, nc*np);
        prep_matrix(this.jac_VT, nc*np, nc*np);
        prep_matrix(this.jac_inv, nc*np, nc*np);
        prep_matrix(this.resid_sens, nc*np, nd);
        prep_vector(this.obj_sens, nc*np);
        prep_matrix(this.dRdq, np, np);
        prep_matrix(this.dqdQ, np, np);
        prep_matrix(this.Lext, np, np);

        ifaces_cpy.length = nifaces;
        foreach (ref iface; this.ifaces_cpy) {
            iface = new FVInterface(0, to!number(0.0));
        }

        cells_cpy.length = nghost+ncells+nghost;
        foreach (ref cell; this.cells_cpy) {
            cell = new FVCell(0, to!number(0.0), to!number(0.0));
        }
    }
}
