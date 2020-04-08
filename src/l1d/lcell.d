// lcell.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
//
// Classes for the faces and cell bodies that make up
// the internal structure of a GasSlug.
//
// PA Jacobs
// 2020-04-09
//
module lcell;

import geom;
import gas;
import gasflow;
import config;


class LFace {
public:
    double x;
    double area;
    
    this()
    {
        // Do nothing.
    }
    
} // end class LFace


class LCell {
public:
    double xmid;
    double volume;
    double vel;
    double L_bar;
    GasState gas;
    double shear_stress;
    double heat_flux;
    double dt_chem;
    double dt_therm;

    this(GasModel gm)
    {
        gas = new GasState(gm);
    }

} // end class LCell
