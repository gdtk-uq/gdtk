// lcell.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
//
// Classes for the faces and cell bodies that make up
// the internal structure of a GasSlug.
//
// PA Jacobs
// 2020-04-09
//
module lcell;

import std.format;

import geom;
import gas;
import gasflow;
import config;


class LFace {
public:
    double x;
    double area;

    double x0;       // initial position for time step
    double[2] dxdt;  // time derivatives, predictor-corrector levels
    this()
    {
        // Do nothing.
    }
    
} // end class LFace


class BadCellException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}


class LCell {
public:
    double xmid;         // position of centre, metres
    double volume;       // m^^3
    double vel;          // bulk velocity, m/s
    double L_bar;        // distance travelled by cell, m
    GasState gas;
    double diam;         // local diameter of tube, m
    double Twall;        // interpolates wall temperature, K
    double K_over_L;     // head-loss loss coefficient per unit length, 1/m
    double shear_stress; // from wall
    double heat_flux;    // from wall
    double dt_chem;      // saved suggested time step, s
    double dt_therm;

    // Conserved quantities.
    double mass;         // in kg
    double moment;       // x-momentum per unit volume, kg/(m^^2.s)
    double energy;       // total energy per unit volume, J/(m^^3)
    // Time derivatives for predictor-corrector time-stepping.
    double[2] dmassdt;
    double[2] dmomdt;
    double[2] dEdt;
    // Record of conserved variables for adaptive time stepping.
    double mass0;
    double moment0;
    double energy0;
    
    this(GasModel gm)
    {
        gas = new GasState(gm);
    }

    @nogc
    void encode_conserved(GasModel gm)
    {
        mass = gas.rho * volume;
        moment = mass * vel;
        energy = gm.internal_energy(gas) + 0.5*vel^^2;
        return;
    }

    @nogc
    void decode_conserved(GasModel gm)
    {
        if (mass <= 0.0 || volume <= 0.0) {
            string msg = "Invalid mass or volume for cell";
            debug {
                msg ~= format(" mass=%g, volume=%g at x=%g", mass, volume, xmid);
            }
            throw new Exception(msg);
        }
        gas.rho = mass/volume;
        vel = moment/mass;
        double ke = 0.5*vel^^2;
        double u_modes_sum = 0.0;
        foreach (u_mode; gas.u_modes) { u_modes_sum += u_mode; }
        gas.u = energy/mass - ke - u_modes_sum;
        gm.update_thermo_from_rhou(gas);
        gm.update_sound_speed(gas);
        gm.update_trans_coeffs(gas);
        return;
    }
} // end class LCell
