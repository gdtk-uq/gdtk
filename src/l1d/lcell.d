// lcell.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
//
// Classes for the faces and cell bodies that make up
// the internal structure of a GasSlug.
//
// PA Jacobs
// 2020-04-09
//
module lcell;

import std.stdio;
import std.math;
import std.format;
import std.algorithm;

import geom;
import gas;
import kinetics;
import gasdyn.gasflow;
import config;


class LFace {
public:
    double x;
    double area;

    double x0;       // initial position for time step
    double[2] dxdt;  // time derivatives, predictor-corrector levels
    double p;        // current pressure from the Riemann problem solution.
    double pLstar, pRstar;  // pressures either side of a closed (hypothetical) valve
    double fopen;    // fraction open for a valve 0.0 <= fopen <= 1.0

    this()
    {
        // Do nothing.
    }

    @nogc
    void record_state()
    {
        x0 = x;
        return;
    }

    @nogc
    void restore_state()
    {
        x = x0;
        return;
    }

    @nogc
    void predictor_step(double dt)
    {
        x = x0 + dxdt[0]*dt;
        return;
    }

    @nogc
    void corrector_step(double dt)
    {
        x = x0 + 0.5*(dxdt[0]+dxdt[1])*dt;
        return;
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
    GasState gas_ref;
    debug { GasState gas_save; } // for reporting, if there is a failure in thermochem update.
    double D;            // local diameter of tube, m
    double L;            // local cell length, m
    double Twall;        // interpolates wall temperature, K
    double vf;           // Viscous factor, nominally 1.0 for full viscous effects
                         // but may be reduced to 0.0 for some sections of tube
                         // in order to get ideal flow behaviour.
    double htcf;         // Heat-transfer factor, nominally 1.0,
                         // may be increased to augment the heat-transfer in
                         // passages that have significantly more surface area
                         // than simple circular pipe.
    double K_over_L;     // head-loss loss coefficient per unit length, 1/m
    double shear_stress; // for reporting via solution file
    double heat_flux;    // for reporting via solution file
    double dt_chem;      // saved suggested time step, s
    double dt_therm;
    //
    // Conserved quantities.
    double mass;         // in kg
    double moment;       // x-momentum per unit volume, kg/(m^^2.s)
    double energy;       // total energy per unit volume, J/(m^^3)
    //
    // Time derivatives for predictor-corrector time-stepping.
    double[2] dmassdt;
    double[2] dmomdt;
    double[2] dEdt;
    double[2] dL_bardt;
    //
    // Record of conserved variables for adaptive time stepping.
    double mass0;
    double moment0;
    double energy0;
    double L_bar0;
    //
    // Source terms accumulators.
    double Q_mass;
    double Q_moment;
    double Q_energy;


    this(GasModel gm)
    {
        gas = GasState(gm);
        gas_ref = GasState(gm);
        debug { gas_save = GasState(gm); }
    }

    @nogc
    void encode_conserved(GasModel gm)
    {
        mass = gas.rho * volume;
        moment = mass * vel;
        energy = mass * (gm.internal_energy(gas) + 0.5*vel^^2);
        return;
    }

    @nogc
    void decode_conserved(GasModel gm)
    {
        if (mass <= 0.0 || volume <= 0.0) {
            string msg = "Invalid mass or volume for cell";
            debug { msg ~= format(" mass=%g, volume=%g at x=%g", mass, volume, xmid); }
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

    @nogc
    void record_state()
    {
        mass0 = mass;
        moment0 = moment;
        energy0 = energy;
        L_bar0 = L_bar;
        return;
    }

    @nogc
    void restore_state(GasModel gm)
    {
        mass = mass0;
        moment = moment0;
        energy = energy0;
        L_bar = L_bar0;
        return;
    }

    @nogc
    void predictor_step(double dt, GasModel gm)
    {
        mass = mass0 + dmassdt[0]*dt;
        moment = moment0 + dmomdt[0]*dt;
        energy = energy0 + dEdt[0]*dt;
        L_bar = L_bar0 + dL_bardt[0]*dt;
        return;
    }

    @nogc
    void corrector_step(double dt, GasModel gm)
    {
        mass = mass0 + 0.5*(dmassdt[0]+dmassdt[1])*dt;
        moment = moment0 + 0.5*(dmomdt[0]+dmomdt[1])*dt;
        energy = energy0 + 0.5*(dEdt[0]+dEdt[1])*dt;
        L_bar = L_bar0 + 0.5*(dL_bardt[0]+dL_bardt[1])*dt;
        return;
    }

    @nogc
    void chemical_increment(double dt, GasModel gmodel, ThermochemicalReactor reactor)
    {
        if (gas.T <= L1dConfig.T_frozen) return;
        double[maxParams] params; // Not using these, but some gas models need them.
        // The gas-dynamic time step may become smaller than dt_chem, dt_therm
        // and that may cause trouble for Rowan's chemistry.
        dt_chem = (dt_chem <= dt) ? dt_chem : dt;
        debug {
            // Keep a copy for reporting, if there is a failure.
            double dt_chem_save = dt_chem;
            gas_save.copy_values_from(gas);
        }
        try {
            reactor(gas, dt, dt_chem, params);
        } catch(ThermochemicalReactorUpdateException err) {
            // It's probably worth one more try but setting dt_chem = -1.0 to give
            // the ODE solver a fresh chance to find a good timestep.
            dt_chem = -1.0;
            try {
                 reactor(gas, dt, dt_chem, params);
            } catch(ThermochemicalReactorUpdateException err) {
                string msg = "chemical_increment() failed.";
                debug {
                    msg ~= format("\n%s", err.msg);
                    msg ~= format("\nThis cell is located at: %s", xmid);
                    msg ~= format("\nThe flow timestep is: %12.6e", dt);
                    msg ~= format("\nThe initial attempted dt_chem is: %12.6e", dt_chem_save);
                    msg ~= format("\nBEFORE the failed update gas state is:\n  %s", gas_save);
                    msg ~= format("\nAFTER the failed update gas state is:\n  %s", gas);
                }
                throw new ThermochemicalReactorUpdateException(msg);
            }
        }
        // The update only changes mass fractions; we need to impose
        // a thermodynamic constraint based on a call to the equation of state.
        try {
            gmodel.update_thermo_from_rhou(gas);
        }
        catch (Exception err) {
            string msg = "chemical_increment() failed update_thermo_from_rhou";
            debug {
                msg ~= format("\n%s", err.msg);
                msg ~= format("\nThis cell is located at: %g", xmid);
                msg ~= "\nThis failure occurred when trying to update the thermo state after";
                msg ~= "\ncomputing the species change due to chemical reactions.";
                msg ~= format("\nGas state:\n  %s", gas);
            }
            throw new GasModelException(msg);
        }
        return;
    } // end chemical_increment()

    @nogc
    void source_terms(int viscous_effects, bool adiabatic, GasModel gmodel)
    {
        // Compute the components of the source vector, Q_*.
        // These quantities are used to include generic flow losses
        // and wall heat-transfer.
        // The influence of these processes will be added directly to
        // the time-derivatives of the bulk-flow conserved quantities.
        //
        // Start with a clean slate.
        Q_mass = 0.0;
        Q_moment = 0.0;
        Q_energy = 0.0;
        shear_stress = 0.0;
        heat_flux = 0.0;
        if (viscous_effects == 0) { return; } // There is nothing to do.
        //
        double abs_vel = fabs(vel);
        double Prandtl = 0.75; // A constant value, for the moment.
        double omega = pow(Prandtl, 0.333); // Recovery factor, assume turbulent.
        double M = abs_vel/gas.a; // Local Mach number.
        double lambda = 1.0 + (gmodel.gamma(gas)-1.0)*0.5*omega*M*M;
        double T_aw = lambda * gas.T; // Adiabatic wall temperature.
        double T_wall_seen = Twall;
	if (adiabatic) {
	    T_wall_seen = T_aw;
	} else {
	    // John Hunter's suggestion for heat transfer, from 1991/2.
            // Heat transfer will tend to increase the wall temperature
	    // to something close to the core temperature of the gas.
            // This was activated with adiabatic==2 in the old code.
            // T_wall_seen = max(gas.T-400.0, Twall);
	}
	// Transport properties based on Eckert reference conditions.
	double T_ref = gas.T + 0.5*(T_wall_seen-gas.T) + 0.22*(T_aw-gas.T);
        gas_ref.copy_values_from(gas);
        gas_ref.T = T_ref;
        gas_ref.rho = gas.rho * gas.T/T_ref;
        gas_ref.p = gas_ref.rho * gmodel.R(gas_ref) * gas_ref.T;
        gmodel.update_thermo_from_pT(gas_ref);
        gmodel.update_trans_coeffs(gas_ref);
        double f = 0.0;
        if (viscous_effects == 2 && L_bar > 0.0) {
            // Friction and heat flux based on a flat plate calculation
	    // No compressibility correction apart from
	    // property evaluation at Eckert reference conditions
	    // and adiabatic wall temperature based on recovery factor.
            // Use distance moved as reference distance for another value of Re.
            double Re_L = gas_ref.rho*L_bar*abs_vel/gas_ref.mu;
	    f = f_flat_plate(Re_L);
	} else {
            // Presumably, viscous_effects == 1 at this point.
	    // Default: friction factor determined from
	    // fully-developed pipe flow correlation.
            // Local Reynolds number based on diameter and reference conditions.
            double Re_D = gas_ref.rho*D*abs_vel/gas_ref.mu;
	    f = f_darcy_weisbach(Re_D)/lambda;
	}
	// Local shear stress, in Pa, computed from the friction factor and
        // modulated by the viscous factor at this location.
        // This is the stress that the wall applies to the gas and it is signed.
	shear_stress = -(0.125*f) * gas.rho * vel*abs_vel * vf;
	// Rate of energy transfer into the cell.
        // Shear stress does no work on the cell because the wall velocity is zero.
        double w_dot = 0.0;
	// Convective heat transfer coefficient.
	double St = (f*0.125) * pow(Prandtl, -0.667);
	double h = gas_ref.rho * gmodel.Cp(gas) * abs_vel * St;
	// Convective heat transfer from the wall into the gas.
        // Note the modulation by the viscous factor and Matt's htc augmentation.
	if (adiabatic) {
	    heat_flux = 0.0;
	} else {
	    heat_flux = h * (T_wall_seen - T_aw) * vf * htcf; // units W/m^^2
	}
        // Shear stress and heat transfer act over the wall-constact surface.
        Q_moment += shear_stress * PI*D*L;
        Q_energy += heat_flux * PI*D*L;
        //
        // Pipe fitting loss is a momentum loss on top of other viscous effects.
        // It is applied to the whole volume of the cell.
	double F_loss = K_over_L * 0.5*gas.rho*vel*abs_vel;
	Q_moment -= volume * F_loss;
        return;
    } // end source_terms()
} // end class LCell

//----------------------------------------------------------------------

// Friction factors for a couple of different models.

@nogc
double f_darcy_weisbach(double Re)
{
    // Darcy-Weisbach friction factor for fully-developed flow in a smooth pipe.
    // Re is based on pipe diameter, D.
    double f;
    if (Re < 10.0) {
	// A reasonable limit for very low speeds.
	f = 6.4;
    } else if (Re < 2000.0) {
	// Laminar regime.
	f = 64.0 / Re;
    } else if (Re < 4000.0) {
	// Transition regime.
	f = 0.032 / pow(Re/2000.0, 0.3187);
    } else {
	// Fully turbulent, smooth wall
	double lgtmp = 1.327359 - 0.9*log10(Re);
	f = 1.14 - 2.0 * lgtmp;
	f = 1.0 / (f * f);
    }
    return f;
} // end f_darcy_weisbach()

@nogc
double f_flat_plate(double Re)
{
    // Flat-plate friction factor with Re based on length along plate.
    // See Holman (1986) Heat Transfer, Eq 5.125 to 5.127
    double f;
    if (Re < 1.0e4) {
	// Somewhat arbitrary lower Re limit to stop
	// possibility of near-infinite friction and heat flux.
	f = 8 * 0.332/sqrt(10000.0);
    } else if (Re < 5.0e5) {
        // Laminar regime
	f = 8 * 0.332/sqrt(Re);
    } else if (Re < 1.0e7) {
        // lower Re turbulent regime.
	f = 8 * 0.0296*pow(Re, -0.2);
    } else {
        // high Re tubulent regime.
	f = 8 * 0.185*pow(log10(Re), -2.584);
    }
    return f;
} // end f_flat_plate()
