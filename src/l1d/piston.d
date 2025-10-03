// piston.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
// 2025-03-19 : allow variable piston mass.
//
module l1d.piston;

import std.conv;
import std.stdio;
import std.file;
import std.string;
import std.json;
import std.format;
import std.algorithm;
import std.math;

import l1d.config;
import l1d.endcondition;
import l1d.lcell;
import l1d.misc;

import gas;
import gasdyn.gasflow;
import geom;
import util.json_helper;


class Piston {
public:
    size_t indx;
    string label = "";
    double mass;  // mass, kg
    double diam;  // diameter, m
    double area;  // m^2
    double L;     // length, m
    double x;     // position, m
    double vel;   // velocity, m/s
    double massf; // current fraction of original mass that remains
    double dmassfdt; // rate of change of mass fraction (negative for mass decay)
    double massfmin; // minimum mass fraction, below which interaction assumes zero mass
    double front_seal_f; // friction factor
    double front_seal_area; // area over which pressure acts
    double back_seal_f;
    double back_seal_area;
    double p_restrain;
    bool is_restrain;
    bool with_brakes;
    bool brakes_on;
    bool full_stop;
    double brakes_friction_force; // in N
    double x_buffer;
    bool on_buffer;
    int ecL_id;
    int ecR_id;
    EndCondition ecL;
    EndCondition ecR;

    bool is_restrain0;
    bool brakes_on0;
    bool on_buffer0;
    double x0;
    double vel0;
    double massf0;
    double[2] dxdt;
    double[2] dvdt;

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct piston[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        label = getJSONstring(jsonData, "label", "");
        mass = getJSONdouble(jsonData, "mass", 0.0);
        diam = getJSONdouble(jsonData, "diameter", 0.0);
        area = 0.25*PI*diam*diam;
        L = getJSONdouble(jsonData, "length", 0.0);
        dmassfdt = getJSONdouble(jsonData, "dmassfdt", 0.0);
        massfmin = getJSONdouble(jsonData, "massfmin", 1.0e-3);
        front_seal_f = getJSONdouble(jsonData, "front_seal_f", 0.0);
        front_seal_area = getJSONdouble(jsonData, "front_seal_area", 0.0);
        back_seal_f = getJSONdouble(jsonData, "back_seal_f", 0.0);
        back_seal_area = getJSONdouble(jsonData, "back_seal_area", 0.0);
        p_restrain = getJSONdouble(jsonData, "p_restrain", 0.0);
        x_buffer = getJSONdouble(jsonData, "x_buffer", 1.0e6);
        with_brakes = getJSONbool(jsonData, "with_brakes", false);
        brakes_friction_force = getJSONdouble(jsonData, "brakes_friction_force", 0.0);
        ecL_id = getJSONint(jsonData, "ecL_id", -1);
        ecR_id = getJSONint(jsonData, "ecR_id", -1);
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Piston[", indx, "]:");
            writeln("  mass= ", mass);
            writeln("  diam= ", diam);
            writeln("  area= ", area);
            writeln("  L= ", L);
            writeln("  dmassfdt= ", dmassfdt);
            writeln("  massfmin= ", massfmin);
            writeln("  front_seal_f= ", front_seal_f);
            writeln("  front_seal_area= ", front_seal_area);
            writeln("  back_seal_f= ", back_seal_f);
            writeln("  back_seal_area= ", back_seal_area);
            writeln("  p_restrain= ", p_restrain);
            writeln("  x_buffer= ", x_buffer);
            writeln("  with_brakes= ", with_brakes);
            writeln("  brakes_friction_force= ", brakes_friction_force);
            writeln("  ecL_id= ", ecL_id);
            writeln("  ecR_id= ", ecR_id);
        }
    } // end Piston constructor

    void read_data(File fp, int tindx)
    {
        string text = fp.readln().chomp();
        while (text.canFind("#")) { text = fp.readln().chomp(); }
        string[] items = text.split();
        int myTindx = to!int(items[0]);
        while (myTindx < tindx) {
            text = fp.readln().chomp();
            items = text.split();
            myTindx = to!int(items[0]);
        }
        // We should be at the line that contains the requested tindx.
        x = to!double(items[1]);
        vel = to!double(items[2]);
        is_restrain = (to!int(items[3]) == 1);
        brakes_on = (to!int(items[4]) == 1);
        on_buffer = (to!int(items[5]) == 1);
        // We introduce the current mass fraction variable 2025-03-19
        // so it may not be always available in older simulation files.
        if (items.length > 6) massf = to!double(items[6]);
    } // end read_data()

    void write_data(File fp, int tindx, bool write_header)
    {
        if (write_header) {
            fp.writeln("# tindx  x  vel  is_restrain  brakes_on  on_buffer  massf");
        }
        fp.writeln(format("%d %e %e %d %d %d %e", tindx, x, vel,
                          ((is_restrain)?1:0), ((brakes_on)?1:0),
                          ((on_buffer)?1:0), massf));
    } // end write_data()

    @nogc @property
    double energy()
    {
        return mass*massf*0.5*vel*vel;
    }

    @nogc
    void record_state()
    {
        is_restrain0 = is_restrain;
        brakes_on0 = brakes_on;
        on_buffer0 = on_buffer;
        x0 = x;
        vel0 = vel;
        massf0 = massf;
        return;
    }

    @nogc
    void restore_state()
    {
        is_restrain = is_restrain0;
        brakes_on = brakes_on0;
        on_buffer = on_buffer0;
        x = x0;
        vel = vel0;
        massf = massf0;
        return;
    }

    void check_for_buffer_interaction(double t)
    {
        if (on_buffer) {
            if (vel < -0.001) {
                on_buffer = false;
            }
        } else { // not (yet) identified as being on the buffer
            if ((x > x_buffer) && (vel > 0.0)) {
                on_buffer = true;
                if (vel > 0.001) {
                    // Only report an impact with significant speed.
                    string msg = format("t=%e Piston[%d] buffer strike with vel=%e\n", t, indx, vel);
                    write(msg);
                    append(L1dConfig.job_name~"/events.txt", msg);
                }
            }
        }
        if (on_buffer) { vel = 0.0; }
        return;
    }

    @nogc
    void time_derivatives(int level)
    {
        // Computation of the time derivatives for the normal piston dynamics.
        // We compute the rate of change for velocity and position of the piston
        // assuming that the usual surface forces are at play.
        //
        // Pressure on each piston face.
        double pL = 0.0;
        if (ecL && ecL.slugL) {
            pL = (ecL.slugL_end == End.L) ?
                ecL.slugL.faces[0].p : ecL.slugL.faces[$-1].p;
        }
        double pR = 0.0;
        if (ecR && ecR.slugR) {
            pR = (ecR.slugR_end == End.L) ?
                ecR.slugR.faces[0].p : ecR.slugR.faces[$-1].p;
        }
        // Pressures drive the piston dynamics.
        if (is_restrain) {
            if (pL > p_restrain) { is_restrain = false; }
        }
        if (is_restrain) {
            dxdt[level] = 0.0;
            dvdt[level] = 0.0;
        }
        // The (signed) pressure force.
        double pressure_force = area*(pL-pR);
        // The magnitude of the friction force from pressurized seals.
        double friction_force = front_seal_f*front_seal_area*pR +
            back_seal_f*back_seal_area*pL;
        // Brakes modelled on those used on T4 shock tunnel.
        if (with_brakes && vel < 0.0) {
            brakes_on = true;
            friction_force += brakes_friction_force;
        } else {
            brakes_on = false;
        }
        full_stop = false;
        if (massf > massfmin) {
            // We have the usual dynamic processes with a non-zero mass piston.
            // Rate of change of velocity is acceleration.
            // Normally, we would use these derivative values in the
            // predictor and corrector updates.
            double effectivemass = mass * massf;
            // If the brakes have stopped the piston, we have an unusual update
            // and don't want to use those derivatives.
            immutable vel_tol = 1.0e-6;
            if (vel > vel_tol) {
                // Moving forward, apply full friction in reverse.
                dvdt[level] = (pressure_force-friction_force)/effectivemass;
            } else if (vel < -vel_tol) {
                // Moving backward, apply full friction forward.
                dvdt[level] = (pressure_force+friction_force)/effectivemass;
            } else {
                // We are effectively stationary.
                if (fabs(pressure_force) > friction_force) {
                    // Pressure force overcomes friction.
                    dvdt[level] = (pressure_force > 0.0) ?
                        (pressure_force-friction_force)/effectivemass :
                        (pressure_force+friction_force)/effectivemass;
                } else {
                    // Friction force dominates; let's remain stationary for now.
                    vel = 0.0;
                    dvdt[level] = 0.0;
                    full_stop = true;
                } // end if sufficient pressure to accelerate
            } // end if vel...
        } else {
            // We pretend that the piston has zero mass
            // and allow the two gas slugs to push directly
            // against each other, with the pressure and velocity
            // determined via a solution of the Riemann problem.
            // We can only do this if we actually have both slugs of gas.
            if (ecL && ecL.slugL && ecR && ecR.slugR) {
                auto slugL = ecL.slugL;
                LCell cL = (ecL.slugL_end == End.L) ? slugL.cells[0] : slugL.cells[$-1];
                auto slugR = ecR.slugR;
                LCell cR = (ecR.slugR_end == End.L) ? slugR.cells[0] : slugR.cells[$-1];
                double p;
                lrivp(cL.gas, cR.gas, cL.vel, cR.vel, slugL.gmodel, slugR.gmodel, vel, p);
                dvdt[level] = 0.0;
            } else {
                assert(0, "FATAL: Piston has effectively zero mass and we do NOT have both slugs!");
            }
        }
        //
        // Rate of change of position is velocity.
        dxdt[level] = vel;
        return;
    } // end time_derivatives()

    @nogc
    void predictor_step(double dt)
    {
        if (is_restrain) {
            x = x0;
            vel = vel0;
        } else {
            x = x0 + dxdt[0]*dt;
            vel = (full_stop) ? 0.0 : vel0 + dvdt[0]*dt;
        }
        return;
    }

    @nogc
    void corrector_step(double dt)
    {
        if (is_restrain) {
            x = x0;
            vel = vel0;
        } else {
            x = x0 + 0.5*(dxdt[0]+dxdt[1])*dt;
            vel = (full_stop) ? 0.0 : vel0 + 0.5*(dvdt[0]+dvdt[1])*dt;
        }
        return;
    }

    @nogc
    void change_mass(double dt)
    {
        // For modelling the secondary diaphragms in an expansion tube,
        // where we assume that the diaphram is initially punched out
        // of its restraint and temporarily acts as a piston.
        if (!is_restrain) {
            massf += dmassfdt * dt;
            massf = fmax(massf, 0.0);
        } else {
            // Do nothing while the piston is still restrained.
        }
        return;
    }
} // end class Piston
