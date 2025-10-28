/**
 * flowstate.d
 * FlowState class for use in the main solver.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module lmr.flowstate;

import std.algorithm;
import std.array;
import std.conv;
import std.format;
import std.json;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.zip;

import gas;
import geom;
import gzip;
import nm.number;
import ntypes.complex;
import util.json_helper;

import lmr.lmrexceptions : LmrException;
import lmr.globalconfig;

@nogc
void into_rotating_frame(ref Vector3 v, ref const(Vector3) pos, double omegaz)
// Velocity vector becomes relative to the rotating frame of the block
// by subtracting the entrainment velocity (-y*omegaz i + x*omegaz j).
{
    v.x += pos.y * omegaz;
    v.y -= pos.x * omegaz;
}

@nogc
void into_nonrotating_frame(ref Vector3 v, ref const(Vector3) pos, double omegaz)
// Velocity vector becomes relative to a nonrotating frame
// by adding the entrainment velocity (-y*omegaz i + x*omegaz j).
{
    v.x -= pos.y * omegaz;
    v.y += pos.x * omegaz;
}


struct FlowState {
public:
    GasState gas;  // gas state
    Vector3 vel;   // flow velocity, m/s
    version(MHD) {
        Vector3 B;     // magnetic field strength
        number psi;    // divergence cleaning parameter
        number divB;   // divergence of the magnetic field
    }
    version(turbulence) {
        number[2] turb; // turbulence primitives
    }
    double[2] electric_field;
    number mu_t;   // turbulence viscosity
    number k_t;    // turbulence thermal-conductivity
    number S;         // shock indicator

    @disable this();

    this(GasModel gm,
         in double p_init,
         in double T_init,
         in double[] T_modes_init,
         in Vector3 vel_init,
         in double[2] turb_init,
         in double[] massf_init=[1.0,],
         in double quality_init=1.0,
         in Vector3 B_init=Vector3(0.0,0.0,0.0),
         in double psi_init=0.0, in double divB_init=1.0,
         in double mu_t_init=0.0, in double k_t_init=0.0,
         in double S_init=0.0)
    {
        gas = GasState(gm, p_init, T_init, T_modes_init, massf_init, quality_init);
        vel = vel_init;
        version(MHD) {
            B = B_init;
            psi = psi_init;
            divB = divB_init;
        }
        version(turbulence) {
            foreach (i; 0 .. turb.length) turb[i] = turb_init[i];
        }
        mu_t = mu_t_init;
        k_t = k_t_init;
        S = S_init;
    }

    this(in FlowState other, GasModel gm)
    {
        gas = GasState(gm);
        gas.copy_values_from(other.gas);
        vel = other.vel;
        version(MHD) {
            B = other.B;
            psi = other.psi;
            divB = other.divB;
        }
        version(turbulence) {
            foreach (i; 0 .. turb.length) turb[i] = other.turb[i];
        }
        mu_t = other.mu_t;
        k_t = other.k_t;
        S = other.S;
    }

    this(in FlowState other)
    {
        gas = GasState(other.gas);
        vel.set(other.vel);
        version(MHD) {
            B.set(other.B);
            psi = other.psi;
            divB = other.divB;
        }
        version(turbulence) {
            turb = other.turb.dup;
        }
        mu_t = other.mu_t;
        k_t = other.k_t;
        S = other.S;
    }

    this(GasModel gm, size_t nturb)
    {
        gas = GasState(gm, 100.0e3, 300.0, [1.0,], 1.0);
        vel.set(0.0,0.0,0.0);
        version(MHD) {
            B.set(0.0,0.0,0.0);
            psi = 0.0;
            divB = 0.0;
        }
        version(turbulence) {
            foreach (i; 0 .. turb.length) turb[i] = 0.0;
        }
        mu_t = 0.0;
        k_t = 0.0;
        S = 0.0;
    }

    this(in JSONValue json_data, GasModel gm)
    {
        double p = getJSONdouble(json_data, "p", 100.0e3);
        double T = getJSONdouble(json_data, "T", 300.0e3);
        double[] T_modes;
        version(multi_T_gas) {
            foreach(i; 0 .. gm.n_modes) { T_modes ~= T; }
            T_modes = getJSONdoublearray(json_data, "T_modes", []);
        }
        double[] massf;
        version(multi_species_gas) {
            massf = getJSONdoublearray(json_data, "massf", [1.0,]);
        }
        double quality = getJSONdouble(json_data, "quality", 1.0);
        gas = GasState(gm, p, T, T_modes, massf, quality);
        double velx = getJSONdouble(json_data, "velx", 0.0);
        double vely = getJSONdouble(json_data, "vely", 0.0);
        double velz = getJSONdouble(json_data, "velz", 0.0);
        vel.set(velx,vely,velz);
        version(MHD) {
            double Bx = getJSONdouble(json_data, "Bx", 0.0);
            double By = getJSONdouble(json_data, "By", 0.0);
            double Bz = getJSONdouble(json_data, "Bz", 0.0);
            B.set(Bx,By,Bz);
            psi = getJSONdouble(json_data, "psi", 0.0);
            divB = getJSONdouble(json_data, "divB", 0.0);
        }
        version(turbulence) {
            double[2] turb_in;
            turb_in = getJSONdoublearray(json_data, "turb", [0.0, 1.0]);
            foreach (i; 0 .. turb.length) turb[i] = turb_in[i];
        }
        mu_t = getJSONdouble(json_data, "mu_t", 0.0);
        k_t = getJSONdouble(json_data, "k_t", 0.0);
        S = getJSONdouble(json_data, "S", 0.0);
    }

    FlowState dup() const
    {
        return FlowState(this);
    }

    @nogc void opAssign(in FlowState other)
    {
        gas = other.gas;
        vel.set(other.vel);
        version(MHD) {
            B.set(other.B);
            psi = other.psi;
            divB = other.divB;
        }
        version(turbulence) {
            debug {
                if (turb.length != other.turb.length) { throw new Error("Incorrect turb length."); }
            }
            foreach (i; 0 .. turb.length) { turb[i] = other.turb[i]; }
        }
        mu_t = other.mu_t;
        k_t = other.k_t;
        S = other.S;
    }

    @nogc
    void copy_values_from(in FlowState other)
    {
        this = other;
    }

    @nogc
    void copy_values_from(in FlowState* other)
    {
        this = *other;
    }

    @nogc
    void copy_average_values_from(in FlowState fs0, in FlowState fs1, double w0=0.5)
    // Avoids memory allocation, it's all in place.
    {
        double w1 = 1.0 - w0;
        gas.copy_average_values_from(fs0.gas, fs1.gas, w0);
        vel.set(w0*fs0.vel.x + w1*fs1.vel.x, w0*fs0.vel.y + w1*fs1.vel.y, w0*fs0.vel.z + w1*fs1.vel.z);
        version(MHD) {
            B.set(w0*fs0.B.x + w1*fs1.B.x, w0*fs0.B.y + w1*fs1.B.y, w0*fs0.B.z + w1*fs1.B.z);
            psi = w0*fs0.psi + w1*fs1.psi;
            divB = w0*fs0.divB + w1*fs1.divB;
        }
        version(turbulence) {
            foreach (i; 0 .. turb.length) { turb[i] =  w0*fs0.turb[i] + w1*fs1.turb[i]; }
        }
        mu_t = w0*fs0.mu_t + w1*fs1.mu_t;
        k_t = w0*fs0.k_t + w1*fs1.k_t;
    } // end copy_average_values_from()

    void copy_average_values_from(in FlowState*[] others, GasModel gm)
    // Note that we must not send the current object in the others list as well.
    // Involves some memory allocation.
    {
        size_t n = others.length;
        if (n == 0) throw new FlowSolverException("Need to average from a nonempty array.");
        GasState*[] gasList;
        // Note that, because we cast away their "const"ness,
        // we need to be honest and not to fiddle with the other gas states.
        foreach(other; others) {
            if ( &this is other ) {
                throw new FlowSolverException("Must not include destination in source list.");
            }
            gasList ~= cast(GasState*)&(other.gas);
        }
        gas.copy_average_values_from(gasList, gm);
        // Accumulate from a clean slate and then divide.
        vel.clear();
        version(MHD) {
            B.clear();
            psi = 0.0;
            divB = 0.0;
        }
        version(turbulence) {
            foreach(i; 0 .. turb.length) turb[i] = 0.0;
        }
        mu_t = 0.0;
        k_t = 0.0;
        S = 0.0; // Shock detector is no longer necessarily an integer flag.
        foreach(other; others) {
            vel.x += other.vel.x; vel.y += other.vel.y; vel.z += other.vel.z;
            version(MHD) {
                B.x += other.B.x; B.y += other.B.y; B.z += other.B.z;
                psi += other.psi;
                divB += other.divB;
            }
            version(turbulence) {
                foreach (i; 0 .. turb.length) turb[i] += other.turb[i];
            }
            mu_t += other.mu_t;
            k_t += other.k_t;
            S += other.S;
        }
        number scale = 1.0/to!number(n);
        vel *= scale;
        version(MHD) {
            B *= scale;
            psi *= scale;
            divB *= scale;
        }
        version(turbulence) {
            foreach (i; 0 .. turb.length) turb[i] *= scale;
        }
        mu_t *= scale;
        k_t *= scale;
        S = S;
    } // end copy_average_values_from()

    string toString() const
    {
        char[] repr;
        repr ~= "FlowState(";
        repr ~= "gas=" ~ to!string(gas);
        repr ~= ", vel=" ~ to!string(vel);
        version(MHD) {
            repr ~= ", B=" ~ to!string(B);
            repr ~= ", psi=" ~ to!string(psi);
            repr ~= ", divB=" ~ to!string(psi);
        }
        version(turbulence) {
            repr ~= ", turb=" ~ to!string(turb);
        }
        repr ~= ", mu_t=" ~ to!string(mu_t);
        repr ~= ", k_t=" ~ to!string(k_t);
        repr ~= ", S=" ~ to!string(S);
        repr ~= ", E=" ~ to!string(electric_field);
        repr ~= ")";
        return to!string(repr);
    }

    string toJSONString() const
    {
        auto writer = appender!string();
        formattedWrite(writer, "{");
        formattedWrite(writer, "\"p\": %.18e", gas.p.re);
        formattedWrite(writer, ", \"T\": %.18e", gas.T.re);
        version(multi_T_gas) {
            // zero or more T_modes
            formattedWrite(writer, ", \"T_modes\": [");
            if (gas.T_modes.length > 0) { formattedWrite(writer, " %.18e", gas.T_modes[0].re); }
            foreach (i; 1 .. gas.T_modes.length) { formattedWrite(writer, ", %.18e", gas.T_modes[i].re); }
            formattedWrite(writer, "]");
        }
        version(multi_species_gas) {
            // one or more mass fractions
            formattedWrite(writer, ", \"massf\": [ %.18e", gas.massf[0].re);
            foreach (i; 1 .. gas.massf.length) {
                formattedWrite(writer, ", %.18e", gas.massf[i].re);
            }
            formattedWrite(writer, "]");
        }
        formattedWrite(writer, ", \"quality\": %.18e", gas.quality.re);
        formattedWrite(writer, ", \"velx\": %.18e", vel.x.re);
        formattedWrite(writer, ", \"vely\": %.18e", vel.y.re);
        formattedWrite(writer, ", \"velz\": %.18e", vel.z.re);
        version(MHD) {
            formattedWrite(writer, ", \"Bx\": %.18e", B.x.re);
            formattedWrite(writer, ", \"By\": %.18e", B.y.re);
            formattedWrite(writer, ", \"Bz\": %.18e", B.z.re);
            formattedWrite(writer, ", \"psi\": %.18e", psi.re);
            formattedWrite(writer, ", \"divB\": %.18e", divB.re);
        }
        version(turbulence) {
            formattedWrite(writer, ", \"turb\": [");
            if (turb.length > 0) { formattedWrite(writer, " %.18e", turb[0].re); }
            foreach (i; 1 .. turb.length) { formattedWrite(writer, ", %.18e", turb[i].re); }
            formattedWrite(writer, "]");
        }
        formattedWrite(writer, ", \"mu_t\": %.18e", mu_t.re);
        formattedWrite(writer, ", \"k_t\": %.18e", k_t.re);
        formattedWrite(writer, ", \"S\": %.18e", S.re);
        formattedWrite(writer, "}");
        return writer.data;
    } // end toJSONString()

    @nogc
    bool check_data(ref Vector3 pos, ref const(LocalConfig) lc) const
    {
        auto flowstate_limits = lc.flowstate_limits;
        bool is_data_valid = gas.check_values(true);
        if (fabs(vel.x) > flowstate_limits.max_velocity ||
            fabs(vel.y) > flowstate_limits.max_velocity ||
            fabs(vel.z) > flowstate_limits.max_velocity) {
            debug { writeln("Velocity too high ", vel); }
            is_data_valid = false;
        }
        if (gas.T < flowstate_limits.min_temp) {
            debug { writeln("Temperature below minimum ", gas.T); }
            is_data_valid = false;
        }
        if (gas.T > flowstate_limits.max_temp) {
            debug { writeln("Temperature above maximum ", gas.T); }
            is_data_valid = false;
        }
        version(turbulence) {
            if (!lc.turb_model.is_valid(flowstate_limits, turb)) {
            is_data_valid = false;
            }
        }
        if (!is_data_valid) {
            debug { writeln("   at position ", pos); }
        }
        return is_data_valid;
    } // end check_data()

    @nogc
    void reorient_vector_quantities(const(double[]) Rmatrix)
    {
        vel.apply_matrix_transform(Rmatrix);
        version(MHD) {
            B.apply_matrix_transform(Rmatrix);
        }
    }

version(complex_numbers) {
    @nogc
    void clear_imaginary_components()
    // When performing the complex-step Frechet derivative in the Newton-Krylov accelerator,
    // the flowstate values accumulate imaginary components, so we have to start with a clean slate, so to speak.
    {
        gas.clear_imaginary_components();
        vel.x.im = 0.0;
        vel.y.im = 0.0;
        vel.z.im = 0.0;
        version(MHD) {
            B.x.im = 0.0;
            B.y.im = 0.0;
            B.z.im = 0.0;
            psi.im = 0.0;
            divB.im = 0.0;
        }
        version(turbulence) {
            foreach (i; 0..turb.length) turb[i].im = 0.0;
        }
        mu_t.im = 0.0;
        k_t.im = 0.0;
    } // end clear_imaginary_components()
} // end version(complex)

} // end class FlowState


class StaticFlowProfile {
    // For use in the classes that implement the InflowBC_StaticProfile boundary condition.
    // GhostCellFlowStateCopyFromStaticProfile, BIE_FlowStateCopyFromStaticProfile
    // There are non-obvious options for the match parameter in the constructor call.
    // See the switch statement in the compute_distance() function for some hints.

public:
    string fileName;
    string posMatch;
    FlowState[] fstate;
    Vector3[] pos;
    size_t[size_t] which_point; // A place to memoize the mapped indices and we find them.
    // Below, we search for the profile point nearest to the initial position.
    // This position will only change for moving-grid simulations and we will not try
    // to deal with that complication.

    // Construct a StaticFlowProfile from a data file.
    //
    // Format will be a header line followed by one sample point per line.
    // The header line will specify the names of the columns and, at a minimum, be:
    // pos.x pos.y p T vel.x vel.y
    // The data lines will have numerical values for the corresponding quantities.
    // This format should be compatible with that provided by the sliced output
    // of the postprocessing programs and also be compatible with gnuplot.
    this (string fileName, string match)
    {
        this.fileName = fileName;
        this.posMatch = match;
        //
        // Some context for the expected flow data.
        alias cfg = GlobalConfig;
        GasModel gmodel = cfg.gmodel_master;
        size_t n_species = gmodel.n_species;
        size_t n_modes = gmodel.n_modes;
        size_t n_turb = cfg.turb_model.nturb;
        //
        string[] expected_names = ["pos.x", "pos.y", "p", "T", "vel.x", "vel.y"];
        if (cfg.dimensions == 3) { expected_names ~= ["pos.z", "vel.z"]; }
        if (cfg.MHD) { expected_names ~= ["B.x", "B.y", "B.z"]; }
        string[] speciesList;
        if (n_species > 1) {
            foreach (i; 0..n_species) { speciesList ~= "massf-" ~ gmodel.species_name(i); }
            expected_names ~= speciesList;
        }
        string[] TmodeList;
        foreach (i; 0..n_modes) { TmodeList ~= "T-" ~ gmodel.energy_mode_name(to!int(i)); }
        expected_names ~= TmodeList;
        string[] turbList;
        foreach (i; 0..n_turb) { turbList ~= "tq-" ~ cfg.turb_model.primitive_variable_name(i); }
        expected_names ~= turbList;
        //
        // Open filename and read all data points.
        auto f = new File(fileName);
        auto range = f.byLine();
        //
        // First line is expected to name the columns with the correct variable names.
        // Assume that it does, and proceed to identify the columns.
        auto line = range.front;
        string txt = to!string(line);
        txt = stripLeft(txt, "# ").stripRight();
        string[] varnames = txt.split();
        size_t[string] column_dict; foreach (i, varname; varnames) { column_dict[varname] = i; }
        range.popFront();
        foreach (name; expected_names) {
            if (!canFind(varnames, name)) {
                string msg = format("Could not find required variable \"%s\" in file: %s", name, fileName);
                msg ~= format("\nFile header consisted of the following variables:\n %s", varnames);
                throw new LmrException(msg);
            }
        }
        // Start picking up the data lines.
        line = range.front;
        while (!line.empty) {
            txt = to!string(line);
            auto mypos = Vector3();
            auto myfs = FlowState(gmodel, n_turb);
            txt = stripLeft(txt).stripRight();
            string[] items = txt.split();
            mypos.x = to!number(items[column_dict["pos.x"]]);
            mypos.y = to!number(items[column_dict["pos.y"]]);
            mypos.z = (cfg.dimensions == 3) ? to!number(items[column_dict["pos.z"]]) : to!number(0.0);
            myfs.gas.p = to!number(items[column_dict["p"]]);
            myfs.gas.T = to!number(items[column_dict["T"]]);
            foreach (i, name; TmodeList) {
                myfs.gas.T_modes[i] = to!number(items[column_dict[name]]);
            }
            if (n_species > 1) {
                foreach (i, name; speciesList) { myfs.gas.massf[i] = to!number(items[column_dict[name]]); }
            } else {
                myfs.gas.massf[0] = 1.0;
            }
            gmodel.update_thermo_from_pT(myfs.gas);
            gmodel.update_sound_speed(myfs.gas);
            gmodel.update_trans_coeffs(myfs.gas);
            foreach (i; 0 .. n_species) { myfs.gas.rho_s[i] = myfs.gas.massf[i] * myfs.gas.rho; }
            myfs.vel.x = to!number(items[column_dict["vel.x"]]);
            myfs.vel.y = to!number(items[column_dict["vel.y"]]);
            myfs.vel.z = (cfg.dimensions == 3) ? to!number(items[column_dict["vel.z"]]) : to!number(0.0);
            foreach (i, name; turbList) { myfs.turb[i] = to!number(items[column_dict[name]]); }
            if (cfg.MHD) {
                myfs.B.x = to!number(items[column_dict["B.x"]]);
                myfs.B.y = to!number(items[column_dict["B.y"]]);
                myfs.B.z = to!number(items[column_dict["B.z"]]);
            } else {
                myfs.B.set(0.0, 0.0, 0.0);
            }
            fstate ~= myfs;
            pos ~= mypos;
            range.popFront();
            line = range.front;
        } // end while
        //
        // The mapping of the nearest profile point to each ghost-cell or interface location
        // will be done as needed, at application time.
        // This way, all of the necessary cell and position data should be valid.
    } // end this()

    @nogc
    double compute_distance(ref const(Vector3) my_pos, ref const(Vector3) other_pos)
    {
        double distance, other_r, my_r, dx, dy, dz, dr;
        switch (posMatch) {
        case "xyz-to-xyz":
            // 2D or 3D, closest match on all components of position.
            // In 2D all z-components are supposed to be zero (and so, not matter).
            dx = my_pos.x.re - other_pos.x.re;
            dy = my_pos.y.re - other_pos.y.re;
            dz = my_pos.z.re - other_pos.z.re;
            distance = sqrt(dx*dx + dy*dy + dz*dz);
            break;
        case "xyA-to-xyA":
            // 2D or 3D; don't care about z-component of position.
            dx = my_pos.x.re - other_pos.x.re;
            dy = my_pos.y.re - other_pos.y.re;
            distance = sqrt(dx^^2 + dy^^2);
            break;
        case "AyA-to-AyA":
            // 2D or 3D; only care about the y-component of position.
            dy = my_pos.y.re - other_pos.y.re;
            distance = fabs(dy);
            break;
        case "xy-to-xR":
            // Starting with a profile from a 2D simulation, map it to
            // a radial profile in a 3D simulation, considering the x-component
            // of the position of the ghost cells when computing distance and
            // picking the nearest point in the profile.
            dx = my_pos.x.re - other_pos.x.re;
            other_r = sqrt(other_pos.y.re^^2 + other_pos.z.re^^2);
            my_r = sqrt(my_pos.y.re^^2 + my_pos.z.re^^2);
            dr = my_r - other_r;
            distance = sqrt(dx*dx + dr*dr);
            break;
        case "Ay-to-AR":
            // Starting with a profile from a 2D simulation, map it to
            // a radial profile in a 3D simulation, ignoring the x-component
            // of the position of the ghost cells when computing distance and
            // picking the nearest point in the profile.
            other_r = sqrt(other_pos.y.re^^2 + other_pos.z.re^^2);
            my_r = sqrt(my_pos.y.re^^2 + my_pos.z.re^^2);
            dr = my_r - other_r;
            distance = fabs(dr);
            break;
        case "2Daxi-to-3D-rotating":
            // Map the y axis of a 2D simulation to the x-y plane of a 3D
            // simulation, ignoring the other axes. (NNG, June 2024)
            other_r = sqrt(other_pos.y.re^^2);
            my_r = sqrt(my_pos.y.re^^2 + my_pos.x.re^^2);
            dr = my_r - other_r;
            distance = fabs(dr);
            break;
        default:
            throw new FlowSolverException("Invalid match option.");
        }
        return distance;
    } // end compute_distance()

    @nogc
    size_t find_nearest_profile_point(ref const(Vector3) my_pos)
    {
        size_t ip = 0; // Start looking here, assuming that there is at least one point.
        double min_distance = compute_distance(my_pos, pos[0]);
        foreach (i; 1 .. pos.length) {
            double new_distance = compute_distance(my_pos, pos[i]);
            if (new_distance < min_distance) { ip = i; min_distance = new_distance; }
        }
        return ip;
    } // end find_nearest_profile_point()

    // not @nogc because of associative array lookup
    FlowState get_flowstate(size_t my_id, ref const(Vector3) my_pos)
    {
        assert(fstate.length > 0, "StaticFlowProfile is empty.");
        if (my_id in which_point) {
            return fstate[which_point[my_id]];
        } else {
            size_t ip = find_nearest_profile_point(my_pos);
            which_point[my_id] = ip;
            return fstate[ip];
        }
    } // end get_flowstate()

    @nogc
    void adjust_velocity(ref FlowState fs, ref const(Vector3) my_pos, double omegaz)
    {
        switch (posMatch) {
        case "xyz-to-xyz": /* 3D, do nothing. */ break;
        case "xyA-to-xyA": /* 3D, do nothing. */ break;
        case "AyA-to-AyA": /* 3D, do nothing. */ break;
        case "xy-to-xR": goto case "Ay-to-AR";
        case "Ay-to-AR":
            // We are assuming that the original 2D simulation had y>0.
            double r = sqrt(my_pos.y.re^^2 + my_pos.z.re^^2);
            double vel_yz = sqrt(fs.vel.y.re^^2 + fs.vel.z.re^^2);
            double vely_sign = (fs.vel.y < 0.0) ? -1.0 : 1.0;
            fs.vel.y = vely_sign * vel_yz * my_pos.y.re / r;
            fs.vel.z = vely_sign * vel_yz * my_pos.z.re / r;
            break;
        case "2Daxi-to-3D-rotating":
            double axi_vely = sqrt(fs.vel.y.re^^2);
            double axi_velx = fs.vel.x.re;

            double cyl_theta = atan2(my_pos.y.re, my_pos.x.re);
            double cyl_vel_x = axi_vely*cos(cyl_theta);
            double cyl_vel_y = axi_vely*sin(cyl_theta);
            double cyl_vel_z = axi_velx;

            fs.vel.x = cyl_vel_x;
            fs.vel.y = cyl_vel_y;
            fs.vel.z = cyl_vel_z;
            break;
        default:
            throw new FlowSolverException("Invalid match option.");
        }
        if (omegaz != 0.0) { into_rotating_frame(fs.vel, my_pos, omegaz); }
    }

    @nogc
    void adjust_velocity(FlowState* fs, ref const(Vector3) my_pos, double omegaz)
    {
        adjust_velocity(*fs, my_pos, omegaz);
    }
} // end class StaticFlowProfile


class TransientFlowProfile {
    // For use in the classes that implement the InflowBC_TransientProfile boundary condition.
    // GhostCellFlowStateCopyFromTransientProfile, BIE_FlowStateCopyFromTransientProfile
    // There are non-obvious options for the match parameter in the constructor call.
    // See the switch statement in the compute_distance() function for some hints.

public:
    string fileName;
    string posMatch;
    double[] times;
    Vector3[] pos;
    FlowState[][] fstate;

    // Construct a TransientFlowProfile from a zip archive.
    //
    // The zip archive contains a JSON file with metadata, as specified below,
    // and one data file for each time instant.
    // As for the StaticFlowProfile, the format for each data file will be a header line
    // followed by one sample point per line.
    // The header line will specify the names of the columns and, at a minimum, be:
    // pos.x pos.y p T vel.x vel.y
    // These names should match the variable names retrieved from the metadata.json file
    // and, although they are redundant, allow the direct use of files generated by
    // the slicing operations of the postprocessing programs.
    // The data lines will have numerical values for the corresponding quantities.
    this (string fileName, string match)
    {
        this.fileName = fileName;
        this.posMatch = match;
        //
        // Some context for the expected flow data.
        alias cfg = GlobalConfig;
        GasModel gmodel = cfg.gmodel_master;
        size_t n_species = gmodel.n_species;
        size_t n_modes = gmodel.n_modes;
        size_t n_turb = cfg.turb_model.nturb;
        //
        string[] expected_names = ["pos.x", "pos.y", "p", "T", "vel.x", "vel.y"];
        if (cfg.dimensions == 3) { expected_names ~= ["pos.z", "vel.z"]; }
        if (cfg.MHD) { expected_names ~= ["B.x", "B.y", "B.z"]; }
        string[] speciesList;
        if (n_species > 1) {
            foreach (i; 0..n_species) { speciesList ~= "massf-" ~ gmodel.species_name(i); }
            expected_names ~= speciesList;
        }
        string[] TmodeList;
        foreach (i; 0..n_modes) { TmodeList ~= "T-" ~ gmodel.energy_mode_name(to!int(i)); }
        expected_names ~= TmodeList;
        string[] turbList;
        foreach (i; 0..n_turb) { turbList ~= "tq-" ~ cfg.turb_model.primitive_variable_name(i); }
        expected_names ~= turbList;
        //
        // Open zip archive and read the metadata and the profiles.
        auto zip = new ZipArchive(read(fileName));
        auto zipMembers = zip.directory;
        foreach (name, member; zipMembers) { zip.expand(member); }
        //
        // The metadata comes as a JSON string that can be parsed for its values.
        //
        string metadataStr = assumeUTF(zipMembers["metadata.json"].expandedData);
        JSONValue metadata = parseJSON(metadataStr);
        //
        int nspecies = getJSONint(metadata, "nspecies", 0);
        if (nspecies != n_species) {
            string msg = text("Expected nspecies:", n_species, " got:", nspecies);
            throw new LmrException(msg);
        }
        int nmodes = getJSONint(metadata, "nmodes", 0);
        if (nmodes != n_modes) {
            string msg = text("Expected nmodes:", n_modes, " got:", nmodes);
            throw new LmrException(msg);
        }
        int nvar = getJSONint(metadata, "nvar", 0);
        string[] varnames; foreach (i; 0..nvar) { varnames ~= ""; }
        varnames = getJSONstringarray(metadata, "varnames", varnames);
        foreach (name; expected_names) {
            if (!canFind(varnames, name)) {
                string msg = format("Could not find required variable \"%s\" in file: %s", name, fileName);
                msg ~= format("\nFile header consisted of the following variables:\n %s", varnames);
                throw new LmrException(msg);
            }
        }
        size_t[string] column_dict; foreach (i, varname; varnames) { column_dict[varname] = i; }
        //
        ntimes = to!size_t(getJSONint(metadata, "ntimes", 0));
        if (ntimes == 0) {
            throw new LmrException("Expected at least one time instant in archive.");
        }
        foreach (k; 0..ntimes) { times ~= 0.0; }
        times = getJSONdoublearray(metadata, "times", times);
        //
        int npoints = getJSONint(metadata, "npoints", 0);
        if (npoints == 0) {
            throw new LmrException("Expected at least one sample point in profile.");
        }
        // Now that we know how much data is coming,
        // we can allocate suitable storage.
        foreach (j; 0..npoints) { pos ~= Vector3(); }
        fstate.length = ntimes;
        foreach (k; 0..ntimes) {
            foreach (j; 0..npoints) { fstate[k] ~= FlowState(gmodel, n_turb); }
        }
        //
        foreach (k; 0..ntimes) {
            string memberName = format("data.%d", k);
            string contentStr = assumeUTF(zipMembers[memberName].expandedData);
            // Format for each profile will be a header line followed by one sample point per line.
            // The header line will specify the names of the columns and, at a minimum, be:
            // pos.x pos.y p T vel.x vel.y
            // The data lines will have numerical values for the corresponding quantities.
            auto range = lineSplitter(contentStr);
            // First line is expected to name the columns with the correct variable names.
            // Just discard it..
            auto line = range.front;
            range.popFront();
            // Pick up the data lines and scan them.
            foreach (j; 0..npoints) {
                string[] items = to!string(range.front).strip().split();
                if (k == 0) {
                    // We use the position data from the first time-instant only.
                    Vector3* mypos = &pos[j];
                    mypos.x = to!number(items[column_dict["pos.x"]]);
                    mypos.y = to!number(items[column_dict["pos.y"]]);
                    mypos.z = (cfg.dimensions == 3) ? to!number(items[column_dict["pos.z"]]) : to!number(0.0);
                }
                FlowState* myfs = &fstate[k][j];
                myfs.gas.p = to!number(items[column_dict["p"]]);
                myfs.gas.T = to!number(items[column_dict["T"]]);
                foreach (i, name; TmodeList) {
                    myfs.gas.T_modes[i] = to!number(items[column_dict[name]]);
                }
                if (n_species > 1) {
                    foreach (i, name; speciesList) { myfs.gas.massf[i] = to!number(items[column_dict[name]]); }
                } else {
                    myfs.gas.massf[0] = 1.0;
                }
                gmodel.update_thermo_from_pT(myfs.gas);
                gmodel.update_sound_speed(myfs.gas);
                gmodel.update_trans_coeffs(myfs.gas);
                foreach (i; 0 .. n_species) { myfs.gas.rho_s[i] = myfs.gas.massf[i] * myfs.gas.rho; }
                myfs.vel.x = to!number(items[column_dict["vel.x"]]);
                myfs.vel.y = to!number(items[column_dict["vel.y"]]);
                myfs.vel.z = (cfg.dimensions == 3) ? to!number(items[column_dict["vel.z"]]) : to!number(0.0);
                foreach (i, name; turbList) { myfs.turb[i] = to!number(items[column_dict[name]]); }
                if (cfg.MHD) {
                    myfs.B.x = to!number(items[column_dict["B.x"]]);
                    myfs.B.y = to!number(items[column_dict["B.y"]]);
                    myfs.B.z = to!number(items[column_dict["B.z"]]);
                } else {
                    myfs.B.set(0.0, 0.0, 0.0);
                }
                range.popFront();
            } // end foreach j
        } // end foreach k
        //
        // The mapping of the nearest profile point to each ghost-cell or interface location
        // will be done as needed, at application time.
        // This way, all of the necessary cell and position data should be valid.
    } // end this()

    @nogc
    double compute_distance(ref const(Vector3) my_pos, ref const(Vector3) other_pos)
    {
        double distance, other_r, my_r, dx, dy, dz, dr;
        switch (posMatch) {
        case "xyz-to-xyz":
            // 2D or 3D, closest match on all components of position.
            // In 2D all z-components are supposed to be zero (and so, not matter).
            dx = my_pos.x.re - other_pos.x.re;
            dy = my_pos.y.re - other_pos.y.re;
            dz = my_pos.z.re - other_pos.z.re;
            distance = sqrt(dx*dx + dy*dy + dz*dz);
            break;
        case "xyA-to-xyA":
            // 2D or 3D; don't care about z-component of position.
            dx = my_pos.x.re - other_pos.x.re;
            dy = my_pos.y.re - other_pos.y.re;
            distance = sqrt(dx^^2 + dy^^2);
            break;
        case "AyA-to-AyA":
            // 2D or 3D; only care about the y-component of position.
            dy = my_pos.y.re - other_pos.y.re;
            distance = fabs(dy);
            break;
        case "xy-to-xR":
            // Starting with a profile from a 2D simulation, map it to
            // a radial profile in a 3D simulation, considering the x-component
            // of the position of the ghost cells when computing distance and
            // picking the nearest point in the profile.
            dx = my_pos.x.re - other_pos.x.re;
            other_r = sqrt(other_pos.y.re^^2 + other_pos.z.re^^2);
            my_r = sqrt(my_pos.y.re^^2 + my_pos.z.re^^2);
            dr = my_r - other_r;
            distance = sqrt(dx*dx + dr*dr);
            break;
        case "Ay-to-AR":
            // Starting with a profile from a 2D simulation, map it to
            // a radial profile in a 3D simulation, ignoring the x-component
            // of the position of the ghost cells when computing distance and
            // picking the nearest point in the profile.
            other_r = sqrt(other_pos.y.re^^2 + other_pos.z.re^^2);
            my_r = sqrt(my_pos.y.re^^2 + my_pos.z.re^^2);
            dr = my_r - other_r;
            distance = fabs(dr);
            break;
        case "2Daxi-to-3D-rotating":
            // Map the y axis of a 2D simulation to the x-y plane of a 3D
            // simulation, ignoring the other axes. (NNG, June 2024)
            other_r = sqrt(other_pos.y.re^^2);
            my_r = sqrt(my_pos.y.re^^2 + my_pos.x.re^^2);
            dr = my_r - other_r;
            distance = fabs(dr);
            break;
        default:
            throw new FlowSolverException("Invalid match option.");
        }
        return distance;
    } // end compute_distance()

    @nogc
    size_t find_nearest_profile_point(ref const(Vector3) my_pos)
    {
        size_t ip = 0; // Start looking here, assuming that there is at least one point.
        double min_distance = compute_distance(my_pos, pos[0]);
        foreach (i; 1 .. pos.length) {
            double new_distance = compute_distance(my_pos, pos[i]);
            if (new_distance < min_distance) { ip = i; min_distance = new_distance; }
        }
        return ip;
    } // end find_nearest_profile_point()

    @nogc
    void set_time_interpolation(double t)
    {
        if (ntimes == 1) {
            // Only one instance; much like a static profile.
            it0 = 0;
            it1 = 0;
            w0 = 1.0;
        } else {
            // Search for the correct pair of instances.
            it0 = ntimes - 2;
            while (it0 > 0 && times[it0] > t) { it0 -= 1; }
            it1 = it0 + 1;
            // Clipped interpolation for the weight.
            if (t < times[it0]) {
                w0 = 1.0;
            } else if (t > times[it1]) {
                w0 = 0.0;
            } else {
                w0 = (times[it1] - t)/(times[it1] - times[it0]);
            }
        }
        return;
    } // end set_time_interpolation()

    // not @nogc because of associative array lookup
    void set_flowstate(FlowState* fs, size_t my_id, ref const(Vector3) my_pos)
    {
        assert(fstate.length > 0, "TransientFlowProfile is empty.");
        size_t ip = (my_id in which_point) ? which_point[my_id] : find_nearest_profile_point(my_pos);
        fs.copy_average_values_from(fstate[it0][ip], fstate[it1][ip], w0);
        return;
    } // end set_flowstate()

    @nogc
    void adjust_velocity(ref FlowState fs, ref const(Vector3) my_pos, double omegaz)
    {
        switch (posMatch) {
        case "xyz-to-xyz": /* 3D, do nothing. */ break;
        case "xyA-to-xyA": /* 3D, do nothing. */ break;
        case "AyA-to-AyA": /* 3D, do nothing. */ break;
        case "xy-to-xR": goto case "Ay-to-AR";
        case "Ay-to-AR":
            // We are assuming that the original 2D simulation had y>0.
            double r = sqrt(my_pos.y.re^^2 + my_pos.z.re^^2);
            double vel_yz = sqrt(fs.vel.y.re^^2 + fs.vel.z.re^^2);
            double vely_sign = (fs.vel.y < 0.0) ? -1.0 : 1.0;
            fs.vel.y = vely_sign * vel_yz * my_pos.y.re / r;
            fs.vel.z = vely_sign * vel_yz * my_pos.z.re / r;
            break;
        case "2Daxi-to-3D-rotating":
            double axi_vely = sqrt(fs.vel.y.re^^2);
            double axi_velx = fs.vel.x.re;

            double cyl_theta = atan2(my_pos.y.re, my_pos.x.re);
            double cyl_vel_x = axi_vely*cos(cyl_theta);
            double cyl_vel_y = axi_vely*sin(cyl_theta);
            double cyl_vel_z = axi_velx;

            fs.vel.x = cyl_vel_x;
            fs.vel.y = cyl_vel_y;
            fs.vel.z = cyl_vel_z;
            break;
        default:
            throw new FlowSolverException("Invalid match option.");
        }
        if (omegaz != 0.0) { into_rotating_frame(fs.vel, my_pos, omegaz); }
    }

    @nogc
    void adjust_velocity(FlowState* fs, ref const(Vector3) my_pos, double omegaz)
    {
        adjust_velocity(*fs, my_pos, omegaz);
    }

private:
    size_t ntimes;
    size_t it0 = 0;
    size_t it1 = 0;
    double w0 = 0.0;
    size_t[size_t] which_point; // A place to memoize the mapped indices and we find them.
    // Above, we search for the profile point nearest to the initial position.
    // This position will only change for moving-grid simulations and we will not try
    // to deal with that complication.
} // end class TransientFlowProfile


class FlowHistory {
    // For use in the classes that implement the InflowBC_Transient boundary condition.
    // GhostCellFlowStateCopyFromHistory, BIE_FlowStateCopyFromHistory

public:
    string fileName;
    FlowState[] fstate;
    double[] times;

    // Construct a FlowHistory object from a file with columns of data values.
    //
    // Format will be a header line followed by one sample time per line.
    // The header line will specify the names of the columns and, at a minimum, be:
    // time p T vel.x vel.y
    // The data lines will have numerical values for the corresponding quantities.
    this (string fileName)
    {
        this.fileName = fileName;
        //
        // Some context for the expected flow data.
        alias cfg = GlobalConfig;
        GasModel gmodel = cfg.gmodel_master;
        size_t n_species = gmodel.n_species;
        size_t n_modes = gmodel.n_modes;
        size_t n_turb = cfg.turb_model.nturb;
        //
        string[] expected_names = ["time", "p", "T", "vel.x", "vel.y"];
        if (cfg.MHD) { expected_names ~= ["B.x", "B.y", "B.z"]; }
        string[] speciesList;
        if (n_species > 1) {
            foreach (i; 0..n_species) { speciesList ~= "massf-" ~ gmodel.species_name(i); }
        }
        string[] TmodeList;
        foreach (i; 0..n_modes) { TmodeList ~= "T-" ~ gmodel.energy_mode_name(to!int(i)); }
        string[] turbList;
        foreach (i; 0..n_turb) { turbList ~= "tq-" ~ cfg.turb_model.primitive_variable_name(i); }
        //
        // Open filename and read it.
        auto f = new File(fileName);
        auto range = f.byLine();
        //
        // First line is expected to name the columns with the correct variable names.
        // Assume that it does, and proceed to identify the columns.
        auto line = range.front;
        string txt = to!string(line);
        txt = stripLeft(txt, "# ").stripRight();
        string[] varnames = txt.split();
        size_t[string] column_dict; foreach (i, varname; varnames) { column_dict[varname] = i; }
        range.popFront();
        foreach (name; expected_names) {
            if (!canFind(varnames, name)) {
                string msg = text("FlowHistory: Did not find ", name, " in variables: ", varnames);
                throw new LmrException(msg);
            }
        }
        // Start picking up the data lines.
        line = range.front;
        while (!line.empty) {
            txt = to!string(line);
            double tme;
            auto myfs = FlowState(gmodel, n_turb);
            txt = stripLeft(txt).stripRight();
            string[] items = txt.split();
            tme = to!double(items[column_dict["time"]]);
            myfs.gas.p = to!number(items[column_dict["p"]]);
            myfs.gas.T = to!number(items[column_dict["T"]]);
            foreach (i, name; TmodeList) {
                myfs.gas.T_modes[i] = to!number(items[column_dict[name]]);
            }
            if (n_species > 1) {
                foreach (i, name; speciesList) { myfs.gas.massf[i] = to!number(items[column_dict[name]]); }
            } else {
                myfs.gas.massf[0] = 1.0;
            }
            gmodel.update_thermo_from_pT(myfs.gas);
            gmodel.update_sound_speed(myfs.gas);
            gmodel.update_trans_coeffs(myfs.gas);
            foreach (i; 0 .. n_species) { myfs.gas.rho_s[i] = myfs.gas.massf[i] * myfs.gas.rho; }
            myfs.vel.x = to!number(items[column_dict["vel.x"]]);
            myfs.vel.y = to!number(items[column_dict["vel.y"]]);
            myfs.vel.z = (cfg.dimensions == 3) ? to!number(items[column_dict["vel.z"]]) : to!number(0.0);
            foreach (i, name; turbList) { myfs.turb[i] = to!number(items[column_dict[name]]); }
            if (cfg.MHD) {
                myfs.B.x = to!number(items[column_dict["B.x"]]);
                myfs.B.y = to!number(items[column_dict["B.y"]]);
                myfs.B.z = to!number(items[column_dict["B.z"]]);
            } else {
                myfs.B.set(0.0, 0.0, 0.0);
            }
            fstate ~= myfs;
            times ~= tme;
            range.popFront();
            line = range.front;
        } // end while
        //
        if (fstate.length < 2) {
            throw new Error("FlowHistory did not have sufficient instances (minimum 2).");
        }
    } // end this()

    @nogc
    void set_flowstate(ref FlowState fs, double t, LocalConfig cfg)
    {
        GasModel gmodel = cfg.gmodel;
        size_t n_species = gmodel.n_species;
        size_t n_modes = gmodel.n_modes;
        size_t n_turb = cfg.turb_model.nturb;
        // Find where we are in history and interpolate the flow state.
        size_t nt = times.length;
        size_t k = 0;
        while ((k < nt-1) && t > times[k+1]) { k++; }
        k = min(k, nt-1);
        if (k < nt-1 && t <= times[$-1]) {
            // Linearly interpolate between states k, k+1
            double frac1 = (t-times[k])/(times[k+1]-times[k]);
            double frac0 = 1.0-frac1;
            fs.vel.x = fstate[k].vel.x*frac0 + fstate[k+1].vel.x*frac1;
            fs.vel.y = fstate[k].vel.y*frac0 + fstate[k+1].vel.y*frac1;
            fs.vel.z = fstate[k].vel.z*frac0 + fstate[k+1].vel.z*frac1;
            fs.gas.p = fstate[k].gas.p*frac0 + fstate[k+1].gas.p*frac1;
            fs.gas.T = fstate[k].gas.T*frac0 + fstate[k+1].gas.T*frac1;
            foreach (i; 0 .. n_species) {
                fs.gas.massf[i] = fstate[k].gas.massf[i]*frac0 + fstate[k+1].gas.massf[i]*frac1;
            }
            foreach (i; 0 .. n_modes) {
                fs.gas.T_modes[i] = fstate[k].gas.T_modes[i]*frac0 + fstate[k+1].gas.T_modes[i]*frac1;
            }
            gmodel.update_thermo_from_pT(fs.gas);
            gmodel.update_sound_speed(fs.gas);
            gmodel.update_trans_coeffs(fs.gas);
            foreach (i; 0 .. n_species) { fs.gas.rho_s[i] = fs.gas.massf[i] * fs.gas.rho; }
            if (cfg.MHD) {
                fs.B.x = fstate[k].B.x*frac0 + fstate[k+1].B.x*frac1;
                fs.B.y = fstate[k].B.y*frac1 + fstate[k+1].B.y*frac1;
                fs.B.z = fstate[k].B.z*frac0 + fstate[k+1].B.z*frac1;
            } else {
                fs.B.set(0.0, 0.0, 0.0);
            }
            version(turbulence) {
                foreach(i; 0 .. n_turb) {
                    fs.turb[i] = fstate[k].turb[i]*frac0 + fstate[k+1].turb[i]*frac1;
                }
            }
            fs.mu_t = 0.0;
            fs.k_t = 0.0;
        } else {
            // Keep condition constant beyond the largest time.
            fs.copy_values_from(fstate[$-1]);
        }
        return;
    } // end set_flowstate()

} // end FlowHistory


class SyntheticFlowState {
    // For use in the classes that implement the InflowBC_Synthetic boundary condition.
    // GhostCellSynthesiseFlowState, BIE_SynthesiseFlowState
    //
    // 2021-07-20 PJ
    // This is a place-holder for Lachlan's synthetic boundary condition.
    // For the moment, just leave the history code in place.

public:
    string fileName;

    this (string fileName)
    {
        this.fileName = fileName;
        // Open filename and read the defining data in JSON format.
        auto gm = GlobalConfig.gmodel_master;
        import std.json;
        import util.json_helper;
        JSONValue jsonData = readJSONfile(fileName);
        //
        // Lachlan, you may decide what you want to do here and below in set_flowstate().
        //
        auto baseFlow = jsonData["base"];
        base_velx = getJSONdouble(baseFlow, "velx", 0.0);
        base_vely = getJSONdouble(baseFlow, "vely", 0.0);
        base_velz = getJSONdouble(baseFlow, "velz", 0.0);
        base_p = getJSONdouble(baseFlow, "p", 1.0e5);
        base_T = getJSONdouble(baseFlow, "T", 300.0);
    } // end this()

    @nogc
    void set_flowstate(ref FlowState fs, double t, ref Vector3 pos, GasModel gm)
    {
        fs.vel.x = base_velx;
        fs.vel.y = base_vely;
        fs.vel.z = base_velz;
        fs.gas.p = base_p;
        fs.gas.T = base_T;
        fs.gas.massf[0] = 1.0;
        foreach (j; 1 .. gm.n_species) { fs.gas.massf[j] = 0.0; }
        foreach (j; 0 .. gm.n_modes) { fs.gas.T_modes[j] = base_T; }
        gm.update_thermo_from_pT(fs.gas);
        foreach (j; 0 .. gm.n_species) { fs.gas.rho_s[j] = fs.gas.massf[j] * fs.gas.rho; }
        return;
    } // end set_flowstate()

private:
    double base_velx;
    double base_vely;
    double base_velz;
    double base_p;
    double base_T;
} // end class SyntheticFlowState


class SourceFlow {
    // For computing flow state increments when sampling flow from a conical nozzle.
    // We assume that the sample points are nearby the nominal distance from the
    // virtual origin of the source flow.
    // See PJ's workbook pages 24-27, 2021-07-29.

public:
    double r;  // Distance from virtual origin till nominal flow condition.
    double p, rho, u, v; // Nominal flow condition.
    double dfdrho, dfdu; // Sensitivities for EOS p = f(rho, u)

    this(GasModel gmodel, ref const(FlowState) fs, double r)
    {
        this.r = r;
        // Keep a copy of the interesting parts of the nominal state.
        p = fs.gas.p.re;
        rho = fs.gas.rho.re;
        u = gmodel.internal_energy(fs.gas).re;
        v = fs.vel.x.re;
        // Derivatives of EOS that are needed when computing increments later.
        GasState gs = GasState(fs.gas);
        gmodel.update_thermo_from_rhou(gs);
        if (fabs(gs.p.re - p)/p > 1.0e-5) {
            throw new Error("Pressure mismatch in gas states that should be equal. What's up?");
        }
        double drho = 0.001 * rho;
        gs.rho += drho;
        gmodel.update_thermo_from_rhou(gs);
        dfdrho = (gs.p.re - p) / drho;
        double du = 0.001 * fabs(u) + 100.0;  // It is possible that u might be zero.
        gs.rho = rho; gs.u += du;
        gmodel.update_thermo_from_rhou(gs);
        dfdu = (gs.p.re - p) / du;
    }

    @nogc
    double[4] get_rho_v_p_u_increments(double dr)
    {
        // Compute flow-state increments for a point at a slightly different distance
        // from the virtual origin.
        double dAonA = 2.0*dr/r;
        double q1 = dfdrho*rho + dfdu*p/rho;
        double denom = rho*v*v - q1;
        double q2 = dAonA/denom;
        double drho = -rho*rho*v*v*q2;
        double dv = v*q1*q2;
        double dp = -rho*v*v*q1*q2;
        double du = -p*v*v*q2;
        return [drho, dv, dp, du];
    }
} // end SourceFlow
