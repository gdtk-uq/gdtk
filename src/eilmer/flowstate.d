/**
 * flowstate.d
 * FlowState class for use in the main solver.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module flowstate;

import std.string;
import std.conv;
import std.json;
import std.array;
import std.format;
import std.stdio;
import json_helper;
import geom;
import gas;
import fvcell;
import sgrid;
import globalconfig;

class FlowState {
public:
    GasState gas;  // gas state
    Vector3 vel;   // flow velocity, m/s
    Vector3 B;     // magnetic field strength
    double tke;    // turbulent kinetic energy 0.5(u'^2+v'^2+w'^2)
    double omega;  // turbulence 'frequency' in k-omega model
    double mu_t;   // turbulence viscosity
    double k_t;    // turbulence thermal-conductivity
    int S;         // shock indicator, value 0 or 1

    this(GasModel gm, in double p_init, in double[] T_init, in Vector3 vel_init,
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in Vector3 B_init=(0.0,0.0,0.0),
	 in double tke_init=0.0, in double omega_init=1.0,
	 in double mu_t_init=0.0, in double k_t_init=0.0,
	 in int S_init=0)
    {
	gas = new GasState(gm, p_init, T_init, massf_init, quality_init);
	vel = vel_init;
	B = B_init;
	tke = tke_init;
	omega = omega_init;
	mu_t = mu_t_init;
	k_t = k_t_init;
	S = S_init;
    }

    this(in FlowState other, GasModel gm)
    {
	gas = new GasState(gm, other.gas.p, other.gas.T, other.gas.massf, other.gas.quality); 
	vel = other.vel;
	B = other.B;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    this(in FlowState other)
    {
	gas = new GasState(to!int(other.gas.massf.length), to!int(other.gas.T.length));
	gas.copy_values_from(other.gas);
	vel.refx = other.vel.x; vel.refy = other.vel.y; vel.refz = other.vel.z;
	B.refx = other.B.x; B.refy = other.B.y; B.refz = other.B.z;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    this(GasModel gm)
    {
	gas = new GasState(gm, 100.0e3, [300.0,], [1.0,], 1.0); 
	vel = Vector3(0.0,0.0,0.0);
	B = Vector3(0.0,0.0,0.0);
	tke = 0.0;
	omega = 1.0;
	mu_t = 0.0;
	k_t = 0.0;
	S = 0;
    }

    this(in JSONValue json_data, GasModel gm)
    {
	double p = getJSONdouble(json_data, "p", 100.0e3);
	double[] T = getJSONdoublearray(json_data, "T", [300.0,]);
	double[] massf = getJSONdoublearray(json_data, "massf", [1.0,]);
	double quality = 1.0;
	gas = new GasState(gm, p, T, massf, quality);
	double velx = getJSONdouble(json_data, "velx", 0.0);
	double vely = getJSONdouble(json_data, "vely", 0.0);
	double velz = getJSONdouble(json_data, "velz", 0.0);
	vel = Vector3(velx,vely,velz);
	double Bx = getJSONdouble(json_data, "Bx", 0.0);
	double By = getJSONdouble(json_data, "By", 0.0);
	double Bz = getJSONdouble(json_data, "Bz", 0.0);
	B = Vector3(Bx,By,Bz);
	tke = getJSONdouble(json_data, "tke", 0.0);
	omega = getJSONdouble(json_data, "omega", 1.0);
	mu_t = getJSONdouble(json_data, "mu_t", 0.0);
	k_t = getJSONdouble(json_data, "k_t", 0.0);
	S = getJSONint(json_data, "S", 0);
    }

    this() {} // makes no sense to define the data in the absence of a model

    FlowState dup() const
    {
	return new FlowState(this);
    }

    @nogc 
    void copy_values_from(in FlowState other)
    {
	gas.copy_values_from(other.gas);
	vel.refx = other.vel.x; vel.refy = other.vel.y; vel.refz = other.vel.z;
	B.refx = other.B.x; B.refy = other.B.y; B.refz = other.B.z;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    @nogc 
    void copy_average_values_from(in FlowState fs0, in FlowState fs1)
    // Avoids memory allocation, it's all in place.
    {
	gas.copy_average_values_from(fs0.gas, fs1.gas);
	vel.refx = 0.5 * (fs0.vel.x + fs1.vel.x);
	vel.refy = 0.5 * (fs0.vel.y + fs1.vel.y);
	vel.refz = 0.5 * (fs0.vel.z + fs1.vel.z);
	B.refx = 0.5 * (fs0.B.x + fs1.B.x);
	B.refy = 0.5 * (fs0.B.y + fs1.B.y);
	B.refz = 0.5 * (fs0.B.z + fs1.B.z);
	tke = 0.5 * (fs0.tke + fs1.tke);
	omega = 0.5 * (fs0.omega + fs1.omega);
	mu_t = 0.5 * (fs0.mu_t + fs1.mu_t);
	k_t = 0.5 * (fs0.k_t + fs1.k_t);
    } // end copy_average_values_from()

    void copy_average_values_from(in FlowState[] others, GasModel gm)
    // Note that we must not send the current object in the others list as well.
    // Involves some memory allocation.
    {
	size_t n = others.length;
	if (n == 0) throw new Error("Need to average from a nonempty array.");
	GasState[] gasList;
	// Note that, because we cast away their "const"ness,
	// we need to be honest and not to fiddle with the other gas states.
	foreach(other; others) {
	    if ( this is other ) {
		throw new Error("Must not include destination in source list.");
	    }
	    gasList ~= cast(GasState)other.gas;
	}
	gas.copy_average_values_from(gasList, gm);
	// Accumulate from a clean slate and then divide.
	vel.refx = 0.0; vel.refy = 0.0; vel.refz = 0.0;
	B.refx = 0.0; B.refy = 0.0; B.refz = 0.0;
	tke = 0.0;
	omega = 0.0;
	mu_t = 0.0;
	k_t = 0.0;
	S = 0; // Remember that shock detector is an integer flag.
	foreach(other; others) {
	    vel.refx += other.vel.x;
	    vel.refy += other.vel.y;
	    vel.refz += other.vel.z;
	    B.refx += other.B.x;
	    B.refx += other.B.x;
	    B.refx += other.B.x;
	    tke += other.tke;
	    omega += other.omega;
	    mu_t += other.mu_t;
	    k_t += other.k_t;
	    S += other.S;
	}
	vel /= n;
	B /= n;
	tke /= n;
	omega /= n;
	mu_t /= n;
	k_t /= n;
	S = (S > 0) ? 1 : 0;
    } // end copy_average_values_from()

    override string toString() const
    {
	char[] repr;
	repr ~= "FlowState(";
	repr ~= "gas=" ~ to!string(gas);
	repr ~= ", vel=" ~ to!string(vel);
	repr ~= ", B=" ~ to!string(B);
	repr ~= ", tke=" ~ to!string(tke);
	repr ~= ", omega=" ~ to!string(omega);
	repr ~= ", mu_t=" ~ to!string(mu_t);
	repr ~= ", k_t=" ~ to!string(k_t);
	repr ~= ", S=" ~ to!string(S);
	repr ~= ")";
	return to!string(repr);
    }

    string toJSONString() const
    {
	auto writer = appender!string();
	formattedWrite(writer, "{");
	formattedWrite(writer, "\"p\": %.12e", gas.p);
	formattedWrite(writer, ", \"T\": [ %.12e", gas.T[0]);
	foreach (i; 1 .. gas.T.length) {
	    formattedWrite(writer, ", %.12e", gas.T[i]);
	}
	formattedWrite(writer, "]");
	formattedWrite(writer, ", \"massf\": [ %.12e", gas.massf[0]);
	foreach (i; 1 .. gas.massf.length) {
	    formattedWrite(writer, ", %.12e", gas.massf[i]);
	}
	formattedWrite(writer, "]");
	// double quality = 1.0;
	formattedWrite(writer, ", \"velx\": %.12e", vel.x);
	formattedWrite(writer, ", \"vely\": %.12e", vel.y);
	formattedWrite(writer, ", \"velz\": %.12e", vel.z);
	formattedWrite(writer, ", \"Bx\": %.12e", B.x);
	formattedWrite(writer, ", \"By\": %.12e", B.y);
	formattedWrite(writer, ", \"Bz\": %.12e", B.z);
	formattedWrite(writer, ", \"tke\": %.12e", tke);
	formattedWrite(writer, ", \"omega\": %.12e", omega);
	formattedWrite(writer, ", \"mu_t\": %.12e", mu_t);
	formattedWrite(writer, ", \"k_t\": %.12e", k_t);
	formattedWrite(writer, ", \"S\": %d", S);
	formattedWrite(writer, "}");
	return writer.data;
    } // end toJSONString()

/+ [TODO]
    double * copy_values_to_buffer(double *buf) const;
    double * copy_values_from_buffer(double *buf);
+/
} // end class FlowState

void write_initial_flow_file(string fileName, ref StructuredGrid grid,
			     in FlowState fs, double t0, GasModel gmodel)
{
    // Numbers of cells derived from numbers of vertices.
    int nic = grid.niv - 1;
    int njc = grid.njv - 1;
    int nkc = grid.nkv - 1;
    if (GlobalConfig.dimensions == 2) nkc = 1;
    //	
    string cell_data_as_string(int i, int j, int k)
    {
	auto p000 = grid[i,j,k];
	auto p100 = grid[i+1,j,k];
	auto p110 = grid[i+1,j+1,k];
	auto p010 = grid[i,j+1,k];
	// [TODO] provide better calculation using geom module.
	// For the moment, it doesn't matter greatly because the solver 
	// will compute it's own approximations
	auto pos = 0.25*(p000 + p100 + p110 + p010);
	auto volume = 0.0;
	if (GlobalConfig.dimensions == 3) {
	    auto p001 = grid[i,j,k+1];
	    auto p101 = grid[i+1,j,k+1];
	    auto p111 = grid[i+1,j+1,k+1];
	    auto p011 = grid[i,j+1,k+1];
	    pos = 0.5*pos + 0.125*(p001 + p101 + p111 + p011);
	}
	// Should match FVCell.write_values_to_string()
	auto writer = appender!string();
	formattedWrite(writer, "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e",
		       pos.x, pos.y, pos.z, volume, fs.gas.rho,
		       fs.vel.x, fs.vel.y, fs.vel.z);
	if (GlobalConfig.MHD)
	    formattedWrite(writer, " %.12e %.12e %.12e", fs.B.x, fs.B.y, fs.B.z); 
	formattedWrite(writer, " %.12e %.12e %.12e", fs.gas.p, fs.gas.a, fs.gas.mu);
	foreach (kvalue; fs.gas.k) formattedWrite(writer, " %.12e", kvalue); 
	int S = 0;  // zero for shock detector
	formattedWrite(writer, " %.12e %.12e %d", fs.mu_t, fs.k_t, S);
	if (GlobalConfig.radiation) {
	    double Q_rad_org = 0.0; double f_rad_org = 0.0; double Q_rE_rad = 0.0;
	    formattedWrite(writer, " %.12e %.12e %.12e", Q_rad_org, f_rad_org, Q_rE_rad);
	}
	formattedWrite(writer, " %.12e %.12e", fs.tke, fs.omega);
	foreach (massfvalue; fs.gas.massf) formattedWrite(writer, " %.12e", massfvalue); 
	double dt_chem = -1.0;
	if (fs.gas.massf.length > 1) formattedWrite(writer, " %.12e", dt_chem); 
	foreach (imode; 0 .. fs.gas.e.length) 
	    formattedWrite(writer, " %.12e %.12e", fs.gas.e[imode], fs.gas.T[imode]); 
	double dt_therm = -1.0;
	if (fs.gas.e.length > 1) formattedWrite(writer, " %.12e", dt_therm); 
	return writer.data;
    } // end  cell_data_as_string()

    // Write the data for the whole structured block.
    auto f = File(fileName, "w");
    f.writefln("%20.12e", t0);
    // Variable list for cell on one line.
    auto writer = appender!string();
    foreach(varname; variable_list_for_cell(gmodel)) {
	formattedWrite(writer, " \"%s\"", varname);
    }
    f.writeln(writer.data);
    // Numbers of cells.
    f.writefln("%d %d %d", nic, njc, nkc);
    // The actual cell data.
    foreach (k; 0 .. nkc) {
	foreach (j; 0 .. njc) {
	    foreach (i; 0 .. nic) {
		f.writefln(" %s", cell_data_as_string(i,j,k));
	    }
	}
    }
    f.close();
    return;
} // end write_initial_flow_file()
