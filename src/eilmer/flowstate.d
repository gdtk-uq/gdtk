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
import std.algorithm;
import std.json;
import std.array;
import std.format;
import std.stdio;
import std.math;
import json_helper;
import gzip;
import geom;
import gas;
import fvcore;
import fvcell;
import sgrid;
import usgrid;
import globalconfig;

class FlowState {
public:
    GasState gas;  // gas state
    Vector3 vel;   // flow velocity, m/s
    Vector3 B;     // magnetic field strength
    double psi;	   // divergence cleaning parameter
    double divB;   // divergence of the magnetic field
    double tke;    // turbulent kinetic energy 0.5(u'^2+v'^2+w'^2)
    double omega;  // turbulence 'frequency' in k-omega model
    double mu_t;   // turbulence viscosity
    double k_t;    // turbulence thermal-conductivity
    int S;         // shock indicator, value 0 or 1

    this(GasModel gm,
	 in double p_init,
	 in double Ttr_init,
	 in double[] T_modes_init,
	 in Vector3 vel_init,
	 in double[] massf_init=[1.0,],
	 in double quality_init=1.0,
	 in Vector3 B_init=Vector3(0.0,0.0,0.0),
	 in double psi_init=0.0, in double divB_init=1.0,
	 in double tke_init=0.0, in double omega_init=1.0,
	 in double mu_t_init=0.0, in double k_t_init=0.0,
	 in int S_init=0)
    {
	gas = new GasState(gm, p_init, Ttr_init, T_modes_init,
			   massf_init, quality_init);
	vel = vel_init;
	B = B_init;
	psi = psi_init;
	divB = divB_init;
	tke = tke_init;
	omega = omega_init;
	mu_t = mu_t_init;
	k_t = k_t_init;
	S = S_init;
    }

    this(in FlowState other, GasModel gm)
    {
	gas = new GasState(gm, other.gas.p, other.gas.Ttr, other.gas.T_modes,
			   other.gas.massf, other.gas.quality); 
	vel = other.vel;
	B = other.B;
	psi = other.psi;
	divB = other.divB;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    this(in FlowState other)
    {
	gas = new GasState(to!int(other.gas.massf.length),
			   to!int(other.gas.T_modes.length));
	gas.copy_values_from(other.gas);
	vel.refx = other.vel.x; vel.refy = other.vel.y; vel.refz = other.vel.z;
	B.refx = other.B.x; B.refy = other.B.y; B.refz = other.B.z;
	psi = other.psi;
	divB = other.divB;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    this(GasModel gm)
    {
	gas = new GasState(gm, 100.0e3, 300.0, [1.0,], 1.0); 
	vel = Vector3(0.0,0.0,0.0);
	B = Vector3(0.0,0.0,0.0);
	psi = 0.0;
	divB = 0.0;
	tke = 0.0;
	omega = 1.0;
	mu_t = 0.0;
	k_t = 0.0;
	S = 0;
    }

    this(in JSONValue json_data, GasModel gm)
    {
	double p = getJSONdouble(json_data, "p", 100.0e3);
	double Ttr = getJSONdouble(json_data, "Ttr", 300.0e3);
	double[] T_modes;
	foreach(i; 0 .. gm.n_modes) { T_modes ~= Ttr; }
	T_modes = getJSONdoublearray(json_data, "T_modes", []);
	double[] massf = getJSONdoublearray(json_data, "massf", [1.0,]);
	double quality = getJSONdouble(json_data, "quality", 1.0);
	gas = new GasState(gm, p, Ttr, T_modes, massf, quality);
	double velx = getJSONdouble(json_data, "velx", 0.0);
	double vely = getJSONdouble(json_data, "vely", 0.0);
	double velz = getJSONdouble(json_data, "velz", 0.0);
	vel = Vector3(velx,vely,velz);
	double Bx = getJSONdouble(json_data, "Bx", 0.0);
	double By = getJSONdouble(json_data, "By", 0.0);
	double Bz = getJSONdouble(json_data, "Bz", 0.0);
	B = Vector3(Bx,By,Bz);
	psi = getJSONdouble(json_data, "psi", 0.0);
	divB = getJSONdouble(json_data, "divB", 0.0);
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
	psi = other.psi;
	divB = other.divB;
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
	psi = 0.5 * (fs0.psi + fs1.psi);
	divB = 0.5 * (fs0.divB + fs1.divB);
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
	if (n == 0) throw new FlowSolverException("Need to average from a nonempty array.");
	GasState[] gasList;
	// Note that, because we cast away their "const"ness,
	// we need to be honest and not to fiddle with the other gas states.
	foreach(other; others) {
	    if ( this is other ) {
		throw new FlowSolverException("Must not include destination in source list.");
	    }
	    gasList ~= cast(GasState)other.gas;
	}
	gas.copy_average_values_from(gasList, gm);
	// Accumulate from a clean slate and then divide.
	vel.refx = 0.0; vel.refy = 0.0; vel.refz = 0.0;
	B.refx = 0.0; B.refy = 0.0; B.refz = 0.0;
	psi = 0.0;
	divB = 0.0;
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
	    psi += other.psi;
	    divB += other.divB;
	    tke += other.tke;
	    omega += other.omega;
	    mu_t += other.mu_t;
	    k_t += other.k_t;
	    S += other.S;
	}
	vel /= n;
	B /= n;
	psi /= n;
	divB /= n;
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
	repr ~= ", psi=" ~ to!string(psi);
	repr ~= ", divB=" ~ to!string(psi);
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
	formattedWrite(writer, "\"p\": %.18e", gas.p);
	formattedWrite(writer, ", \"Ttr\": %.18e", gas.Ttr);
	// zero or more T_modes
	formattedWrite(writer, ", \"T_modes\": [");
	if (gas.T_modes.length > 0) { formattedWrite(writer, " %.18e", gas.T_modes[0]); }
	foreach (i; 1 .. gas.T_modes.length) { formattedWrite(writer, ", %.18e", gas.T_modes[i]); }
	// one or more mass fractions
	formattedWrite(writer, "]");
	formattedWrite(writer, ", \"massf\": [ %.18e", gas.massf[0]);
	foreach (i; 1 .. gas.massf.length) {
	    formattedWrite(writer, ", %.18e", gas.massf[i]);
	}
	formattedWrite(writer, "]");
	formattedWrite(writer, ", \"quality\": %.18e", gas.quality);
	formattedWrite(writer, ", \"velx\": %.18e", vel.x);
	formattedWrite(writer, ", \"vely\": %.18e", vel.y);
	formattedWrite(writer, ", \"velz\": %.18e", vel.z);
	formattedWrite(writer, ", \"Bx\": %.18e", B.x);
	formattedWrite(writer, ", \"By\": %.18e", B.y);
	formattedWrite(writer, ", \"Bz\": %.18e", B.z);
	formattedWrite(writer, ", \"psi\": %.18e", psi);
	formattedWrite(writer, ", \"divB\": %.18e", divB);
	formattedWrite(writer, ", \"tke\": %.18e", tke);
	formattedWrite(writer, ", \"omega\": %.18e", omega);
	formattedWrite(writer, ", \"mu_t\": %.18e", mu_t);
	formattedWrite(writer, ", \"k_t\": %.18e", k_t);
	formattedWrite(writer, ", \"S\": %d", S);
	formattedWrite(writer, "}");
	return writer.data;
    } // end toJSONString()


    bool check_data(ref Vector3 pos, ref const(FlowStateLimits) flowstate_limits) const
    {
	bool is_data_valid = gas.check_values(true);
	if (fabs(vel.x) > flowstate_limits.max_velocity ||
	    fabs(vel.y) > flowstate_limits.max_velocity ||
	    fabs(vel.z) > flowstate_limits.max_velocity) {
	    writeln("Velocity too high ", vel);
	    is_data_valid = false;
	}
	if (!isFinite(tke)) {
	    writeln("Turbulence KE invalid number ", tke);
	    is_data_valid = false;
	}
	if (tke < flowstate_limits.min_tke) {
	    writeln("Turbulence KE below minimum ", tke);
	    is_data_valid = false;
	}
	if (tke > flowstate_limits.max_tke) {
	    writeln("Turbulence KE above maximum ", tke);
	    is_data_valid = false;
	}
	if (gas.Ttr < flowstate_limits.min_temp) {
	    writeln("Temperature below minimum ", gas.Ttr);
	    is_data_valid = false;
	}
	if (gas.Ttr > flowstate_limits.max_temp) {
	    writeln("Temperature above maximum ", gas.Ttr);
	    is_data_valid = false;
	}
	if (!isFinite(omega)) {
	    writeln("Turbulence frequency invalid number ", omega);
	    is_data_valid = false;
	}
	if (omega <= 0.0) {
	    writeln("Turbulence frequency nonpositive ", omega);
	    is_data_valid = false;
	}
	if (!is_data_valid) {
	    writeln("   at position ", pos);
	}
	return is_data_valid;
    } // end check_data()

    @nogc
    void reorient_vector_quantities(const(double[]) Rmatrix)
    {
	vel.apply_matrix_transform(Rmatrix);
	B.apply_matrix_transform(Rmatrix);
    }
} // end class FlowState
 
void write_initial_flow_file(string fileName, ref StructuredGrid grid,
			     in FlowState fs, double t0, GasModel gmodel)
// Keep in sync with SBlock.write_solution.
{
    // Numbers of cells derived from numbers of vertices in grid.
    auto nicell = grid.niv - 1;
    auto njcell = grid.njv - 1;
    auto nkcell = grid.nkv - 1;
    if (GlobalConfig.dimensions == 2) nkcell = 1;
    //	
    // Write the data for the whole structured block.
    auto outfile = new GzipOut(fileName);
    auto writer = appender!string();
    formattedWrite(writer, "structured_grid_flow 1.0\n");
    formattedWrite(writer, "label: %s\n", grid.label);
    formattedWrite(writer, "sim_time: %.18e\n", t0);
    auto variable_list = variable_list_for_cell(gmodel, GlobalConfig.include_quality, 
						GlobalConfig.MHD, GlobalConfig.divergence_cleaning,
						GlobalConfig.radiation);
    formattedWrite(writer, "variables: %d\n", variable_list.length);
    // Variable list for cell on one line.
    foreach(varname; variable_list) {
	formattedWrite(writer, " \"%s\"", varname);
    }
    formattedWrite(writer, "\n");
    // Numbers of cells
    formattedWrite(writer, "dimensions: %d\n", GlobalConfig.dimensions);
    formattedWrite(writer, "nicell: %d\n", nicell);
    formattedWrite(writer, "njcell: %d\n", njcell);
    formattedWrite(writer, "nkcell: %d\n", nkcell);
    outfile.compress(writer.data);
    // The actual cell data.
    foreach (k; 0 .. nkcell) {
	foreach (j; 0 .. njcell) {
	    foreach (i; 0 .. nicell) {
		Vector3 p000 = *grid[i,j,k];
		Vector3 p100 = *grid[i+1,j,k];
		Vector3 p110 = *grid[i+1,j+1,k];
		Vector3 p010 = *grid[i,j+1,k];
		// [TODO] provide better calculation using geom module.
		// For the moment, it doesn't matter greatly because the solver 
		// will compute it's own approximations
		auto pos = 0.25*(p000 + p100 + p110 + p010);
		double volume = 0.0; 
		if (GlobalConfig.dimensions == 3) {
		    Vector3 p001 = *grid[i,j,k+1];
		    Vector3 p101 = *grid[i+1,j,k+1];
		    Vector3 p111 = *grid[i+1,j+1,k+1];
		    Vector3 p011 = *grid[i,j+1,k+1];
		    pos = 0.5*pos + 0.125*(p001 + p101 + p111 + p011);
		}
		outfile.compress(" " ~ cell_data_as_string(pos, volume, fs,
							   0.0, 0.0, 0.0, -1.0, -1.0,
							   GlobalConfig.include_quality,
							   GlobalConfig.MHD,
							   GlobalConfig.divergence_cleaning,
							   GlobalConfig.radiation) ~ "\n");
	    }
	}
    }
    outfile.finish();
    return;
} // end write_initial_flow_file() StructuredGrid version
 
void write_initial_flow_file(string fileName, ref UnstructuredGrid grid,
			     in FlowState fs, double t0, GasModel gmodel)
// Keep in sync with UBlock.write_solution.
{
    // Numbers of cells derived from numbers of vertices in grid.
    auto ncells = grid.ncells;
    //	
    // Write the data for the whole unstructured block.
    auto outfile = new GzipOut(fileName);
    auto writer = appender!string();
    formattedWrite(writer, "unstructured_grid_flow 1.0\n");
    formattedWrite(writer, "label: %s\n", grid.label);
    formattedWrite(writer, "sim_time: %.18e\n", t0);
    auto variable_list = variable_list_for_cell(gmodel, GlobalConfig.include_quality,
						GlobalConfig.MHD, GlobalConfig.divergence_cleaning,
						GlobalConfig.radiation);
    formattedWrite(writer, "variables: %d\n", variable_list.length);
    // Variable list for cell on one line.
    foreach(varname; variable_list) {
	formattedWrite(writer, " \"%s\"", varname);
    }
    formattedWrite(writer, "\n");
    // Numbers of cells
    formattedWrite(writer, "dimensions: %d\n", GlobalConfig.dimensions);
    formattedWrite(writer, "ncells: %d\n", ncells);
    outfile.compress(writer.data);
    // The actual cell data.
    foreach (i; 0 .. ncells) {
	Vector3 pos = Vector3(0.0, 0.0, 0.0);
	foreach (id; grid.cells[i].vtx_id_list) { pos += grid.vertices[id]; }
	pos /= grid.cells[i].vtx_id_list.length;
	double volume = 0.0; 
	outfile.compress(" " ~ cell_data_as_string(pos, volume, fs,
						   0.0, 0.0, 0.0, -1.0, -1.0,
						   GlobalConfig.include_quality,
						   GlobalConfig.MHD,
						   GlobalConfig.divergence_cleaning,
						   GlobalConfig.radiation) ~ "\n");
    }
    outfile.finish();
    return;
} // end write_initial_flow_file() UnstructuredGrid version


class FlowProfile {
    // For use in the classes that implement the InflowBC_StaticProfile boundary condition.
    // GhostCellFlowStateCopyFromProfile, BIE_FlowStateCopyFromProfile
    
public:
    string fileName;
    string posMatch;
    FlowState[] fstate;
    Vector3[] pos;
    size_t[size_t] which_point; // A place to memoize the mapped indices and we find them.
    // Below, we search for the profile point nearest to the initial position.
    // This position will only change for moving-grid simulations and we will not try
    // to deal with that complication.

    this (string fileName, string match)
    {
	this.fileName = fileName;
	this.posMatch = match;
	// Open filename and read all data points.
	// Format will be sample point data as per the postprocessor.
	auto f = new File(fileName);
	auto range = f.byLine();
	auto line = range.front;
	int npoints = 0;
	while (!line.empty) {
	    string text = to!string(line);
	    if (!canFind(text, "#")) {
		// Assume that we have a line of data rather than variable names.
		fstate ~= new FlowState(GlobalConfig.gmodel_master);
		pos ~= Vector3();
		double volume, Q_rad_org, f_rad_org, Q_rE_rad, dt_chem, dt_therm;
		scan_cell_data_from_string(text, pos[$-1], volume, fstate[$-1],
					   Q_rad_org, f_rad_org, Q_rE_rad, dt_chem, dt_therm,
					   GlobalConfig.include_quality, GlobalConfig.MHD,
					   GlobalConfig.divergence_cleaning, GlobalConfig.radiation);
		npoints += 1;
	    }
	    range.popFront();
	    line = range.front;
	} // end while
	// writefln("FlowProfile: file=\"%s\", match=\"%s\", npoints=%d", fileName, match, npoints);
	//
	// The mapping of the nearest profile point to each ghost-cell or interface location
	// will be done as needed, at application time.
	// This way, all of the necessary cell and position data should be valid.
    } // end this()

    double compute_distance(ref const(Vector3) my_pos, ref const(Vector3) other_pos)
    {
	double distance, other_r, my_r, dx, dy, dz, dr;
	switch (posMatch) {
	case "xyz-to-xyz": goto default;
	case "xyA-to-xyA":
	    dx = my_pos.x - other_pos.x;
	    dy = my_pos.y - other_pos.y;
	    distance = sqrt(dx^^2 + dy^^2);
	    break;
	case "AyA-to-AyA":
	    dy = my_pos.y - other_pos.y;
	    distance = fabs(dy);
	    break;
	case "xy-to-xR":
	    dx = my_pos.x - other_pos.x;
	    other_r = sqrt(other_pos.y^^2 + other_pos.z^^2);
	    my_r = sqrt(my_pos.y^^2 + my_pos.z^^2);
	    dr = my_r - other_r;
	    distance = sqrt(dx*dx + dr*dr);
	    break;
	case "Ay-to-AR":
	    other_r = sqrt(other_pos.y^^2 + other_pos.z^^2);
	    my_r = sqrt(my_pos.y^^2 + my_pos.z^^2);
	    dr = my_r - other_r;
	    distance = fabs(dr);
	    break;
	default:
	    // 3D, closest match on all position components.
	    dx = my_pos.x - other_pos.x;
	    dy = my_pos.y - other_pos.y;
	    dz = my_pos.z - other_pos.z;
	    distance = sqrt(dx*dx + dy*dy + dz*dz);
	    break; 
	}
	return distance;
    } // end compute_distance()
    
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

    FlowState get_flowstate(size_t my_id, ref const(Vector3) my_pos)
    {
	assert(fstate.length > 0, "FlowProfile is empty.");
	if (my_id in which_point) {
	    return fstate[which_point[my_id]];
	} else {
	    size_t ip = find_nearest_profile_point(my_pos);
	    which_point[my_id] = ip;
	    return fstate[ip];
	}
    } // end get_flowstate()

    void adjust_velocity(ref FlowState fs, ref const(Vector3) my_pos)
    {
	switch (posMatch) {
	case "xyz-to-xyz": break;
	case "xyA-to-xyA": break;
	case "AyA-to-AyA": break;
	case "xy-to-xR": goto case "Ay-to-AR";
	case "Ay-to-AR":
	    // We are assuming that the original 2D simulation had y>0.
	    double r = sqrt(my_pos.y^^2 + my_pos.z^^2);
	    double vel_yz = sqrt(fs.vel.y^^2 + fs.vel.z);
	    double vely_sign = (fs.vel.y < 0.0) ? -1.0 : 1.0;
	    fs.vel.refy = vely_sign * vel_yz * my_pos.y / r;
	    fs.vel.refz = vely_sign * vel_yz * my_pos.z / r;
	    break;
	default: break; // 3D, do nothing. 
	}
    }
} // end class FlowProfile
