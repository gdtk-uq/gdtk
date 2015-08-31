/**
 * fvcell.d
 * Finite-volume cell class for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcell;

import std.conv;
import std.string;
import std.array;
import std.format;
import std.stdio;
import std.math;
import geom;
import gas;
import kinetics;
import fvcore;
import flowstate;
import conservedquantities;
import fvvertex;
import fvinterface;
import globalconfig;

// The following two functions are used at compile time.
// They generate a string formula for averaging some quantity
// over the set of cell vertices. 
string avg_over_vtx_2D(string name)
{
    return "0.25 * (vtx[0]." ~ name ~ " + vtx[1]." ~ name ~
	" + vtx[2]." ~ name ~ " + vtx[3]." ~ name ~ ")";
}
string avg_over_vtx_3D(string name)
{
    return "0.125 * (vtx[0]." ~ name ~ " + vtx[1]." ~ name ~
	" + vtx[2]." ~ name ~ " + vtx[3]." ~ name ~ 
	" + vtx[4]." ~ name ~ " + vtx[5]." ~ name ~ 
	" + vtx[6]." ~ name ~ " + vtx[7]." ~ name ~ ")";
}

class FVCell {
public:
    uint id;  // allows us to work out where, in the block, the cell is
    bool fr_reactions_allowed; // if true, will call chemical_increment (also thermal_increment)
    double dt_chem; // acceptable time step for finite-rate chemistry
    double dt_therm; // acceptable time step for thermal relaxation
    bool in_turbulent_zone; // if true, we will keep the turbulence viscosity
    double base_qdot; // base-level of heat addition to cell, W/m**3
    // Geometry
    Vector3[] pos; // Centre x,y,z-coordinates for time-levels, m,m,m
    double[] volume; // Cell volume for time-levels (per unit depth or radian in 2D), m**3
    double[] areaxy; // (x,y)-plane area for time-levels, m**2
    double iLength; // length in the i-index direction
    double jLength; // length in the j-index direction
    double kLength; // length in the k-index direction
    double L_min;   // minimum length scale for cell
    double distance_to_nearest_wall; // for turbulence model correction.
    double half_cell_width_at_wall;  // ditto
    FVCell cell_at_nearest_wall;   // ditto
    // Connections
    FVInterface[] iface;  // references to defining interfaces of cell
    FVVertex[] vtx;  // references to vertices for quad (2D) and hexahedral (3D) cells
    // Flow
    FlowState fs; // Flow properties
    ConservedQuantities[] U;  // Conserved flow quantities for the update stages.
    ConservedQuantities[] dUdt; // Time derivatives for the update stages.
    ConservedQuantities Q; // source (or production) terms
    // Terms for loose-coupling of radiation.
    double Q_rad_org;
    double f_rad_org;
    double Q_rE_rad; // Rate of energy addition to cell via radiation.
    double Q_rE_rad_save; // Presently, the radiation source term is calculated
                          // at the first update stage.  We need to retain that
                          // value for all of the update stages.
    // Data for computing residuals.
    double rho_at_start_of_step, rE_at_start_of_step;
    // [TODO] implicit variables

private:
    LocalConfig myConfig;

public:
    this(LocalConfig myConfig, int id_init=0)
    {
	this.myConfig = myConfig;
	auto gmodel = myConfig.gmodel;
	int n_species = gmodel.n_species;
	int n_modes = gmodel.n_modes;
	id = id_init;
	pos.length = n_time_levels;
	volume.length = n_time_levels;
	areaxy.length = n_time_levels;
	fs = new FlowState(gmodel, 100.0e3, [300.0,], Vector3(0.0,0.0,0.0));
	foreach(i; 0 .. n_time_levels) {
	    U ~= new ConservedQuantities(n_species, n_modes);
	    dUdt ~= new ConservedQuantities(n_species, n_modes);
	}
	Q = new ConservedQuantities(n_species, n_modes);
    }

    @nogc
    void copy_values_from(FVCell other, int type_of_copy)
    {
	switch ( type_of_copy ) {
	case CopyDataOption.minimal_flow:
	    fs.copy_values_from(other.fs);
	    break;
	case CopyDataOption.all_flow:
	    fs.copy_values_from(other.fs);
	    Q.copy_values_from(other.Q);
	    foreach(i; 0 .. n_time_levels) {
		U[i].copy_values_from(other.U[i]);
		dUdt[i].copy_values_from(other.dUdt[i]);
	    }
	    break;
	case CopyDataOption.grid:
	    foreach(i; 0 .. n_time_levels) {
		pos[i].refx = other.pos[i].x;
		pos[i].refy = other.pos[i].y;
		pos[i].refz = other.pos[i].z;
		volume[i] = other.volume[i];
		areaxy[i] = other.areaxy[i];
	    }
	    iLength = other.iLength;
	    jLength = other.jLength;
	    kLength = other.kLength;
	    L_min = other.L_min;
	    break;
	case CopyDataOption.cell_lengths_only:
	    iLength = other.iLength;
	    jLength = other.jLength;
	    kLength = other.kLength;
	    L_min = other.L_min;
	    break;
	case CopyDataOption.all: 
	default:
	    // [TODO] really need to think about what needs to be copied...
	    id = other.id;
	    myConfig = other.myConfig;
	    foreach(i; 0 .. n_time_levels) {
		pos[i].refx = other.pos[i].x;
		pos[i].refy = other.pos[i].y;
		pos[i].refz = other.pos[i].z;
		volume[i] = other.volume[i];
		areaxy[i] = other.areaxy[i];
	    }
	    iLength = other.iLength;
	    jLength = other.jLength;
	    kLength = other.kLength;
	    L_min = other.L_min;
	    fs.copy_values_from(other.fs);
	    Q.copy_values_from(other.Q);
	    foreach(i; 0 .. n_time_levels) {
		U[i].copy_values_from(other.U[i]);
		dUdt[i].copy_values_from(other.dUdt[i]);
	    }
	} // end switch
    }

    @nogc
    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
	pos[to_level] = pos[from_level];
	volume[to_level] = volume[from_level];
	areaxy[to_level] = areaxy[from_level];
	// When working over all cells in a block, the following copies
	// will no doubt do some doubled-up work, but it should be otherwise benign.
	foreach(ref face; iface) {
	    if (face) face.copy_grid_level_to_level(from_level, to_level);
	}
	foreach(ref v; vtx) {
	    if (v) v.copy_grid_level_to_level(from_level, to_level);
	}
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FVCell(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", pos=" ~ to!string(pos);
	repr ~= ", volume=" ~ to!string(volume);
	repr ~= ", areaxy=" ~ to!string(areaxy);
	repr ~= ", iLength=" ~ to!string(iLength);
	repr ~= ", jLength=" ~ to!string(jLength);
	repr ~= ", kLength=" ~ to!string(kLength);
	repr ~= ", L_min=" ~ to!string(L_min);
	repr ~= ", dt_chem=" ~ to!string(dt_chem);
	repr ~= ", dt_therm=" ~ to!string(dt_therm);
	repr ~= ", in_turbulent_zone=" ~ to!string(in_turbulent_zone);
	repr ~= ", fr_reactions_allowed=" ~ to!string(fr_reactions_allowed);
	repr ~= ", fs=" ~ to!string(fs);
	repr ~= ", U=" ~ to!string(U);
	repr ~= ", dUdt=" ~ to!string(dUdt);
	repr ~= ")";
	return to!string(repr);
    }

    bool point_is_inside(in Vector3 p, int dimensions, int gtl) const
    // Returns true if the point p is inside or on the cell surface.
    {
	if ( dimensions == 2 ) {
	    // In 2 dimensions,
	    // we split the x,y-plane into half-planes and check which side p is on.
	    double xA = vtx[1].pos[gtl].x; double yA = vtx[1].pos[gtl].y;
	    double xB = vtx[1].pos[gtl].x; double yB = vtx[2].pos[gtl].y;
	    double xC = vtx[3].pos[gtl].x; double yC = vtx[3].pos[gtl].y;
	    double xD = vtx[0].pos[gtl].x; double yD = vtx[0].pos[gtl].y;
	    // Now, check to see if the specified point is on the
	    // left of (or on) each bloundary line AB, BC, CD and DA.
	    if ((p.x - xB) * (yA - yB) >= (p.y - yB) * (xA - xB) &&
		(p.x - xC) * (yB - yC) >= (p.y - yC) * (xB - xC) &&
		(p.x - xD) * (yC - yD) >= (p.y - yD) * (xC - xD) &&
		(p.x - xA) * (yD - yA) >= (p.y - yA) * (xD - xA)) {
		return true;
	    } else {
		return false;
	    }
	} else {
	    // In 3 dimensions,
	    // the test consists of dividing the 6 cell faces into triangular facets
	    // with outwardly-facing normals and then computing the volumes of the
	    // tetrahedra formed by these facets and the sample point p.
	    // If any of the tetrahedra volumes are positive
	    // (i.e. p is on the positive side of a facet) and we assume a convex cell,
	    // it means that the point is outside the cell and we may say so
	    // without further testing.

	    // North
	    if ( tetrahedron_volume(vtx[2].pos[gtl], vtx[3].pos[gtl], vtx[7].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[7].pos[gtl], vtx[6].pos[gtl], vtx[2].pos[gtl], p) > 0.0 ) return false;
	    // East
	    if ( tetrahedron_volume(vtx[1].pos[gtl], vtx[2].pos[gtl], vtx[6].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[6].pos[gtl], vtx[5].pos[gtl], vtx[1].pos[gtl], p) > 0.0 ) return false;
	    // South
	    if ( tetrahedron_volume(vtx[0].pos[gtl], vtx[1].pos[gtl], vtx[5].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[5].pos[gtl], vtx[4].pos[gtl], vtx[0].pos[gtl], p) > 0.0 ) return false;
	    // West
	    if ( tetrahedron_volume(vtx[3].pos[gtl], vtx[0].pos[gtl], vtx[4].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[4].pos[gtl], vtx[7].pos[gtl], vtx[3].pos[gtl], p) > 0.0 ) return false;
	    // Bottom
	    if ( tetrahedron_volume(vtx[1].pos[gtl], vtx[0].pos[gtl], vtx[3].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[3].pos[gtl], vtx[2].pos[gtl], vtx[1].pos[gtl], p) > 0.0 ) return false;
	    // Top
	    if ( tetrahedron_volume(vtx[4].pos[gtl], vtx[5].pos[gtl], vtx[6].pos[gtl], p) > 0.0 ) return false;
	    if ( tetrahedron_volume(vtx[6].pos[gtl], vtx[7].pos[gtl], vtx[4].pos[gtl], p) > 0.0 ) return false;
	    // If we arrive here, we haven't determined that the point is outside...
	    return true;
	} // end dimensions != 2
    } // end point_is_inside()

    @nogc
    void copy_values_to_buffer(ref double[] buf, int type_of_copy, int gtl) const
    {
	assert(false, "[TODO] FVCell.copy_values_to_buffer() not yet implemented");
    }

    @nogc
    void copy_values_from_buffer(in double buf, int type_of_copy, int gtl) 
    {
	assert(false, "[TODO] FVCell.copy_values_from_buffer() not yet implemented");
    }

    void replace_flow_data_with_average(in FVCell[] others) 
    {
	auto gmodel = myConfig.gmodel;
	size_t n = others.length;
	if (n == 0) throw new Error("Need to average from a nonempty array.");
	FlowState[] fsList;
	// We need to be honest and not to fiddle with the other gas states.
	foreach(other; others) {
	    if ( this is other ) throw new Error("Must not include destination in source list.");
	    fsList ~= cast(FlowState)other.fs;
	}
	fs.copy_average_values_from(fsList, gmodel);
	// Accumulate from a clean slate and then divide.
	Q_rE_rad = 0.0;
	foreach(other; others) {
	    Q_rE_rad += other.Q_rE_rad;
	}
	Q_rE_rad /= n;
    }

    void scan_values_from_string(string buffer)
    // Note that the position data is read into grid_time_level 0.
    {
	auto gm = myConfig.gmodel;
	auto items = split(buffer);
	pos[0].refx = to!double(items.front); items.popFront();
	pos[0].refy = to!double(items.front); items.popFront();
	pos[0].refz = to!double(items.front); items.popFront();
	volume[0] = to!double(items.front); items.popFront();
	fs.gas.rho = to!double(items.front); items.popFront();
	fs.vel.refx = to!double(items.front); items.popFront();
	fs.vel.refy = to!double(items.front); items.popFront();
	fs.vel.refz = to!double(items.front); items.popFront();
	if ( myConfig.MHD ) {
	    fs.B.refx = to!double(items.front); items.popFront();
	    fs.B.refy = to!double(items.front); items.popFront();
	    fs.B.refz = to!double(items.front); items.popFront();
	}
	fs.gas.p = to!double(items.front); items.popFront();
	fs.gas.a = to!double(items.front); items.popFront();
	fs.gas.mu = to!double(items.front); items.popFront();
	foreach(i; 0 .. gm.n_modes) {
	    fs.gas.k[i] = to!double(items.front); items.popFront();
	}
	fs.mu_t = to!double(items.front); items.popFront();
	fs.k_t = to!double(items.front); items.popFront();
	fs.S = to!int(items.front); items.popFront();
	if ( myConfig.radiation ) {
	    Q_rad_org = to!double(items.front); items.popFront();
	    f_rad_org = to!double(items.front); items.popFront();
	    Q_rE_rad = to!double(items.front); items.popFront();
	} else {
	    Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
	}
	fs.tke = to!double(items.front); items.popFront();
	fs.omega = to!double(items.front); items.popFront();
	foreach(i; 0 .. gm.n_species) {
	    fs.gas.massf[i] = to!double(items.front); items.popFront();
	}
	if ( gm.n_species > 1 ) {
	    dt_chem = to!double(items.front); items.popFront();
	}
	foreach(i; 0 .. gm.n_modes) {
	    fs.gas.e[i] = to!double(items.front); items.popFront();
	    fs.gas.T[i] = to!double(items.front); items.popFront();
	}
	if ( gm.n_modes > 1 ) {
	    dt_therm = to!double(items.front); items.popFront(); 
	}
    }

    string write_values_to_string() const
    {
	// Should match cell_data_as_string() in flowstate.d
	auto writer = appender!string();
	formattedWrite(writer, "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e",
		       pos[0].x, pos[0].y, pos[0].z, volume[0], fs.gas.rho,
		       fs.vel.x, fs.vel.y, fs.vel.z);
	if ( myConfig.MHD ) 
	    formattedWrite(writer, " %.12e %.12e %.12e", fs.B.x, fs.B.y, fs.B.z); 
	formattedWrite(writer, " %.12e %.12e %.12e", fs.gas.p, fs.gas.a, fs.gas.mu);
	foreach(i; 0 .. fs.gas.k.length) formattedWrite(writer, " %.12e", fs.gas.k[i]); 
	formattedWrite(writer, " %.12e %.12e %d", fs.mu_t, fs.k_t, fs.S);
	if ( myConfig.radiation ) 
	    formattedWrite(writer, " %.12e %.12e %.12e", Q_rad_org, f_rad_org, Q_rE_rad); 
	formattedWrite(writer, " %.12e %.12e", fs.tke, fs.omega);
	foreach(i; 0 .. fs.gas.massf.length) formattedWrite(writer, " %.12e", fs.gas.massf[i]); 
	if ( fs.gas.massf.length > 1 ) formattedWrite(writer, " %.12e", dt_chem); 
	foreach(i; 0 .. fs.gas.e.length)
	    formattedWrite(writer, " %.12e %.12e", fs.gas.e[i], fs.gas.T[i]); 
	if ( fs.gas.e.length > 1 ) formattedWrite(writer, " %.12e", dt_therm);
	return writer.data;
    }

    @nogc
    void encode_conserved(int gtl, int ftl, double omegaz)
    // gtl = grid time level
    // ftl = flow time level
    {
	ConservedQuantities myU = U[ftl];
	bool with_k_omega = (myConfig.turbulence_model == TurbulenceModel.k_omega);

	myU.mass = fs.gas.rho;
	// X-, Y- and Z-momentum per unit volume.
	myU.momentum.refx = fs.gas.rho * fs.vel.x;
	myU.momentum.refy = fs.gas.rho * fs.vel.y;
	myU.momentum.refz = fs.gas.rho * fs.vel.z;
	// Magnetic field
	myU.B.refx = fs.B.x;
	myU.B.refy = fs.B.y;
	myU.B.refz = fs.B.z;
	// Total Energy / unit volume = density
	// (specific internal energy + kinetic energy/unit mass).
	double e = 0.0; foreach(elem; fs.gas.e) e += elem;
	double ke = 0.5 * (fs.vel.x * fs.vel.x + fs.vel.y * fs.vel.y + fs.vel.z * fs.vel.z);
	if (with_k_omega) {
	    myU.tke = fs.gas.rho * fs.tke;
	    myU.omega = fs.gas.rho * fs.omega;
	    myU.total_energy = fs.gas.rho * (e + ke + fs.tke);
	} else {
	    myU.tke = 0.0;
	    myU.omega = fs.gas.rho * 1.0;
	    myU.total_energy = fs.gas.rho * (e + ke);
	}
	if (myConfig.MHD) {
	    double me = 0.5 * (fs.B.x * fs.B.x + fs.B.y * fs.B.y + fs.B.z * fs.B.z);
	    myU.total_energy += me;
	}
	// Species densities: mass of species is per unit volume.
	foreach(isp; 0 .. myU.massf.length) {
	    myU.massf[isp] = fs.gas.rho * fs.gas.massf[isp];
	}
	// Individual energies: energy in mode per unit volume
	foreach(imode; 0 .. myU.energies.length) {
	    myU.energies[imode] = fs.gas.rho * fs.gas.e[imode];
	}
    
	if (omegaz != 0.0) {
	    // Rotating frame.
	    // Finally, we adjust the total energy to make rothalpy.
	    // We do this last because the gas models don't know anything
	    // about rotating frames and we don't want to mess their
	    // energy calculations around.
	    double rho = fs.gas.rho;
	    double x = pos[gtl].x;
	    double y = pos[gtl].y;
	    double rsq = x*x + y*y;
	    // The conserved quantity is rothalpy. I = E - (u**2)/2
	    // where rotating frame velocity  u = omegaz * r.
	    myU.total_energy -= rho * 0.5 * omegaz * omegaz * rsq;
	}
	assert(U[ftl].mass > 0.0, "invalid density in conserved quantities vector" ~
	       " when leaving FVCell.encode_conserved().");
	return;
    } // end encode_conserved()

    void decode_conserved(int gtl, int ftl, double omegaz) 
    {
	auto gmodel = myConfig.gmodel;
	ConservedQuantities myU = U[ftl];
	bool with_k_omega = (myConfig.turbulence_model == TurbulenceModel.k_omega);
	double e, ke, dinv, rE, me;
	// Mass / unit volume = Density
	if (!(myU.mass > 0.0)) {
	    writeln("FVCell.decode_conserved(): Density invalid in conserved quantities.");
	    writeln("  id= ", id, " x= ", pos[gtl].x, " y= ", pos[gtl].y, " z= ", pos[gtl].z);
	    writeln("  gas= ", fs.gas);
	    writeln("  U= ", myU);
	    throw new Error("Bad cell, give up.");
	}
	double rho = myU.mass;
	fs.gas.rho = rho; // This is limited to nonnegative and finite values.
	dinv = 1.0 / rho;
	if (omegaz != 0.0) {
	    // Rotating frame.
	    // The conserved quantity is rothalpy so we need to convert
	    // back to enthalpy to do the rest of the decode.
	    double x = pos[gtl].x;
	    double y = pos[gtl].y;
	    double rsq = x*x + y*y;
	    rE = myU.total_energy + rho * 0.5 * omegaz * omegaz * rsq;
	} else {
	    // Non-rotating frame.
	    rE = myU.total_energy;
	}
	// Velocities from momenta.
	fs.vel.refx = myU.momentum.x * dinv;
	fs.vel.refy = myU.momentum.y * dinv;
	fs.vel.refz = myU.momentum.z * dinv;
	// Magnetic field
	fs.B.refx = myU.B.x;
	fs.B.refy = myU.B.y;
	fs.B.refz = myU.B.z;
	// Specific internal energy from total energy per unit volume.
	ke = 0.5 * (fs.vel.x * fs.vel.x + fs.vel.y * fs.vel.y + fs.vel.z * fs.vel.z);
	if ( myConfig.MHD ) {
	    me = 0.5*(fs.B.x*fs.B.x + fs.B.y*fs.B.y + fs.B.z*fs.B.z);
	} else {
	    me = 0.0;
	}
	if (with_k_omega) {
	    fs.tke = myU.tke * dinv;
	    fs.omega = myU.omega * dinv;
	    e = (rE - myU.tke - me) * dinv - ke;
	} else {
	    fs.tke = 0.0;
	    fs.omega = 1.0;
	    e = (rE - me) * dinv - ke;
	}
	foreach(isp; 0 .. gmodel.n_species) fs.gas.massf[isp] = myU.massf[isp] * dinv; 
	if (gmodel.n_species > 1) scale_mass_fractions(fs.gas.massf);
	foreach(imode; 0 .. gmodel.n_modes) fs.gas.e[imode] = myU.energies[imode] * dinv; 
	// We can recompute e[0] from total energy and component
	// modes NOT in translation.
	if (gmodel.n_modes > 1) {
	    double e_tmp = 0.0;
	    foreach(imode; 1 .. gmodel.n_modes) e_tmp += fs.gas.e[imode];
	    fs.gas.e[0] = e - e_tmp;
	} else {
	    fs.gas.e[0] = e;
	}
	// Fill out the other variables: P, T, a, and viscous transport coefficients.
	gmodel.update_thermo_from_rhoe(fs.gas);
	gmodel.update_sound_speed(fs.gas);
	if (myConfig.viscous) gmodel.update_trans_coeffs(fs.gas);
	// if (myConfig.diffusion) gmodel.update_diff_coeffs(fs.gas);
	return;
    } // end decode_conserved()

    bool check_flow_data() const
    {
	bool is_data_valid = fs.gas.check_values(true);
	const double MAXVEL = 30000.0;
	if (fabs(fs.vel.x) > MAXVEL || fabs(fs.vel.y) > MAXVEL || fabs(fs.vel.z) > MAXVEL) {
	    writeln("Velocity too high ", fs.vel.x, " ", fs.vel.y, " ", fs.vel.z);
	    is_data_valid = false;
	}
	if ( !isFinite(fs.tke) ) {
	    writeln("Turbulence KE invalid number ", fs.tke);
	    is_data_valid = false;
	}
	if ( fs.tke < 0.0 ) {
	    writeln("Turbulence KE negative ", fs.tke);
	    is_data_valid = false;
	}
	if ( !isFinite(fs.omega) ) {
	    writeln("Turbulence frequency invalid number ", fs.omega);
	    is_data_valid = false;
	}
	if ( fs.omega <= 0.0 ) {
	    writeln("Turbulence frequency nonpositive ", fs.omega);
	    is_data_valid = false;
	}
	if ( !is_data_valid ) {
	    writeln("cell pos=", pos[0]);
	    writeln(fs);
	    writeln("----------------------------------------------------------");
	}
	return is_data_valid;
    } // end check_flow_data()

    @nogc
    void time_derivatives(int gtl, int ftl, bool with_k_omega) 
    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    // gtl : (grid-time-level) flow derivatives are evaluated at this grid level
    // ftl : (flow-time-level) specifies where computed derivatives are to be stored.
    //       0: Start of stage-1 update.
    //       1: End of stage-1.
    //       2: End of stage-2.
    {
	FVInterface IFn = iface[Face.north];
	FVInterface IFe = iface[Face.east];
	FVInterface IFs = iface[Face.south];
	FVInterface IFw = iface[Face.west];
	FVInterface IFt, IFb;
	auto dimensions = myConfig.dimensions;
	if (dimensions == 3) {
	    IFt = iface[Face.top];
	    IFb = iface[Face.bottom];
	}
	// Cell volume (inverted).
	double vol_inv = 1.0 / volume[gtl];
	double integral;
    
	// Time-derivative for Mass/unit volume.
	// Note that the unit normals for the interfaces are oriented
	// such that the unit normals for the east, north and top faces
	// are outward and the unit normals for the south, west and
	// bottom faces are inward.
	integral = -IFe.F.mass * IFe.area[gtl] - IFn.F.mass * IFn.area[gtl]
	    + IFw.F.mass * IFw.area[gtl] + IFs.F.mass * IFs.area[gtl];
	if (dimensions == 3)
	    integral += IFb.F.mass * IFb.area[gtl] - IFt.F.mass * IFt.area[gtl];
	dUdt[ftl].mass = vol_inv * integral + Q.mass;

	// Time-derivative for X-Momentum/unit volume.
	integral = -IFe.F.momentum.x * IFe.area[gtl] - IFn.F.momentum.x * IFn.area[gtl]
	    + IFw.F.momentum.x * IFw.area[gtl] + IFs.F.momentum.x * IFs.area[gtl];
	if (dimensions == 3)
	    integral += IFb.F.momentum.x * IFb.area[gtl] - IFt.F.momentum.x * IFt.area[gtl];
	dUdt[ftl].momentum.refx = vol_inv * integral + Q.momentum.x;
	// Time-derivative for Y-Momentum/unit volume.
	integral = -IFe.F.momentum.y * IFe.area[gtl] - IFn.F.momentum.y * IFn.area[gtl]
	    + IFw.F.momentum.y * IFw.area[gtl] + IFs.F.momentum.y * IFs.area[gtl];
	if (dimensions == 3)
	    integral += IFb.F.momentum.y * IFb.area[gtl] - IFt.F.momentum.y * IFt.area[gtl];
	dUdt[ftl].momentum.refy = vol_inv * integral + Q.momentum.y;
    
	// we require the z-momentum for MHD even in 2D
	if ((dimensions == 3) || ( myConfig.MHD )) {
	    // Time-derivative for Z-Momentum/unit volume.
	    integral = -IFe.F.momentum.z * IFe.area[gtl] - IFn.F.momentum.z * IFn.area[gtl]
		+ IFw.F.momentum.z * IFw.area[gtl] + IFs.F.momentum.z * IFs.area[gtl];
	}
	if (dimensions == 3) {
	    integral += IFb.F.momentum.z * IFb.area[gtl] - IFt.F.momentum.z * IFt.area[gtl];
	}
	if ((dimensions == 3) || ( myConfig.MHD )) {
	    dUdt[ftl].momentum.refz = vol_inv * integral + Q.momentum.z;
	} else {
	    dUdt[ftl].momentum.refz = 0.0;
	}
    
	if ( myConfig.MHD ) {
	    // Time-derivative for X-Magnetic Field/unit volume.
	    integral = -IFe.F.B.x * IFe.area[gtl] - IFn.F.B.x * IFn.area[gtl]
		+ IFw.F.B.x * IFw.area[gtl] + IFs.F.B.x * IFs.area[gtl];
	    if (dimensions == 3)
		integral += IFb.F.B.x * IFb.area[gtl] - IFt.F.B.x * IFt.area[gtl];
	    dUdt[ftl].B.refx = vol_inv * integral + Q.B.x;
	    // Time-derivative for Y-Magnetic Field/unit volume.
	    integral = -IFe.F.B.y * IFe.area[gtl] - IFn.F.B.y * IFn.area[gtl]
		+ IFw.F.B.y * IFw.area[gtl] + IFs.F.B.y * IFs.area[gtl];
	    if (dimensions == 3)
		integral += IFb.F.B.y * IFb.area[gtl] - IFt.F.B.y * IFt.area[gtl];
	    dUdt[ftl].B.refy = vol_inv * integral + Q.B.y;
	    // Time-derivative for Z-Magnetic Field/unit volume.
	    integral = -IFe.F.B.z * IFe.area[gtl] - IFn.F.B.z * IFn.area[gtl]
		+ IFw.F.B.z * IFw.area[gtl] + IFs.F.B.z * IFs.area[gtl];
	    if (dimensions == 3) {
		integral += IFb.F.B.z * IFb.area[gtl] - IFt.F.B.z * IFt.area[gtl];
	    }
	    dUdt[ftl].B.refz = vol_inv * integral + Q.B.z;
	} else {
	    dUdt[ftl].B.refx = 0.0;
	    dUdt[ftl].B.refy = 0.0;
	    dUdt[ftl].B.refz = 0.0;
	}

	// Time-derivative for Total Energy/unit volume.
	integral = -IFe.F.total_energy * IFe.area[gtl] - IFn.F.total_energy * IFn.area[gtl]
	    + IFw.F.total_energy * IFw.area[gtl] + IFs.F.total_energy * IFs.area[gtl];
	if ( dimensions == 3 )
	    integral += IFb.F.total_energy * IFb.area[gtl] - IFt.F.total_energy * IFt.area[gtl];
	dUdt[ftl].total_energy = vol_inv * integral + Q.total_energy;
    
	if (with_k_omega) {
	    integral = -IFe.F.tke * IFe.area[gtl] - IFn.F.tke * IFn.area[gtl]
		+ IFw.F.tke * IFw.area[gtl] + IFs.F.tke * IFs.area[gtl];
	    if (dimensions == 3)
		integral += IFb.F.tke * IFb.area[gtl] - IFt.F.tke * IFt.area[gtl];
	    dUdt[ftl].tke = vol_inv * integral + Q.tke;
	
	    integral = -IFe.F.omega * IFe.area[gtl] - IFn.F.omega * IFn.area[gtl]
		+ IFw.F.omega * IFw.area[gtl] + IFs.F.omega * IFs.area[gtl];
	    if (dimensions == 3)
		integral += IFb.F.omega * IFb.area[gtl] - IFt.F.omega * IFt.area[gtl];
	    dUdt[ftl].omega = vol_inv * integral + Q.omega;
	} else {
	    dUdt[ftl].tke = 0.0;
	    dUdt[ftl].omega = 0.0;
	}
	// Time-derivative for individual species.
	// The conserved quantity is the mass per unit
	// volume of species isp and
	// the fluxes are mass/unit-time/unit-area.
	// Units of DmassfDt are 1/sec.
	foreach(isp; 0 .. IFe.F.massf.length) {
	    integral =
		-IFe.F.massf[isp] * IFe.area[gtl]
		- IFn.F.massf[isp] * IFn.area[gtl]
		+ IFw.F.massf[isp] * IFw.area[gtl]
		+ IFs.F.massf[isp] * IFs.area[gtl];
	    if ( dimensions == 3 )
		integral += IFb.F.massf[isp] * IFb.area[gtl] - IFt.F.massf[isp] * IFt.area[gtl];
	    dUdt[ftl].massf[isp] = vol_inv * integral + Q.massf[isp];
	}
	// Individual energies.
	// We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
	// Instead we get this from the conservation of total energy
	foreach(imode; 1 .. IFe.F.energies.length) {
	    integral =
		-IFe.F.energies[imode] * IFe.area[gtl]
		- IFn.F.energies[imode] * IFn.area[gtl]
		+ IFw.F.energies[imode] * IFw.area[gtl]
		+ IFs.F.energies[imode] * IFs.area[gtl];
	    if (dimensions == 3)
		integral += IFb.F.energies[imode] * IFb.area[gtl] - IFt.F.energies[imode] * IFt.area[gtl];
	    dUdt[ftl].energies[imode] = vol_inv * integral + Q.energies[imode];
	}
    } // end time_derivatives()

    @nogc
    void stage_1_update_for_flow_on_fixed_grid(double dt, bool force_euler, bool with_k_omega) 
    {
	ConservedQuantities dUdt0 = dUdt[0];
	ConservedQuantities U0 = U[0];
	ConservedQuantities U1 = U[1];
	double gamma_1 = 1.0; // for normal Predictor-Corrector or Euler update.
	// In some parts of the code (viscous updates, k-omega updates)
	// we use this function as an Euler update even when the main
	// gasdynamic_update_scheme is of higher order.
	// force_euler is set true for these situations.
	if (!force_euler) {
	    final switch (myConfig.gasdynamic_update_scheme) {
	    case GasdynamicUpdate.euler:
	    case GasdynamicUpdate.pc: gamma_1 = 1.0; break;
	    case GasdynamicUpdate.midpoint: gamma_1 = 0.5; break;
	    case GasdynamicUpdate.classic_rk3: gamma_1 = 0.5; break;
	    case GasdynamicUpdate.tvd_rk3: gamma_1 = 1.0; break;
	    case GasdynamicUpdate.denman_rk3: gamma_1 = 8.0/15.0; break;
	    }
	}
	U1.mass = U0.mass + dt * gamma_1 * dUdt0.mass;
	// Side note: 
	// It would be convenient (codewise) for the updates of these Vector3 quantities to
	// be done with the Vector3 arithmetic operators but I suspect that the implementation
	// of those oerators is such that a whole lot of Vector3 temporaries would be created.
	U1.momentum.refx = U0.momentum.x + dt * gamma_1 * dUdt0.momentum.x;
	U1.momentum.refy = U0.momentum.y + dt * gamma_1 * dUdt0.momentum.y;
	U1.momentum.refz = U0.momentum.z + dt * gamma_1 * dUdt0.momentum.z;
	if (myConfig.MHD) {
	    // Magnetic field
	    U1.B.refx = U0.B.x + dt * gamma_1 * dUdt0.B.x;
	    U1.B.refy = U0.B.y + dt * gamma_1 * dUdt0.B.y;
	    U1.B.refz = U0.B.z + dt * gamma_1 * dUdt0.B.z;
	}
	U1.total_energy = U0.total_energy + dt * gamma_1 * dUdt0.total_energy;
	if (with_k_omega) {
	    U1.tke = U0.tke + dt * gamma_1 * dUdt0.tke;
	    U1.tke = fmax(U1.tke, 0.0);
	    U1.omega = U0.omega + dt * gamma_1 * dUdt0.omega;
	    U1.omega = fmax(U1.omega, U0.mass);
	    // ...assuming a minimum value of 1.0 for omega
	    // It may occur (near steps in the wall) that a large flux of romega
	    // through one of the cell interfaces causes romega within the cell
	    // to drop rapidly.
	    // The large values of omega come from Menter's near-wall correction that may be
	    // applied outside the control of this finite-volume core code.
	    // These large values of omega will be convected along the wall and,
	    // if they are convected past a corner with a strong expansion,
	    // there will be an unreasonably-large flux out of the cell.
	} else {
	    U1.tke = U0.tke;
	    U1.omega = U0.omega;
	}
	foreach(isp; 0 .. U1.massf.length) {
	    U1.massf[isp] = U0.massf[isp] + dt * gamma_1 * dUdt0.massf[isp];
	}
	// We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
	// Instead we get this from the conservation of total energy
	foreach(imode; 1 .. U1.energies.length) {
	    U1.energies[imode] = U0.energies[imode] + dt * gamma_1 * dUdt0.energies[imode];
	}
	return;
    } // end stage_1_update_for_flow_on_fixed_grid()

    @nogc
    void stage_2_update_for_flow_on_fixed_grid(double dt, bool with_k_omega) 
    {
	ConservedQuantities dUdt0 = dUdt[0];
	ConservedQuantities dUdt1 = dUdt[1];
	ConservedQuantities U_old = U[0];
	if (myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) U_old = U[1];
	ConservedQuantities U2 = U[2];
	double gamma_1 = 0.5; // Presume predictor-corrector.
	double gamma_2 = 0.5;
	final switch (myConfig.gasdynamic_update_scheme) {
	case GasdynamicUpdate.euler: assert(false, "Euler update has no second stage.");
	case GasdynamicUpdate.pc: gamma_1 = 0.5, gamma_2 = 0.5; break;
	case GasdynamicUpdate.midpoint: gamma_1 = 0.0; gamma_2 = 1.0; break;
	case GasdynamicUpdate.classic_rk3: gamma_1 = -1.0; gamma_2 = 2.0; break;
	case GasdynamicUpdate.tvd_rk3: gamma_1 = 0.25; gamma_2 = 0.25; break;
	case GasdynamicUpdate.denman_rk3: gamma_1 = -17.0/60.0; gamma_2 = 5.0/12.0; break;
	}
	U2.mass = U_old.mass + dt * (gamma_1 * dUdt0.mass + gamma_2 * dUdt1.mass);
	U2.momentum.refx = U_old.momentum.x + dt * (gamma_1 * dUdt0.momentum.x + gamma_2 * dUdt1.momentum.x);
	U2.momentum.refy = U_old.momentum.y + dt * (gamma_1 * dUdt0.momentum.y + gamma_2 * dUdt1.momentum.y);
	U2.momentum.refz = U_old.momentum.z + dt * (gamma_1 * dUdt0.momentum.z + gamma_2 * dUdt1.momentum.z);
	if (myConfig.MHD) {
	    // Magnetic field
	    U2.B.refx = U_old.B.x + dt * (gamma_1 * dUdt0.B.x + gamma_2 * dUdt1.B.x);
	    U2.B.refy = U_old.B.y + dt * (gamma_1 * dUdt0.B.y + gamma_2 * dUdt1.B.y);
	    U2.B.refz = U_old.B.z + dt * (gamma_1 * dUdt0.B.z + gamma_2 * dUdt1.B.z);
	}
	U2.total_energy = U_old.total_energy + 
	    dt * (gamma_1 * dUdt0.total_energy + gamma_2 * dUdt1.total_energy);
	if ( with_k_omega ) {
	    U2.tke = U_old.tke + dt * (gamma_1 * dUdt0.tke + gamma_2 * dUdt1.tke);
	    U2.tke = fmax(U2.tke, 0.0);
	    U2.omega = U_old.omega + dt * (gamma_1 * dUdt0.omega + gamma_2 * dUdt1.omega);
	    U2.omega = fmax(U2.omega, U_old.mass);
	} else {
	    U2.tke = U_old.tke;
	    U2.omega = U_old.omega;
	}
	foreach(isp; 0 .. U2.massf.length) {
	    U2.massf[isp] = U_old.massf[isp] + dt * (gamma_1 * dUdt0.massf[isp] + gamma_2 * dUdt1.massf[isp]);
	}
	// We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
	// Instead we get this from the conservation of total energy
	foreach(imode; 1 .. U2.energies.length) {
	    U2.energies[imode] = U_old.energies[imode] + 
		dt * (gamma_1 * dUdt0.energies[imode] + gamma_2 * dUdt1.energies[imode]);
	}
	return;
    } // end stage_2_update_for_flow_on_fixed_grid()

    @nogc
    void stage_3_update_for_flow_on_fixed_grid(double dt, bool with_k_omega) 
    {
	ConservedQuantities dUdt0 = dUdt[0];
	ConservedQuantities dUdt1 = dUdt[1];
	ConservedQuantities dUdt2 = dUdt[2];
	ConservedQuantities U_old = U[0];
	if (myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) U_old = U[2];
	ConservedQuantities U3 = U[3];
	double gamma_1 = 1.0/6.0; // presume TVD_RK3 scheme.
	double gamma_2 = 1.0/6.0;
	double gamma_3 = 4.0/6.0;
	final switch (myConfig.gasdynamic_update_scheme) {
	case GasdynamicUpdate.euler:
	case GasdynamicUpdate.pc:
	case GasdynamicUpdate.midpoint:
	    assert(false, "Euler PC and midpoint updates have no second stage.");
	case GasdynamicUpdate.classic_rk3: gamma_1 = 1.0/6.0; gamma_2 = 4.0/6.0; gamma_3 = 1.0/6.0; break;
	case GasdynamicUpdate.tvd_rk3: gamma_1 = 1.0/6.0; gamma_2 = 1.0/6.0; gamma_3 = 4.0/6.0; break;
	    // FIX-ME: Check that we have Andrew Denman's scheme ported correctly.
	case GasdynamicUpdate.denman_rk3: gamma_1 = 0.0; gamma_2 = -5.0/12.0; gamma_3 = 3.0/4.0; break;
	}
	U3.mass = U_old.mass + dt * (gamma_1*dUdt0.mass + gamma_2*dUdt1.mass + gamma_3*dUdt2.mass);
	U3.momentum.refx = U_old.momentum.x +
	    dt * (gamma_1*dUdt0.momentum.x + gamma_2*dUdt1.momentum.x + gamma_3*dUdt2.momentum.x);
	U3.momentum.refy = U_old.momentum.y +
	    dt * (gamma_1*dUdt0.momentum.y + gamma_2*dUdt1.momentum.y + gamma_3*dUdt2.momentum.y);
	U3.momentum.refz = U_old.momentum.z + 
	    dt * (gamma_1*dUdt0.momentum.z + gamma_2*dUdt1.momentum.z + gamma_3*dUdt2.momentum.z);
	if (myConfig.MHD) {
	    // Magnetic field
	    U3.B.refx = U_old.B.x + dt * (gamma_1*dUdt0.B.x + gamma_2*dUdt1.B.x + gamma_3*dUdt2.B.x);
	    U3.B.refy = U_old.B.y + dt * (gamma_1*dUdt0.B.y + gamma_2*dUdt1.B.y + gamma_3*dUdt2.B.y);
	    U3.B.refz = U_old.B.z + dt * (gamma_1*dUdt0.B.z + gamma_2*dUdt1.B.z + gamma_3*dUdt2.B.z);
	}
	U3.total_energy = U_old.total_energy + 
	    dt * (gamma_1*dUdt0.total_energy + gamma_2*dUdt1.total_energy + gamma_3*dUdt2.total_energy);
	if (with_k_omega) {
	    U3.tke = U_old.tke + dt * (gamma_1*dUdt0.tke + gamma_2*dUdt1.tke + gamma_3*dUdt2.tke);
	    U3.tke = fmax(U3.tke, 0.0);
	    U3.omega = U_old.omega + dt * (gamma_1*dUdt0.omega + gamma_2*dUdt1.omega + gamma_3*dUdt2.omega);
	    U3.omega = fmax(U3.omega, U_old.mass);
	} else {
	    U3.tke = U_old.tke;
	    U3.omega = U_old.omega;
	}
	foreach(isp; 0 .. U3.massf.length) {
	    U3.massf[isp] = U_old.massf[isp] +
		dt * (gamma_1*dUdt0.massf[isp] + gamma_2*dUdt1.massf[isp] + gamma_3*dUdt2.massf[isp]);
	}
	// We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
	// Instead we get this from the conservation of total energy
	foreach(imode; 1 .. U3.energies.length) {
	    U3.energies[imode] = U_old.energies[imode] +
		dt * (gamma_1*dUdt0.energies[imode] + gamma_2*dUdt1.energies[imode] +
		      gamma_3*dUdt2.energies[imode]);
	}
	return;
    } // end stage_3_update_for_flow_on_fixed_grid()

    @nogc
    void stage_1_update_for_flow_on_moving_grid(double dt, bool with_k_omega) 
    {
	ConservedQuantities dUdt0 = dUdt[0];
	ConservedQuantities U0 = U[0];
	ConservedQuantities U1 = U[1];
	double gamma_1 = 1.0;
	double vr = volume[0] / volume[1];

	U1.mass = vr * (U0.mass + dt * gamma_1 * dUdt0.mass);
	U1.momentum.refx = vr * (U0.momentum.x + dt * gamma_1 * dUdt0.momentum.x);
	U1.momentum.refy = vr * (U0.momentum.y + dt * gamma_1 * dUdt0.momentum.y);
	U1.momentum.refz = vr * (U0.momentum.z + dt * gamma_1 * dUdt0.momentum.z);
	if (myConfig.MHD) {
	    // Magnetic field
	    U1.B.refx = vr * (U0.B.x + dt * gamma_1 * dUdt0.B.x);
	    U1.B.refy = vr * (U0.B.y + dt * gamma_1 * dUdt0.B.y);
	    U1.B.refz = vr * (U0.B.z + dt * gamma_1 * dUdt0.B.z);
	}
	U1.total_energy = vr * (U0.total_energy + dt * gamma_1 * dUdt0.total_energy);
	if (with_k_omega) {
	    U1.tke = vr * (U0.tke + dt * gamma_1 * dUdt0.tke);
	    U1.tke = fmax(U1.tke, 0.0);
	    U1.omega = vr * (U0.omega + dt * gamma_1 * dUdt0.omega);
	    U1.omega = fmax(U1.omega, U0.mass);
	} else {
	    U1.tke = U0.tke;
	    U1.omega = U0.omega;
	}
	foreach(isp; 0 .. U1.massf.length) {
	    U1.massf[isp] = vr * (U0.massf[isp] + dt * gamma_1 * dUdt0.massf[isp]);
	}
	// We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
	// Instead we get this from the conservation of total energy
	foreach(imode; 1 .. U1.energies.length) {
	    U1.energies[imode] = vr * (U0.energies[imode] + dt * gamma_1 * dUdt0.energies[imode]);
	}
	assert(false, "[TODO] not yet ready for use");
    } // end stage_1_update_for_flow_on_moving_grid()

    @nogc
    void stage_2_update_for_flow_on_moving_grid(double dt, bool with_k_omega) 
    {
	ConservedQuantities dUdt0 = dUdt[0];
	ConservedQuantities dUdt1 = dUdt[1];
	ConservedQuantities U0 = U[0];
	// ConservedQuantities U1 = U[1];
	ConservedQuantities U2 = U[2];
	double gamma_2 = 0.5;
	double gamma_1 = 0.5;
	double v_old = volume[0];
	double vol_inv = 1.0 / volume[2];
	gamma_1 *= volume[0]; gamma_2 *= volume[1]; // Roll-in the volumes for convenience below. 
    
	U2.mass = vol_inv * (v_old * U0.mass + dt * (gamma_1 * dUdt0.mass + gamma_2 * dUdt1.mass));
	U2.momentum.refx = vol_inv * (v_old * U0.momentum.x + 
				      dt * (gamma_1 * dUdt0.momentum.x + gamma_2 * dUdt1.momentum.x));
	U2.momentum.refy = vol_inv * (v_old * U0.momentum.y + 
				      dt * (gamma_1 * dUdt0.momentum.y + gamma_2 * dUdt1.momentum.y));
	U2.momentum.refz = vol_inv * (v_old * U0.momentum.z + 
				      dt * (gamma_1 * dUdt0.momentum.z + gamma_2 * dUdt1.momentum.z));
	if ( myConfig.MHD ) {
	    // Magnetic field
	    U2.B.refx = vol_inv * (v_old * U0.B.x + dt * (gamma_1 * dUdt0.B.x + gamma_2 * dUdt1.B.x));
	    U2.B.refy = vol_inv * (v_old * U0.B.y + dt * (gamma_1 * dUdt0.B.y + gamma_2 * dUdt1.B.y));
	    U2.B.refz = vol_inv * (v_old * U0.B.z + dt * (gamma_1 * dUdt0.B.z + gamma_2 * dUdt1.B.z));
	}
	U2.total_energy = vol_inv * (v_old * U0.total_energy + 
				     dt * (gamma_1 * dUdt0.total_energy + gamma_2 * dUdt1.total_energy));
	if ( with_k_omega ) {
	    U2.tke = vol_inv * (v_old * U0.tke + dt * (gamma_1 * dUdt0.tke + gamma_2 * dUdt1.tke));
	    U2.tke = fmax(U2.tke, 0.0);
	    U2.omega = vol_inv * (v_old * U0.omega + dt * (gamma_1 * dUdt0.omega + gamma_2 * dUdt1.omega));
	    U2.omega = fmax(U2.omega, U0.mass);
	} else {
	    U2.tke = vol_inv * (v_old * U0.tke);
	    U2.omega = vol_inv * (v_old * U0.omega);
	}
	foreach(isp; 0 .. U2.massf.length) {
	    U2.massf[isp] = vol_inv * (v_old * U0.massf[isp] +
				       dt * (gamma_1 * dUdt0.massf[isp] + 
					     gamma_2 * dUdt1.massf[isp]));
	}
	// We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
	// Instead we get this from the conservation of total energy
	foreach(imode; 1 .. U2.energies.length) {
	    U2.energies[imode] = vol_inv * (v_old * U0.energies[imode] +
					    dt * (gamma_1 * dUdt0.energies[imode] + 
						  gamma_2 * dUdt1.energies[imode]));
	}
	assert(false, "[TODO] not yet ready for use");
    } // end stage_2_update_for_flow_on_moving_grid()

    void chemical_increment(double dt, double T_frozen) 
    // Use the finite-rate chemistry module to update the species fractions
    // and the other thermochemical properties.
    {
	if (!fr_reactions_allowed || fs.gas.T[0] <= T_frozen) return;
	double T_save = fs.gas.T[0];
	if (myConfig.ignition_zone_active) {
	    // When active, replace gas temperature with an effective ignition temperature
	    foreach(zone; myConfig.ignition_zones) {
		if ( zone.is_inside(pos[0], myConfig.dimensions) ) fs.gas.T[0] = zone.Tig; 
	    }
	}
	try {
	    myConfig.reaction_update.update_state(fs.gas, dt, dt_chem, myConfig.gmodel);
	    if (myConfig.ignition_zone_active) {
		// Restore actual gas temperature
		fs.gas.T[0] = T_save;
	    }
	} catch(Exception err) {
	    writefln("catch %s", err.msg);
	    writeln("The chemical_increment() failed for cell: ", id);
	    writeln("The gas state after the failed update is:");
	    writefln("fs.gas %s", fs.gas);
	}

	// The update only changes mass fractions; we need to impose
	// a thermodynamic constraint based on a call to the equation of state.
	myConfig.gmodel.update_thermo_from_rhoe(fs.gas);

	// If we are doing a viscous sim, we'll need to ensure
	// viscous properties are up-to-date
	if (myConfig.viscous) myConfig.gmodel.update_trans_coeffs(fs.gas);
	// [TODO] if ( myConfig.diffusion ) myConfig.gmodel.update_diffusion_coeffs(fs.gas);

	// Finally, we have to manually update the conservation quantities
	// for the gas-dynamics time integration.
	// Species densities: mass of species isp per unit volume.
	foreach(isp; 0 .. fs.gas.massf.length)
	    U[0].massf[isp] = fs.gas.rho * fs.gas.massf[isp];

    } // end chemical_increment()

    void thermal_increment(double dt, double T_frozen_energy, GasModel gmodel) 
    // Use the nonequilibrium multi-Temperature module to update the
    // energy values and the other thermochemical properties.
    // We are assuming that this is done after a successful gas-dynamic update
    // and that the current conserved quantities are held in U[0].
    {
	if ( !fr_reactions_allowed || fs.gas.T[0] <= T_frozen_energy ) return;
	// [TODO] auto eeupdate = myConfig.energy_exchange_update_scheme;
	// [TODO] eeupdate.update_state(fs.gas, dt, dt_therm, gmodel);
	// The update only changes modal energies, we need to impose
	// a thermodynamic constraint based on a call to the equation
	// of state.
	gmodel.update_thermo_from_rhoe(fs.gas);
	// If we are doing a viscous sim, we'll need to ensure
	// viscous properties are up-to-date
	if ( myConfig.viscous ) gmodel.update_trans_coeffs(fs.gas);
	// [TODO] if ( myConfig.diffusion ) gmodel.update_diff_coeffs(fs.gas);
	// Finally, we have to manually update the conservation quantities
	// for the gas-dynamics time integration.
	// Independent energies energy: Joules per unit volume.
	foreach(imode; 0 .. U[0].energies.length) {
	    U[0].energies[imode] = fs.gas.rho * fs.gas.e[imode];
	}
	assert(false, "[TODO] not yet ready for use");
    } // end thermal_increment()

    double signal_frequency()
    {
	auto gmodel = myConfig.gmodel;
	bool with_k_omega = (myConfig.turbulence_model == TurbulenceModel.k_omega && 
			     !myConfig.separate_update_for_k_omega_source);
	double signal;
	double un_N, un_E, un_T, u_mag;
	double Bn_N = 0.0;
	double Bn_E = 0.0;
	double Bn_T = 0.0;
	double B_mag = 0.0;
	double ca2 = 0.0;
	double cfast = 0.0;
	double gam_eff;
	int statusf;
	//
	// Get the local normal velocities by rotating the local frame of reference.
	// Also, compute the velocity magnitude and recall the minimum length.
	un_N = fabs(fs.vel.dot(iface[Face.north].n));
	un_E = fabs(fs.vel.dot(iface[Face.east].n));
	if (myConfig.dimensions == 3) {
	    un_T = fabs(fs.vel.dot(iface[Face.top].n));
	    u_mag = sqrt(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y + fs.vel.z*fs.vel.z);
	} else {
	    un_T = 0.0;
	    u_mag = sqrt(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y);
	}
	if (myConfig.MHD) {
	    Bn_N = fabs(fs.B.dot(iface[Face.north].n));
	    Bn_E = fabs(fs.B.dot(iface[Face.east].n));
	    if (myConfig.dimensions == 3) {
		Bn_T = fabs(fs.B.dot(iface[Face.top].n));
	    }
	    u_mag = sqrt(fs.vel.x * fs.vel.x + fs.vel.y * fs.vel.y + fs.vel.z * fs.vel.z);
	    B_mag = sqrt(fs.B.x * fs.B.x + fs.B.y * fs.B.y + fs.B.z * fs.B.z);
	}
	// Check the INVISCID time step limit first,
	// then add a component to ensure viscous stability.
	if (myConfig.stringent_cfl) {
	    // Make the worst case.
	    if (myConfig.MHD) {
		ca2 = B_mag*B_mag / fs.gas.rho;
		cfast = sqrt(ca2 + fs.gas.a * fs.gas.a);
		signal = (u_mag + cfast) / L_min;
	    } else {
		// Hydrodynamics only
		signal = (u_mag + fs.gas.a) / L_min;
	    }
	} else {
	    // Standard signal speeds along each face.
	    double signalN, signalE, signalT;
	    if (myConfig.MHD) {
		double catang2_N, catang2_E, cfast_N, cfast_E;
		ca2 = B_mag * B_mag / fs.gas.rho;
		ca2 = ca2 + fs.gas.a * fs.gas.a;
		catang2_N = Bn_N * Bn_N / fs.gas.rho;
		cfast_N = 0.5 * ( ca2 + sqrt( ca2*ca2 - 4.0 * (fs.gas.a * fs.gas.a * catang2_N) ) );
		cfast_N = sqrt(cfast_N);
		catang2_E = Bn_E * Bn_E / fs.gas.rho;
		cfast_E = 0.5 * ( ca2 + sqrt( ca2*ca2 - 4.0 * (fs.gas.a * fs.gas.a * catang2_E) ) );
		cfast_E = sqrt(cfast_E);
		if (myConfig.dimensions == 3) {
		    double catang2_T, cfast_T;
		    catang2_T = Bn_T * Bn_T / fs.gas.rho;
		    cfast_T = 0.5 * ( ca2 + sqrt( ca2*ca2 - 4.0 * (fs.gas.a * fs.gas.a * catang2_T) ) );
		    cfast_T = sqrt(cfast_T);
		    signalN = (un_N + cfast_N) / jLength;
		    signal = signalN;
		    signalE = (un_E + cfast_E) / iLength;
		    if ( signalE > signal ) signal = signalE;
		    signalT = (un_T + cfast_T) / kLength;
		    if ( signalT > signal ) signal = signalT;
		} else {
		    signalN = (un_N + cfast) / jLength;
		    signalE = (un_E + cfast) / iLength;
		    signal = fmax(signalN, signalE);
		}
	    } else if (myConfig.dimensions == 3) {
		// eilmer -- 3D cells
		signalN = (un_N + fs.gas.a) / jLength;
		signal = signalN;
		signalE = (un_E + fs.gas.a) / iLength;
		if ( signalE > signal ) signal = signalE;
		signalT = (un_T + fs.gas.a) / kLength;
		if ( signalT > signal ) signal = signalT;
	    } else {
		// mbcns2 -- 2D cells
		// The velocity normal to the north face is assumed to run
		// along the length of the east face.
		signalN = (un_N + fs.gas.a) / jLength;
		signalE = (un_E + fs.gas.a) / iLength;
		signal = fmax(signalN, signalE);
	    }
	}
	if (myConfig.viscous && fs.gas.mu > 10.0e-23) {
	    // Factor for the viscous time limit.
	    // This factor is not included if viscosity is zero.
	    // See Swanson, Turkel and White (1991)
	    gam_eff = gmodel.gamma(fs.gas);
	    // Need to sum conductivities for TNE
	    double k_total = 0.0;
	    foreach(i; 0 .. fs.gas.k.length) k_total += fs.gas.k[i];
	    double Prandtl = fs.gas.mu * gmodel.Cp(fs.gas) / k_total;
	    if (myConfig.dimensions == 3) {
		signal += 4.0 * myConfig.viscous_factor * (fs.gas.mu + fs.mu_t)
		    * gam_eff / (Prandtl * fs.gas.rho)
		    * (1.0/(iLength*iLength) + 1.0/(jLength*jLength) + 1.0/(kLength*kLength))
		    * myConfig.viscous_signal_factor;
	    } else {
		signal += 4.0 * myConfig.viscous_factor * (fs.gas.mu + fs.mu_t) 
		    * gam_eff / (Prandtl * fs.gas.rho)
		    * (1.0/(iLength*iLength) + 1.0/(jLength*jLength))
		    * myConfig.viscous_signal_factor;
	    }
	}
	if (with_k_omega) {
	    if (fs.omega > signal) signal = fs.omega;
	}
	return signal;
    } // end signal_frequency()

    @nogc
    void turbulence_viscosity_zero() 
    {
	fs.mu_t = 0.0;
	fs.k_t = 0.0;
    }

    @nogc
    void turbulence_viscosity_zero_if_not_in_zone() 
    {
	if ( in_turbulent_zone ) {
	    /* Do nothing, leaving the turbulence quantities as set. */
	} else {
	    /* Presume this part of the flow is laminar; clear turbulence quantities. */
	    fs.mu_t = 0.0;
	    fs.k_t = 0.0;
	}
    }

    @nogc
    void turbulence_viscosity_limit(double factor) 
    // Limit the turbulent viscosity to reasonable values relative to
    // the local molecular viscosity.
    // In shock started flows, we seem to get crazy values on the
    // starting shock structure and the simulations do not progress.
    {
	fs.mu_t = fmin(fs.mu_t, factor * fs.gas.mu);
	fs.k_t = fmin(fs.k_t, factor * fs.gas.k[0]); // ASSUMPTION re k[0]
    }

    @nogc
    void turbulence_viscosity_factor(double factor) 
    // Scale the turbulent viscosity to model effects
    // such as not-fully-developed turbulence that might be expected
    // in short-duration transient flows.
    {
	fs.mu_t *= factor;
	fs.k_t *= factor;
    }

    void turbulence_viscosity_k_omega() 
    {
	auto gmodel = myConfig.gmodel;
	if ( myConfig.turbulence_model != TurbulenceModel.k_omega ) {
	    // [TODO] may have to do something better if another turbulence model is active.
	    fs.mu_t = 0.0;
	    fs.k_t = 0.0;
	    return;
	}
	double dudx, dudy, dvdx, dvdy;
	double S_bar_squared;
	double C_lim = 0.875;
	double beta_star = 0.09;
	if ( myConfig.dimensions == 2 ) {
	    // 2D cartesian or 2D axisymmetric
	    double avg2D(int i, int j)() 
		if ( is(typeof(vtx[0].grad_vel[i][j]) == double) )
	    {
		return 0.25 * (vtx[0].grad_vel[i][j] + vtx[1].grad_vel[i][j] + 
			       vtx[2].grad_vel[i][j] + vtx[3].grad_vel[i][j]);
	    }
	    dudx = avg2D!(0,0)(); dudy = avg2D!(0,1)();
	    dvdx = avg2D!(1,0)(); dvdy = avg2D!(1,1)();
	    if ( myConfig.axisymmetric ) {
		// 2D axisymmetric
		double v_over_y = fs.vel.y / pos[0].y;
		S_bar_squared = dudx*dudx + dvdy*dvdy + v_over_y*v_over_y
		    - 1.0/3.0 * (dudx + dvdy + v_over_y)
		    * (dudx + dvdy + v_over_y)
		    + 0.5 * (dudy + dvdx) * (dudy + dvdx) ;
	    } else {
		// 2D cartesian
		S_bar_squared = dudx*dudx + dvdy*dvdy
		    - 1.0/3.0 * (dudx + dvdy) * (dudx + dvdy)
		    + 0.5 * (dudy + dvdx) * (dudy + dvdx);
	    }
	} else {
	    // 3D cartesian
	    double dudz, dvdz, dwdx, dwdy, dwdz;
	    double avg3D(int i, int j)() 
		if ( is(typeof(vtx[0].grad_vel[i][j]) == double) )
	    {
		return 0.125 * (vtx[0].grad_vel[i][j] + vtx[1].grad_vel[i][j] + 
				vtx[2].grad_vel[i][j] + vtx[3].grad_vel[i][j] +
				vtx[4].grad_vel[i][j] + vtx[5].grad_vel[i][j] + 
				vtx[6].grad_vel[i][j] + vtx[7].grad_vel[i][j]);
	    }
	    dudx = avg3D!(0,0)(); dudy = avg3D!(0,1)(); dudz = avg3D!(0,2)();
	    dvdx = avg3D!(1,0)(); dvdy = avg3D!(1,1)(); dvdz = avg3D!(1,2)();
	    dwdx = avg3D!(2,0)(); dwdy = avg3D!(2,1)(); dwdz = avg3D!(2,2)();
	    S_bar_squared =  dudx*dudx + dvdy*dvdy + dwdz*dwdz
		- 1.0/3.0*(dudx + dvdy + dwdz)*(dudx + dvdy + dwdz)
		+ 0.5 * (dudy + dvdx) * (dudy + dvdx)
		+ 0.5 * (dudz + dwdx) * (dudz + dwdx)
		+ 0.5 * (dvdz + dwdy) * (dvdz + dwdy);
	}
	S_bar_squared = fmax(0.0, S_bar_squared);
	double omega_t = fmax(fs.omega, C_lim*sqrt(2.0*S_bar_squared/beta_star));
	fs.mu_t = fs.gas.rho * fs.tke / omega_t;
	double Pr_t = myConfig.turbulence_prandtl_number;
	fs.k_t = gmodel.Cp(fs.gas) * fs.mu_t / Pr_t;
    } // end turbulence_viscosity_k_omega()

    void update_k_omega_properties(double dt) 
    {
	// Do not update k_omega properties if we are in laminar block
	if ( !in_turbulent_zone ) return;

	double DrtkeDt_perturbTke, DromegaDt_perturbTke;
	double DrtkeDt_perturbOmega, DromegaDt_perturbOmega;
	double DGkDzetak, DGkDzetaw, DGwDzetak, DGwDzetaw;
	double DfkDk, DfkDw, DfwDk, DfwDw;
	double Gk, Gw;
	double delta_rtke, delta_romega;
	double tke, omega;
	double tke_current, omega_current;
	double tke_updated, omega_updated;
	double DrtkeDt_current, DromegaDt_current;
	double DrtkeDt_updated, DromegaDt_updated;
	double perturbFactor = 1.01;  // Perturbation factor for perturbation
	// analysis to get derivatives
	double tol = 1.0e-6;          // Tolerance for the Newton-solve loop

	// Encode conserved quantities for cell.
	U[0].tke = fs.gas.rho * fs.tke;
	U[0].omega = fs.gas.rho * fs.omega;

	// Start of implicit updating scheme.
	tke_current = fs.tke; omega_current = fs.omega;  // Current values of tke and omega
	tke_updated = fs.tke; omega_updated = fs.omega;  // First guess of updated values

	// Work out values of Drtke_current and DromegaDt_current.
	this.k_omega_time_derivatives(DrtkeDt_current, DromegaDt_current, tke_current, omega_current);

	// Implicit updating scheme.
	// A for-loop is used to limit the Newton-solve to 20 steps
	// just in case convergence does not occur.
	for ( int i = 1; i <= 20; ++i ) {
	    // Work out unperturbed values of Drtke_updated and DromegaDt_updated.
	    this.k_omega_time_derivatives(DrtkeDt_updated, DromegaDt_updated, tke_updated, omega_updated);
	    // Perturb tke and obtain perturbed values of DrtkeDt_updated and DromegaDt_updated.
	    tke = perturbFactor * tke_updated; omega = omega_updated;
	    this.k_omega_time_derivatives(DrtkeDt_perturbTke, DromegaDt_perturbTke, tke, omega);
	    // Perturb omega and obtain perturbed values of DrtkeDt_updated and DromegaDt_updated.
	    tke = tke_updated; omega = perturbFactor * omega_updated;
	    this.k_omega_time_derivatives(DrtkeDt_perturbOmega, DromegaDt_perturbOmega, tke, omega);
	    // Compute derivatives from perturb values.
	    // FIX-ME : Dividing by tke and omega (instead of rtke and romega) seems to work (gives
	    //          same results as explicit update scheme), but we will keep this note here for
	    //          future reference..
	    DfkDk = (DrtkeDt_perturbTke - DrtkeDt_updated) / ((perturbFactor - 1.0) * tke_updated);
	    DfkDw = (DrtkeDt_perturbOmega - DrtkeDt_updated) / ((perturbFactor - 1.0) * omega_updated);
	    DfwDk = (DromegaDt_perturbTke - DromegaDt_updated) / ((perturbFactor - 1.0) * tke_updated);
	    DfwDw = (DromegaDt_perturbOmega - DromegaDt_updated) / ((perturbFactor - 1.0) * omega_updated);
	    // Compute components in matrix A of Ax=B problem.
	    DGkDzetak = -1.0 + 0.5 * dt * DfkDk;
	    DGkDzetaw = 0.5 * dt * DfkDw;
	    DGwDzetak = 0.5 * dt * DfwDk;
	    DGwDzetaw = -1.0 + 0.5 * dt * DfwDw;
	    // Compute vector B of Ax=B problem.
	    Gk = fs.gas.rho * tke_updated - fs.gas.rho * tke_current -
		0.5 * dt * (DrtkeDt_updated + DrtkeDt_current);
	    Gw = fs.gas.rho * omega_updated - fs.gas.rho * omega_current -
		0.5 * dt * (DromegaDt_updated + DromegaDt_current);
	    // Solve Ax=B algebraically.
	    delta_rtke = (DGkDzetaw * Gw - DGwDzetaw * Gk) /
		(DGwDzetak * DGkDzetaw - DGwDzetaw * DGkDzetak);
	    delta_romega = (Gk - DGkDzetak * delta_rtke) / DGkDzetaw;
	    // Assign the updated tke and omega values if delta_rtke and
	    // delta_romega are both smaller than the given tolerance, and
	    // then break out from the Newton-solve loop.
	    if (fabs(delta_rtke) <= tol && fabs(delta_romega) <= tol) {
		fs.tke = tke_updated;
		fs.omega = omega_updated;
		break;
	    } else {
		// Compute next estimates for rtke and romega from
		// delta_rtke and delta_romega.
		if (delta_rtke + fs.gas.rho * tke_updated < 0.0) {
		    // Don't let rtke go negative.
		    U[0].tke = fs.gas.rho * tke_updated;
		} else {
		    // Next estimate for rtke.
		    U[0].tke = delta_rtke + fs.gas.rho * tke_updated;
		}
		if (delta_romega + fs.gas.rho * omega_updated < 0.0) {
		    // Don't let romega go negative.
		    U[0].omega = fs.gas.rho * omega_updated;
		} else {
		    // Next estimate for romega.
		    U[0].omega = delta_romega + fs.gas.rho * omega_updated;
		}
		// Decode for the next step of the Newton-solve loop
		tke_updated = U[0].tke / fs.gas.rho;
		omega_updated = U[0].omega / fs.gas.rho;
	    }
	} // End of Newton-solve loop for implicit update scheme
    } // end update_k_omega_properties()

    @nogc
    void k_omega_time_derivatives(ref double Q_rtke, ref double Q_romega, double tke, double omega) 
    // Compute k-omega source terms.
    //
    // Production and Dissipation expressions for turbulence kinetic energy
    // and turbulence frequency (or pseudo-vorticity). Based on Wilcox's 2006 model.
    //
    // Jan 2007: Initial implementation (Jan-Pieter Nap, PJ)
    // Dec 2008: Implementation of the 3D terms (W Chan)
    // Jan 2011: Minor modification to allow for implicit updating of tke and omega (W Chan)
    //           All "fs->tke" and "fs->omega" instances are replaced with tke and omega.
    // Jul 2014: Port to D by PJ
    {
	if ( myConfig.turbulence_model != TurbulenceModel.k_omega ) {
	    // [TODO] may need to do something better is another turbulence model is active.
	    Q_rtke = 0.0;
	    Q_romega = 0.0;
	    return;
	}
	double dudx, dudy, dvdx, dvdy;
	double dtkedx, dtkedy, domegadx, domegady;
	double alpha = 0.52;
	double beta_0 = 0.0708;
	double beta;
	double beta_star = 0.09;
	double P_K, D_K, P_W, D_W;
	double cross_diff;
	double sigma_d = 0.0;
	double WWS, X_w, f_beta;
	if ( myConfig.dimensions == 2 ) {
	    // 2D cartesian or 2D axisymmetric
	    // The following compile-time function is more complicated 
	    // than the resulting 2D code and actually take up just as much space, 
	    // however, it's practice for the main event (3D code, below).
	    dudx = mixin(avg_over_vtx_2D("grad_vel[0][0]"));
	    dudy = mixin(avg_over_vtx_2D("grad_vel[0][1]"));
	    dvdx = mixin(avg_over_vtx_2D("grad_vel[1][0]"));
	    dvdy = mixin(avg_over_vtx_2D("grad_vel[1][1]"));
	    dtkedx = mixin(avg_over_vtx_2D("grad_tke.x"));
	    dtkedy = mixin(avg_over_vtx_2D("grad_tke.y"));
	    domegadx = mixin(avg_over_vtx_2D("grad_omega.x"));
	    domegady = mixin(avg_over_vtx_2D("grad_omega.y"));
	    if ( myConfig.axisymmetric ) {
		// 2D axisymmetric
		double v_over_y = fs.vel.y / pos[0].y;
		// JP.Nap correction from 03-May-2007 (-v_over_y in parentheses)
		// P_K -= 0.6667 * mu_t * v_over_y * (dudx+dvdy-v_over_y);
		// Wilson Chan correction to JP Nap's version (13 Dec 2008)
		P_K = 2.0 * fs.mu_t * (dudx*dudx + dvdy*dvdy)
		    + fs.mu_t * (dudy + dvdx) * (dudy + dvdx)
		    - 2.0/3.0 * fs.mu_t * (dudx + dvdy + v_over_y)
		    * (dudx + dvdy + v_over_y)
		    + 2.0 * fs.mu_t * (v_over_y) * (v_over_y)
		    - 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy + v_over_y);
		WWS = 0.25 * (dvdx - dudy) * (dvdx - dudy) * v_over_y ;
	    } else {
		// 2D cartesian
		P_K = 1.3333 * fs.mu_t * (dudx*dudx - dudx*dvdy + dvdy*dvdy)
		    + fs.mu_t * (dudy + dvdx) * (dudy + dvdx)
		    - 0.66667 * fs.gas.rho * tke * (dudx + dvdy);
		WWS = 0.0 ;
	    }
	    cross_diff = dtkedx * domegadx + dtkedy * domegady ;
	} else {
	    // 3D cartesian
	    double dudz, dvdz, dwdx, dwdy, dwdz;
	    double dtkedz, domegadz;
	    dudx = mixin(avg_over_vtx_3D("grad_vel[0][0]"));
	    dudy = mixin(avg_over_vtx_3D("grad_vel[0][1]"));
	    dudz = mixin(avg_over_vtx_3D("grad_vel[0][2]"));
	    dvdx = mixin(avg_over_vtx_3D("grad_vel[1][0]"));
	    dvdy = mixin(avg_over_vtx_3D("grad_vel[1][1]"));
	    dvdz = mixin(avg_over_vtx_3D("grad_vel[1][2]"));
	    dwdx = mixin(avg_over_vtx_3D("grad_vel[2][0]"));
	    dwdy = mixin(avg_over_vtx_3D("grad_vel[2][1]"));
	    dwdz = mixin(avg_over_vtx_3D("grad_vel[2][2]"));
	    dtkedx = mixin(avg_over_vtx_3D("grad_tke.x"));
	    dtkedy = mixin(avg_over_vtx_3D("grad_tke.y"));
	    dtkedz = mixin(avg_over_vtx_3D("grad_tke.z"));
	    domegadx = mixin(avg_over_vtx_3D("grad_omega.x"));
	    domegady = mixin(avg_over_vtx_3D("grad_omega.y"));
	    domegadz = mixin(avg_over_vtx_3D("grad_omega.z"));
	    P_K = 2.0 * fs.mu_t * (dudx*dudx + dvdy*dvdy + dwdz*dwdz)
		- 2.0/3.0 * fs.mu_t * (dudx + dvdy + dwdz) * (dudx + dvdy + dwdz)
		- 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy + dwdz)
		+ fs.mu_t * (dudy + dvdx) * (dudy + dvdx)
		+ fs.mu_t * (dudz + dwdx) * (dudz + dwdx)
		+ fs.mu_t * (dvdz + dwdy) * (dvdz + dwdy) ;
	    cross_diff = dtkedx * domegadx + dtkedy * domegady + dtkedz * domegadz ;
	    WWS = 0.25 * (dudy - dvdx) * (dudy - dvdx) * dwdz
		+ 0.25 * (dudz - dwdx) * (dudz - dwdx) * dvdy
		+ 0.25 * (dvdz - dwdy) * (dvdz - dwdy) * dudx
		+ 0.25 * (dudy - dvdx) * (dvdz - dwdy) * (dwdx + dudz)
		+ 0.25 * (dudz - dwdx) * (dwdy - dvdz) * (dudy + dvdx)
		+ 0.25 * (dvdx - dudy) * (dudz - dwdx) * (dwdy + dvdx) ;
	}

	D_K = beta_star * fs.gas.rho * tke * omega;
    
	// Apply a limit to the tke production as suggested by Jeff White, November 2007.
	const double P_OVER_D_LIMIT = 25.0;
	P_K = fmin(P_K, P_OVER_D_LIMIT*D_K);

	if ( cross_diff > 0 ) sigma_d = 0.125;
	P_W = alpha * omega / fmax(tke, small_tke) * P_K +
	    sigma_d * fs.gas.rho / fmax(omega, small_omega) * cross_diff;

	X_w = fabs(WWS / pow(beta_star*omega, 3)) ;
	f_beta = (1.0 + 85.0 * X_w) / (1.0 + 100.0 * X_w) ;
	beta = beta_0 * f_beta;
	D_W = beta * fs.gas.rho * omega * omega;

	Q_rtke = P_K - D_K;
	Q_romega = P_W - D_W;
    } // end k_omega_time_derivatives()

    @nogc
    void clear_source_vector() 
    // When doing the gasdynamic update stages, the source vector values
    // are accumulated for the inviscid and then viscous terms, so we 
    // have to start with a clean slate, so to speak.
    {
	Q.mass = 0.0;
	Q.momentum.refx = 0.0; Q.momentum.refy = 0.0; Q.momentum.refz = 0.0;
	Q.B.refx = 0.0; Q.B.refy = 0.0; Q.B.refz = 0.0;
	Q.total_energy = 0.0;
	Q.tke = 0.0;
	Q.omega = 0.0;
	foreach(ref elem; Q.massf) elem = 0.0;
	foreach(ref elem; Q.energies) elem = 0.0;
	Q_rE_rad = 0.0;
    } // end clear_source_vector()

    @nogc
    void add_inviscid_source_vector(int gtl, double omegaz=0.0) 
    // Add the components of the source vector, Q, for inviscid flow.
    //
    // Currently, the axisymmetric equations include the
    // pressure contribution to the y-momentum equation
    // here rather than in the boundary fluxes.
    // By default, assume 2D-planar, or 3D-Cartesian flow.
    {
	if (omegaz != 0.0) {
	    // Rotating frame.
	    double rho = fs.gas.rho;
	    double x = pos[gtl].x;
	    double y = pos[gtl].y;
	    double wx = fs.vel.x;
	    double wy = fs.vel.y;
	    // Coriolis and centrifugal forces contribute to momenta.
	    Q.momentum.refx += rho * (omegaz*omegaz*x + 2.0*omegaz*wy);
	    Q.momentum.refy += rho * (omegaz*omegaz*y - 2.0*omegaz*wx);
	    // There is no contribution to the energy equation in the rotating frame
	    // because it is implicit in the use of rothalpy as the conserved quantity.
	}
	if (myConfig.axisymmetric) {
	    // For axisymmetric flow:
	    // pressure contribution from the Front and Back (radial) interfaces.
	    Q.momentum.refy += fs.gas.p * areaxy[gtl] / volume[gtl];
	}
	// Species production (other than chemistry).
	// For the chemistry, see chemical_increment().
	// Individual energies (other than energy exchange)
	// For the energy exchange, see thermal_increment()
	// Radiation can potentially be removed from both the electronic and
	// total energy source terms.
	if (myConfig.radiation) {
	    // Radiative source term should be already calculated
	    // Add value to total energy
	    // FIX-ME: - assuming electronic mode is the last in the vector of energies
	    //         - what about Q_renergies[0]?
	    Q.total_energy += Q_rE_rad;
	    Q.energies.back() += Q_rE_rad;
	}
	return;
    } // end add_inviscid_source_vector()

    @nogc
    void add_viscous_source_vector(bool with_k_omega) 
    {
	if (myConfig.axisymmetric) {
	    // For viscous, axisymmetric flow:
	    double v_over_y = fs.vel.y / pos[0].y;
	    double dudx = mixin(avg_over_vtx_2D("grad_vel[0][0]"));
	    double dvdy = mixin(avg_over_vtx_2D("grad_vel[1][1]"));
	    double mu = 
		0.25 * (iface[Face.east].fs.gas.mu + iface[Face.west].fs.gas.mu +
			iface[Face.north].fs.gas.mu + iface[Face.south].fs.gas.mu) +
		0.25 * (iface[Face.east].fs.mu_t + iface[Face.west].fs.mu_t +
			iface[Face.north].fs.mu_t + iface[Face.south].fs.mu_t);
	    mu *= myConfig.viscous_factor;
	    double lmbda = -2.0/3.0 * mu;
	    double tau_00 = 2.0 * mu * v_over_y + lmbda * (dudx + dvdy + v_over_y);
	    // Y-Momentum; viscous stress contribution from the front and Back interfaces.
	    // Note that these quantities are approximated at the
	    // mid-point of the cell face and so should never be
	    // singular -- at least I hope that this is so.
	    Q.momentum.refy -= tau_00 * areaxy[0] / volume[0];
	} // end if ( myConfig.axisymmetric )

	if (with_k_omega) {
	    double Q_tke = 0.0; double Q_omega = 0.0;
	    if ( in_turbulent_zone ) {
		this.k_omega_time_derivatives(Q_tke, Q_omega, fs.tke, fs.omega);
	    }
	    Q.tke += Q_tke; Q.omega += Q_omega;
	}

	if (myConfig.electric_field_work) {
	    // Work done on electrons due to electric field induced by charge separation
	    // on scales less than the Debye length
	    // FIXME: Only consistent with ambipolar diffusion. Currently this is up to
	    //        the user to enforce.
	    double udivpe = 0.0;
	    if ( myConfig.dimensions == 2 ) {
		// Estimate electron pressure gradient as average of all vertices
		double dpedx = mixin(avg_over_vtx_2D("grad_pe.x"));
		double dpedy = mixin(avg_over_vtx_2D("grad_pe.y"));
		// Approximation for work done on electrons: u dot div(pe)
		udivpe = fs.vel.x * dpedx + fs.vel.y * dpedy;
	    } else {
		// Estimate electron pressure gradient as average of all vertices
		double dpedx = mixin(avg_over_vtx_3D("grad_pe.x"));
		double dpedy = mixin(avg_over_vtx_3D("grad_pe.y"));
		double dpedz = mixin(avg_over_vtx_3D("grad_pe.z"));
		// Approximation for work done on electrons: u dot div(pe)
		udivpe = fs.vel.x * dpedx + fs.vel.y * dpedy + fs.vel.z * dpedz;
	    }
	    // [TODO] FIXME: Assuming the free electron energy is included in the last mode
	    Q.energies.back() += udivpe * myConfig.diffusion_factor;
	} // end if ( myConfig.electric_field_work )
	return;
    } // end add_viscous_source_vector()

    double calculate_wall_Reynolds_number(int which_boundary, GasModel gmodel)
    {
	FVInterface IFace = iface[which_boundary];
	gmodel.update_thermo_from_rhoT(IFace.fs.gas); // Note that we adjust IFace here.
	double a_wall = IFace.fs.gas.a;
	double cell_width = 0.0;
	if ( which_boundary == Face.east || which_boundary == Face.west )
	    cell_width = iLength;
	else if ( which_boundary == Face.north || which_boundary == Face.south )
	    cell_width = jLength;
	else if ( which_boundary == Face.top || which_boundary == Face.bottom )
	    cell_width = kLength;
	double Re_wall = IFace.fs.gas.rho * a_wall * cell_width / IFace.fs.gas.mu;
	return Re_wall;
    } // end calculate_wall_Reynolds_number()

    void store_rad_scaling_params() 
    // Store parameters for (re-)scaling of radiative source term.
    // Simple rho x T**4 scaling seems to be adequate.
    {
	// 1. Store the freshly computed radiative flux as the 'original'
	Q_rad_org = Q_rE_rad;
	// 2. Compute the scaling factor based on local gas properties
	// NOTE: - The idea is that f_rad_org is proportional to actual value
	//       - Assuming that the last temperature is the electronic temperature
	double T = fs.gas.T.back();
	if ( Q_rad_org <= 0.0 ) {
	    // This cell is a net emitter
	    f_rad_org = fs.gas.rho * pow(T, 4);
	} else if ( Q_rad_org > 0.0 ) {
	    // This cell is a net absorber
	    f_rad_org = fs.gas.rho / pow(T, 4);
	}
    } // end store_rad_scaling_params()

    void rescale_Q_rE_rad() 
    {
	// 1. Compute the current scaling factor based on local gas properties
	double T = fs.gas.T[0];
	double f_rad_new = 1.0;
	if ( Q_rad_org <= 0.0 ) {
	    // This cell is a net emitter
	    f_rad_new = fs.gas.rho * pow(T, 4);
	}
	else if ( Q_rad_org > 0.0 ) {
	    // This cell is a net absorber
	    f_rad_new = fs.gas.rho / pow(T, 4);
	}
	// 2. (Re-)scale the original source term
	Q_rE_rad = ( f_rad_new / f_rad_org ) * Q_rad_org;
    } // end rescale_Q_rE_rad()

    void reset_Q_rad_to_zero() 
    {
	Q_rE_rad = 0.0;
    } // end reset_Q_rad_to_zero()

    double rad_scaling_ratio() 
    {
	// 1. Compute the current scaling factor based on local gas properties
	double T = fs.gas.T[0];
	double f_rad = 1.0;
	if ( Q_rE_rad <= 0.0 ) {
	    // This cell is a net emitter
	    f_rad = fs.gas.rho * pow(T, 4);
	}
	else if ( Q_rE_rad > 0.0 ) {
	    // This cell is a net absorber
	    f_rad = fs.gas.rho / pow(T, 4);
	}
	return fabs( f_rad - f_rad_org ) / f_rad_org;
    } // end rad_scaling_ratio()

} // end class FVCell


int number_of_values_in_cell_copy(int type_of_copy)
// This function must match the copy-to/from-buffer methods above.
// The buffers are used for communication between worker processes.
{
    return 0; // [TODO] something sensible, eventually.
} // end number_of_values_in_cell_copy()

string[] variable_list_for_cell(GasModel gmodel)
{
    // This function needs to be kept consistent with functions
    // FVCell.write_values_to_string, FVCell.scan_values_from_string
    // (found above) and with the corresponding Python functions
    // write_cell_data and variable_list_for_cell
    // that may be found in app/eilmer3/source/e3_flow.py.
    string[] list;
    list ~= ["pos.x", "pos.y", "pos.z", "volume"];
    list ~= ["rho", "vel.x", "vel.y", "vel.z"];
    if ( GlobalConfig.MHD ) list ~= ["B.x", "B.y", "B.z"];
    list ~= ["p", "a", "mu"];
    foreach(i; 0 .. gmodel.n_modes) list ~= "k[" ~ to!string(i) ~ "]";
    list ~= ["mu_t", "k_t", "S"];
    if ( GlobalConfig.radiation ) list ~= ["Q_rad_org", "f_rad_org", "Q_rE_rad"];
    list ~= ["tke", "omega"];
    foreach(i; 0 .. gmodel.n_species) {
	auto name = cast(char[]) gmodel.species_name(i);
	name = tr(name, " \t", "--", "s"); // Replace internal whitespace with dashes.
	list ~= ["massf[" ~ to!string(i) ~ "]-" ~ to!string(name)];
    }
    if ( gmodel.n_species > 1 ) list ~= ["dt_chem"];
    foreach(i; 0 .. gmodel.n_modes) list ~= ["e[" ~ to!string(i) ~ "]", "T[" ~ to!string(i) ~ "]"];
    if ( gmodel.n_modes > 1 ) list ~= ["dt_therm"];
    return list;
} // end variable_list_for_cell()
