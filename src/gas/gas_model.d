/**
 * gasmodel.d
 * Storage arrangement for the data defining a gas state,
 * interface description of the gas model functionality and
 * utilities to create specific gas model objects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22, first cut, exploring the options.
 */

module gas.gas_model;

import std.conv;
import std.math;
import std.stdio;

import util.msg_service;
import gas.physical_constants;

immutable double SMALL_MOLE_FRACTION = 1.0e-15;
immutable double MIN_MASS_FRACTION = 1.0e-30;
immutable double MIN_MOLES = 1.0e-30;
immutable double T_MIN = 20.0; 

class GasModel {
public:
    @nogc @property uint n_species() const { return _n_species; }
    @nogc @property uint n_modes() const { return _n_modes; }
    @property ref double[] mol_masses() { return _mol_masses; }
    final string species_name(int i) const { return _species_names[i]; }

    // Methods to be overridden.
    abstract void update_thermo_from_pT(GasState Q);
    abstract void update_thermo_from_rhoe(GasState Q);
    abstract void update_thermo_from_rhoT(GasState Q);
    abstract void update_thermo_from_rhop(GasState Q);
    abstract void update_thermo_from_ps(GasState Q, double s);
    abstract void update_thermo_from_hs(GasState Q, double h, double s);
    abstract void update_sound_speed(GasState Q);
    abstract void update_trans_coeffs(GasState Q);
    // const void update_diff_coeffs(ref GasState Q) {}

    // Methods to be overridden.
    abstract double dedT_const_v(in GasState Q); 
    abstract double dhdT_const_p(in GasState Q); 
    abstract double dpdrho_const_T(in GasState Q); 
    abstract double gas_constant(in GasState Q);
    abstract double internal_energy(in GasState Q);
    abstract double enthalpy(in GasState Q);
    abstract double entropy(in GasState Q);
    
    final double Cv(in GasState Q) { return dedT_const_v(Q); }
    final double Cp(in GasState Q) { return dhdT_const_p(Q); }
    final double R(in GasState Q)  { return gas_constant(Q); }
    final double gamma(in GasState Q) { return Cp(Q)/Cv(Q); }
    final double molecular_mass(in GasState Q) 
    in {
	assert(Q.massf.length == _mol_masses.length, brokenPreCondition("Inconsistent array lengths."));
    }
    body {
	return mixture_molecular_mass(Q.massf, _mol_masses);
    }
    final void massf2molef(in GasState Q, double[] molef) 
    in {
	assert(Q.massf.length == molef.length, brokenPreCondition("Inconsistent array lengths."));
    }
    body {
	gas.gas_model.massf2molef(Q.massf, _mol_masses, molef);
    }
    final void molef2massf(in double[] molef, GasState Q) 
    in {
	assert(Q.massf.length == molef.length, brokenPreCondition("Inconsistent array lengths."));
    }
    body {
	gas.gas_model.molef2massf(molef, _mol_masses, Q.massf);
    }
    final void massf2conc(in GasState Q, double[] conc) 
    in {
	assert(Q.massf.length == conc.length, brokenPreCondition("Inconsistent array lenghts."));
    }
    body {
	foreach ( i; 0.._n_species ) {
	    conc[i] = Q.massf[i]*Q.rho / _mol_masses[i];
	    if ( conc[i] < MIN_MOLES ) conc[i] = 0.0;
	}
    }
    final void conc2massf(in double[] conc, GasState Q) 
    in {
	assert(Q.massf.length == conc.length, brokenPreCondition("Inconsisten array lenghts."));
    }
    body {
	foreach ( i; 0.._n_species ) {
	    Q.massf[i] = conc[i]*_mol_masses[i] / Q.rho;
	    if ( Q.massf[i] < MIN_MASS_FRACTION ) Q.massf[i] = 0.0;
	}
    }

protected:
    // These data need to be properly initialized by the derived class.
    uint _n_species;
    uint _n_modes;
    string[] _species_names;
    double[] _mol_masses;
}


class GasState {
public:
    /// Thermodynamic properties.
    double rho;  /// density, kg/m**3
    double p;    /// presure, Pa
    double p_e;  /// electron pressure
    double a;    /// sound speed, m/s
    double[] e;  /// specific internal energies, J/kg
    double[] T;  /// temperatures, K
    /// Transport properties
    double mu;   /// viscosity, Pa.s
    double[] k;  /// thermal conductivities, W/(m.k)
    // double[][] D_AB; /// binary diffusion coefficients
    double sigma;    /// electrical conductivity, S/m
    /// Composition
    double[] massf;  /// species mass fractions
    double quality;  /// vapour quality

    this(uint n_species, uint n_modes)
    {
	massf.length = n_species;
	e.length = n_modes;
	T.length = n_modes;
	k.length = n_modes;
    }

    this(GasModel gm, in double p_init, in double[] T_init, 
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in double sigma_init=0.0)
    {
	p = p_init;
	p_e = p_init;
	T.length = gm.n_modes;
	foreach(i; 0 .. gm.n_modes) {
	    try
		{ T[i] = T_init[i]; }
	    catch (Exception e)
		// We assume that at least 1 temperature arrived.
		{ T[i] = T_init[0]; }
	}
	e.length = gm.n_modes;
	k.length = gm.n_modes;
	massf.length = gm.n_species;
	foreach(i; 0 .. gm.n_species) {
	    if (i < massf_init.length) { 
		massf[i] = massf_init[i];
	    } else {
		massf[i] = 0.0;
	    }
	}
	quality = quality_init;
	sigma = sigma_init;
	// Now, evaluate the rest of the properties using the gas model.
	gm.update_thermo_from_pT(this);
	gm.update_sound_speed(this);
	gm.update_trans_coeffs(this);
    }

    this(GasModel gm, in double p_init, in double T_init, 
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in double sigma_init=0.0)
    {
	double[] Tlocal;
	Tlocal.length = gm.n_modes;
	foreach(ref Tmode; Tlocal) {
	    Tmode = T_init;
	}
	this(gm, p_init, Tlocal, massf_init, quality_init, sigma_init);
    }

    this(in GasState other) 
    {
	rho = other.rho;
	p = other.p;
	p_e = other.p_e;
	a = other.a;
	e = other.e.dup;
	T = other.T.dup;
	mu = other.mu;
	k = other.k.dup;
	// D_AB
	sigma = other.sigma;
	massf = other.massf.dup;
	quality = other.quality;
    }

    // Postblit constructor for struct
    // Not needed for class
    // this(this)
    // {
    // 	massf = massf.dup;
    // 	e = e.dup;
    // 	T = T.dup;
    // 	k = k.dup;
    // }

    @nogc void copy_values_from(ref const(GasState) other) 
    {
	rho = other.rho;
	p = other.p;
	p_e = other.p_e;
	a = other.a;
	foreach (i; 0 .. e.length) e[i] = other.e[i];
	foreach (i; 0 .. T.length) T[i] = other.T[i];
	mu = other.mu;
	foreach (i; 0 .. k.length) k[i] = other.k[i];
	// D_AB
	sigma = other.sigma;
	foreach (i; 0 .. massf.length) massf[i] = other.massf[i];
	quality = other.quality;
    }

    @nogc void copy_average_values_from(ref const(GasState) gs0, ref const(GasState) gs1) 
    // Avoids memory allocation, it's all in place.
    {
	rho = 0.5 * (gs0.rho + gs1.rho);
	p = 0.5 * (gs0.p + gs1.p);
	p_e = 0.5 * (gs0.p_e + gs1.p_e);
	a = 0.5 * (gs0.a + gs1.a);
	foreach(i; 0 .. e.length) e[i] = 0.5 * (gs0.e[i] + gs1.e[i]);
	foreach(i; 0 .. T.length) T[i] = 0.5 * (gs0.T[i] + gs1.T[i]);
	mu = 0.5 * (gs0.mu + gs1.mu);
	foreach(i; 0 .. k.length) k[i] = 0.5 * (gs0.k[i] + gs1.k[i]);
	// D_AB
	sigma = 0.5 * (gs0.sigma * gs1.sigma);
	foreach(i; 0 .. massf.length) massf[i] = 0.5 * (gs0.massf[i] * gs1.massf[i]);
	quality = 0.5 * (gs0.quality + gs1.quality);
    }

    void copy_average_values_from(in GasState[] others, GasModel gm) 
    // Note that we must not send the current object in the others list as well.
    {
	size_t n = others.length;
	if (n == 0) throw new Error("Need to average from a nonempty array.");
	foreach(other; others) {
	    if ( this is other ) throw new Error("Must not include destination in source list.");
	}
	// Accumulate from a clean slate and then divide.
	p = 0.0;
	p_e = 0.0;
	foreach(ref elem; T) elem = 0.0;
	sigma = 0.0;
	foreach(ref elem; massf) elem = 0.0;
	quality = 0.0;
	foreach(other; others) {
	    p += other.p;
	    p_e += other.p_e;
	    foreach(i; 0 .. T.length) T[i] += other.T[i];
	    sigma += other.sigma;
	    foreach(i; 0 .. massf.length) massf[i] += other.massf[i];
	    quality += other.quality;
	}
	p /= n;
	p_e /= n;
	foreach(ref elem; T) elem /= n;
	sigma /= n;
	foreach(ref elem; massf) elem /= n;
	quality /= n;
	// Now, evaluate the rest of the properties using the gas model.
	gm.update_thermo_from_pT(this);
	gm.update_sound_speed(this);
	gm.update_trans_coeffs(this);
    } // end copy_average_values_from()

    @nogc bool check_values(bool print_message=true) const
    // Have commented out the print statements to ensure @nogc.
    {
	double RHOMIN = 0.0;
	double TMIN = 1.0;
	bool is_data_valid = true;
	if (!(isFinite(rho)) || rho < 1.01 * RHOMIN) {
	    // if (print_message) writeln("Density invalid: ", rho);
	    is_data_valid = false;
	}
	auto nmodes = e.length;
	foreach(imode; 0 .. nmodes) {
	    if (!isFinite(T[imode]) || T[imode] < 1.01 * TMIN) {
		// if (print_message) writeln("Temperature[", imode, "] invalid: ", T[imode]);
		is_data_valid = false;
	    }
	    if ( !isFinite(e[imode]) ) {
		// if (print_message) writeln("Energy[", imode, "] invalid: ", e[imode]);
		is_data_valid = false;
	    }
	}
	if (!isFinite(p)) {
	    // if (print_message) writeln("Total pressure invalid: ", p);
	    is_data_valid = false;
	}
	if (!isFinite(p_e)) {
	    // if (print_message) writeln("Electron pressure invalid: ", p_e);
	    is_data_valid = false;
	}
	if (!isFinite(a)) {
	    // if (print_message) writeln("Sound speed invalid: ", a);
	    is_data_valid = false;
	}
	double f_sum = 0.0; foreach(elem; massf) f_sum += elem;
	if (f_sum < 0.99 || f_sum > 1.01 || !isFinite(f_sum)) {
	    // if (print_message) writeln("Mass fraction sum bad: ", f_sum);
	    is_data_valid = false;
	}
	return is_data_valid;
    } // end check_values()

    override string toString() const
    {
	char[] repr;
	repr ~= "GasState(";
	repr ~= "rho=" ~ to!string(rho);
	repr ~= ", p=" ~ to!string(p);
	repr ~= ", p_e=" ~ to!string(p_e);
	repr ~= ", a=" ~ to!string(a);
	repr ~= ", T=" ~ to!string(T);
	repr ~= ", e=" ~ to!string(e);
	repr ~= ", mu=" ~ to!string(mu);
	repr ~= ", k=" ~ to!string(k);
	repr ~= ", massf=" ~ to!string(massf);
	repr ~= ", quality=" ~ to!string(quality);
	repr ~= ", sigma=" ~ to!string(sigma);
	repr ~= ")";
	return to!string(repr);
    }

/+ [TODO]
    double * copy_values_to_buffer(double *buf) const;
    double * copy_values_from_buffer(double *buf);
+/
} // end class GasState


@nogc void scale_mass_fractions(ref double[] massf, double tolerance=0.0)
{
    auto my_nsp = massf.length;
    if (my_nsp == 1) {
	// Single species, always expect massf[0]==1.0, so we can take a short-cut.
	assert(fabs(massf[0] - 1.0) < 0.1, "Single species mass fraction far from 1.0");
	massf[0] = 1.0;
    } else {
	// Multiple species, do the full job.
	double massf_sum = 0.0;
	foreach(isp; 0 .. my_nsp) {
	    massf[isp] = massf[isp] >= 0.0 ? massf[isp] : 0.0;
	    massf_sum += massf[isp];
	}
	assert(fabs(massf_sum - 1.0) < 0.1, "Sum of species mass fractions far from 1.0");
	if ( fabs(massf_sum - 1.0) > tolerance ) {
	    foreach(isp; 0 .. my_nsp) massf[isp] /= massf_sum;
	}
    }
    return;
} // end scale_mass_fractions()

@nogc pure double mass_average(in GasState Q, in double[] phi)
in {
    assert(Q.massf.length == phi.length, "Inconsistent array lengths.");
}
body {
    double result = 0.0;
    foreach ( i; 0..Q.massf.length ) result += Q.massf[i] * phi[i];
    return result;
}

@nogc pure double mole_average(in double[] molef, in double[] phi)
in {
    assert(molef.length == phi.length, "Inconsistent array lengths.");
}
body {
    double result = 0.0;
    foreach ( i; 0..molef.length ) result += molef[i] * phi[i];
    return result;
}

@nogc pure double mixture_molecular_mass(in double[] massf, in double[] mol_masses)
in {
    assert(massf.length == mol_masses.length, "Inconsistent array lengths.");
}
body {
    double M_inv = 0.0;
    foreach (i; 0 .. massf.length) M_inv += massf[i] / mol_masses[i];
    return 1.0/M_inv;
}

@nogc void massf2molef(in double[] massf, in double[] mol_masses, double[] molef)
in {
    assert(massf.length == mol_masses.length, "Inconsistent array lengths.");
    assert(massf.length == molef.length, "Inconsistent array lengths.");
}
body {
    double mixMolMass = mixture_molecular_mass(massf, mol_masses);
    foreach ( i; 0..massf.length ) molef[i] = massf[i] * mixMolMass / mol_masses[i];
}

@nogc void molef2massf(in double[] molef, in double[] mol_masses, double[] massf)
in {
    assert(massf.length == mol_masses.length, "Inconsistent array lengths.");
    assert(massf.length == molef.length, "Inconsistent array lengths.");
}
body {
    double mixMolMass = mole_average(molef, mol_masses);
    foreach ( i; 0..massf.length ) massf[i] = molef[i] * mol_masses[i] / mixMolMass;
}



unittest {
    auto gd = new GasState(2, 1);
    gd.massf[0] = 0.8;
    gd.massf[1] = 0.2;
    double[] phi = [9.0, 16.0];
    assert(approxEqual(10.4, mass_average(gd, phi)));
}
