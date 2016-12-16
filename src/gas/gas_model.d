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
import std.string;
import util.lua;
import util.lua_service;

import util.msg_service;
import gas.physical_constants;

immutable double SMALL_MOLE_FRACTION = 1.0e-15;
immutable double MIN_MASS_FRACTION = 1.0e-30;
immutable double MIN_MOLES = 1.0e-30;
immutable double T_MIN = 20.0; 
immutable double MASSF_ERROR_TOL = 1.0e-6;

class GasModel {
public:
    @nogc @property uint n_species() const { return _n_species; }
    @nogc @property uint n_modes() const { return _n_modes; }
    @property ref double[] mol_masses() { return _mol_masses; }
    final string species_name(int i) const { return _species_names[i]; }
    final int species_index(string spName) const { return _species_indices.get(spName, -1); }

    void create_species_reverse_lookup()
    {
	foreach ( int isp; 0 .. _n_species ) {
	    _species_indices[_species_names[isp]] = isp;
	}
    }
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
    double enthalpy(in GasState Q, int isp)
    {
	// For the single-species gases, provide a default implementation
	// but we need to be careful to override this for multi-component gases.
	return enthalpy(Q);
    }
    abstract double entropy(in GasState Q);
    double entropy(in GasState Q, int isp)
    {
	// For the single-species gases, provide a default implementation
	// but we need to be carefule to override this for multi-component gases.
	return entropy(Q);
    }

    double gibbs_free_energy(GasState Q, int isp)
    {
	double h = enthalpy(Q, isp);
	double s = entropy(Q, isp);
	double g = h - Q.Ttr*s;
	return g;
    }
    
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
	assert(Q.massf.length == conc.length, brokenPreCondition("Inconsistent array lengths."));
    }
    body {
	foreach ( i; 0.._n_species ) {
	    conc[i] = Q.massf[i]*Q.rho / _mol_masses[i];
	    if ( conc[i] < MIN_MOLES ) conc[i] = 0.0;
	}
    }
    final void conc2massf(in double[] conc, GasState Q) 
    in {
	assert(Q.massf.length == conc.length, brokenPreCondition("Inconsisten array lengths."));
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
    int[string] _species_indices;
    double[] _mol_masses;
}
class GasState {
public:
    /// Thermodynamic properties.
    double rho;  /// density, kg/m**3
    double p;    /// presure, Pa
    double p_e;  /// electron pressure
    double a;    /// sound speed, m/s
    // For a gas in thermal equilibrium, all of the internal energies
    // are bundled together into u and are represented by a single
    // temperature Ttr.
    double Ttr;  /// thermal temperature, K
    double u;    /// specific thermal energy, J/kg
    // For a gas in thermal nonequilibrium, the internal energies are
    // stored unbundled, with u being the trans-rotational thermal energy.
    // The array length will be determined by the specific model and,
    // to get the total internal energy,
    // the gas-dynamics part of the code will need to sum the array elements.
    double[] e_modes;  /// specific internal energies for thermal nonequilibrium model, J/kg
    double[] T_modes;  /// temperatures for internal energies for thermal nonequilibrium, K
    /// Transport properties
    double mu;   /// viscosity, Pa.s
    double kth;  /// thermal conductivity for a single temperature gas, W/(m.K)
    double[] k_modes;  /// thermal conductivities for the nonequilibrium model, W/(m.K)
    double sigma;    /// electrical conductivity, S/m
    /// Composition
    double[] massf;  /// species mass fractions
    double quality;  /// vapour quality

    this(uint n_species, uint n_modes)
    {
	massf.length = n_species;
	e_modes.length = n_modes;
	T_modes.length = n_modes;
	k_modes.length = n_modes;
    }

    this(GasModel gm)
    {
	this(gm.n_species, gm.n_modes);
    }

    this(GasModel gm, in double p_init, in double Ttr_init, in double[] T_modes_init,
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in double sigma_init=0.0)
    {
	p = p_init;
	p_e = p_init;
	Ttr = Ttr_init;
	T_modes.length = gm.n_modes;
	foreach(i; 0 .. gm.n_modes) { T_modes[i] = T_modes_init[i]; }
	e_modes.length = gm.n_modes;
	k_modes.length = gm.n_modes;
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
	double[] T_modes;
	T_modes.length = gm.n_modes;
	foreach(ref Tmode; T_modes) { Tmode = T_init; }
	this(gm, p_init, T_init, T_modes, massf_init, quality_init, sigma_init);
    }

    this(in GasState other) 
    {
	rho = other.rho;
	p = other.p;
	p_e = other.p_e;
	Ttr = other.Ttr;
	u = other.u;
	a = other.a;
	e_modes = other.e_modes.dup;
	T_modes = other.T_modes.dup;
	mu = other.mu;
	kth = other.kth;
	k_modes = other.k_modes.dup;
	sigma = other.sigma;
	massf = other.massf.dup;
	quality = other.quality;
    }

    @nogc void copy_values_from(ref const(GasState) other) 
    {
	rho = other.rho;
	p = other.p;
	Ttr = other.Ttr;
	u = other.u;
	p_e = other.p_e;
	a = other.a;
	foreach (i; 0 .. e_modes.length) { e_modes[i] = other.e_modes[i]; }
	foreach (i; 0 .. T_modes.length) { T_modes[i] = other.T_modes[i]; }
	mu = other.mu;
	kth = other.kth;
	foreach (i; 0 .. k_modes.length) { k_modes[i] = other.k_modes[i]; }
	sigma = other.sigma;
	foreach (i; 0 .. massf.length) { massf[i] = other.massf[i]; }
	quality = other.quality;
    }

    @nogc void copy_average_values_from(ref const(GasState) gs0, ref const(GasState) gs1) 
    // Avoids memory allocation, it's all in place.
    {
	rho = 0.5 * (gs0.rho + gs1.rho);
	p = 0.5 * (gs0.p + gs1.p);
	Ttr = 0.5 * (gs0.Ttr + gs1.Ttr);
	u = 0.5 * (gs0.u + gs1.u);
	p_e = 0.5 * (gs0.p_e + gs1.p_e);
	a = 0.5 * (gs0.a + gs1.a);
	foreach(i; 0 .. e_modes.length) { e_modes[i] = 0.5 * (gs0.e_modes[i] + gs1.e_modes[i]); }
	foreach(i; 0 .. T_modes.length) { T_modes[i] = 0.5 * (gs0.T_modes[i] + gs1.T_modes[i]); }
	mu = 0.5 * (gs0.mu + gs1.mu);
	kth = 0.5 * (gs0.kth + gs1.kth);
	foreach(i; 0 .. k_modes.length) { k_modes[i] = 0.5 * (gs0.k_modes[i] + gs1.k_modes[i]); }
	sigma = 0.5 * (gs0.sigma + gs1.sigma);
	foreach(i; 0 .. massf.length) massf[i] = 0.5 * (gs0.massf[i] + gs1.massf[i]);
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
	Ttr = 0.0;
	u = 0.0;
	p_e = 0.0;
	foreach(ref elem; T_modes) { elem = 0.0; }
	sigma = 0.0;
	foreach(ref elem; massf) { elem = 0.0; }
	quality = 0.0;
	foreach(other; others) {
	    p += other.p;
	    Ttr += other.Ttr;
	    u += other.u;
	    p_e += other.p_e;
	    foreach(i; 0 .. T_modes.length) { T_modes[i] += other.T_modes[i]; }
	    sigma += other.sigma;
	    foreach(i; 0 .. massf.length) { massf[i] += other.massf[i]; }
	    quality += other.quality;
	}
	p /= n;
	Ttr /= n;
	u /= n;
	p_e /= n;
	foreach(ref elem; T_modes) { elem /= n; }
	sigma /= n;
	foreach(ref elem; massf) { elem /= n; }
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
	if (!isFinite(Ttr) || Ttr < 1.01 * TMIN) {
	    // if (print_message) writeln("Transrotational temperature invalid: ", Ttr);
	    is_data_valid = false;
	}
	auto nmodes = e_modes.length;
	foreach(imode; 0 .. nmodes) {
	    if (!isFinite(T_modes[imode]) || T_modes[imode] < 1.01 * TMIN) {
		// if (print_message) writeln("Temperature[", imode, "] invalid: ", T_modes[imode]);
		is_data_valid = false;
	    }
	    if ( !isFinite(e_modes[imode]) ) {
		// if (print_message) writeln("Energy[", imode, "] invalid: ", e_modes[imode]);
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
	repr ~= ", Ttr=" ~ to!string(Ttr);
	repr ~= ", u=" ~ to!string(u);
	repr ~= ", p_e=" ~ to!string(p_e);
	repr ~= ", a=" ~ to!string(a);
	repr ~= ", T_modes=" ~ to!string(T_modes);
	repr ~= ", e_modes=" ~ to!string(e_modes);
	repr ~= ", mu=" ~ to!string(mu);
	repr ~= ", kth=" ~ to!string(kth);
	repr ~= ", k_modes=" ~ to!string(k_modes);
	repr ~= ", massf=" ~ to!string(massf);
	repr ~= ", quality=" ~ to!string(quality);
	repr ~= ", sigma=" ~ to!string(sigma);
	repr ~= ")";
	return to!string(repr);
    }
} // end class GasState


@nogc void scale_mass_fractions(ref double[] massf, double tolerance=0.0,
				double assert_error_tolerance=0.1)
{
    auto my_nsp = massf.length;
    if (my_nsp == 1) {
	// Single species, always expect massf[0]==1.0, so we can take a short-cut.
	assert(fabs(massf[0] - 1.0) < assert_error_tolerance,
	       "Single species mass fraction far from 1.0");
	massf[0] = 1.0;
    } else {
	// Multiple species, do the full job.
	double massf_sum = 0.0;
	foreach(isp; 0 .. my_nsp) {
	    massf[isp] = massf[isp] >= 0.0 ? massf[isp] : 0.0;
	    massf_sum += massf[isp];
	}
	assert(fabs(massf_sum - 1.0) < assert_error_tolerance,
	       "Sum of species mass fractions far from 1.0");
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


/* The following functions:
   update_thermo_state_pT(), update_thermo_state_rhoT(), update_thermo_state_rhop() 
   are for updating the thermo state from when  the gas model does not  have a method
   for doing so with those variables, but does have a defined method for
   update_thermo_from_rhoe(). (e.g. in the UniformLUT class)
   A guess is made for rho & e and that guess is iterated using the  Newton-Raphson 
   method.
   A GasModel object is a function parameter so that the update method from rho,e for
   that gas model can be called.
   The funcions:
   update_thermo_state_ps(), and update_thermo_state_hs() actually iterate on the update
   update_thermo_from_pT(), as  called from the gas model (though the p,T update can be
   defined as the function defined here that itself iterates on the update method for
   rho,e)
*/

immutable MAX_RELATIVE_STEP = 0.1;
immutable MAX_STEPS = 20;

void update_thermo_state_pT(GasModel gmodel, GasState Q)
{
    double drho, rho_old, rho_new, e_old, e_new, de;
    double drho_sign, de_sign;
    double Cv_eff, R_eff, T_old;
    double fp_old, fT_old, fp_new, fT_new;
    double dfp_drho, dfT_drho, dfp_de, dfT_de, det;
    int converged, count;

    double p_given = Q.p;
    double T_given = Q.Ttr;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
 
    double fT_tol = 1.0e-6 * T_given;
    double fp_tol = 1.0e-6 * p_given;
    double fp_tol_fail = 0.02 * p_given;
    double fT_tol_fail = 0.02 * T_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
    Q.rho = 1.0; // kg/m**3 
    Q.u = 2.0e5; // J/kg 
    gmodel.update_thermo_from_rhoe(Q);
    
    T_old = Q.Ttr;
    R_eff = Q.p / (Q.rho * T_old);
    de = 0.01 * Q.u;
    Q.u += de;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 1 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

    Cv_eff = de / (Q.Ttr - T_old);
    // Now, get a better guess for the appropriate density and
    // internal energy.
    e_old = Q.u + (T_given - Q.Ttr) * Cv_eff;
    rho_old = p_given / (R_eff * T_given);

    // Evaluate state variables using this guess.
    Q.rho = rho_old;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 2 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

    fp_old = p_given - Q.p;
    fT_old = T_given - Q.Ttr;
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fp_old) < fp_tol) && (fabs(fT_old) < fT_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	// Perturb first dimension to get derivatives.
	rho_new = rho_old * 1.001;
	e_new = e_old;
	Q.rho = rho_new;
	Q.u = e_new;
	try { gmodel.update_thermo_from_rhoe(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s failed at call A in %s\n", count, __FUNCTION__); 
	    msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	fp_new = p_given - Q.p;
	fT_new = T_given - Q.Ttr;
	dfp_drho = (fp_new - fp_old) / (rho_new - rho_old);
	dfT_drho = (fT_new - fT_old) / (rho_new - rho_old);
	// Perturb other dimension to get derivatives.
	rho_new = rho_old;
	e_new = e_old * 1.001;
	Q.rho = rho_new;
	Q.u = e_new;

	try { gmodel.update_thermo_from_rhoe(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s failed at call B in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}

	fp_new = p_given - Q.p;
	fT_new = T_given - Q.Ttr;
	dfp_de = (fp_new - fp_old) / (e_new - e_old);
	dfT_de = (fT_new - fT_old) / (e_new - e_old);
	det = dfp_drho * dfT_de - dfT_drho * dfp_de;
       	if( fabs(det) < 1.0e-12 ) {
	    string msg;
	    msg ~= format("Error in function %s\n", __FUNCTION__);
	    msg ~= format("    Nearly zero determinant, det = ", det);
	    throw new Exception(msg);
	}
	drho = (-dfT_de * fp_old + dfp_de * fT_old) / det;
	de = (dfT_drho * fp_old - dfp_drho * fT_old) / det;
	if( fabs(drho) > MAX_RELATIVE_STEP * rho_old ) {
	    // move a little toward the goal 
	    drho_sign = (drho > 0.0 ? 1.0 : -1.0);
	    drho = drho_sign * MAX_RELATIVE_STEP * rho_old;
	} 
	if( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
	    // move a little toward the goal
	    de_sign = (de > 0.0 ? 1.0 : -1.0);
	    de = de_sign * MAX_RELATIVE_STEP * e_old;
	} 
	rho_old += drho;
	e_old += de;
	// Make sure of consistent thermo state.
	Q.rho = rho_old;
	Q.u = e_old;
	try { gmodel.update_thermo_from_rhoe(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s failed in %s\n", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	// Prepare for next iteration.
	fp_old = p_given - Q.p;
	fT_old = T_given - Q.Ttr;
	converged = (fabs(fp_old) < fp_tol) && (fabs(fT_old) < fT_tol);
	++count;
    } // end while 

    if ( count >= MAX_STEPS ) {
	string msg;
	msg ~= format("Warning in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations did not converge.\n");
	msg ~= format("    fp_old = %g, fT_old = %g\n", fp_old, fT_old);
	msg ~= format("    p_given = %.10s, T_given, %.5s\n", p_given, T_given); 
	msg ~= "  Supplied Q:" ~ Q.toString();
	writeln(msg);

    }

    if( (fabs(fp_old) > fp_tol_fail) || (fabs(fT_old) > fT_tol_fail) ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations failed badly.\n");
	msg ~= format("    p_given = %.10s, T_given, %.5s\n", p_given, T_given); 
	msg ~= format("    fp_old = %g, fT_old = %g\n", fp_old, fT_old);
	msg ~= "  Supplied Q:" ~ Q.toString();
	throw new Exception(msg);
    }
}
   


void update_thermo_state_rhoT(GasModel gmodel, GasState Q)  
{
    // This method can be applied to single-species models only
    double e_old, e_new, de, tmp, de_sign;
    double Cv_eff, T_old;
    double dfT_de, fT_old, fT_new;
    int converged, count;

    double rho_given = Q.rho;
    double T_given = Q.Ttr;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fT_tol = 1.0e-6 * T_given;
    double fT_tol_fail = 0.02 * T_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
       
    Q.rho = rho_given; // kg/m**3 
    Q.u = 2.0e5; // J/kg 
    gmodel.update_thermo_from_rhoe(Q);

    T_old = Q.Ttr;
    de = 0.01 * Q.u;
    Q.u += de;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 0 failed in %s", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }
   
   
    Cv_eff = de / (Q.Ttr - T_old);
    // Now, get a better guess for the appropriate density and internal energy.
    e_old = Q.u + (T_given - Q.Ttr) * Cv_eff;
    // Evaluate state variables using this guess.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
	string msg;
	msg ~= format("Starting guess at iteration 1 failed in %s", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }
   
    fT_old = T_given - Q.Ttr;
    // Perturb to get derivative.
    e_new = e_old * 1.001;
    Q.rho = rho_given;
    Q.u = e_new;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
	string msg;
	msg ~= format("Starting guess at iteration 2 failed in %s", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

    fT_new = T_given - Q.Ttr;
    dfT_de = (fT_new - fT_old) / (e_new - e_old);

    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fT_new) < fabs(fT_old) ) {
	tmp = fT_new; fT_new = fT_old; fT_old = tmp;
	tmp = e_new; e_new = e_old; e_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fT_old) < fT_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	de = -fT_old / dfT_de;
	if ( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
	    // move a little toward the goal 
	    de_sign = (de > 0.0 ? 1.0 : -1.0);
	    de = de_sign * MAX_RELATIVE_STEP * fabs(e_old);
	} 
	e_new = e_old + de;
	Q.rho = rho_given;
	Q.u = e_new;
	try { gmodel.update_thermo_from_rhoe(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s failed in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	fT_new = T_given - Q.Ttr;
	dfT_de = (fT_new - fT_old) / (e_new - e_old);
	// Prepare for the next iteration.
	++count;
	fT_old = fT_new;
	e_old = e_new;
	converged = fabs(fT_old) < fT_tol;
    }   // end while 
    // Ensure that we have the current data for all EOS variables.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
	string msg;
	msg ~= format("Function %s failed after finishing iterations", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

     if ( count >= MAX_STEPS ) {
	string msg;
	msg ~= format("Warning in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations did not converge.\n");
	msg ~= format("    fT_old = %g\n", fT_old);
	msg ~= format("    rho_given = %.5s, T_given, %.5s\n", rho_given, T_given); 
        msg ~= "  Supplied Q:" ~ Q.toString;
	writeln(msg);

    }
    if ( fabs(fT_old) > fT_tol_fail ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations failed badly.\n");
	msg ~= format("    rho_given = %.5s, T_given, %.5s\n", rho_given, T_given); 
	msg ~= "  Supplied Q:" ~ Q.toString();
	throw new Exception(msg);
    }  
}

void update_thermo_state_rhop(GasModel gmodel, GasState Q)
{
    double e_old, e_new, de, dedp, tmp, de_sign;
    double p_old;
    double dfp_de, fp_old, fp_new;
    int converged, count;

    double rho_given = Q.rho;
    double p_given = Q.p;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fp_tol = 1.0e-6 * p_given;
    double fp_tol_fail = 0.02 * p_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
    Q.rho = rho_given; // kg/m**3
    Q.u = 2.0e5; // J/kg 
    gmodel.update_thermo_from_rhoe(Q);
    p_old = Q.p;
    de = 0.01 * Q.u;
    Q.u += de;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 0 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

    dedp = de / (Q.p - p_old);
    // Now, get a better guess for the appropriate internal energy.
    e_old = Q.u + (p_given - Q.p) * dedp;
    //     printf( "Initial guess e_old= %g dedp= %g\n", e_old, dedp );
    // Evaluate state variables using this guess.
    Q.rho = rho_given;
    Q.u = e_old;


    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 1 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

    fp_old = p_given - Q.p;
    // Perturb to get derivative.
    e_new = e_old * 1.001;
    Q.rho = rho_given;
    Q.u = e_new;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 2 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

    fp_new = p_given - Q.p;
    dfp_de = (fp_new - fp_old) / (e_new - e_old);

    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fp_new) < fabs(fp_old) ) {
	tmp = fp_new; fp_new = fp_old; fp_old = tmp;
	tmp = e_new; e_new = e_old; e_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fp_old) < fp_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	de = -fp_old / dfp_de;
	if ( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
	    // move a little toward the goal
	    de_sign = (de > 0.0 ? 1.0 : -1.0);
	    de = de_sign * MAX_RELATIVE_STEP * fabs(e_old);
	} 
	e_new = e_old + de;
	Q.rho = rho_given;
	Q.u = e_new;

	try { gmodel.update_thermo_from_rhoe(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s failed in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}

	fp_new = p_given - Q.p;
	dfp_de = (fp_new - fp_old) / (e_new - e_old);
	// Prepare for next iteration.
	++count;
	fp_old = fp_new;
	e_old = e_new;
	converged = fabs(fp_old) < fp_tol;
    }   // end while 
    // Ensure that we have the current data for all EOS variables.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhoe(Q); }
    catch (Exception caughtException) {
	string msg;
	msg ~= format("Function %s failed after finishing iterations", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_rhoe() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }

      if ( count >= MAX_STEPS ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations did not converge.\n");
	msg ~= format("    fp_old = %g, e_old = %g\n", fp_old, e_old);
	msg ~= format("    rho_given = %.5s, p_given, %.8s\n", rho_given, p_given); 
        msg ~= "  Supplied Q:" ~ Q.toString;
	writeln(msg);
    }

    if ( fabs(fp_old) > fp_tol_fail ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations failed badly.\n");
	msg ~= format("    rho_given = %.5s, T_given, %.8s\n", rho_given, p_given); 
	msg ~= "  Supplied Q:" ~ Q.toString();
	throw new Exception(msg);
    }
}


void update_thermo_state_ps(GasModel gmodel, GasState Q, double s) 
{
    double T_old, T_new, dT, tmp, dT_sign;
    double dfs_dT, fs_old, fs_new;
    int converged, count;

    double s_given = s;
    double p_given = Q.p;
   
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fs_tol = 1.0e-6 * s_given;
    double fs_tol_fail = 0.02 * s_given;

    // Guess the thermo state assuming that T is a good guess.
    T_old = Q.Ttr;
    try { gmodel.update_thermo_from_pT(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 0 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }
    ////**** Need to check this is the correct method - is called 2 more times*****/////
    double s_old = gmodel.entropy(Q);   
    fs_old = s_given - s_old;
    // Perturb T to get a derivative estimate
    T_new = T_old * 1.001;
    Q.Ttr = T_new;

    try { gmodel.update_thermo_from_pT(Q); }
    catch (Exception caughtException) {
  	string msg;
	msg ~= format("Starting guess at iteration 1 failed in %s\n", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }
    double s_new = gmodel.entropy(Q);
    fs_new = s_given - s_new;
    dfs_dT = (fs_new - fs_old)/(T_new - T_old);
    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fs_new) < fabs(fs_old) ) {
	tmp = fs_new; fs_new = fs_old; fs_old = tmp;
	tmp = s_new; s_new = s_old; s_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fs_old) < fs_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	dT = -fs_old / dfs_dT;
	if ( fabs(dT) > MAX_RELATIVE_STEP * T_old ) {
	    // move a little toward the goal
	    dT_sign = (dT > 0.0 ? 1.0 : -1.0);
	    dT = dT_sign * MAX_RELATIVE_STEP * fabs(T_old);
	} 
	T_new = T_old + dT;
	Q.Ttr = T_new;
	try { gmodel.update_thermo_from_pT(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s failed in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	s_new = gmodel.entropy(Q);
	fs_new = s_given - s_new;
	dfs_dT = (fs_new - fs_old) / (T_new - T_old);
	// Prepare for next iteration.
	++count;
	fs_old = fs_new;
	T_old = T_new;
	converged = (fabs(fs_old) < fs_tol);
    }   // end while 
    // Ensure that we have the current data for all EOS variables.
    Q.Ttr = T_old;


    try { gmodel.update_thermo_from_pT(Q); }
    catch (Exception caughtException) {
	string msg;
	msg ~= format("Function %s failed after finishing iterations", __FUNCTION__);
	msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	msg ~= to!string(caughtException);
	throw new Exception(msg);
    }
    if ( count >= MAX_STEPS ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations did not converge.\n");
	msg ~= format("    fs_old = %g\n", fs_old);
	msg ~= format("    p_given = %.8s, s_given, %.5s\n", p_given, s_given); 
        msg ~= "  Supplied Q:" ~ Q.toString;
	writeln(msg);
    }

    if ( fabs(fs_old) > fs_tol_fail ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations failed badly.\n");
	msg ~= format("    p_given = %.8s, s_given, %.5s\n", p_given, s_given); 
	msg ~= "  Supplied Q:" ~ Q.toString();
	throw new Exception(msg);
    }
}

void update_thermo_state_hs(GasModel gmodel, GasState Q, double h, double s)
{
double dp, p_old, p_new, T_old, T_new, dT;
    double dp_sign, dT_sign;
    double fh_old, fs_old, fh_new, fs_new;
    double dfh_dp, dfs_dp, dfh_dT, dfs_dT, det;
    int converged, count;

    double h_given = h;
    double s_given = s;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fh_tol = 1.0e-6 * h_given;
    double fs_tol = 1.0e-6 * s_given;
    double fh_tol_fail = 0.02 * h_given;
    double fs_tol_fail = 0.02 * s_given;

    // Use current gas state as guess
    p_old = Q.p;
    T_old = Q.Ttr;
    double h_new = gmodel.enthalpy(Q);
    double s_new = gmodel.entropy(Q);
    fh_old = h_given - h_new;
    fs_old = s_given - s_new;

    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fh_old) < fh_tol) && (fabs(fs_old) < fs_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	// Perturb first dimension to get derivatives.
	p_new = p_old * 1.001;
	T_new = T_old;
	Q.p = p_new;
	Q.Ttr = T_new;
	try { gmodel.update_thermo_from_pT(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s at call A failed in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	h_new = gmodel.enthalpy(Q);
	s_new = gmodel.entropy(Q);
	fh_new = h_given - h_new;
	fs_new = s_given - s_new;
	dfh_dp = (fh_new - fh_old) / (p_new - p_old);
	dfs_dp = (fs_new - fs_old) / (p_new - p_old);
	// Perturb other dimension to get derivatives.
	p_new = p_old;
	T_new = T_old * 1.001;
	Q.p = p_new;
	Q.Ttr = T_new;
	try { gmodel.update_thermo_from_pT(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s at call B failed in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	h_new = gmodel.enthalpy(Q);
	s_new = gmodel.entropy(Q);
	fh_new = h_given - h_new;
	fs_new = s_given - s_new;
	dfh_dT = (fh_new - fh_old) / (T_new - T_old);
	dfs_dT = (fs_new - fs_old) / (T_new - T_old);

	det = dfh_dp * dfs_dT - dfs_dp * dfh_dT;
      
       	if( fabs(det) < 1.0e-12 ) {
	    string msg;
	    msg ~= format("Error in function %s\n", __FUNCTION__);
	    msg ~= format("    Nearly zero determinant, det = ", det);
	    throw new Exception(msg);
	}
	dp = (-dfs_dT * fh_old + dfh_dT * fs_old) / det;
	dT = (dfs_dp * fh_old - dfh_dp * fs_old) / det;
	if( fabs(dp) > MAX_RELATIVE_STEP * p_old ) {
	    // move a little toward the goal 
	    dp_sign = (dp > 0.0 ? 1.0 : -1.0);
	    dp = dp_sign * MAX_RELATIVE_STEP * p_old;
	} 
	if( fabs(dT) > MAX_RELATIVE_STEP * T_old ) {
	    // move a little toward the goal
	    dT_sign = (dT > 0.0 ? 1.0 : -1.0);
	    dT = dT_sign * MAX_RELATIVE_STEP * T_old;
	} 
	p_old += dp;
	T_old += dT;
	// Make sure of consistent thermo state.
	Q.p = p_old;
	Q.Ttr = T_old;
	try { gmodel.update_thermo_from_pT(Q); }
	catch (Exception caughtException) {
	    string msg;
	    msg ~= format("Iteration %s at call C failed in %", count, __FUNCTION__);
	    msg ~= format("Excpetion message from update_thermo_from_pT() was:\n\n");
	    msg ~= to!string(caughtException);
	    throw new Exception(msg);
	}
	h_new = gmodel.enthalpy(Q);
	s_new = gmodel.entropy(Q);
	// Prepare for next iteration.
	fh_old = h_given - h_new;
	fs_old = s_given - s_new;
	converged = (fabs(fh_old) < fh_tol) && (fabs(fs_old) < fs_tol);
	++count;
    } // end while 

    if ( count >= MAX_STEPS ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations did not converge.\n");
	msg ~= format("    fh_old = %g, fs_old = %g\n", fh_old, fs_old);
	msg ~= format("    h_given = %.10s, h_given, %.5s\n", h_given, s_given); 
	msg ~= "  Supplied Q:" ~ Q.toString();
	writeln(msg);
    }

    if( (fabs(fh_old) > fh_tol_fail) || (fabs(fs_old) > fs_tol_fail) ) {
	string msg;
	msg ~= format("Error in function: %s:\n", __FUNCTION__);
	msg ~= format("    Iterations failed badly.\n");
	msg ~= format("    h_given = %.10s, h_given, %.5s\n", h_given, s_given); 	
	msg ~= "  Supplied Q:" ~ Q.toString();
	throw new Exception(msg);
    }

} // end update_thermo_state_hs()


// Utility function to construct specific gas models needs to know about
// all of the gas-model modules.
import gas.ideal_gas;
import gas.therm_perf_gas;
import gas.very_viscous_air;
import gas.co2gas;
import gas.co2gas_sw;
import gas.sf6virial;
import gas.uniform_lut;
import gas.adaptive_lut_CEA;
import core.stdc.stdlib : exit;


GasModel init_gas_model(in string file_name="gas-model.lua")
/**
 * Get the instructions for setting up the GasModel object from a Lua file.
 * The first item in the file should be a model name which we use to select 
 * the specific GasModel class.
 * The constructor for each specific gas model will know how to pick out the
 * specific parameters of interest.
 * As new GasModel classes are added to the collection, just 
 * add a new case to the switch statement below.
 */
{
    lua_State* L;
   
    try { 
        L = init_lua_State(file_name);
    } catch (Exception e) {
        writeln("ERROR: in function init_gas_model() in gas_model.d");
        writeln("ERROR: There was a problem parsing the input file: ", file_name);
	writeln("ERROR: There could be a Lua syntax error OR the file might not exist.");
	writeln("ERROR: Quitting at this point.");
 	exit(1);
    }
    string gas_model_name;
    try {
    	gas_model_name = getString(L, LUA_GLOBALSINDEX, "model");
    } catch (Exception e) {
        writeln("ERROR: in function init_gas_model() in gas_model.d");
        writeln("ERROR: There was a problem reading the 'model' name" );
	writeln("ERROR: in the gas model input Lua file.");
	writeln("ERROR: Quitting at this point.");
        exit(1);
    }
    GasModel gm;
    switch ( gas_model_name ) {
    case "IdealGas":
	gm = new IdealGas(L);
	break;
    case "ThermallyPerfectGas":
	gm = new ThermallyPerfectGas(L);
	break;
    case "VeryViscousAir":
	gm = new VeryViscousAir(L);
	break;
    case "CO2Gas":
	gm = new CO2Gas(L);
	break;
    case "CO2GasSW":
	gm = new CO2GasSW(L);
	break;
    case "SF6Virial":
	gm = new SF6Virial(L);
	break;
    case "look-up table":  
	gm = new  UniformLUT(L);
	break;
    case "CEA adaptive look-up table":
	gm = new AdaptiveLUT(L);
	break;
    default:
	string errMsg = format("The gas model '%s' is not available.", gas_model_name);
	throw new Error(errMsg);
    }
    lua_close(L);
    return gm;
}

//---------------------------------------------------------------------------------

version(gas_model_test) {
    int main() {
	// Methods for testing gas state class
	auto gd = new GasState(2, 1);
	gd.massf[0] = 0.8;
	gd.massf[1] = 0.2;
	double[] phi = [9.0, 16.0];
	assert(approxEqual(10.4, mass_average(gd, phi), 1.0e-6));
	
	// Iterative methods test using idealgas single species model
	// These assume IdealGas class is working properly
	GasModel gm;
	try {
	    gm = init_gas_model("sample-data/ideal-air-gas-model.lua");
	}
	catch (Exception e) {
	    writeln(e.msg);
	    string msg;
	    msg ~= "Test of iterative methods in gas_model.d require file:";
	    msg ~= " ideal-air-gas-model.lua in directory: gas/sample_data";
	    throw new Exception(msg);
	}

	gd = new GasState(gm, 100.0e3, 300.0);
	assert(approxEqual(gm.R(gd), 287.086, 1.0e-4), "gas constant");
	assert(gm.n_modes == 0, "number of energy modes");
	assert(gm.n_species == 1, "number of species");
	assert(approxEqual(gd.p, 1.0e5), "pressure");
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), "static temperature");
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), "massf[0]");

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	assert(approxEqual(gd.rho, 1.16109, 1.0e-4), "density");
	assert(approxEqual(gd.u, 215314.0, 1.0e-4), "internal energy");
	assert(approxEqual(gd.a, 347.241, 1.0e-4), "sound speed");
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 1.84691e-05, 1.0e-6), "viscosity");
	assert(approxEqual(gd.kth, 0.0262449, 1.0e-6), "conductivity");

	// Select arbitrary energy and density and establish a set of 
	// variables that are thermodynamically consistent
	double e_given = 1.0e7;
	double rho_given = 2.0;
	auto Q = new GasState(gm);
	Q.u = e_given;
	Q.rho = rho_given;
	gm.update_thermo_from_rhoe(Q);
	double p_given = Q.p;
	double T_given = Q.Ttr;
	
	// Initialise the same state from the different property combinations
	// Test pT iterative update
	Q.p = p_given;
	Q.Ttr = T_given;
     	update_thermo_state_pT(gm, Q); 
	// Determine correct entropy/enthalpy for updates that use them
	double s_given = gm.entropy(Q); 
	double h_given = gm.enthalpy(Q);
	assert(approxEqual(Q.rho, rho_given, 1.0e-6),  failedUnitTest());
	assert(approxEqual(Q.u, e_given, 1.0e-6), failedUnitTest());
	// Test rhoT iterative update
	Q.rho = rho_given;
	Q.Ttr = T_given;
	update_thermo_state_rhoT(gm, Q);
	assert(approxEqual(Q.u, e_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.p, p_given, 1.0e-6),  failedUnitTest());
	// Test rhop iterative update
	Q.rho = rho_given;
	Q.p = p_given;
	assert(approxEqual(Q.Ttr, T_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.u, e_given, 1.0e-6), failedUnitTest());
	// Test  ps iterative update
  	Q.p = p_given;
	update_thermo_state_ps(gm, Q, s_given);	
	assert(approxEqual(Q.Ttr, T_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.u, e_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.rho, rho_given, 1.0e-6), failedUnitTest());
	// Test hs iterative update
	assert(approxEqual(Q.Ttr, T_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.u, e_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.rho, rho_given, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.p, p_given, 1.0e-6), failedUnitTest());

	// Try 
	return 0;
    }
}

