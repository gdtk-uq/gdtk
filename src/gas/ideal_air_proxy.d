/**
 * ideal_air_proxy.d
 * Ideal-air gas model, implemented in fortran, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2016-12-27: initial cut, to explore mixed-language binding.
 */

module gas.ideal_air_proxy;

import gas.gas_model;
import gas.physical_constants;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import core.stdc.stdlib : exit;

extern(C) {
    void iaf_init();
    int iaf_n_species();
    int iaf_n_modes();
    double iaf_mol_mass(int i);
    void iaf_update_thermo_from_pT(double *p, double *Ttr, double *rho, double *u, double *massf);
    void iaf_update_thermo_from_rhou(double *p, double *Ttr, double *rho, double *u, double *massf);
    void iaf_update_thermo_from_rhoT(double *p, double *Ttr, double *rho, double *u, double *massf);
    void iaf_update_thermo_from_rhop(double *p, double *Ttr, double *rho, double *u, double *massf);
    void iaf_update_thermo_from_ps(double *p, double *Ttr, double *rho, double *u, double *massf,
				   double *s);
    void iaf_update_thermo_from_hs(double *p, double *Ttr, double *rho, double *u, double *massf,
				   double *h, double *s);
    void iaf_update_sound_speed(double *p, double *Ttr, double *rho, double *u, double *massf,
				double *a);
    void iaf_update_trans_coeffs(double *p, double *Ttr, double *rho, double *u, double *massf,
				 double *mu, double *k);
    double iaf_get_Cv();
    double iaf_get_Cp();
    double iaf_get_Rgas();
    double iaf_entropy(double *p, double *Ttr);
}

class IdealAirProxy: GasModel {
public:

    this() {
	// This proxy class delegates most tasks to the Fortran module.
	iaf_init();
	// but a few things are set in the D-domain, as well.
	_n_species = iaf_n_species();
	assert(_n_species == 1, "oops, wrong n_species");
	_n_modes = iaf_n_modes();
	assert(_n_modes == 0, "oops, wrong n_modes");
	_species_names.length = 1;
	_species_names[0] = "air";
	_mol_masses.length = 1;
	_mol_masses[0] = iaf_mol_mass(1); // Fortran array index starts at 1
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "IdealAirProxy =(";
	repr ~= "name=\"" ~ _species_names[0] ~"\"";
	repr ~= ", Mmass=" ~ to!string(_mol_masses[0]);
	// Should delegate the following to the Fortran domain
	// when we work out how to send strings. 
	// repr ~= ", gamma=" ~ to!string(_gamma);
	// repr ~= ", s1=" ~ to!string(_s1);
	// repr ~= ", T1=" ~ to!string(_T1);
	// repr ~= ", p1=" ~ to!string(_p1);
	// repr ~= ", constPrandtl=" ~ to!string(_constPrandtl);
	// repr ~= ", Prandtl=" ~ to!string(_Prandtl);
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	iaf_update_thermo_from_pT(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr);
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
	iaf_update_thermo_from_rhou(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr);
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	iaf_update_thermo_from_rhoT(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr);
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	iaf_update_thermo_from_rhop(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr);
    }
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	iaf_update_thermo_from_ps(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr, &s);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	iaf_update_thermo_from_hs(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr, &h, &s);
    }
    override void update_sound_speed(GasState Q) const
    {
	iaf_update_sound_speed(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr, &(Q.a));
    }
    override void update_trans_coeffs(GasState Q)
    {
	iaf_update_trans_coeffs(&(Q.p), &(Q.Ttr), &(Q.rho), &(Q.u), Q.massf.ptr, &(Q.mu), &(Q.k));
    }
    override double dudT_const_v(in GasState Q) const
    {
	return iaf_get_Cv(); // May need something more general for a more complex gas.
    }
    override double dhdT_const_p(in GasState Q) const
    {
	return iaf_get_Cp();
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	double R = iaf_get_Rgas();
	return R*Q.Ttr;
    }
    override double gas_constant(in GasState Q) const
    {
	return iaf_get_Rgas();
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.u;
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.u + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	double p = Q.p;
	double Ttr = Q.Ttr;
	return iaf_entropy(&p, &Ttr);
    }
} // end class IdealAirProxy

version(ideal_air_proxy_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	auto gm = new IdealAirProxy();
	auto gd = new GasState(1, 0);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	gd.massf[0] = 1.0;
	assert(approxEqual(gm.R(gd), 287.086, 1.0e-4), failedUnitTest());
	assert(gm.n_modes == 0, failedUnitTest());
	assert(gm.n_species == 1, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	assert(approxEqual(gd.rho, 1.16109, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.u, 215314.0, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.a, 347.241, 1.0e-4), failedUnitTest());
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 1.84691e-05, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0262449, 1.0e-6), failedUnitTest());

	return 0;
    }
}
