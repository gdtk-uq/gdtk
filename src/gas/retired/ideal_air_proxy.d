/**
 * ideal_air_proxy.d
 * Ideal-air gas model, implemented in fortran, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version:
 * 2016-12-27: initial cut, to explore mixed-language binding.
 * 2018-06-03: Adapted to accept complex numbers but not pass them through.
 *             We assume that the Fortran side is only coded for double numbers.
 */

module gas.ideal_air_proxy;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;

extern(C) {
    void iaf_init();
    @nogc int iaf_n_species();
    @nogc int iaf_n_modes();
    @nogc double iaf_mol_mass(int i);
    @nogc void iaf_update_thermo_from_pT(double *p, double *T, double *rho, double *u, double *massf);
    @nogc void iaf_update_thermo_from_rhou(double *p, double *T, double *rho, double *u, double *massf);
    @nogc void iaf_update_thermo_from_rhoT(double *p, double *T, double *rho, double *u, double *massf);
    @nogc void iaf_update_thermo_from_rhop(double *p, double *T, double *rho, double *u, double *massf);
    @nogc void iaf_update_thermo_from_ps(double *p, double *T, double *rho, double *u, double *massf,
                                         double *s);
    @nogc void iaf_update_thermo_from_hs(double *p, double *T, double *rho, double *u, double *massf,
                                         double *h, double *s);
    @nogc void iaf_update_sound_speed(double *p, double *T, double *rho, double *u, double *massf,
                                      double *a);
    @nogc void iaf_update_trans_coeffs(double *p, double *T, double *rho, double *u, double *massf,
                                       double *mu, double *k);
    @nogc double iaf_get_Cv();
    @nogc double iaf_get_Cp();
    @nogc double iaf_get_Rgas();
    @nogc double iaf_entropy(double *p, double *T);
}

class IdealAirProxy: GasModel {
public:

    this() {
        type_str = "IdealAirProxy";
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
        version(complex_numbers) {
            throw new Error("Do not use with complex numbers.");
        }
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

    version(complex_numbers) {
        override void update_thermo_from_pT(ref GasState Q) const
        {
            double p = Q.p.re;
            double T = Q.T.re;
            double rho = Q.rho.re;
            double u = Q.u.re;
            double[1] massf;
            iaf_update_thermo_from_pT(&p, &T, &rho, &u, massf.ptr);
            Q.p = to!number(p);
            Q.T = to!number(T);
            Q.rho = to!number(rho);
            Q.u = to!number(u);
            Q.massf[0] = to!number(massf[0]);
        }
        override void update_thermo_from_rhou(ref GasState Q) const
        {
            // [TODO] iaf_update_thermo_from_rhou(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_rhoT(ref GasState Q) const
        {
            // [TODO] iaf_update_thermo_from_rhoT(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_rhop(ref GasState Q) const
        {
            // [TODO] iaf_update_thermo_from_rhop(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_ps(ref GasState Q, number s) const
        {
            // [TODO] iaf_update_thermo_from_ps(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &s);
        }
        override void update_thermo_from_hs(ref GasState Q, number h, number s) const
        {
            // [TODO] iaf_update_thermo_from_hs(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &h, &s);
        }
        override void update_sound_speed(ref GasState Q) const
        {
            // [TODO] iaf_update_sound_speed(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &(Q.a));
        }
        override void update_trans_coeffs(ref GasState Q)
        {
            // [TODO] iaf_update_trans_coeffs(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &(Q.mu), &(Q.k));
        }
        override number dudT_const_v(in GasState Q) const
        {
            return to!number(iaf_get_Cv()); // May need something more general for a more complex gas.
        }
        override number dhdT_const_p(in GasState Q) const
        {
            return to!number(iaf_get_Cp());
        }
        override number dpdrho_const_T(in GasState Q) const
        {
            double R = iaf_get_Rgas();
            return R*Q.T;
        }
        override number gas_constant(in GasState Q) const
        {
            return to!number(iaf_get_Rgas());
        }
        override number internal_energy(in GasState Q) const
        {
            return Q.u;
        }
        override number enthalpy(in GasState Q) const
        {
            return Q.u + Q.p/Q.rho;
        }
        override number entropy(in GasState Q) const
        {
            double p = Q.p.re;
            double T = Q.T.re;
            return to!number(iaf_entropy(&p, &T));
        }
    } else {
        override void update_thermo_from_pT(ref GasState Q) const
        {
            iaf_update_thermo_from_pT(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_rhou(ref GasState Q) const
        {
            iaf_update_thermo_from_rhou(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_rhoT(ref GasState Q) const
        {
            iaf_update_thermo_from_rhoT(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_rhop(ref GasState Q) const
        {
            iaf_update_thermo_from_rhop(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr);
        }
        override void update_thermo_from_ps(ref GasState Q, double s) const
        {
            iaf_update_thermo_from_ps(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &s);
        }
        override void update_thermo_from_hs(ref GasState Q, double h, double s) const
        {
            iaf_update_thermo_from_hs(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &h, &s);
        }
        override void update_sound_speed(ref GasState Q) const
        {
            iaf_update_sound_speed(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &(Q.a));
        }
        override void update_trans_coeffs(ref GasState Q)
        {
            iaf_update_trans_coeffs(&(Q.p), &(Q.T), &(Q.rho), &(Q.u), Q.massf.ptr, &(Q.mu), &(Q.k));
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
            return R*Q.T;
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
            double T = Q.T;
            return iaf_entropy(&p, &T);
        }
    }
} // end class IdealAirProxy

version(ideal_air_proxy_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        auto gm = new IdealAirProxy();
        auto gd = GasState(1, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 1.0;
        assert(isClose(gm.R(gd), 287.086, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        assert(isClose(gd.rho, 1.16109, 1.0e-4), failedUnitTest());
        assert(isClose(gd.u, 215314.0, 1.0e-4), failedUnitTest());
        assert(isClose(gd.a, 347.241, 1.0e-4), failedUnitTest());
        gm.update_trans_coeffs(gd);
        assert(isClose(gd.mu, 1.84691e-05, 1.0e-3), failedUnitTest());
        assert(isClose(gd.k, 0.0262449, 1.0e-6), failedUnitTest());

        return 0;
    }
}
