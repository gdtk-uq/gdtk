/**
 * two_temperature_reacting_argon.d
 *
 * Two-temperature reacting argon based off of:
 * "Quasi-One-Dimensional, Nonequilibrium Gas Dynamics of Partially Ionised Two-Temperature Argon"
 * Martin I. Hoffert and Hwachii Lien
 * 
 * The Physics of Fluids 10, 1769 (1967); doi 10.1063/1.1762356
 *
 *
 * Authors: Daniel Smith and Rory Kelly
 * Version: 19-July-2017: initial cut.
 */

module gas.two_temperature_reacting_argon;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;

// First, the basic gas model.

class TwoTemperatureReactingArgon: GasModel {
public:

    this(lua_State *L) {
        // Some parameters are fixed and some come from the gas model file.
        _is_plasma = true;
        _n_species = 3;
        _n_modes = 1;
        _species_names.length = 3;
        _species_names[0] = "Ar";
        _species_names[1] = "Ar_plus";
        _species_names[2] = "e_minus";
        lua_getglobal(L, "TwoTemperatureReactingArgon");
        // [TODO] test that we actually have the table as item -1
        // Now, pull out the remaining numeric value parameters.
        _Rgas = 208.0;
        _theta_ion = 183100.0;
        _theta_A1star = 135300.0;
        _ion_tol = getDouble(L, -1, "ion_tol");
        _mol_masses.length = 3;
        _mol_masses[0] = 39.948e-3; // Units are kg/mol
        _mol_masses[2] = 5.485799e-7; // Units are kg/mol
        _mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
        _Cp = 520.0;
        _Cv = 312.0;
        lua_pop(L, 1); // dispose of the table
        // Compute derived parameters
        create_species_reverse_lookup();
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= "TwoTemperatureReactingArgon =(";
        repr ~= "species=[\"Ar\", \"Ar_plus\", \"e_minus\"]";
        repr ~= ", Mmass=[" ~ to!string(_mol_masses[0]);
        repr ~= "," ~ to!string(_mol_masses[1]);
        repr ~= "," ~ to!string(_mol_masses[2]) ~ "]";
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
        number alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
        if (Q.T <= 0.0 || Q.p <= 0.0) {
            string msg = "Temperature and/or pressure was negative for update_thermo_from_pT."; 
            throw new GasModelException(msg);
        }
        Q.rho = Q.p/(_Rgas*(Q.T + alpha*Q.T_modes[0]));
        Q.u = 3.0/2.0*_Rgas*Q.T;
        if (alpha<=_ion_tol) {
            Q.T_modes[0] = Q.T;
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        } else {
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        }
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
        number alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
        if (Q.u <= 0.0 || Q.rho <= 0.0) {
            string msg = "Internal energy and/or density was negative for update_thermo_from_rhou."; 
            throw new GasModelException(msg);
        }
        Q.T = 2.0/3.0*Q.u/_Rgas;
        if (alpha <= _ion_tol) {
            Q.T_modes[0] = Q.T;
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        } else {
            Q.T_modes[0] = (Q.u_modes[0]/alpha-_Rgas*_theta_ion)*2.0/3.0/_Rgas;
        }
        Q.p = Q.rho*_Rgas*(Q.T+alpha*Q.T_modes[0]);                                     // Q.rho*_Rgas*Q.T;
    }

    override void update_thermo_from_rhoT(GasState Q) const
    {
        number alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
        if (Q.T <= 0.0 || Q.rho <= 0.0) {
            string msg = "Temperature and/or density was negative for update_thermo_from_rhoT."; 
            throw new GasModelException(msg);
        }
        Q.p = Q.rho*_Rgas*(Q.T+alpha*Q.T_modes[0]);     //Q.rho*_Rgas*Q.T;
        Q.u = 3.0/2.0*_Rgas*Q.T;
        if (alpha <= _ion_tol) {
            Q.T_modes[0] = Q.T;
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        } else {
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        }
    }

    override void update_thermo_from_rhop(GasState Q) const
    {
        number alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
        if (Q.p <= 0.0 || Q.rho <= 0.0) {
            string msg = "Pressure and/or density was negative for update_thermo_from_rhop."; 
            throw new GasModelException(msg);
        }
        Q.T = Q.p/Q.rho/_Rgas - alpha*Q.T_modes[0];
        // Assume Q.T_modes[0] is set independently, and correct.
        Q.u = 3.0/2.0*_Rgas*Q.T;
        if (alpha <= _ion_tol) {
            Q.T_modes[0] = Q.T;
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        } else {
            Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
        }
    }
    override void update_thermo_from_ps(GasState Q, number s) const
    {
        throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureReactingArgon.");
    }
    override void update_thermo_from_hs(GasState Q, number h, number s) const
    {
        throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureReactingArgon.");
    }
    override void update_sound_speed(GasState Q) const
    {
        if (Q.T <= 0.0) {
            string msg = "Temperature was negative for update_sound_speed."; 
            throw new GasModelException(msg);
        }
        number _gamma = dhdT_const_p(Q)/dudT_const_v(Q);
        Q.a = sqrt(_gamma*_Rgas*Q.T);
        //[TODO] update the _Cv and _Cp properties to be dependant on alpha...
    }
    override void update_trans_coeffs(GasState Q)
    {
        // TODO (DANIEL) - Should add these in such that viscous effects can be modelled
        Q.mu = 0.0;
        Q.k = 0.0;
        Q.k_modes[0] = 0.0;
    }
    override number dudT_const_v(in GasState Q) const
    {
        //number alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
        //return 3.0/2.0*_Rgas*(1+alpha) + alpha-_Rgas*alpha*(1-alpha)/(2-alpha)*pow((3.0/2.0*Q.T+alpha*_theta_ion)/Q.T,2);
        return to!number(312.0);
    }
    override number dhdT_const_p(in GasState Q) const
    {
        //number alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
        //return 5.0/2.0*_Rgas*(1+alpha) + _Rgas/2*alpha*(1-pow(alpha,2))*pow((5.0/2.0*Q.T+alpha*_theta_ion)/Q.T,2);
        return to!number(520.0);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        return _Rgas*Q.T; //TODO (Daniel) Check this
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(_Rgas);
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
        throw new GasModelException("entropy not implemented in TwoTemperatureReactingArgon.");
    }

private:
    // Thermodynamic constants
    double _Rgas; // J/kg/K
    double _theta_ion;
    double _theta_A1star;
    double _Cp;
    double _Cv;
    double _ion_tol;
} // end class

//// Unit test of the basic gas model...

version(two_temperature_reacting_argon_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : approxEqual;
    int main() {
        writeln("Beginning the unit test...");
        writeln("Testing the gas state functions...");
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/two-temperature-reacting-argon-model.lua");
        auto gm = new TwoTemperatureReactingArgon(L);
        lua_close(L);
        auto gd = new GasState(3, 1);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.T_modes[0] = 300;
        gd.massf[0] = 1.0; gd.massf[1] = 0.0; gd.massf[2] = 0.0;

        assert(approxEqual(gm.R(gd), 208.0, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 1, failedUnitTest());
        assert(gm.n_species == 3, failedUnitTest());
        assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.massf[1], 0.0, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.massf[2], 0.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        number my_rho = 1.0e5 / (208.0 * 300.0);
        assert(approxEqual(gd.rho, my_rho, 1.0e-4), failedUnitTest());

        number my_Cv = gm.dudT_const_v(gd);
        number my_u = my_Cv*300.0; 
        assert(approxEqual(gd.u, my_u, 1.0e-3), failedUnitTest());

        number my_Cp = gm.dhdT_const_p(gd);
        number my_a = sqrt(my_Cp/my_Cv*208.0*300.0);
        assert(approxEqual(gd.a, my_a, 1.0e-3), failedUnitTest());

        gm.update_trans_coeffs(gd);
        assert(approxEqual(gd.mu, 0.0, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.k, 0.0, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.k, 0.0, 1.0e-6), failedUnitTest());

        return 0;
    }
}
