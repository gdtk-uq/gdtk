/**
 * two_temperature_argon_plus_ideal.d
 *
 * Composite gas model of a generic ideal gas plus the ionizing argon gas.
 * We will compute and carry the electron energy in the context
 * of the argon-model species, only, ignoring the other ideal-gas species.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2020-02-17.
 */

module gas.two_temperature_argon_plus_ideal;

import gas.gas_model;
import gas.gas_state;
import gas.ideal_gas;
import gas.two_temperature_reacting_argon;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;
import nm.bracketing;
import nm.brent;

enum Species {ideal=0, Ar, Ar_plus, e_minus}

class TwoTemperatureArgonPlusIdealGas: GasModel {
public:

    this(lua_State *L) {
        type_str = "TwoTemperatureArgonPlusIdealGas";
        _n_species = 4;
        _species_names.length = 4;
        _species_names[Species.ideal] = "ideal";
        _species_names[Species.Ar] = "Ar";
        _species_names[Species.Ar_plus] = "Ar+";
        _species_names[Species.e_minus] = "e-";
        _n_modes = 1; // for storage of electron energy
        //
        // Bring table to TOS
        lua_getglobal(L, "TwoTemperatureArgonPlusIdealGas");
        // There are just two file-names to be read from the composite-gas-model file.
        string ideal_file = getString(L, -1, "ideal_file");
        string argon_file = getString(L, -1, "argon_file");
        doLuaFile(L, ideal_file);
        ideal_gas = new IdealGas(L);
        doLuaFile(L, argon_file);
        argon_gas = new TwoTemperatureReactingArgon(L);
        //
        _mol_masses.length = 4;
        _mol_masses[Species.ideal] = ideal_gas._mol_masses[0];
        _mol_masses[Species.Ar] = argon_gas._mol_masses[0];
        _mol_masses[Species.Ar_plus] = argon_gas._mol_masses[1];
        _mol_masses[Species.e_minus] = argon_gas._mol_masses[2];
        _e_mass_over_ion_mass = _mol_masses[Species.e_minus] / _mol_masses[Species.Ar_plus];
        create_species_reverse_lookup();
        Q_ideal = GasState(ideal_gas);
        Q_argon = GasState(argon_gas);
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "TwoTemperatureArgonPlusIdealGas =(";
        repr ~= "ideal_gas=" ~ to!string(ideal_gas);
        repr ~= ", argon_gas=" ~ to!string(argon_gas);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        if (Q.T <= 0.0 || Q.p <= 0.0 || Q.T_modes[0] <= 0.0) {
            string msg = "Temperature and/or pressure was negative for update_thermo_from_pT.";
            throw new GasModelException(msg);
        }
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        // As an initial guess, let's assume that the gas constants are equal.
        if (with_ideal) {
            Q_ideal.p = Q.massf[0]*Q.p; Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_pT(Q_ideal);
        }
        if (with_argon) {
            Q_argon.p = argon_massf*Q.p; Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_pT(Q_argon);
        }
        if (with_argon && with_ideal) {
            // Need to determine the density at which
            // the partial pressures sum to the mixture pressure.
            number p_error(number rho)
            {
                Q_ideal.rho = Q.massf[0]*rho;
                ideal_gas.update_thermo_from_rhoT(Q_ideal);
                Q_argon.rho = argon_massf*rho;
                argon_gas.update_thermo_from_rhoT(Q_argon);
                return (Q_argon.p + Q_ideal.p) - Q.p;
            }
            // Initial guess for mixture density.
            number rho2 = 1.05*(Q_argon.rho + Q_ideal.rho);
            number rho1 = rho2*0.9;
            bracket!(p_error, number)(rho1, rho2);
            Q.rho = solve!(p_error, number)(rho1, rho2, 1.0e-6);
            // Now that we have density, go back and get internal energy.
            Q_ideal.rho = Q.massf[0]*Q.rho;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
            Q_argon.rho = argon_massf*Q.rho;
            argon_gas.update_thermo_from_rhoT(Q_argon);
            // Internal energy is a weighted average.
            Q.u = Q.massf[0]*Q_ideal.u + argon_massf*Q_argon.u;
            Q.u_modes[0] = Q_argon.u_modes[0]; // Just carry the electron energy, unscaled.
        } else if (with_argon) {
            // argon_massf very close to 1.0
            Q.rho = Q_argon.rho;
            Q.u = Q_argon.u;
            Q.u_modes[0] = Q_argon.u_modes[0];
        } else {
            Q.rho = Q_ideal.rho;
            Q.u = Q_ideal.u;
            Q.u_modes[0] = to!number(0.0);
        }
    }
    override void update_thermo_from_rhou(ref GasState Q)
    {
        if (Q.u <= 0.0 || Q.rho <= 0.0 || Q.u_modes[0] < 0.0) {
            string msg = "Internal energy and/or density was negative for update_thermo_from_rhou.";
            throw new GasModelException(msg);
        }
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon && with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.u_modes[0] = Q.u_modes[0]; // unscaled electron energy
            Q_argon.T_modes[0] = Q.T_modes[0]; // will remain fixed
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            // Need to determine the temperature at which the thermal energies sum
            // to the mixture thermal energy.
            // Freeze the electron energy and temperature.
            // Let's assume that the form of the internal energies is similar,
            // such that we can just form an average T from the known energy.
            // Although this might be a poor approximation, it is only used
            // to get the initial guess for temperature.
            Q_ideal.u = Q.u;
            ideal_gas.update_thermo_from_rhou(Q_ideal);
            Q_argon.u = Q.u;
            argon_gas.update_thermo_from_rhou(Q_argon);
            // Initial guess for Temperature
            number T = Q.massf[0]*Q_ideal.T + argon_massf*Q_argon.T;
            number u_error(number T)
            {
                Q_argon.T = T; Q_ideal.T = T;
                argon_gas.update_thermo_from_rhoT(Q_argon);
                ideal_gas.update_thermo_from_rhoT(Q_ideal);
                return (Q.massf[0]*Q_ideal.u + argon_massf*Q_argon.u) - Q.u;
            }
            number T2 = T*1.1; number T1 = T*0.9;
            bracket!(u_error, number)(T1, T2);
            T = solve!(u_error, number)(T1, T2, 1.0e-6);
            Q_ideal.T = T; Q_argon.T = T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
            argon_gas.update_thermo_from_rhoT(Q_argon);
            // Partial pressures just sum to the mixture pressure.
            Q.p = Q_ideal.p + Q_argon.p;
            Q.T = T;
        } else if (with_argon) {
            // argon_massf very close to 1.0
            Q_argon.rho = Q.rho*argon_massf;
            Q_argon.u = Q.u;
            Q_argon.u_modes[0] = Q.u_modes[0];
            Q_argon.T_modes[0] = Q.T_modes[0]; // will remain fixed
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhou(Q_argon);
            Q.p = Q_argon.p;
            Q.T = Q_argon.T;
        } else {
            // ideal massf very close to 1.0
            Q_ideal.rho = Q.rho*Q.massf[0];
            Q_ideal.u = Q.u;
            ideal_gas.update_thermo_from_rhou(Q_ideal);
            Q.p = Q_ideal.p;
            Q.T = Q_ideal.T;
        }
    }
    override void update_thermo_from_rhoT(ref GasState Q)
    {
        if (Q.T <= 0.0 || Q.rho <= 0.0 || Q.T_modes[0] <= 0.0) {
            string msg = "Temperature and/or density was negative for update_thermo_from_rhoT.";
            throw new GasModelException(msg);
        }
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.rho = argon_massf * Q.rho;
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.T = Q.T;
            Q_ideal.rho = Q.massf[0] * Q.rho;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
        }
        if (with_argon && with_ideal) {
            Q.p = Q_argon.p + Q_ideal.p; // partial pressures add
            Q.u = Q.massf[0]*Q_ideal.u + argon_massf*Q_argon.u;
            Q.u_modes[0] = Q_argon.u_modes[0];
        } else if (with_argon) {
            // argon_massf very close to 1.0
            Q.p = Q_argon.p;
            Q.u = Q_argon.u;
            Q.u_modes[0] = Q_argon.u_modes[0];
        } else {
            Q.p = Q_ideal.p;
            Q.u = Q_ideal.u;
            Q.u_modes[0] = to!number(0.0);
        }
    }
    override void update_thermo_from_rhop(ref GasState Q)
    {
        if (Q.p <= 0.0 || Q.rho <= 0.0) {
            string msg = "Pressure and/or density was negative for update_thermo_from_rhop.";
            throw new GasModelException(msg);
        }
        // Assume Q.T_modes[0] remains fixed.
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon && with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            // Need to determine the temperature at which the partial pressures sum
            // to the mixture pressure.
            // First guess is just ideal gas
            number T = Q.p/ideal_gas.gas_constant(Q_ideal)/Q.rho;
            number p_error(number T)
            {
                Q_argon.T = T;
                Q_argon.T_modes[0] = Q.T_modes[0];
                Q_ideal.T = T;
                argon_gas.update_thermo_from_rhoT(Q_argon);
                ideal_gas.update_thermo_from_rhoT(Q_ideal);
                return (Q_argon.p + Q_ideal.p) - Q.p;
            }
            number T2 = T*1.1; number T1 = T*0.9;
            bracket!(p_error, number)(T1, T2);
            T = solve!(p_error, number)(T1, T2, 1.0e-6);
            Q_argon.T = T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_ideal.T = T;
            argon_gas.update_thermo_from_rhoT(Q_argon);
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
            Q.T = T;
            Q.u = Q.massf[0]*Q_ideal.u + argon_massf*Q_argon.u;
            Q.u_modes[0] = Q_argon.u_modes[0];
        } else if (with_argon) {
            // argon_massf very close to 1.0
            Q_argon.rho = Q.rho; Q_argon.p = Q.p;
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhop(Q_argon);
            Q.T = Q_argon.T;
            Q.u = Q_argon.u;
            Q.u_modes[0] = Q_argon.u_modes[0];
        } else {
            Q_ideal.rho = Q.rho; Q_ideal.p = Q.p;
            ideal_gas.update_thermo_from_rhop(Q_ideal);
            Q.T = Q_ideal.T;
            Q.u = Q_ideal.u;
            Q.u_modes[0] = to!number(0.0);
        }
    }

    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon && with_ideal) {
            throw new GasModelException("update_thermo_from_ps() not implemented for a mix");
        } else if (with_argon) {
            throw new GasModelException("update_thermo_from_ps() not implemented for argon");
        } else {
            Q_ideal.p = Q.p;
            ideal_gas.update_thermo_from_ps(Q_ideal, s);
            Q.T = Q_ideal.T; Q.T_modes[0] = Q_ideal.T;
            Q.u = Q_ideal.u; Q.u_modes[0] = 0.0;
            Q.rho = Q_ideal.rho;
        }
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon && with_ideal) {
            throw new Error("update_thermo_from_hs() not implemented for a mix");
        } else if (with_argon) {
            throw new Error("update_thermo_from_hs() not implemented for argon");
        } else {
            ideal_gas.update_thermo_from_hs(Q_ideal, h, s);
            Q.p = Q_ideal.p;
            Q.T = Q_ideal.T; Q.T_modes[0] = Q_ideal.T;
            Q.u = Q_ideal.u; Q.u_modes[0] = 0.0;
            Q.rho = Q_ideal.rho;
        }
    }
    override void update_sound_speed(ref GasState Q)
    {
        if (Q.T <= 0.0) {
            string msg = "Temperature was negative for update_sound_speed.";
            throw new GasModelException(msg);
        }
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
            argon_gas.update_sound_speed(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
            ideal_gas.update_sound_speed(Q_ideal);
        }
        if (with_argon && with_ideal) {
            Q.a = Q.massf[0]*Q_ideal.a + argon_massf*Q_argon.a;
        } else if (with_argon) {
            Q.a = Q_argon.a;
        } else {
            Q.a = Q_ideal.a;
        }
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
            argon_gas.update_trans_coeffs(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
            ideal_gas.update_trans_coeffs(Q_ideal);
        }
        if (with_argon && with_ideal) {
            Q.mu = Q.massf[0]*Q_ideal.mu + argon_massf*Q_argon.mu;
            Q.k = Q.massf[0]*Q_ideal.k + argon_massf*Q_argon.k;
        } else if (with_argon) {
            Q.mu = Q_argon.mu;
            Q.k = Q_argon.k;
        } else {
            Q.mu = Q_ideal.mu;
            Q.k = Q_ideal.k;
        }
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
        throw new Exception("not implemented");
    }
    */
    override number dudT_const_v(in GasState Q)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
        }
        number result;
        if (with_argon && with_ideal) {
            result = Q.massf[0]*ideal_gas.dudT_const_v(Q_ideal) +
                argon_massf*argon_gas.dudT_const_v(Q_argon);
        } else if (with_argon) {
            result = argon_gas.dudT_const_v(Q_argon);
        } else {
            result = ideal_gas.dudT_const_v(Q_ideal);
        }
        return result;
    }
    override number dhdT_const_p(in GasState Q)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
        }
        number result;
        if (with_argon && with_ideal) {
            result = Q.massf[0]*ideal_gas.dhdT_const_p(Q_ideal) +
                argon_massf*argon_gas.dhdT_const_p(Q_argon);
        } else if (with_argon) {
            result = argon_gas.dhdT_const_p(Q_argon);
        } else {
            result = ideal_gas.dhdT_const_p(Q_ideal);
        }
        return result;
    }
    override number dpdrho_const_T(in GasState Q)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
        }
        number result;
        if (with_argon && with_ideal) {
            result = Q.massf[0]*ideal_gas.dpdrho_const_T(Q_ideal) +
                argon_massf*argon_gas.dpdrho_const_T(Q_argon);
        } else if (with_argon) {
            result = argon_gas.dpdrho_const_T(Q_argon);
        } else {
            result = ideal_gas.dpdrho_const_T(Q_ideal);
        }
        return result;
    }
    override number gas_constant(in GasState Q)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            Q_argon.rho = argon_massf*Q.rho;
            Q_argon.T = Q.T;
            Q_argon.T_modes[0] = Q.T_modes[0];
            Q_argon.massf[0] = Q.massf[1]/argon_massf;
            Q_argon.massf[1] = Q.massf[2]/argon_massf;
            Q_argon.massf[2] = Q.massf[3]/argon_massf;
            argon_gas.update_thermo_from_rhoT(Q_argon);
        }
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
        }
        number result;
        if (with_argon && with_ideal) {
            result = Q.massf[0]*ideal_gas.gas_constant(Q_ideal) +
                argon_massf*argon_gas.gas_constant(Q_argon);
        } else if (with_argon) {
            result = argon_gas.gas_constant(Q_argon);
        } else {
            result = ideal_gas.gas_constant(Q_ideal);
        }
        return result;
    }
    override @nogc number internal_energy(in GasState Q)
    {
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        return Q.u + argon_massf*Q.u_modes[0];
    }
    override number enthalpy(in GasState Q)
    {
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        return Q.u + argon_massf*Q.u_modes[0] + Q.p/Q.rho;
    }
    override number entropy(in GasState Q)
    {
        bool with_ideal = Q.massf[0] > massf_tiny;
        number argon_massf = Q.massf[1]+Q.massf[2]+Q.massf[3];
        bool with_argon = argon_massf > massf_tiny;
        if (with_argon) {
            throw new GasModelException("entropy not implemented in reacting argon.");
        }
        number result;
        if (with_ideal) {
            Q_ideal.rho = Q.massf[0]*Q.rho;
            Q_ideal.T = Q.T;
            ideal_gas.update_thermo_from_rhoT(Q_ideal);
            result = ideal_gas.entropy(Q_ideal);
        }
        return result;
    }

    override void balance_charge(ref GasState Q) const
    {
        Q.massf[Species.e_minus] = Q.massf[Species.Ar_plus] * _e_mass_over_ion_mass;
    }

    @nogc number ionisation_fraction_from_mass_fractions(const(GasState) Q) const
    {
        // Note that this is the ionization fraction within the argon species, only.
        number ions = Q.massf[Species.Ar_plus] / _mol_masses[Species.Ar_plus];
        number atoms = Q.massf[Species.Ar] / _mol_masses[Species.Ar];
        return ions/(ions+atoms);
    }

public:
    // This gas model is composed to two other gas models and the kinetics
    // module will want access to these components.
    GasModel argon_gas;
    GasModel ideal_gas;
    immutable double massf_tiny = 1.0e-6;
private
    GasState Q_argon, Q_ideal;
    double _e_mass_over_ion_mass; // set once mol masses are set in constructor
} // end class TwoTemperatureArgonPlusIdealGas

version(two_temperature_argon_plus_ideal_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        // [TODO] Make a proper test that doesn't assume just ideal air is in the mix.
        lua_State* L = init_lua_State();
        doLuaFile(L, "two-temperature-argon-plus-ideal-air-gas-model.lua");
        auto gm = new TwoTemperatureArgonPlusIdealGas(L);
        lua_close(L);
        auto gd = GasState(4, 1);

        // (1) Try pure ideal air.
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[Species.ideal] = 1.0;
        gd.massf[Species.Ar] = 0.0;
        gd.massf[Species.Ar_plus] = 0.0;
        gd.massf[Species.e_minus] = 0.0;
        assert(approxEqualNumbers(gm.R(gd), to!number(287.086), 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 1, failedUnitTest());
        assert(gm.n_species == 4, failedUnitTest());
        assert(approxEqualNumbers(gd.p, to!number(1.0e5), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[0], to!number(1.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[1], to!number(0.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[2], to!number(0.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[3], to!number(0.0), 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        writeln("rho=", gd.rho);
        assert(approxEqualNumbers(gd.rho, to!number(1.16102), 1.0e-4), failedUnitTest());
        writeln("u=", gd.u);
        assert(approxEqualNumbers(gd.u, to!number(215327.0), 1.0e-4), failedUnitTest());
        writeln("a=", gd.a);
        assert(approxEqualNumbers(gd.a, to!number(347.251), 1.0e-4), failedUnitTest());
        gm.update_trans_coeffs(gd);
        writeln("mu=", gd.mu);
        assert(approxEqualNumbers(gd.mu, to!number(1.84691e-05), 1.0e-6), failedUnitTest());
        writeln("k=", gd.k);
        assert(approxEqualNumbers(gd.k, to!number(0.0262449), 1.0e-6), failedUnitTest());

        gm.update_thermo_from_rhou(gd);
        writeln("same condition: p=", gd.p);
        assert(approxEqualNumbers(gd.p, to!number(1.0e5), 1.0e-6), failedUnitTest());
        writeln("T=", gd.T);
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), failedUnitTest());

        gd.u *= 1.2;
        gm.update_thermo_from_rhou(gd);
        writeln("increase u: p=", gd.p);
        assert(approxEqualNumbers(gd.p, to!number(1.2e5), 1.0e-6), failedUnitTest());
        writeln("T=", gd.T);
        assert(approxEqualNumbers(gd.T, to!number(360.0), 1.0e-6), failedUnitTest());

        gd.p /= 1.2;
        gm.update_thermo_from_rhop(gd);
        writeln("decrease p: u=", gd.u);
        assert(approxEqualNumbers(gd.u, to!number(215327.0), 1.0e-4), failedUnitTest());
        writeln("T=", gd.T);
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = gd.u; // copy unperturbed value, but we don't really need it
            double h = 1.0e-20;
            gd.T += complex(0.0,h);
            gm.update_thermo_from_rhoT(gd);
            double myCv = gd.u.im/h;
            assert(isClose(myCv, gm.dudT_const_v(gd).re), failedUnitTest());
        }

        // (2) Try pure argon.
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.T_modes[0] = 300.0;
        gd.massf[Species.ideal] = 0.0;
        gd.massf[Species.Ar] = 1.0;
        gd.massf[Species.Ar_plus] = 0.0;
        gd.massf[Species.e_minus] = 0.0;

        assert(isClose(gm.R(gd), 208.0, 1.0e-4), failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[Species.ideal], 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[Species.Ar], 1.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[Species.Ar_plus], 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[Species.e_minus], 0.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        number my_rho = 1.0e5 / (208.0 * 300.0);
        assert(isClose(gd.rho, my_rho, 1.0e-4), failedUnitTest());

        number my_Cv = gm.dudT_const_v(gd);
        number my_u = my_Cv*300.0;
        assert(isClose(gd.u, my_u, 1.0e-3), failedUnitTest());

        number my_Cp = gm.dhdT_const_p(gd);
        number alpha = gm.ionisation_fraction_from_mass_fractions(gd);
        number my_a = sqrt(5.0/3.0*208.0*(gd.T + alpha*gd.T_modes[0]));
        assert(isClose(gd.a, my_a, 1.0e-3), failedUnitTest());

        gm.update_trans_coeffs(gd);
        assert(isClose(gd.mu, 22.912e-6, 1.0e-3), failedUnitTest());
        assert(isClose(gd.k, 0.0178625, 1.0e-3), failedUnitTest());

        return 0;
    }
}
