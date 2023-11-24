/**
 * Authors: Rowan G.
 * Date: 2019-09-22
 *
 * This model for two-temperature dissociating nitrogen
 * began life as a subset of the two-temperature air module.
 */

module gas.two_temperature_dissociating_nitrogen;

import std.math;
import std.conv;
import std.stdio;
import std.string;

import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.thermo.perf_gas_mix_eos;
import gas.thermo.cea_thermo_curves;


immutable double T_REF = 298.15; // K

// The following symbols are for indexing into the thermo-coefficient database.
enum Species {N2=0, N}

// For table parameters see end of file.
// They are declared in the static this() function.

class TwoTemperatureDissociatingNitrogen : GasModel {
public:

    this(lua_State* L)
    {
        type_str = "TwoTemperatureDissociatingNitrogen";
        _n_species = 2;
        _is_plasma = false;
        _n_modes = 1;
        _energy_mode_names.length = 1;
        _energy_mode_names[0] = "vibroelectronic";
        create_energy_mode_reverse_lookup();

        _species_names.length = 2;
        _species_names[Species.N2] = "N2";
        _species_names[Species.N] = "N";
        create_species_reverse_lookup();

        _mol_masses.length = 2;
	_R.length = _n_species;
        _molef.length = _n_species;
        _particleMass.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _mol_masses[isp] = getDouble(L, -1, "M");
            _R[isp] = R_universal/_mol_masses[isp];
            _particleMass[isp] = _mol_masses[isp]/Avogadro_number;
            _particleMass[isp] *= 1000.0; // kg -> g
            lua_getfield(L, -1, "thermoCoeffs");
            _curves ~= new CEAThermoCurve(L, _R[isp]);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        _pgMixEOS = new PerfectGasMixEOS(_R, false, -1, -1);

	/* Heat of formation is value of enthalpy at 298.15K */
        _del_hf.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _del_hf[isp] = _curves[isp].eval_h(to!number(298.15));
        }

        _Cp_tr_rot.length = _n_species;
        _Cp_tr_rot[Species.N2] = (7./2.)*_R[Species.N2];
        _Cp_tr_rot[Species.N] = (5./2.)*_R[Species.N];

        // Setup storage of parameters for collision integrals.
	initialiseParameters();
        _A_11.length = _n_species;
        _B_11.length = _n_species;
        _C_11.length = _n_species;
        _D_11.length = _n_species;
        _Delta_11.length = _n_species;
        _alpha.length = _n_species;
        _A_22.length = _n_species;
        _B_22.length = _n_species;
        _C_22.length = _n_species;
        _D_22.length = _n_species;
        _Delta_22.length = _n_species;
        _mu.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _A_11[isp].length = isp+1;
            _B_11[isp].length = isp+1;
            _C_11[isp].length = isp+1;
            _D_11[isp].length = isp+1;
            _A_22[isp].length = isp+1;
            _B_22[isp].length = isp+1;
            _C_22[isp].length = isp+1;
            _D_22[isp].length = isp+1;
            _mu[isp].length = isp+1;
            // The following are NOT ragged arrays, unlike above.
            _Delta_11[isp].length = _n_species;
            _Delta_22[isp].length = _n_species;
            _alpha[isp].length = _n_species;
            foreach (jsp; 0 .. isp+1) {
                string key = _species_names[isp] ~ ":" ~ _species_names[jsp];
                if (!(key in A_11)) {
                    // Just reverse the order, eg. N2:O2 --> O2:N2
                    key = _species_names[jsp] ~ ":" ~ _species_names[isp];
                }
                _A_11[isp][jsp] = A_11[key];
                _B_11[isp][jsp] = B_11[key];
                _C_11[isp][jsp] = C_11[key];
                _D_11[isp][jsp] = D_11[key];
                _A_22[isp][jsp] = A_22[key];
                _B_22[isp][jsp] = B_22[key];
                _C_22[isp][jsp] = C_22[key];
                _D_22[isp][jsp] = D_22[key];
                double M_isp = mol_masses[isp];
                double M_jsp = mol_masses[jsp];
                _mu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                _mu[isp][jsp] *= 1000.0; // convert kg/mole --> g/mole
                double M_ratio = M_isp/M_jsp;
                double numer = (1.0 - M_ratio)*(0.45 - 2.54*M_ratio);
                double denom = (1.0 + M_ratio)^^2;
                _alpha[isp][jsp] = 1.0 + numer/denom;
                _alpha[jsp][isp] = _alpha[isp][jsp];
            }
        }
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "TwoTemperatureDissociatingNitrogen";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        _pgMixEOS.update_density(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_rhou(ref GasState Q)
    {
        // We can compute T by direct inversion since the Cp in
        // in translation and rotation are fully excited,
        // and, as such, constant.
        number sumA = 0.0;
        number sumB = 0.0;
        foreach (isp; 0 .. _n_species) {
            sumA += Q.massf[isp]*(_Cp_tr_rot[isp]*T_REF - _del_hf[isp]);
            sumB += Q.massf[isp]*(_Cp_tr_rot[isp] - _R[isp]);
        }
        Q.T = (Q.u + sumA)/sumB;
        // Next, we can compute Q.T_modes by iteration.
        // We'll use a Newton method since the function
        // should vary smoothly at the polynomial breaks.
        Q.T_modes[0] = vibTemperature(Q);
        // Now we can compute pressure from the perfect gas
        // equation of state.
        _pgMixEOS.update_pressure(Q);
    }

    override void update_thermo_from_rhoT(ref GasState Q)
    {
        _pgMixEOS.update_pressure(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_rhop(ref GasState Q)
    {
        // In this function, we assume that T_modes is set correctly
        // in addition to density and pressure.
        _pgMixEOS.update_temperature(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureDissociatingNitrogen.");
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureDissociatingNitrogen.");
    }

    override void update_sound_speed(ref GasState Q)
    {
        // We compute the frozen sound speed based on an effective gamma
        number R = gas_constant(Q);
        Q.a = sqrt(gamma(Q)*R*Q.T);
    }

    override void update_trans_coeffs(ref GasState Q)
    {
        massf2molef(Q, _molef);
        // Computation of transport coefficients via collision integrals.
        // Equations follow those in Gupta et al. (1990)
        double kB = Boltzmann_constant;
        number T = Q.T;
        number mylogT = log(Q.T);
        foreach (isp; 0 .. _n_species) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_22[isp][jsp]*(mylogT)^^2 + _B_22[isp][jsp]*mylogT + _C_22[isp][jsp];
                number pi_Omega_22 = exp(_D_22[isp][jsp])*pow(T, expnt);
                _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_22;
                _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
            }
        }
        number sumA = 0.0;
        number sumB;
        foreach (isp; 0 .. n_species) {
            sumB = 0.0;
            foreach (jsp; 0 .. n_species) {
                sumB += _molef[jsp]*_Delta_22[isp][jsp];
            }
            sumA += _particleMass[isp]*_molef[isp]/sumB;
        }
        Q.mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)

        // k = k_tr + k_rot
        sumA = 0.0;
        foreach (isp; 0 .. _n_species) {
            sumB = 0.0;
            foreach (jsp; 0 .. n_species) {
                sumB += _alpha[isp][jsp]*_molef[jsp]*_Delta_22[isp][jsp];
            }
            sumA += _molef[isp]/sumB;
        }
        double kB_erg = 1.38066e-16; // erg/K
        number k_tr = 2.3901e-8*(15./4.)*kB_erg*sumA;
        k_tr *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        foreach (isp; 0 .. _n_species) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_11[isp][jsp]*(mylogT)^^2 + _B_11[isp][jsp]*mylogT + _C_11[isp][jsp];
                number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(T, expnt);
                _Delta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_11;
                _Delta_11[jsp][isp] = _Delta_11[isp][jsp];
            }
        }
        number k_rot = 0.0;
        number k_vib = 0.0;
        // Only N2 has a contribution to k_vib and k_rot
        sumB = 0.0;
        foreach (jsp; 0 .. _n_species) {
            sumB += _molef[jsp]*_Delta_11[Species.N2][jsp];
        }
        k_rot += _molef[Species.N2]/sumB;
        number Cp_vib = vibSpecHeatConstV(Q.T_modes[0], Species.N2);
        k_vib += (Cp_vib*_mol_masses[Species.N2]/R_universal)*_molef[Species.N2]/sumB;

        k_rot *= 2.3901e-8*kB_erg;
        k_rot *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        Q.k = k_tr + k_rot;

        k_vib *= 2.3901e-8*kB_erg;
        k_vib *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        Q.k_modes[0] = k_vib;
    }

    override number dudT_const_v(in GasState Q)
    {
        number Cv = 0.0;
        number Cv_tr_rot, Cv_vib;
        foreach (isp; 0 .. _n_species) {
            Cv_tr_rot = transRotSpecHeatConstV(isp);
            Cv_vib = vibSpecHeatConstV(Q.T_modes[0], isp);
            Cv += Q.massf[isp] * (Cv_tr_rot + Cv_vib);
        }
        return Cv;
    }
    override number dhdT_const_p(in GasState Q)
    {
        // Using the fact that internal structure specific heats
        // are equal, that is, Cp_vib = Cv_vib
        number Cp = 0.0;
        number Cp_vib;
        foreach (isp; 0 .. _n_species) {
            Cp_vib = vibSpecHeatConstV(Q.T_modes[0], isp);
            Cp += Q.massf[isp] * (_Cp_tr_rot[isp] + Cp_vib);
        }

        return Cp;
    }
    override number dpdrho_const_T(in GasState Q)
    {
        number R = gas_constant(Q);
        return R * Q.T;
    }
    override number gas_constant(in GasState Q)
    {
        return  mass_average(Q, _R);
    }
    override number internal_energy(in GasState Q)
    {
        return Q.u + Q.u_modes[0];
    }
    override number enthalpy(in GasState Q)
    {
        number e = transRotEnergy(Q) + vibEnergy(Q, Q.T_modes[0]);
        number h = e + Q.p/Q.rho;
        return h;
    }

    override number enthalpy(in GasState Q, int isp)
    {
	number h_tr_rot = _Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp];
	number h_at_Tve = _curves[isp].eval_h(Q.T_modes[0]);
        number h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Q.T_modes[0] - T_REF) - _del_hf[isp];
	return h_tr_rot + h_ve;
    }

    override number enthalpyPerSpeciesInMode(in GasState Q, int isp, int imode)
    {
        return vibEnergy(Q.T_modes[imode], isp);
    }
    override number entropy(in GasState Q)
    {
	number s = 0.0;
        foreach ( isp; 0.._n_species ) {
            s += Q.massf[isp]*(_curves[isp].eval_s(Q.T) - _R[isp]*log(Q.p/P_atm));
        }
        return s;
    }
    override number entropy(in GasState Q, int isp)
    {
        return _curves[isp].eval_s(Q.T) - _R[isp]*log(Q.p/P_atm);
    }

    override void balance_charge(ref GasState Q) const
    {
        // Do nothing: no charge to balance.
        return;
    }

    @nogc number vibEnergy(number Tve, int isp)
    {
        number h_at_Tve = _curves[isp].eval_h(Tve);
        number h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Tve - T_REF) - _del_hf[isp];
        return h_ve;
    }

    @nogc number vibEnergy(in GasState Q, number Tve)
    {
        number e_ve = 0.0;
        foreach (isp; 0 .. _n_species) {
            e_ve += Q.massf[isp] * vibEnergy(Tve, isp);
        }
        return e_ve;
    }

private:
    PerfectGasMixEOS _pgMixEOS;
    double _R_U_cal = 1.987; // cal/(mole.K)
    number[] _molef; // will be getting mole-fractions from outside, so may be complex
    double[] _particleMass;
    double[] _R;
    CEAThermoCurve[] _curves;
    number[] _del_hf;
    double[] _Cp_tr_rot;
    number[7] _A; // working storage of coefficients
    number[][] _A_11, _B_11, _C_11, _D_11, _Delta_11, _alpha;
    number[][] _A_22, _B_22, _C_22, _D_22, _Delta_22, _mu;

    @nogc number transRotEnergy(in GasState Q)
    {
        number e_tr_rot = 0.0;
        foreach (isp; 0 .. _n_species) {
            number h_tr_rot = _Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp];
            e_tr_rot += Q.massf[isp]*(h_tr_rot - _R[isp]*Q.T);
        }
        return e_tr_rot;
    }

    @nogc number vibSpecHeatConstV(number Tve, int isp)
    {
	return _curves[isp].eval_Cp(Tve) - _Cp_tr_rot[isp];
    }

    @nogc number vibSpecHeatConstV(in GasState Q, number Tve)
    {
        number Cv_vib = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv_vib += Q.massf[isp] * vibSpecHeatConstV(Tve, isp);
        }
        return Cv_vib;
    }

    @nogc number transRotSpecHeatConstV(int isp)
    {
        return to!number(_Cp_tr_rot[isp] - _R[isp]);
    }

    @nogc number transRotSpecHeatConstV(in GasState Q)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv += Q.massf[isp]*transRotSpecHeatConstV(isp);
        }
        return Cv;
    }

    @nogc number vibTemperature(in GasState Q)
    {
        int MAX_ITERATIONS = 20;
        // We'll keep adjusting our temperature estimate
        // until it is less than TOL.
        double TOL = 1.0e-6;

        // Take the supplied T_modes[0] as the initial guess.
        number T_guess = Q.T_modes[0];
        number f_guess = vibEnergy(Q, T_guess) - Q.u_modes[0];
        // Before iterating, check if the supplied guess is
        // good enough. Define good enough as 1/100th of a Joule.
        double E_TOL = 0.01;
        if (fabs(f_guess) < E_TOL) {
            // Given temperature is good enough.
            return Q.T_modes[0];
        }

        // Begin iterating.
        int count = 0;
        number Cv, dT;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            Cv = vibSpecHeatConstV(Q, T_guess);
            dT = -f_guess/Cv;
            T_guess += dT;
            if (fabs(dT) < TOL) {
                break;
            }
            f_guess = vibEnergy(Q, T_guess) - Q.u_modes[0];
            count++;
        }

        if (count == MAX_ITERATIONS) {
            string msg = "The 'vibTemperature' function failed to converge.\n";
            debug {
                msg ~= format("The final value for Tvib was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= Q.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }

        return T_guess;
    }

}

version(two_temperature_dissociating_nitrogen_test) {
    int main()
    {
        return 0;
    }
}

double[string] A_11, B_11, C_11, D_11;
double[string] A_22, B_22, C_22, D_22;

void initialiseParameters()
{
    // Parameters for collision integrals
    // Collision cross-section Omega_11
    A_11["N2:N2"]   =  0.0;    B_11["N2:N2"]   = -0.0112; C_11["N2:N2"]   =  -0.1182; D_11["N2:N2"]   =    4.8464;
    A_11["N:N2"]    =  0.0;    B_11["N:N2"]    = -0.0194; C_11["N:N2"]    =   0.0119; D_11["N:N2"]    =    4.1055;
    A_11["N:N"]     =  0.0;    B_11["N:N"]     = -0.0033; C_11["N:N"]     =  -0.0572; D_11["N:N"]     =    5.0452;

    // Collision cross-section Omega_22
    A_22["N2:N2"]   =  0.0;    B_22["N2:N2"]   = -0.0203; C_22["N2:N2"]   =   0.0683; D_22["N2:N2"]   =   4.0900;
    A_22["N:N2"]    =  0.0;    B_22["N:N2"]    = -0.0190; C_22["N:N2"]    =   0.0239; D_22["N:N2"]    =   4.1782;
    A_22["N:N"]     =  0.0;    B_22["N:N"]     = -0.0118; C_22["N:N"]     =  -0.0960; D_22["N:N"]     =   4.3252;
}
