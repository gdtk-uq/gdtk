/**
 * Authors: Rowan G.
 * Date: 2017-12-05
 *
 */

module gas.two_temperature_air;

import std.math;
import std.conv;
import std.stdio;
import std.string;
import std.algorithm : canFind;

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
static string[] molecularSpeciesNames = ["N2", "O2", "NO", "N2+", "O2+", "NO+"];

// The following symbols are for indexing into the thermo-coefficient database.
enum Species {N=0, O, N2, O2, NO, Nplus, Oplus, N2plus, O2plus, NOplus, eminus}
@nogc Species getSpeciesId(string name)
{
    switch (name) {
    case "N": return Species.N;
    case "O": return Species.O;
    case "N2": return Species.N2;
    case "O2": return Species.O2;
    case "NO": return Species.NO;
    case "N+": return Species.Nplus;
    case "O+": return Species.Oplus;
    case "N2+": return Species.N2plus;
    case "O2+": return Species.O2plus;
    case "NO+": return Species.NOplus;
    case "e-": return Species.eminus;
    default: throw new Error("invalid species name.");
    } // end switch
} // end getSpeciesId

// For table parameters see end of file.
// They are declared in the static this() function.

class TwoTemperatureAir : GasModel {
public:
    int[] molecularSpecies;
    this(lua_State* L)
    {
        type_str = "TwoTemperatureAir";
        getArrayOfStrings(L, "species", _species_names);
        _n_species = to!uint(_species_names.length);
        _n_modes = 1;
        _energy_mode_names.length = 1;
        _energy_mode_names[0] = "vibroelectronic";
        create_energy_mode_reverse_lookup();
        string model;
        if (_n_species == 5) {
            model = "5-species";
        }
        else if (_n_species == 7) {
            model = "7-species";
        }
        else if (_n_species == 11) {
            model = "11-species";
        }
        else {
            throw new Error("");
        }
        // For each of the cases, 5-, 7- and 11-species, we know which species
        // must be supplied, however, we don't know what order the user will
        // hand them to us. So we do a little dance below in which we check
        // that the species handed to us are correct and then set the
        // appropriate properties.
        switch (model) {
        case "5-species":
            _is_plasma = false;
            _mol_masses.length = 5;
            bool[string] validSpecies = ["N":true, "O":true, "N2":true, "O2":true, "NO":true];
            foreach (isp, sp; _species_names) {
                if (sp in validSpecies) {
                    validSpecies.remove(sp);
                }
                else {
                    string errMsg = "The species you supplied is not part of the 5-species 2-T air model,\n";
                    errMsg ~= "or you have supplied a duplicate species.\n";
                    errMsg ~= format("The error occurred for supplied species: %s\n", sp);
                    errMsg ~= "The valid species for the 5-species 2-T air model are:\n";
                    errMsg ~= "   'N', 'O', 'N2', 'O2', 'NO'\n";
                    throw new Error(errMsg);
                }
            }
            create_species_reverse_lookup();
            break;
        case "7-species":
            _is_plasma = true;
            _mol_masses.length = 7;
            bool[string] validSpecies = ["N":true, "O":true, "N2":true, "O2":true, "NO":true,
                                         "NO+":true, "e-":true];
            foreach (isp, sp; _species_names) {
                if (sp in validSpecies) {
                    validSpecies.remove(sp);
                }
                else {
                    string errMsg = "The species you supplied is not part of the 7-species 2-T air model,\n";
                    errMsg ~= "or you have supplied a duplicate species.\n";
                    errMsg ~= format("The error occurred for supplied species: %s\n", sp);
                    errMsg ~= "The valid species for the 7-species 2-T air model are:\n";
                    errMsg ~= "   'N', 'O', 'N2', 'O2', 'NO', 'NO+', 'e-'\n";
                    throw new Error(errMsg);
                }
            }
            if (!(_species_names[$-1] == "e-" || _species_names[$-1] == "eminus")) {
                throw new Error("Electrons should be the last species.\n");
            }
            create_species_reverse_lookup();
            break;
        case "11-species":
            _is_plasma = true;
            _mol_masses.length = 11;
            bool[string] validSpecies = ["N":true, "O":true, "N2":true, "O2":true, "NO":true,
                                         "NO+":true, "N+":true, "O+":true, "N2+":true, "O2+":true, "e-":true];
            foreach (isp, sp; _species_names) {
                if (sp in validSpecies) {
                    validSpecies.remove(sp);
                }
                else {
                    string errMsg = "The species you supplied is not part of the 11-species 2-T air model,\n";
                    errMsg ~= "or you have supplied a duplicate species.\n";
                    errMsg ~= format("The error occurred for supplied species: %s\n", sp);
                    errMsg ~= "The valid species for the 7-species 2-T air model are:\n";
                    errMsg ~= "   'N', 'O', 'N2', 'O2', 'NO', 'NO+', 'N+', 'O+', 'N2+', 'O2+', 'e-'\n";
                    throw new Error(errMsg);
                }
            }
            if (!(_species_names[$-1] == "e-" || _species_names[$-1] == "eminus")) {
                throw new Error("Electrons should be the last species.\n");
            }
            create_species_reverse_lookup();
            break;
        default:
            string errMsg = format("The model name '%s' is not a valid selection for the two-temperature air model.", model);
            throw new Error(errMsg);
        }

        _species_ids.length = _n_species;
        foreach (isp, sp; _species_names) { _species_ids[isp] = getSpeciesId(sp); }

        _R.length = _n_species;
        _molef.length = _n_species;
        _s.length = _n_species;
        _particleMass.length = _n_species;
        _charge.length = n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _mol_masses[isp] = getDouble(L, -1, "M");
            _R[isp] = R_universal/_mol_masses[isp];
            _charge[isp] = getInt(L, -1, "charge");
            _particleMass[isp] = _mol_masses[isp]/Avogadro_number;
            _particleMass[isp] *= 1000.0; // kg -> g
            lua_getfield(L, -1, "thermoCoeffs");
            _curves ~= new CEAThermoCurve(L, _R[isp]);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        if (model == "5-species") {
            _pgMixEOS = new PerfectGasMixEOS(_R, false, -1, -1);
        }
        else {
            _pgMixEOS = new PerfectGasMixEOS(_R, true, Species.eminus, 0);
        }

        /* Heat of formation is value of enthalpy at 298.15K */
        _del_hf.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _del_hf[isp] = _curves[isp].eval_h(to!number(298.15));
        }

        _Cp_tr_rot.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.eminus) {
                // The electron translation is governed by T_ve,
                // so it has no energy contribution in the T_tr mode.
                _Cp_tr_rot[isp] = 0.0;
                continue;
            }
            _Cp_tr_rot[isp] = (5./2.)*_R[isp];
            if (canFind(molecularSpeciesNames, _species_names[isp])) {
                _Cp_tr_rot[isp] += _R[isp];
                molecularSpecies ~= isp;
            }
        }

        // Setup storage of parameters for collision integrals.
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
                double M_isp = _mol_masses[isp];
                double M_jsp = _mol_masses[jsp];
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
        repr ~= "TwoTemperatureAir";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        _pgMixEOS.update_density(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibElecEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_rhou(ref GasState Q)
    {
        // We can compute T by direct inversion since the Cp in
        // in translation and rotation are fully excited,
        // and, as such, constant.
        number sumA = 0.0;
        number sumB = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.eminus) continue;
            sumA += Q.massf[isp]*(_Cp_tr_rot[isp]*T_REF - _del_hf[isp]);
            sumB += Q.massf[isp]*(_Cp_tr_rot[isp] - _R[isp]);
        }
        Q.T = (Q.u + sumA)/sumB;
        // Next, we can compute Q.T_modes by iteration.
        // We'll use a Newton method since the function
        // should vary smoothly at the polynomial breaks.
        Q.T_modes[0] = vibElecTemperature(Q);
        // Now we can compute pressure from the perfect gas
        // equation of state.
        _pgMixEOS.update_pressure(Q);
    }

    override void update_thermo_from_rhoT(ref GasState Q)
    {
        _pgMixEOS.update_pressure(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibElecEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_rhop(ref GasState Q)
    {
        // In this function, we assume that T_modes is set correctly
        // in addition to density and pressure.
        _pgMixEOS.update_temperature(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibElecEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureAir.");
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureAir.");
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
        number Te = Q.T_modes[0];
        number mylogTe = log(Te);
        if ( _n_species == 5 ) {
            foreach (isp; 0 .. _n_species) {
                foreach (jsp; 0 .. isp+1) {
                    number expnt = _A_22[isp][jsp]*(mylogT)^^2 + _B_22[isp][jsp]*mylogT + _C_22[isp][jsp];
                    number pi_Omega_22 = exp(_D_22[isp][jsp])*pow(T, expnt);
                    _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_22;
                    _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
                }
            }
        }
        else {
            // electron present.
            // Compute heavy-particle collision integrals with heavy-particle translation temperature
            // and those involving electron with electron temperature.
            foreach (isp; 0 .. _n_species-1) {
                foreach (jsp; 0 .. isp+1) {
                    number expnt = _A_22[isp][jsp]*(mylogT)^^2 + _B_22[isp][jsp]*mylogT + _C_22[isp][jsp];
                    number pi_Omega_22 = exp(_D_22[isp][jsp])*pow(T, expnt);
                    _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_22;
                    _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
                }
            }
            auto isp = _n_species - 1;
            foreach (jsp; 0 .. _n_species) {
                number expnt = _A_22[isp][jsp]*(mylogTe)^^2 + _B_22[isp][jsp]*mylogTe + _C_22[isp][jsp];
                number pi_Omega_22 = exp(_D_22[isp][jsp])*pow(Te, expnt);
                _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*Te))*pi_Omega_22;
                _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
            }
        }

        // Compute mixture viscosity.
        number sumA = 0.0;
        number sumB;
        if (_n_species == 5) {
            foreach (isp; 0 .. _n_species) {
                sumB = 0.0;
                foreach (jsp; 0 .. _n_species) {
                    sumB += _molef[jsp]*_Delta_22[isp][jsp];
                }
                sumA += _particleMass[isp]*_molef[isp]/sumB;
            }
            Q.mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)
        }
        else {
            // electron present.
            // An additional term is required in the mixture viscosity.
            foreach (isp; 0 .. _n_species-1) {
                sumB = 0.0;
                foreach (jsp; 0 .. _n_species-1) {
                    sumB += _molef[jsp]*_Delta_22[isp][jsp];
                }
                sumB += _molef[_n_species-1]*_Delta_22[isp][_n_species-1];
                sumA += _particleMass[isp]*_molef[isp]/sumB;
            }
            sumB = 0.0;
            foreach (jsp; 0 .. _n_species) {
                sumB += _molef[jsp]*_Delta_22[_n_species-1][jsp];
            }
            sumA += _particleMass[_n_species-1]*_molef[_n_species-1]/sumB;
            Q.mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)
        }

        // Compute component thermal conductivities
        // k in transrotational = k_tr + k_rot
        // k in vibroelectronic = k_ve + k_E
        // 1. k_tr
        sumA = 0.0;
        if (_n_species == 5) {
            foreach (isp; 0 .. _n_species) {
                sumB = 0.0;
                foreach (jsp; 0 .. n_species) {
                    sumB += _alpha[isp][jsp]*_molef[jsp]*_Delta_22[isp][jsp];
                }
                sumA += _molef[isp]/sumB;
            }
        }
        else {
            // electron present.
            foreach (isp; 0 .. _n_species-1) {
                sumB = 0.0;
                foreach (jsp; 0 .. n_species-1) {
                    sumB += _alpha[isp][jsp]*_molef[jsp]*_Delta_22[isp][jsp];
                }
                sumB += 3.54*_molef[_n_species-1]*_Delta_22[isp][_n_species-1];
                sumA += _molef[isp]/sumB;
            }
        }
        double kB_erg = 1.38066e-16; // erg/K
        number k_tr = 2.3901e-8*(15./4.)*kB_erg*sumA;
        k_tr *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        // 2. k_rot
        if (_n_species == 5) {
            foreach (isp; 0 .. _n_species) {
                foreach (jsp; 0 .. isp+1) {
                    number expnt = _A_11[isp][jsp]*(mylogT)^^2 + _B_11[isp][jsp]*mylogT + _C_11[isp][jsp];
                    number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(T, expnt);
                    _Delta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_11;
                    _Delta_11[jsp][isp] = _Delta_11[isp][jsp];
                }
            }
        }
        else {
            // electron present.
            // Compute heavy-particle collision integrals with heavy-particle translation temperature
            // and those involving electron with electron temperature.
            foreach (isp; 0 .. _n_species-1) {
                foreach (jsp; 0 .. isp+1) {
                    number expnt = _A_11[isp][jsp]*(mylogT)^^2 + _B_11[isp][jsp]*mylogT + _C_11[isp][jsp];
                    number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(T, expnt);
                    _Delta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_11;
                    _Delta_11[jsp][isp] = _Delta_11[isp][jsp];
                }
            }
            auto isp = _n_species - 1;
            foreach (jsp; 0 .. _n_species) {
                number expnt = _A_11[isp][jsp]*(mylogTe)^^2 + _B_11[isp][jsp]*mylogTe + _C_11[isp][jsp];
                number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(Te, expnt);
                _Delta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*Te))*pi_Omega_11;
                _Delta_11[jsp][isp] = _Delta_11[isp][jsp];
            }
        }

        number k_rot = 0.0;
        number k_vib_e = 0.0;
        foreach (isp; molecularSpecies) {
            sumB = 0.0;
            if (_n_species == 5) {
                foreach (jsp; 0 .. _n_species) {
                    sumB += _molef[jsp]*_Delta_11[isp][jsp];
                }
            }
            else {
                foreach (jsp; 0 .. _n_species-1) {
                    sumB += _molef[jsp]*_Delta_11[isp][jsp];
                }
                sumB += _molef[_n_species-1]*_Delta_11[isp][_n_species-1];
            }
            k_rot += _molef[isp]/sumB;
            number Cp_vib = vibElecSpecHeatConstV(Q.T_modes[0], isp);
            k_vib_e += (Cp_vib*_mol_masses[isp]/R_universal)*_molef[isp]/sumB;
        }
        k_rot *= 2.3901e-8*kB_erg;
        k_rot *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        Q.k = k_tr + k_rot;

        k_vib_e *= 2.3901e-8*kB_erg;
        k_vib_e *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        number k_E = 0.0;
        if (_n_species > 5) {
            // electron present.
            sumB = 0.0;
            foreach (jsp; 0 .. _n_species-1) {
                sumB += 1.45*_molef[jsp]*_Delta_22[_n_species-1][jsp];
            }
            sumB += _molef[_n_species-1]*_Delta_22[_n_species-1][_n_species-1];
            k_E = _molef[_n_species-1]/sumB;
            k_E *= 2.3901e-8*(15./4.)*kB_erg;
            k_E *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        }
        Q.k_modes[0] = k_vib_e + k_E;
    }

    override number dudT_const_v(in GasState Q)
    {
        number Cv = 0.0;
        number Cv_tr_rot, Cv_vib;
        foreach (isp; 0 .. _n_species) {
            Cv_tr_rot = transRotSpecHeatConstV(isp);
            Cv_vib = vibElecSpecHeatConstV(Q.T_modes[0], isp);
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
            Cp_vib = vibElecSpecHeatConstV(Q.T_modes[0], isp);
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
        number e = transRotEnergy(Q) + vibElecEnergy(Q, Q.T_modes[0]);
        number h = e + Q.p/Q.rho;
        return h;
    }
     override number enthalpy(in GasState Q, int isp)
    {
        number h_tr_rot = _Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp];
        number h_v_e = vibElecEnergy(Q.T_modes[0], isp);
        return h_tr_rot + h_v_e;
    }
    override number enthalpyPerSpeciesInMode(in GasState Q, int isp, int imode)
    {
        return vibElecEnergy(Q.T_modes[imode], isp);
    }
    override number entropy(in GasState Q)
    {
        foreach ( isp; 0.._n_species ) {
            _s[isp] = _curves[isp].eval_s(Q.T) - _R[isp]*log(Q.p/P_atm);
        }
        return mass_average(Q, _s);
    }
    override number entropy(in GasState Q, int isp)
    {
        return _curves[isp].eval_s(Q.T);
    }

    override void balance_charge(ref GasState Q) const
    {
        if (_is_plasma) {
            throw new Error("[FIX-ME] Not yet implemented.");
        }
    }

    @nogc number vibElecEnergy(number Tve, int isp)
    {
        // The electron possess only energy in translation and this is its contribution
        // in the vibroelectronic mode
        if (_species_ids[isp] == Species.eminus)
            return (3./2.)*_R[isp]*Tve;
        number h_at_Tve = _curves[isp].eval_h(Tve);
        number h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Tve - T_REF) - _del_hf[isp];
        return h_ve;
    }

    @nogc number vibElecEnergy(in GasState Q, number Tve)
    {
        number e_ve = 0.0;
        foreach (isp; 0 .. _n_species) {
            e_ve += Q.massf[isp] * vibElecEnergy(Tve, isp);
        }
        return e_ve;
    }
    @nogc
    override void binary_diffusion_coefficients(ref const(GasState) Q, ref number[][] D)
    {
        debug{ assert(D.length==_n_species); }

        number T = Q.T;
        number mylogT = log(Q.T);
        number Te = Q.T_modes[0];
        number mylogTe = log(Te);
        double kB_erg = 1.38066e-16; // erg/K
        number p_erg = Q.p*10; // kg/m/s2 -> g/cm/s2;
        number Dij;

        // Weighted collision cross section Delta(1) from Gupta et al. (1990) eqn. 34
        size_t ilim = (_n_species==5) ? _n_species : _n_species-1;
        foreach (isp; 0 .. ilim) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_11[isp][jsp]*(mylogT)^^2 + _B_11[isp][jsp]*mylogT + _C_11[isp][jsp];
                number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(T, expnt);
                number Delta_11ij = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_11;

                // Binary Diffusion coefficient from eqn. 42a
                Dij = kB_erg*T/p_erg/Delta_11ij; // bin
                Dij *= 1e-4; // cm^2/sec -> m^2/sec
                D[isp][jsp] = Dij;
                D[jsp][isp] = Dij;
            }
        }
        // If electrons present:
        if (_n_species!=5) {
            auto isp = _n_species - 1;
            foreach (jsp; 0 .. _n_species) {
                number expnt = _A_11[isp][jsp]*(mylogTe)^^2 + _B_11[isp][jsp]*mylogTe + _C_11[isp][jsp];
                number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(Te, expnt);
                number Delta_11ij = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*Te))*pi_Omega_11;

                // Binary Diffusion coefficient from eqn. 42a
                Dij = kB_erg*Te/p_erg/Delta_11ij;
                Dij *= 1e-4; // cm^2/sec -> m^2/sec
                D[isp][jsp] = Dij;
                D[jsp][isp] = Dij;
            }
        }
    }

private:
    PerfectGasMixEOS _pgMixEOS;
    double _R_U_cal = 1.987; // cal/(mole.K)
    number[] _molef; // will be getting mole-fractions from outside, so may be complex
    number[] _s;
    double[] _particleMass;
    double[] _R;
    CEAThermoCurve[] _curves;
    number[] _del_hf;
    double[] _Cp_tr_rot;
    number[][] _A_11, _B_11, _C_11, _D_11, _Delta_11, _alpha;
    number[][] _A_22, _B_22, _C_22, _D_22, _Delta_22, _mu;
    Species[] _species_ids;

    @nogc number transRotEnergy(in GasState Q)
    {
        number e_tr_rot = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.eminus) continue;
            number h_tr_rot = _Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp];
            e_tr_rot += Q.massf[isp]*(h_tr_rot - _R[isp]*Q.T);
        }
        return e_tr_rot;
    }

    @nogc number vibElecSpecHeatConstV(number Tve, int isp)
    {
        if (_species_ids[isp] == Species.eminus) {
            number result = (3./2.)*_R[isp];
            return result;
        }
        else
            return _curves[isp].eval_Cp(Tve) - _Cp_tr_rot[isp];
    }

    @nogc number vibElecSpecHeatConstV(in GasState Q, number Tve)
    {
        number Cv_vib = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv_vib += Q.massf[isp] * vibElecSpecHeatConstV(Tve, isp);
        }
        return Cv_vib;
    }

    @nogc number transRotSpecHeatConstV(int isp)
    {
        if (_species_ids[isp] == Species.eminus)
            return to!number(0.0);
        else
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

    @nogc number vibElecTemperature(in GasState Q)
    {
        int MAX_ITERATIONS = 20;
        // We'll keep adjusting our temperature estimate
        // until it is less than TOL.
        double TOL = 1.0e-6;

        // Take the supplied T_modes[0] as the initial guess.
        number T_guess = Q.T_modes[0];
        number f_guess = vibElecEnergy(Q, T_guess) - Q.u_modes[0];
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
            Cv = vibElecSpecHeatConstV(Q, T_guess);
            dT = -f_guess/Cv;
            T_guess += dT;
            if (fabs(dT) < TOL) {
                break;
            }
            f_guess = vibElecEnergy(Q, T_guess) - Q.u_modes[0];
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

version(two_temperature_air_test) {
    int main()
    {
        /*
        auto gm = new TwoTemperatureAir("5-species", ["N2", "O2", "N", "O", "NO"]);
        auto Q = GasState(5, 1);

        Q.p = 666.0;
        Q.T = 293;
        Q.T_modes[0] = 293;
        Q.massf = [0.78, 0.22, 0.0, 0.0, 0.0];
        gm.update_thermo_from_pT(Q);

        writeln(Q);

        gm.update_thermo_from_rhou(Q);

        writeln(Q);
        */
        return 0;
    }
}

version(two_temp_air_ci_dump)
{
    void main()
    {
        import std.stdio;
        string fname = "gupta_etal_1990_CI_data.lua";
        auto f = File(fname, "w");
        f.writeln("-- Collision integral parameters for 11-species air.");
        f.writeln("-- Source:");
        f.writeln("--    Gupta, Yos, Thompson and Lee (1990)");
        f.writeln("--    A Review of Reaction Rates and Thermodynamics and Transport Properties");
        f.writeln("--    for an 11-Species Air Model for Chemical and Thermal Nonequilibrium Calculations to 30 000 K");
        f.writeln("");
        f.writeln("cis = {}");

        // species order matches Gupta tables.
        string[] species = ["N2", "O2", "N", "O", "NO", "NO+", "e-", "N+", "O+", "N2+", "O2+"];

        foreach (isp; 0 .. species.length) {
            foreach (jsp; 0 .. isp+1) {
                string key = species[isp] ~ ":" ~ species[jsp];
                if (!(key in A_11)) {
                    // Just reverse the order, eg. N2:O2 --> O2:N2
                    key = species[jsp] ~ ":" ~ species[isp];
                }
                f.writefln("cis['%s'] = {", key);
                f.writefln("  pi_Omega_11 = {A= % 6.4f, B= % 6.4f, C= % 6.4f, D= % 6.4f},",
                           A_11[key], B_11[key], C_11[key], D_11[key]);
                f.writefln("  pi_Omega_22 = {A= % 6.4f, B= % 6.4f, C= % 6.4f, D= % 6.4f}",
                           A_22[key], B_22[key], C_22[key], D_22[key]);
                f.writeln("}");
            }
        }

        f.close();
    }
}

static double[string] A_11, B_11, C_11, D_11;
static double[string] A_22, B_22, C_22, D_22;

static this()
{
    // Parameters for collision integrals
    // Collision cross-section Omega_11
    A_11["N2:N2"]   =  0.0;    B_11["N2:N2"]   = -0.0112; C_11["N2:N2"]   =  -0.1182; D_11["N2:N2"]   =    4.8464;
    A_11["O2:N2"]   =  0.0;    B_11["O2:N2"]   = -0.0465; C_11["O2:N2"]   =   0.5729; D_11["O2:N2"]   =    1.6185;
    A_11["O2:O2"]   =  0.0;    B_11["O2:O2"]   = -0.0410; C_11["O2:O2"]   =   0.4977; D_11["O2:O2"]   =    1.8302;
    A_11["N:N2"]    =  0.0;    B_11["N:N2"]    = -0.0194; C_11["N:N2"]    =   0.0119; D_11["N:N2"]    =    4.1055;
    A_11["N:O2"]    =  0.0;    B_11["N:O2"]    = -0.0179; C_11["N:O2"]    =   0.0152; D_11["N:O2"]    =    3.9996;
    A_11["N:N"]     =  0.0;    B_11["N:N"]     = -0.0033; C_11["N:N"]     =  -0.0572; D_11["N:N"]     =    5.0452;
    A_11["O:N2"]    =  0.0;    B_11["O:N2"]    = -0.0139; C_11["O:N2"]    =  -0.0825; D_11["O:N2"]    =    4.5785;
    A_11["O:O2"]    =  0.0;    B_11["O:O2"]    = -0.0226; C_11["O:O2"]    =   0.1300; D_11["O:O2"]    =    3.3363;
    A_11["O:N"]     =  0.0;    B_11["O:N"]     =  0.0048; C_11["O:N"]     =  -0.4195; D_11["O:N"]     =    5.7774;
    A_11["O:O"]     =  0.0;    B_11["O:O"]     = -0.0034; C_11["O:O"]     =  -0.0572; D_11["O:O"]     =    4.9901;
    A_11["NO:N2"]   =  0.0;    B_11["NO:N2"]   = -0.0291; C_11["NO:N2"]   =   0.2324; D_11["NO:N2"]   =    3.2082;
    A_11["NO:O2"]   =  0.0;    B_11["NO:O2"]   = -0.0438; C_11["NO:O2"]   =   0.5352; D_11["NO:O2"]   =    1.7252;
    A_11["NO:N"]    =  0.0;    B_11["NO:N"]    = -0.0185; C_11["NO:N"]    =   0.0118; D_11["NO:N"]    =    4.0590;
    A_11["NO:O"]    =  0.0;    B_11["NO:O"]    = -0.0179; C_11["NO:O"]    =   0.0152; D_11["NO:O"]    =    3.9996;
    A_11["NO:NO"]   =  0.0;    B_11["NO:NO"]   = -0.0364; C_11["NO:NO"]   =   0.3825; D_11["NO:NO"]   =    2.4718;
    A_11["NO+:N2"]  =  0.0;    B_11["NO+:N2"]  =     0.0; C_11["NO+:N2"]  =  -0.4000; D_11["NO+:N2"]  =    6.8543;
    A_11["NO+:O2"]  =  0.0;    B_11["NO+:O2"]  =     0.0; C_11["NO+:O2"]  =  -0.4000; D_11["NO+:O2"]  =    6.8543;
    A_11["NO+:N"]   =  0.0;    B_11["NO+:N"]   =     0.0; C_11["NO+:N"]   =  -0.4000; D_11["NO+:N"]   =    6.8543;
    A_11["NO+:O"]   =  0.0;    B_11["NO+:O"]   =     0.0; C_11["NO+:O"]   =  -0.4000; D_11["NO+:O"]   =    6.8543;
    A_11["NO+:NO"]  =  0.0;    B_11["NO+:NO"]  = -0.0047; C_11["NO+:NO"]  =  -0.0551; D_11["NO+:NO"]  =    4.8737;
    A_11["NO+:NO+"] =  0.0;    B_11["NO+:NO+"] =     0.0; C_11["NO+:NO+"] =  -2.0000; D_11["NO+:NO+"] =   23.8237;
    A_11["e-:N2"]   =  0.1147; B_11["e-:N2"]   = -2.8945; C_11["e-:N2"]   =  24.5080; D_11["e-:N2"]   =  -67.3691;
    A_11["e-:O2"]   =  0.0241; B_11["e-:O2"]   = -0.3467; C_11["e-:O2"]   =   1.3887; D_11["e-:O2"]   =   -0.0110;
    A_11["e-:N"]    =  0.0;    B_11["e-:N"]    =     0.0; C_11["e-:N"]    =      0.0; D_11["e-:N"]    =    1.6094;
    A_11["e-:O"]    =  0.0164; B_11["e-:O"]    = -0.2431; C_11["e-:O"]    =   1.1231; D_11["e-:O"]    =   -1.5561;
    A_11["e-:NO"]   = -0.2202; B_11["e-:NO"]   =  5.2265; C_11["e-:NO"]   = -40.5659; D_11["e-:NO"]   =  104.7126;
    A_11["e-:NO+"]  =  0.0;    B_11["e-:NO+"]  =     0.0; C_11["e-:NO+"]  =  -2.0000; D_11["e-:NO+"]  =   23.8237;
    A_11["e-:e-"]   =  0.0;    B_11["e-:e-"]   =     0.0; C_11["e-:e-"]   =  -2.0000; D_11["e-:e-"]   =   23.8237;
    A_11["N+:N2"]   =  0.0;    B_11["N+:N2"]   =     0.0; C_11["N+:N2"]   =  -0.4000; D_11["N+:N2"]   =    6.8543;
    A_11["N+:O2"]   =  0.0;    B_11["N+:O2"]   =     0.0; C_11["N+:O2"]   =  -0.4000; D_11["N+:O2"]   =    6.8543;
    A_11["N+:N"]    =  0.0;    B_11["N+:N"]    = -0.0033; C_11["N+:N"]    =  -0.0572; D_11["N+:N"]    =    5.0452;
    A_11["N+:O"]    =  0.0;    B_11["N+:O"]    =     0.0; C_11["N+:O"]    =  -0.4000; D_11["N+:O"]    =    6.8543;
    A_11["N+:NO"]   =  0.0;    B_11["N+:NO"]   =     0.0; C_11["N+:NO"]   =  -0.4000; D_11["N+:NO"]   =    6.8543;
    A_11["N+:NO+"]  =  0.0;    B_11["N+:NO+"]  =     0.0; C_11["N+:NO+"]  =  -2.0000; D_11["N+:NO+"]  =   23.8237;
    A_11["N+:e-"]   =  0.0;    B_11["N+:e-"]   =     0.0; C_11["N+:e-"]   =  -2.0000; D_11["N+:e-"]   =   23.8237;
    A_11["N+:N+"]   =  0.0;    B_11["N+:N+"]   =     0.0; C_11["N+:N+"]   =  -2.0000; D_11["N+:N+"]   =   23.8237;
    A_11["O+:N2"]   =  0.0;    B_11["O+:N2"]   =     0.0; C_11["O+:N2"]   =  -0.4000; D_11["O+:N2"]   =    6.8543;
    A_11["O+:O2"]   =  0.0;    B_11["O+:O2"]   =     0.0; C_11["O+:O2"]   =  -0.4000; D_11["O+:O2"]   =    6.8543;
    A_11["O+:N"]    =  0.0;    B_11["O+:N"]    =     0.0; C_11["O+:N"]    =  -0.4000; D_11["O+:N"]    =    6.8543;
    A_11["O+:O"]    =  0.0;    B_11["O+:O"]    = -0.0034; C_11["O+:O"]    =  -0.0572; D_11["O+:O"]    =    4.9901;
    A_11["O+:NO"]   =  0.0;    B_11["O+:NO"]   =     0.0; C_11["O+:NO"]   =  -0.4000; D_11["O+:NO"]   =    6.8543;
    A_11["O+:NO+"]  =  0.0;    B_11["O+:NO+"]  =     0.0; C_11["O+:NO+"]  =  -2.0000; D_11["O+:NO+"]  =   23.8237;
    A_11["O+:e-"]   =  0.0;    B_11["O+:e-"]   =     0.0; C_11["O+:e-"]   =  -2.0000; D_11["O+:e-"]   =   23.8237;
    A_11["O+:N+"]   =  0.0;    B_11["O+:N+"]   =     0.0; C_11["O+:N+"]   =  -2.0000; D_11["O+:N+"]   =   23.8237;
    A_11["O+:O+"]   =  0.0;    B_11["O+:O+"]   =     0.0; C_11["O+:O+"]   =  -2.0000; D_11["O+:O+"]   =   23.8237;
    A_11["N2+:N2"]  =  0.0;    B_11["N2+:N2"]  =     0.0; C_11["N2+:N2"]  =  -0.4000; D_11["N2+:N2"]  =    6.8543;
    A_11["N2+:O2"]  =  0.0;    B_11["N2+:O2"]  =     0.0; C_11["N2+:O2"]  =  -0.4000; D_11["N2+:O2"]  =    6.8543;
    A_11["N2+:N"]   =  0.0;    B_11["N2+:N"]   =     0.0; C_11["N2+:N"]   =  -0.4000; D_11["N2+:N"]   =    6.8543;
    A_11["N2+:O"]   =  0.0;    B_11["N2+:O"]   =     0.0; C_11["N2+:O"]   =  -0.4000; D_11["N2+:O"]   =    6.8543;
    A_11["N2+:NO"]  =  0.0;    B_11["N2+:NO"]  =     0.0; C_11["N2+:NO"]  =  -0.4000; D_11["N2+:NO"]  =    6.8543;
    A_11["N2+:NO+"] =  0.0;    B_11["N2+:NO+"] =     0.0; C_11["N2+:NO+"] =  -2.0000; D_11["N2+:NO+"] =   23.8237;
    A_11["N2+:e-"]  =  0.0;    B_11["N2+:e-"]  =     0.0; C_11["N2+:e-"]  =  -2.0000; D_11["N2+:e-"]  =   23.8237;
    A_11["N2+:N+"]  =  0.0;    B_11["N2+:N+"]  =     0.0; C_11["N2+:N+"]  =  -2.0000; D_11["N2+:N+"]  =   23.8237;
    A_11["N2+:O+"]  =  0.0;    B_11["N2+:O+"]  =     0.0; C_11["N2+:O+"]  =  -2.0000; D_11["N2+:O+"]  =   23.8237;
    A_11["N2+:N2+"] =  0.0;    B_11["N2+:N2+"] =     0.0; C_11["N2+:N2+"] =  -2.0000; D_11["N2+:N2+"] =   23.8237;
    A_11["O2+:N2"]  =  0.0;    B_11["O2+:N2"]  =     0.0; C_11["O2+:N2"]  =  -0.4000; D_11["O2+:N2"]  =    6.8543;
    A_11["O2+:O2"]  =  0.0;    B_11["O2+:O2"]  =     0.0; C_11["O2+:O2"]  =  -0.4000; D_11["O2+:O2"]  =    6.8543;
    A_11["O2+:N"]   =  0.0;    B_11["O2+:N"]   =     0.0; C_11["O2+:N"]   =  -0.4000; D_11["O2+:N"]   =    6.8543;
    A_11["O2+:O"]   =  0.0;    B_11["O2+:O"]   =     0.0; C_11["O2+:O"]   =  -0.4000; D_11["O2+:O"]   =    6.8543;
    A_11["O2+:NO"]  =  0.0;    B_11["O2+:NO"]  =     0.0; C_11["O2+:NO"]  =  -0.4000; D_11["O2+:NO"]  =    6.8543;
    A_11["O2+:NO+"] =  0.0;    B_11["O2+:NO+"] =     0.0; C_11["O2+:NO+"] =  -2.0000; D_11["O2+:NO+"] =   23.8237;
    A_11["O2+:e-"]  =  0.0;    B_11["O2+:e-"]  =     0.0; C_11["O2+:e-"]  =  -2.0000; D_11["O2+:e-"]  =   23.8237;
    A_11["O2+:N+"]  =  0.0;    B_11["O2+:N+"]  =     0.0; C_11["O2+:N+"]  =  -2.0000; D_11["O2+:N+"]  =   23.8237;
    A_11["O2+:O+"]  =  0.0;    B_11["O2+:O+"]  =     0.0; C_11["O2+:O+"]  =  -2.0000; D_11["O2+:O+"]  =   23.8237;
    A_11["O2+:N2+"] =  0.0;    B_11["O2+:N2+"] =     0.0; C_11["O2+:N2+"] =  -2.0000; D_11["O2+:N2+"] =   23.8237;
    A_11["O2+:O2+"] =  0.0;    B_11["O2+:O2+"] =     0.0; C_11["O2+:O2+"] =  -2.0000; D_11["O2+:O2+"] =   23.8237;

    // Collision cross-section Omega_22
    A_22["N2:N2"]   =  0.0;    B_22["N2:N2"]   = -0.0203; C_22["N2:N2"]   =   0.0683; D_22["N2:N2"]   =   4.0900;
    A_22["O2:N2"]   =  0.0;    B_22["O2:N2"]   = -0.0558; C_22["O2:N2"]   =   0.7590; D_22["O2:N2"]   =   0.8955;
    A_22["O2:O2"]   =  0.0;    B_22["O2:O2"]   = -0.0485; C_22["O2:O2"]   =   0.6475; D_22["O2:O2"]   =   1.2607;
    A_22["N:N2"]    =  0.0;    B_22["N:N2"]    = -0.0190; C_22["N:N2"]    =   0.0239; D_22["N:N2"]    =   4.1782;
    A_22["N:O2"]    =  0.0;    B_22["N:O2"]    = -0.0203; C_22["N:O2"]    =   0.0703; D_22["N:O2"]    =   3.8818;
    A_22["N:N"]     =  0.0;    B_22["N:N"]     = -0.0118; C_22["N:N"]     =  -0.0960; D_22["N:N"]     =   4.3252;
    A_22["O:N2"]    =  0.0;    B_22["O:N2"]    = -0.0169; C_22["O:N2"]    =  -0.0143; D_22["O:N2"]    =   4.4195;
    A_22["O:O2"]    =  0.0;    B_22["O:O2"]    = -0.0247; C_22["O:O2"]    =   0.1783; D_22["O:O2"]    =   3.2517;
    A_22["O:N"]     =  0.0;    B_22["O:N"]     =  0.0065; C_22["O:N"]     =  -0.4467; D_22["O:N"]     =   6.0426;
    A_22["O:O"]     =  0.0;    B_22["O:O"]     = -0.0207; C_22["O:O"]     =   0.0780; D_22["O:O"]     =   3.5658;
    A_22["NO:N2"]   =  0.0;    B_22["NO:N2"]   = -0.0385; C_22["NO:N2"]   =   0.4226; D_22["NO:N2"]   =   2.4507;
    A_22["NO:O2"]   =  0.0;    B_22["NO:O2"]   = -0.0522; C_22["NO:O2"]   =   0.7045; D_22["NO:O2"]   =   1.0738;
    A_22["NO:N"]    =  0.0;    B_22["NO:N"]    = -0.0196; C_22["NO:N"]    =   0.0478; D_22["NO:N"]    =   4.0321;
    A_22["NO:O"]    =  0.0;    B_22["NO:O"]    = -0.0203; C_22["NO:O"]    =   0.0703; D_22["NO:O"]    =   3.8818;
    A_22["NO:NO"]   =  0.0;    B_22["NO:NO"]   = -0.0453; C_22["NO:NO"]   =   0.5624; D_22["NO:NO"]   =   1.7669;
    A_22["NO+:N2"]  =  0.0;    B_22["NO+:N2"]  =     0.0; C_22["NO+:N2"]  =  -0.4000; D_22["NO+:N2"]  =   6.7760;
    A_22["NO+:O2"]  =  0.0;    B_22["NO+:O2"]  =     0.0; C_22["NO+:O2"]  =  -0.4000; D_22["NO+:O2"]  =   6.7760;
    A_22["NO+:N"]   =  0.0;    B_22["NO+:N"]   =     0.0; C_22["NO+:N"]   =  -0.4000; D_22["NO+:N"]   =   6.7760;
    A_22["NO+:O"]   =  0.0;    B_22["NO+:O"]   =     0.0; C_22["NO+:O"]   =  -0.4000; D_22["NO+:O"]   =   6.7760;
    A_22["NO+:NO"]  =  0.0;    B_22["NO+:NO"]  =     0.0; C_22["NO+:NO"]  =  -0.4000; D_22["NO+:NO"]  =   6.7760;
    A_22["NO+:NO+"] =  0.0;    B_22["NO+:NO+"] =     0.0; C_22["NO+:NO+"] =  -2.0000; D_22["NO+:NO+"] =  24.3602;
    A_22["e-:N2"]   =  0.1147; B_22["e-:N2"]   = -2.8945; C_22["e-:N2"]   =  24.5080; D_22["e-:N2"]   = -67.3691;
    A_22["e-:O2"]   =  0.0241; B_22["e-:O2"]   = -0.3467; C_22["e-:O2"]   =   1.3887; D_22["e-:O2"]   =  -0.0110; // NOTE, low T range
    A_22["e-:N"]    =  0.0;    B_22["e-:N"]    =     0.0; C_22["e-:N"]    =      0.0; D_22["e-:N"]    =   1.6094;
    A_22["e-:O"]    =  0.0164; B_22["e-:O"]    = -0.2431; C_22["e-:O"]    =   1.1231; D_22["e-:O"]    =  -1.5561; // NOTE, low T range
    A_22["e-:NO"]   = -0.2202; B_22["e-:NO"]   =  5.2265; C_22["e-:NO"]   = -40.5659; D_22["e-:NO"]   = 104.7126; // NOTE, low T range
    A_22["e-:NO+"]  =  0.0;    B_22["e-:NO+"]  =     0.0; C_22["e-:NO+"]  =  -2.0000; D_22["e-:NO+"]  =  24.3061;
    A_22["e-:e-"]   =  0.0;    B_22["e-:e-"]   =     0.0; C_22["e-:e-"]   =  -2.0000; D_22["e-:e-"]   =  24.3061;
    A_22["N+:N2"]   =  0.0;    B_22["N+:N2"]   =     0.0; C_22["N+:N2"]   =  -0.4000; D_22["N+:N2"]   =   6.7760;
    A_22["N+:O2"]   =  0.0;    B_22["N+:O2"]   =     0.0; C_22["N+:O2"]   =  -0.4000; D_22["N+:O2"]   =   6.7760;
    A_22["N+:N"]    =  0.0;    B_22["N+:N"]    =     0.0; C_22["N+:N"]    =  -0.4146; D_22["N+:N"]    =   6.9078;
    A_22["N+:O"]    =  0.0;    B_22["N+:O"]    =     0.0; C_22["N+:O"]    =  -0.4000; D_22["N+:O"]    =   6.7760;
    A_22["N+:NO"]   =  0.0;    B_22["N+:NO"]   =     0.0; C_22["N+:NO"]   =  -0.4000; D_22["N+:NO"]   =   6.7760;
    A_22["N+:NO+"]  =  0.0;    B_22["N+:NO+"]  =     0.0; C_22["N+:NO+"]  =  -2.0000; D_22["N+:NO+"]  =  24.3602;
    A_22["N+:e-"]   =  0.0;    B_22["N+:e-"]   =     0.0; C_22["N+:e-"]   =  -2.0000; D_22["N+:e-"]   =  24.3601;
    A_22["N+:N+"]   =  0.0;    B_22["N+:N+"]   =     0.0; C_22["N+:N+"]   =  -2.0000; D_22["N+:N+"]   =  24.3602;
    A_22["O+:N2"]   =  0.0;    B_22["O+:N2"]   =     0.0; C_22["O+:N2"]   =  -0.4000; D_22["O+:N2"]   =   6.7760;
    A_22["O+:O2"]   =  0.0;    B_22["O+:O2"]   =     0.0; C_22["O+:O2"]   =  -0.4000; D_22["O+:O2"]   =   6.7760;
    A_22["O+:N"]    =  0.0;    B_22["O+:N"]    =     0.0; C_22["O+:N"]    =  -0.4000; D_22["O+:N"]    =   6.7760;
    A_22["O+:O"]    =  0.0;    B_22["O+:O"]    =     0.0; C_22["O+:O"]    =  -0.4235; D_22["O+:O"]    =   6.7787;
    A_22["O+:NO"]   =  0.0;    B_22["O+:NO"]   =     0.0; C_22["O+:NO"]   =  -0.4000; D_22["O+:NO"]   =   6.7760;
    A_22["O+:NO+"]  =  0.0;    B_22["O+:NO+"]  =     0.0; C_22["O+:NO+"]  =  -2.0000; D_22["O+:NO+"]  =  24.3602;
    A_22["O+:e-"]   =  0.0;    B_22["O+:e-"]   =     0.0; C_22["O+:e-"]   =  -2.0000; D_22["O+:e-"]   =  24.3601;
    A_22["O+:N+"]   =  0.0;    B_22["O+:N+"]   =     0.0; C_22["O+:N+"]   =  -2.0000; D_22["O+:N+"]   =  24.3602;
    A_22["O+:O+"]   =  0.0;    B_22["O+:O+"]   =     0.0; C_22["O+:O+"]   =  -2.0000; D_22["O+:O+"]   =  24.3602;
    A_22["N2+:N2"]  =  0.0;    B_22["N2+:N2"]  =     0.0; C_22["N2+:N2"]  =  -0.4000; D_22["N2+:N2"]   =   6.7760;
    A_22["N2+:O2"]  =  0.0;    B_22["N2+:O2"]  =     0.0; C_22["N2+:O2"]  =  -0.4000; D_22["N2+:O2"]   =   6.7760;
    A_22["N2+:N"]   =  0.0;    B_22["N2+:N"]   =     0.0; C_22["N2+:N"]   =  -0.4000; D_22["N2+:N"]    =   6.7760;
    A_22["N2+:O"]   =  0.0;    B_22["N2+:O"]   =     0.0; C_22["N2+:O"]   =  -0.4000; D_22["N2+:O"]    =   6.7760;
    A_22["N2+:NO"]  =  0.0;    B_22["N2+:NO"]  =     0.0; C_22["N2+:NO"]  =  -0.4000; D_22["N2+:NO"]   =   6.7760;
    A_22["N2+:NO+"] =  0.0;    B_22["N2+:NO+"] =     0.0; C_22["N2+:NO+"] =  -2.0000; D_22["N2+:NO+"]  =  24.3602;
    A_22["N2+:e-"]  =  0.0;    B_22["N2+:e-"]  =     0.0; C_22["N2+:e-"]  =  -2.0000; D_22["N2+:e-"]   =  24.3601;
    A_22["N2+:N+"]  =  0.0;    B_22["N2+:N+"]  =     0.0; C_22["N2+:N+"]  =  -2.0000; D_22["N2+:N+"]   =  24.3602;
    A_22["N2+:O+"]  =  0.0;    B_22["N2+:O+"]  =     0.0; C_22["N2+:O+"]  =  -2.0000; D_22["N2+:O+"]   =  24.3602;
    A_22["N2+:N2+"] =  0.0;    B_22["N2+:N2+"] =     0.0; C_22["N2+:N2+"] =  -2.0000; D_22["N2+:N2+"]  =  24.3602;
    A_22["O2+:N2"]  =  0.0;    B_22["O2+:N2"]  =     0.0; C_22["O2+:N2"]  =  -0.4000; D_22["O2+:N2"]   =   6.7760;
    A_22["O2+:O2"]  =  0.0;    B_22["O2+:O2"]  =     0.0; C_22["O2+:O2"]  =  -0.4000; D_22["O2+:O2"]   =   6.7760;
    A_22["O2+:N"]   =  0.0;    B_22["O2+:N"]   =     0.0; C_22["O2+:N"]   =  -0.4000; D_22["O2+:N"]    =   6.7760;
    A_22["O2+:O"]   =  0.0;    B_22["O2+:O"]   =     0.0; C_22["O2+:O"]   =  -0.4000; D_22["O2+:O"]    =   6.7760;
    A_22["O2+:NO"]  =  0.0;    B_22["O2+:NO"]  =     0.0; C_22["O2+:NO"]  =  -0.4000; D_22["O2+:NO"]   =   6.7760;
    A_22["O2+:NO+"] =  0.0;    B_22["O2+:NO+"] =     0.0; C_22["O2+:NO+"] =  -2.0000; D_22["O2+:NO+"]  =  24.3602;
    A_22["O2+:e-"]  =  0.0;    B_22["O2+:e-"]  =     0.0; C_22["O2+:e-"]  =  -2.0000; D_22["O2+:e-"]   =  24.3601;
    A_22["O2+:N+"]  =  0.0;    B_22["O2+:N+"]  =     0.0; C_22["O2+:N+"]  =  -2.0000; D_22["O2+:N+"]   =  24.3602;
    A_22["O2+:O+"]  =  0.0;    B_22["O2+:O+"]  =     0.0; C_22["O2+:O+"]  =  -2.0000; D_22["O2+:O+"]   =  24.3602;
    A_22["O2+:N2+"] =  0.0;    B_22["O2+:N2+"] =     0.0; C_22["O2+:N2+"] =  -2.0000; D_22["O2+:N2+"]  =  24.3602;
    A_22["O2+:O2+"] =  0.0;    B_22["O2+:O2+"] =     0.0; C_22["O2+:O2+"] =  -2.0000; D_22["O2+:O2+"]  =  24.3602;

}
