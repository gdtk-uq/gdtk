/**
 * Two-temperature hydrogen, helium gas model for Gas-Giant entry simulations.
 * Authors: Yu Liu, RG and PJ.
 * Date: 2018-09-05 -- 2018-10-15
 * Adapted from Rowan's two-temperature air gas model,
 * using data from:
 *
 */

module gas.two_temperature_gasgiant;

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


immutable double T_REF = 298.15; // K
static string[] molecularSpeciesNames = ["H2"];

// The following symbols are for indexing into the thermo-coefficient database.
// They also set the order for the list of species in this gas model.
enum Species {H2=0, H, Hplus, eminus, He, Heplus}

// For table parameters see end of file.
// They are declared in the static this() function.

class TwoTemperatureGasGiant : GasModel {
public:
    int[] molecularSpecies;

    this()
    {
        type_str = "TwoTemperatureGasGiant";
        _n_species = 6;
        _n_modes = 1;
        _species_names.length = _n_species;
        // Species_names
        _species_names[Species.H2] = "H2";
        _species_names[Species.H] = "H";
        _species_names[Species.Hplus] = "H+";
        _species_names[Species.eminus] = "e-";
        _species_names[Species.He] = "He";
        _species_names[Species.Heplus] = "He+";
        // Molar masses in kg/mole
        _mol_masses.length = _n_species;
        _mol_masses[Species.H2] = 2.01588e-3;
        _mol_masses[Species.H] = 1.00794e-3;
        _mol_masses[Species.eminus] = 5.48579909070e-7;
        _mol_masses[Species.Hplus] = _mol_masses[Species.H] - _mol_masses[Species.eminus];
        _mol_masses[Species.He] = 4.002602e-3;
        _mol_masses[Species.Heplus] = _mol_masses[Species.He] - _mol_masses[Species.eminus];

        _molef.length = _n_species;

        // Particle mass in g
        _particleMass.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _particleMass[isp] = mol_masses[isp]/Avogadro_number;
            _particleMass[isp] *= 1000.0; // kg -> g
        }

        _R.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _R[isp] = R_universal/_mol_masses[isp];//J/(kg.K)
        }
        _pgMixEOS = new PerfectGasMixEOS(_R, true, Species.eminus, 0);
        // enthalpy of formation
        _del_hf.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _del_hf[isp] = del_hf[isp];
        }
        // Cp_tr_rot at lower temperature (When vibrational energy not excited)
        _Cp_tr_rot.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            if (isp == Species.eminus) {
                _Cp_tr_rot[isp] = 0.0;
                continue;
            }
            _Cp_tr_rot[isp] = (5./2.)*_R[isp];
            if (canFind(molecularSpeciesNames, _species_names[isp])) {
                // Plus rotational contribution for molecule, H2
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
                    // Just reverse the order, eg. H2:H --> H:H2
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
    } // end constructor this()

    override string toString() const
    {
        return "TwoTemperatureGasGiant()";
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        _pgMixEOS.update_density(Q); // rho = P/Rg/T
        Q.u = transRotEnergy(Q);     // Cp(T-Tref)+h_ref-RT
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
        	if (isp == Species.eminus) continue;
            sumA += Q.massf[isp]*(_Cp_tr_rot[isp]*T_REF - _del_hf[isp]); //del_hf: h reference
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
        throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureGasGiant.");
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureGasGiant.");
    }

    override void update_sound_speed(ref GasState Q) //not used
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
        foreach (isp; 0 .. _n_species) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_22[isp][jsp]*(log(T))^^2 + _B_22[isp][jsp]*log(T) + _C_22[isp][jsp];
                number pi_Omega_22 = _D_22[isp][jsp]*pow(T, expnt); //Gupta Page 28
                _Delta_22[isp][jsp] = (16./5.0)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!number(PI)*_R_U_cal*T))*pi_Omega_22;//Gupta Page 23
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
            sumA += _particleMass[isp]*_molef[isp]/sumB;//Gupta Page 24
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
        number kB_erg = 1.38066e-16; // erg/K
        number k_tr = 2.3901e-8*(15./4.)*kB_erg*sumA;
        k_tr *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        foreach (isp; 0 .. _n_species) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_11[isp][jsp]*(log(T))^^2 + _B_11[isp][jsp]*log(T) + _C_11[isp][jsp];
                number pi_Omega_11 = _D_11[isp][jsp]*pow(T, expnt);
                _Delta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!number(PI)*_R_U_cal*T))*pi_Omega_11;
                _Delta_11[jsp][isp] = _Delta_11[isp][jsp];
            }
        }
        number k_rot = 0.0;
        number k_vib = 0.0;
        foreach (isp; molecularSpecies) {
            sumB = 0.0;
            foreach (jsp; 0 .. _n_species) {
                sumB += _molef[jsp]*_Delta_11[isp][jsp];
            }
            k_rot += _molef[isp]/sumB;
            number Cp_vib = vibElecSpecHeatConstV(Q.T_modes[0], isp);
            k_vib += (Cp_vib*_mol_masses[isp]/R_universal)*_molef[isp]/sumB;
        }
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
        number R = gas_constant(Q);
        number h = e + R*Q.T;
        return h;
    }
    override number entropy(in GasState Q)
    {
        throw new GasModelException("entropy not implemented in TwoTemperatureNitrogen.");
    }

    @nogc
    number vibElecEnergy(number Tve, int isp) // Vibrational energy (h) = Total energy - Referenced energy- Rotational Energy
    {
        number h_ve;
        if (isp == Species.H2 && Tve.re < 500.0) {
            // Below 500.0 K, the energy vibrational and electronic energy
            // is essentially 0.
            // We short-circuit the computation here because if we go too much
            // lower in temperature, we run into a temperature region where the
            // rotational mode of hydrogen is not fully excited. So, our assumption
            // of a constant Cp_tr_rot no longer holds at low temperatures.
            h_ve = 0.0;
            return h_ve;
        }
        number h_at_Tve = enthalpyFromCurveFits(Tve, isp);
        if (isp == Species.eminus) {
            return h_at_Tve;
        }
        else {
            h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Tve - T_REF) - _del_hf[isp];
        }
        return h_ve;
    }

private:
    PerfectGasMixEOS _pgMixEOS;
    double _R_U_cal = 1.987; // cal/(mole.K)
    number[] _molef; // will be getting mole-fractions from outside, so may be complex
    double[] _particleMass;
    double[] _R;
    double[] _del_hf;
    double[] _Cp_tr_rot;
    number[9] _A; // working storage of coefficients
    number[][] _A_11, _B_11, _C_11, _D_11, _Delta_11, _alpha;
    number[][] _A_22, _B_22, _C_22, _D_22, _Delta_22, _mu;


    /**
     * In order to get smooth variations in thermodynamic properties
     * at the edges of the polynomial breaks, Gupta et al. suggest
     * taking a linear average of the coefficient values near the
     * the polynomial breaks. They describe this process in words
     * on p.14 of their report, and give a Fortran implementation in
     * Appendix C (pp 40 & 41).
     */
    @nogc
    void determineCoefficients(number T, int isp)
    {
        if (T < 800.0) {
            foreach(i; 0 .. _A.length) { _A[i] = thermoCoeffs[isp][0][i]; }
        }
        if (T >= 800.0 && T <= 1200.0) {
            number wB = (1./400.0)*(T - 800.0);
            number wA = 1.0 - wB;
            foreach(i; 0 .. _A.length) { _A[i] = wA*thermoCoeffs[isp][0][i] + wB*thermoCoeffs[isp][1][i]; }
        }
        if (T > 1200.0 && T < 5500.0) {
            foreach(i; 0 .. _A.length) { _A[i] = thermoCoeffs[isp][1][i]; }
        }
        if (T >= 5500.0 && T <= 6500.0) {
            number wB = (1./1000.0)*(T - 5500.0);
            number wA = 1.0 - wB;
            foreach(i; 0 .. _A.length) { _A[i] = wA*thermoCoeffs[isp][1][i] + wB*thermoCoeffs[isp][2][i]; }
        }
        if ( T > 6500.0) {
            foreach(i; 0 .. _A.length) { _A[i] = thermoCoeffs[isp][2][i]; }
        }
    }

    @nogc
    number CpFromCurveFits(number T, int isp)
    {
        // Assume that Cp is constant off the edges of the curve fits.
        // For T < 200.0, this is a reasonable assumption to make.
        // For T > 20000.0, that assumption might be questionable.
        if (T < T_low[isp] ) { T = T_low[isp]; }
        if (T > 20000.0) { T = 20000.0; }
        // For all other cases, use supplied temperature
        determineCoefficients(T, isp);
        number Cp = (_A[0]/T/T + _A[1]/T + _A[2] + _A[3]*T + _A[4]*T*T+_A[5]*T*T*T+_A[6]*T*T*T*T);
        Cp *= (R_universal/_mol_masses[isp]);
        return Cp;
    }

    @nogc
    number enthalpyFromCurveFits(number T, int isp)
    {
        if (T < T_low[isp]) {
            number Cp = CpFromCurveFits(T_low[isp], isp);
            number h = Cp*(T - T_low[isp]) + enthalpyFromCurveFits(to!number(T_low[isp]), isp);
            return h;
        }
        if ( T > 20000.0) {
            number Cp = CpFromCurveFits(to!number(20000.0), isp);
            number h = Cp*(T - 20000.0) + enthalpyFromCurveFits(to!number(20000.0), isp);
            return h;
        }
        // For all other, determine coefficients and compute specific enthalpy.
        determineCoefficients(T, isp);
        number h = (-1)*_A[0]/T/T + _A[1]*log(T)/T + _A[2] + _A[3]*T/2 + _A[4]*T*T/3 + _A[5]*T*T*T/4+_A[6]*T*T*T*T/5+_A[7]/T;
        h *= (R_universal*T/_mol_masses[isp]);
        return h;
    }

    @nogc
    number vibElecEnergy(in GasState Q, number Tve)
    {
        number e_ve = 0.0;
        foreach (isp; 0 .. _n_species) {
            e_ve += Q.massf[isp] * vibElecEnergy(Tve, isp);
        }
        return e_ve;
    }

    @nogc
    number transRotEnergy(in GasState Q)
    {
        // We made a note that Cp is not constant at low temperature for hydrogen.
        // However, when integrating to get enthalpy, the error introduced by
        // assuming constant Cp is very small. So we use it here for simplicity.
        number e_tr_rot = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (isp == Species.eminus) continue;
            number h_tr_rot = _Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp]; // Cp and del_hf here should be mass based not mole
            e_tr_rot += Q.massf[isp]*(h_tr_rot - _R[isp]*Q.T); //calorically perfect gas
        }
        return e_tr_rot;
    }

    @nogc
    number vibElecSpecHeatConstV(number Tve, int isp)
    {
        return CpFromCurveFits(Tve, isp) - _Cp_tr_rot[isp];
    }

    @nogc
    number vibElecSpecHeatConstV(in GasState Q, number Tve)
    {
        number Cv_vib = 0.0;
        foreach (isp; molecularSpecies) {
            Cv_vib += Q.massf[isp] * vibElecSpecHeatConstV(Tve, isp);
        }
        return Cv_vib;
    }

    @nogc
    number transRotSpecHeatConstV(int isp)
    {
        return to!number(_Cp_tr_rot[isp] - _R[isp]);
    }

    @nogc
    number transRotSpecHeatConstV(in GasState Q)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv += Q.massf[isp]*transRotSpecHeatConstV(isp);
        }
        return Cv;
    }

    @nogc
    number vibElecTemperature(in GasState Q)
    {
        number T_guess;
        // We have made the statement that we consider energy in vibElec
        // mode to be 0.0 at 500.0 K or less. This means that when
        // u_modes is 0.0, we cannot distinguish any vibroelectronic
        // temperature below this.
        // So we'll just set the vibroelectronic temperature to 500.0 K.
        if (Q.u_modes[0] <= 0.0) {
            T_guess = 500.0;
            return T_guess;
        }

    	//Let us assume the density mass fractions and densities are known,
    	//now make the vibrational temperature consistent with this data.
        int MAX_ITERATIONS = 20;
        // We'll keep adjusting our temperature estimate
        // until it is less than TOL.
        double TOL = 1.0e-6;

        // Take the supplied T_modes[0] as the initial guess.
        T_guess = Q.T_modes[0];
        number f_guess = vibElecEnergy(Q, T_guess) - Q.u_modes[0];
        // number f_guess = vibEnergy(Q, T_guess) - Q.u_modes[0];
        // Before iterating, check if the supplied guess is
        // good enough. Define good enough as 1/100th of a Joule.
        double E_TOL = 0.01;
        if (fabs(f_guess) < E_TOL) {
            // Given temperature is good enough.
            return T_guess;
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
            // vibEnergy = Curve_fitted_h - Cp_tran*deltaT + h_ref
            // u_modes = vibEnergy(Q, Q.T_modes[0])
            f_guess = vibElecEnergy(Q, T_guess) - Q.u_modes[0];
            count++;
        }

        if (count == MAX_ITERATIONS) {
            string msg = "The 'vibElecTemperature' function failed to converge.\n";
            debug {
                msg ~= format("The final value for T_ve was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= Q.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }

        return T_guess;
    }

}

static double[6] del_hf;
static number[6] T_low;
static double[9][3][6] thermoCoeffs;
static double[string] A_11, B_11, C_11, D_11;
static double[string] A_22, B_22, C_22, D_22;

static this()
{
    del_hf[Species.H2]  = 0.0;
    del_hf[Species.H]   = 217.999/1.00794e-3*1e3;
    del_hf[Species.Hplus]  = 1536.246/(1.00794e-3 - 5.48579909070e-7)*1e3;
    del_hf[Species.eminus]  = 0.0;
    del_hf[Species.He]  = 0.0;
    del_hf[Species.Heplus]  = 2378.521/(4.002602e-3 - 5.48579909070e-7)*1e3;

    T_low[Species.H2] = 200.0;
    T_low[Species.H] = 200.0;
    T_low[Species.Hplus] = 298.15;
    T_low[Species.eminus] = 298.15;
    T_low[Species.He] = 200.0;
    T_low[Species.Heplus] = 298.15;



//Ref: NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species

    thermoCoeffs[Species.H2] =
        [
         [  4.078323210e+04, -8.009186040e+02, 8.214702010e+00, -1.269714457e-02,  1.753605076e-05, -1.202860270e-08,  3.368093490e-12, 2.682484665e+03, -3.043788844e+01], // 200 -- 1000 K
         [  5.608128010e+05, -8.371504740e+02, 2.975364532e+00,  1.252249124e-03, -3.740716190e-07,  5.936625200e-11, -3.606994100e-15, 5.339824410e+03, -2.202774769e+00], // 1000 -- 6000 K
         [  4.966884120e+08, -3.147547149e+05, 7.984121880e+01, -8.414789210e-03,  4.753248350e-07, -1.371873492e-11,  1.605461756e-16, 2.488433516e+06, -6.695728110e+02], // 6000 -- 20000 K
         ];

    thermoCoeffs[Species.H] =
        [
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, 2.547370801e+04, -4.466828530e-01], // 200 -- 1000 K
         [  6.078774250e+01, -1.819354417e-01, 2.500211817e+00, -1.226512864e-07, 3.732876330e-11, -5.687744560e-15, 3.410210197e-19, 2.547486398e+04, -4.481917770e-01], // 1000 -- 6000 K
         [  2.173757694e+08, -1.312035403e+05, 3.399174200e+01, -3.813999680e-03, 2.432854837e-07, -7.694275540e-12, 9.644105630e-17, 1.067638086e+06, -2.742301051e+02], // 6000 -- 20000 K
         ];

    thermoCoeffs[Species.Hplus] =
        [
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, 1.840214877e+05, -1.140646644e+00], // 298.15 -- 1000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, 1.840214877e+05, -1.140646644e+00], // 1000 -- 6000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, 1.840214877e+05, -1.140646644e+00], // 6000 -- 20000 K
         ];

    thermoCoeffs[Species.eminus] =
        [
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, -7.453750000e+02, -1.172081224e+01], // 298.15 -- 1000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, -7.453750000e+02, -1.172081224e+01], // 1000 -- 6000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, -7.453750000e+02, -1.172081224e+01], // 6000 -- 20000 K
         ];

    thermoCoeffs[Species.He] =
        [
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, -7.453750000e+02,  9.287239740e-01], // 200 -- 1000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00, -7.453750000e+02,  9.287239740e-01], // 1000 -- 6000 K
         [  3.396845420e+06, -2.194037652e+03, 3.080231878e+00, -8.068957550e-05, 6.252784910e-09, -2.574990067e-13, 4.429960218e-18,  1.650518960e+04, -4.048814390e+00], // 6000 -- 20000 K
         ];
    thermoCoeffs[Species.Heplus] =
        [
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00,  2.853233739e+05,  1.621665557e+00], // 298.15 -- 1000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00,  2.853233739e+05,  1.621665557e+00], // 1000 -- 6000 K
         [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,  0.000000000e+00, 0.000000000e+00,  0.000000000e+00, 0.000000000e+00,  2.853233739e+05,  1.621665557e+00], // 6000 -- 20000 K
         ];

    // Parameters for collision integrals
    // Collision cross-section Omega_11
    A_11["H2:H2"]   =  -5.62e-3;     B_11["H2:H2"]   =  9.77e-2;   C_11["H2:H2"]   =  -7.95e-1;   D_11["H2:H2"]   =    2.18e+2;
    A_11["H:H2"]    =  -1.47e-3;     B_11["H:H2"]    = -2.79e-2;   C_11["H:H2"]    =   2.27e-1;   D_11["H:H2"]    =    1.50e+1;
    A_11["H:H"]     =  -1.88e-2;     B_11["H:H"]     =  3.59e-1;   C_11["H:H"]     =  -2.51;      D_11["H:H"]     =    1.15e+4;
    A_11["H+:H2"]   =  -1.80e-2;     B_11["H+:H2"]   =  3.22e-1;   C_11["H+:H2"]   =  -2.17;      D_11["H+:H2"]   =    1.07e+4;
    A_11["H+:H"]    =  -9.97e-4;     B_11["H+:H"]    =  2.26e-2;   C_11["H+:H"]    =  -3.19e-1;   D_11["H+:H"]    =    5.04e+2;
    A_11["H+:H+"]   =   0.0;         B_11["H+:H+"]   = -0.0194;    C_11["H+:H+"]   =   0.0119;    D_11["H+:H+"]   =    4.1055;  // Unknown
    A_11["e-:H2"]   =  -2.27e-2;     B_11["e-:H2"]   =  4.64e-1;   C_11["e-:H2"]   =  -2.90;      D_11["e-:H2"]   =    2.71e+3;
    A_11["e-:H"]    =  -6.18e-3;     B_11["e-:H"]    =  6.82e-2;   C_11["e-:H"]    =  -2.01e-1;   D_11["e-:H"]    =    4.47e+1;
    A_11["e-:H+"]   =   0.0;         B_11["e-:H+"]   = -0.0226;    C_11["e-:H+"]   =   0.1300;    D_11["e-:H+"]   =    3.3363;  // Unknown
    A_11["e-:e-"]   =   0.0;         B_11["e-:e-"]   = -0.0139;    C_11["e-:e-"]   =  -0.0825;    D_11["e-:e-"]   =    4.5785;
    A_11["He:H2"]   =  -3.49e-3;     B_11["He:H2"]   =  5.10e-2;   C_11["He:H2"]   =  -4.60e-1;   D_11["He:H2"]   =    8.39e+1;
    A_11["He:H"]    =  -6.82e-3;     B_11["He:H"]    =  1.06e-1;   C_11["He:H"]    =  -8.08e-1;   D_11["He:H"]    =    1.66e+2;
    A_11["He:H+"]   =  -2.05e-2;     B_11["He:H+"]   =  3.91e-1;   C_11["He:H+"]   =  -2.91;      D_11["He:H+"]   =    7.60e+4;
    A_11["He:e-"]   =  -8.24e-3;     B_11["He:e-"]   =  1.68e-1;   C_11["He:e-"]   =  -1.06;      D_11["He:e-"]   =    4.65e+1;
    A_11["He:He"]   =  -3.58e-3;     B_11["He:He"]   =  5.68e-2;   C_11["He:He"]   =  -5.14e-1;   D_11["He:He"]   =    7.29e+1;
    A_11["He+:H2"]  =   2.91e-4;     B_11["He+:H2"]  =  1.54e-2;   C_11["He+:H2"]  =  -7.92e-1;   D_11["He+:H2"]  =    2.86e+3;
    A_11["He+:H"]   =  -1.78e-2;     B_11["He+:H"]   =  2.82e-1;   C_11["He+:H"]   =  -1.64;      D_11["He+:H"]   =    9.30e+2;
    A_11["He+:H+"]  =  -2.05e-2;     B_11["He+:H+"]  =  3.91e-1;   C_11["He+:H+"]  =  -2.91;      D_11["He+:H+"]  =    7.60e+4; //Unknown
    A_11["He+:e-"]  =  -8.24e-3;     B_11["He+:e-"]  =  1.68e-1;   C_11["He+:e-"]  =  -1.06;      D_11["He+:e-"]  =    4.65e+1; //Unknown
    A_11["He+:He"]  =   8.24e-4;     B_11["He+:He"]  = -2.19e-2;   C_11["He+:He"]  =   3.89e-2;   D_11["He+:He"]  =    1.15e+2;
    A_11["He+:He+"] =  -3.58e-3;     B_11["He+:He+"] =  5.68e-2;   C_11["He+:He+"] =  -5.14e-1;   D_11["He+:He+"] =    7.29e+1;  //Unknown

    // Collision cross-section Omega_22
    A_22["H2:H2"]   =  -5.46e-3;    B_22["H2:H2"]   =  9.60e-2; C_22["H2:H2"]   =  -7.67e-1; D_22["H2:H2"]   =    2.21e+2;
    A_22["H:H2"]    =  -3.31e-3;    B_22["H:H2"]    =  8.68e-3; C_22["H:H2"]    =   2.65e-2; D_22["H:H2"]    =    2.40E+1;
    A_22["H:H"]     =  -1.99e-2;    B_22["H:H"]     =  4.13e-1; C_22["H:H"]     =  -3.10;    D_22["H:H"]     =    7.59e+4;
    A_22["H+:H2"]   =  -8.37e-3;    B_22["H+:H2"]   =  1.49e-1; C_22["H+:H2"]   =  -1.17;    D_22["H+:H2"]   =    1.43e+3;
    A_22["H+:H"]    =  -2.39e-2;    B_22["H+:H"]    =  4.70e-1; C_22["H+:H"]    =  -3.31;    D_22["H+:H"]    =    2.02e+5;
    A_22["H+:H+"]   =   0.0;        B_22["H+:H+"]   = -0.0194;  C_22["H+:H+"]   =   0.0119;  D_22["H+:H+"]   =    2.40e+1; // Unknown
    A_22["e-:H2"]   =  -1.95e-2;    B_22["e-:H2"]   =  3.89e-1; C_22["e-:H2"]   =  -2.30;    D_22["e-:H2"]   =    3.72e+2;
    A_22["e-:H"]    =  -1.90e-3;    B_22["e-:H"]    = -3.28e-2; C_22["e-:H"]    =   5.05e-1; D_22["e-:H"]    =    9.10;
    A_22["e-:H+"]   =   0.0;        B_22["e-:H+"]   = -0.0226;  C_22["e-:H+"]   =   0.1300;  D_22["e-:H+"]   =    3.3363;  //Unknown
    A_22["e-:e-"]   =   0.0;        B_22["e-:e-"]   =  0.0;     C_22["e-:e-"]   =  -2.0000;  D_22["e-:e-"]   =    24.3061;
    A_22["He:H2"]   =  -3.46e-3;    B_22["He:H2"]   =  5.24e-2; C_22["He:H2"]   =  -4.62e-1; D_22["He:H2"]   =    9.41e+1;
    A_22["He:H"]    =  -7.67e-3;    B_22["He:H"]    =  1.28e-1; C_22["He:H"]    =  -9.53e-1; D_22["He:H"]    =    2.66e+2;
    A_22["He:H+"]   =  -2.73e-2;    B_22["He:H+"]   =  5.60e-1; C_22["He:H+"]   =  -4.18;    D_22["He:H+"]   =    1.26e+6;
    A_22["He:e-"]   =  -6.49e-3;    B_22["He:e-"]   =  1.20e-1; C_22["He:e-"]   =  -6.91e-1; D_22["He:e-"]   =    1.95e+1;
    A_22["He:He"]   =  -5.25e-3;    B_22["He:He"]   =  9.66e-2; C_22["He:He"]   =  -8.03e-1; D_22["He:He"]   =    1.65e+2;
    A_22["He+:H2"]  =  -1.15e-3;    B_22["He+:H2"]  =  3.91e-2; C_22["He+:H2"]  =  -8.74e-1; D_22["He+:H2"]  =    2.81e+3;
    A_22["He+:H"]   =  -1.62e-2;    B_22["He+:H"]   =  2.83e-1; C_22["He+:H"]   =  -1.86;    D_22["He+:H"]   =    1.95e+3;
    A_22["He+:H+"]  =  -2.05e-2;    B_22["He+:H+"]  =  3.91e-1; C_22["He+:H+"]  =  -2.91;    D_22["He+:H+"]  =    7.60e+4; //Unknown
    A_22["He+:e-"]  =  -8.24e-3;    B_22["He+:e-"]  =  1.68e-1; C_22["He+:e-"]  =  -1.06;    D_22["He+:e-"]  =    4.65e+1; //Unknown
    A_22["He+:He"]  =  -8.77e-3;    B_22["He+:He"]  =  1.50e-1; C_22["He+:He"]  =  -1.09;    D_22["He+:He"]  =    7.53e+2;
    A_22["He+:He+"] =  -3.58e-3;    B_22["He+:He+"] =  5.68e-2; C_22["He+:He+"] =  -5.14e-1; D_22["He+:He+"] =    7.29e+1;  //Unknown

}


version(two_temperature_gasgiant_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    int main() {
        // Comment out lines that write to console
        // to avoid triggering a test failure.
        // writeln("Beginning the unit test...");
        // writeln("Testing the gas state functions...");
        auto gm = new TwoTemperatureGasGiant();
        auto gd = GasState(6, 1);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.T_modes[0] = 1000.0;
        gd.massf[0] = 1.0; gd.massf[1] = 0.0; gd.massf[2] = 0.0; gd.massf[3] = 0.0; gd.massf[4] = 0.0, gd.massf[5] = 0.0;
        //assert(isClose(gm.R(gd), 4124.506, 1.0e-4), failedUnitTest());
        //assert(gm.n_modes == 1, failedUnitTest());
        //assert(gm.n_species == 5, failedUnitTest());
        //assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        //assert(isClose(gd.T, 310.0, 1.0e-6), failedUnitTest());
        //assert(isClose(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());
        //assert(isClose(gd.massf[1], 0.0, 1.0e-6), failedUnitTest());
        //assert(isClose(gd.massf[2], 0.0, 1.0e-6), failedUnitTest());
        //assert(isClose(gd.massf[3], 0.0, 1.0e-6), failedUnitTest());
        //assert(isClose(gd.massf[4], 0.0, 1.0e-6), failedUnitTest());
        /*
        writeln("before update_thermo_from_pT");
        writeln(gd);
        gm.update_thermo_from_pT(gd);
        writeln("after update_thermo_from_pT");
        writeln(gd);
        gd.T_modes[0] = 2000.0;
        writeln("before update_thermo_from_rhou");
        writeln(gd);
        gm.update_thermo_from_rhou(gd);
        writeln("after update_thermo_from_rhou");
        writeln(gd);
        */
        /*
		//writeln("Cv");
		//writeln(gm. vibSpecHeatConstV(dg, 600);
		//writeln("T_modes");
		//writeln(gd);
        gm.update_trans_coeffs(gd);
        writeln("after update_trans_coeffs");
		////assert(isClose(Cp, 169739.4, 1.0e-1), failedUnitTest());
        writeln(gd);
		//gm.update_thermo_from_rhoT(gd);
		//writeln(gd);
		//gm.update_thermo_from_rhop(gd);
		//gm.dudT_const_v(gd);
		//writeln(gm.dudT_const_v(gd));
        writeln("enthalpy");
        writeln(gm.enthalpy(gd));
        */
        return 0;
    }
}
