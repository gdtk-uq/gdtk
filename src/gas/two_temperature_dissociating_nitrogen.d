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
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.thermo.perf_gas_mix_eos;


immutable double T_REF = 298.15; // K

// The following symbols are for indexing into the thermo-coefficient database.
enum Species {N2=0, N}

// For table parameters see end of file.
// They are declared in the static this() function.

class TwoTemperatureDissociatingNitrogen : GasModel {
public:

    this(lua_State* L)
    {
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
        _mol_masses[Species.N2] = __mol_masses[Species.N2];
        _mol_masses[Species.N] = __mol_masses[Species.N];


        _molef.length = _n_species;
        _particleMass.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _particleMass[isp] = mol_masses[isp]/Avogadro_number;
            _particleMass[isp] *= 1000.0; // kg -> g
        }

        _R.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _R[isp] = R_universal/_mol_masses[isp];
        }

        _pgMixEOS = new PerfectGasMixEOS(_R, false, -1, -1);

        _del_hf.length = _n_species;
        _del_hf[Species.N2] = del_hf[Species.N2];
        _del_hf[Species.N] = del_hf[Species.N];

        _Cp_tr_rot.length = _n_species;
        _Cp_tr_rot[Species.N2] = (7./2.)*_R[Species.N2];
        _Cp_tr_rot[Species.N] = (5./2.)*_R[Species.N];
        
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

    override void update_thermo_from_pT(GasState Q)
    {
        _pgMixEOS.update_density(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_rhou(GasState Q)
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

    override void update_thermo_from_rhoT(GasState Q)
    {
        _pgMixEOS.update_pressure(Q);
        Q.u = transRotEnergy(Q); 
        Q.u_modes[0] = vibEnergy(Q, Q.T_modes[0]);
    }
    
    override void update_thermo_from_rhop(GasState Q)
    {
        // In this function, we assume that T_modes is set correctly
        // in addition to density and pressure.
        _pgMixEOS.update_temperature(Q);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibEnergy(Q, Q.T_modes[0]);
    }

    override void update_thermo_from_ps(GasState Q, number s)
    {
        throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureDissociatingNitrogen.");
    }

    override void update_thermo_from_hs(GasState Q, number h, number s)
    {
        throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureDissociatingNitrogen.");
    }

    override void update_sound_speed(GasState Q)
    {
        // We compute the frozen sound speed based on an effective gamma
        number R = gas_constant(Q);
        Q.a = sqrt(gamma(Q)*R*Q.T);
    }

    override void update_trans_coeffs(GasState Q)
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
        number Cv_tr, Cv_tr_rot, Cv_vib;
        // N2 contribution
        Cv_tr_rot = transRotSpecHeatConstV(Species.N2);
        Cv_vib = vibSpecHeatConstV(Q.T_modes[0], Species.N2);
        Cv += Q.massf[Species.N2] * (Cv_tr_rot + Cv_vib);
        // N contribution
        Cv_tr = transRotSpecHeatConstV(Species.N);
        Cv += Q.massf[Species.N2] * Cv_tr_rot;

        return Cv;
    }
    override number dhdT_const_p(in GasState Q)
    {
        // Using the fact that internal structure specific heats
        // are equal, that is, Cp_vib = Cv_vib
        number Cp = 0.0;
        number Cp_vib;
        // N2 contribution
        Cp_vib = vibSpecHeatConstV(Q.T_modes[0], Species.N2);
        Cp += Q.massf[Species.N2] * (_Cp_tr_rot[Species.N2] + Cp_vib);
        // N contribution
        Cp += Q.massf[Species.N] * _Cp_tr_rot[Species.N];

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
    override number entropy(in GasState Q)
    {
        throw new GasModelException("entropy not implemented in TwoTemperatureNitrogen.");
    }

    override void balance_charge(GasState Q) const
    {
        // Do nothing: no charge to balance.
        return;
    }

    @nogc number vibEnergy(number Tve, int isp)
    {
        number h_at_Tve = enthalpyFromCurveFits(Tve, isp);
        number h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Tve - T_REF) - _del_hf[isp];
        return h_ve;
    }

    @nogc number vibEnergy(in GasState Q, number Tve)
    {
        return Q.massf[Species.N2] * vibEnergy(Tve, Species.N2);
    }
    
private:
    PerfectGasMixEOS _pgMixEOS;
    double _R_U_cal = 1.987; // cal/(mole.K)
    number[] _molef; // will be getting mole-fractions from outside, so may be complex
    double[] _particleMass;
    double[] _R;
    double[] _del_hf;
    double[] _Cp_tr_rot;
    number[7] _A; // working storage of coefficients
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
    @nogc void determineCoefficients(number T, int isp)
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
        if (T > 6500.0 && T < 14500.0) {
            foreach(i; 0 .. _A.length) { _A[i] = thermoCoeffs[isp][2][i]; }
        }
        if (T >= 14500.0 && T <= 15500.0) {
            number wB = (1./1000.0)*(T - 14500.0);
            number wA = 1.0 - wB;
            foreach(i; 0 .. _A.length) { _A[i] = wA*thermoCoeffs[isp][2][i] + wB*thermoCoeffs[isp][3][i]; }
        }
        if (T > 15500.0 && T < 24500.0) {
            foreach(i; 0 .. _A.length) { _A[i] = thermoCoeffs[isp][3][i]; }
        }
        if (T >= 24500.0 && T <= 25500.0) {
            number wB = (1./1000.0)*(T - 24500.0);
            number wA = 1.0 - wB;
            foreach(i; 0 .. _A.length) { _A[i] = wA*thermoCoeffs[isp][3][i] + wB*thermoCoeffs[isp][4][i]; }
        }
        if ( T > 25500.0) {
            foreach(i; 0 .. _A.length) { _A[i] = thermoCoeffs[isp][4][i]; }
        }
    }

    @nogc number CpFromCurveFits(number T, int isp)
    {
        /* Assume that Cp is constant off the edges of the curve fits.
         * For T < 300.0, this is a reasonable assumption to make.
         * For T > 30000.0, that assumption might be questionable.
         */
        if (T < 300.0) {
            determineCoefficients(to!number(300.0), isp);
            T = 300.0;
        }
        if (T > 30000.0) {
            determineCoefficients(to!number(30000.0), isp);
            T = 30000.0;
        }
        // For all other cases, use supplied temperature
        determineCoefficients(T, isp);
        number T2 = T*T;
        number T3 = T2*T;
        number T4 = T3*T;
        number Cp = (_A[0] + _A[1]*T + _A[2]*T2 + _A[3]*T3 + _A[4]*T4);
        Cp *= (R_universal/_mol_masses[isp]);
        return Cp;
    }

    @nogc number enthalpyFromCurveFits(number T, int isp)
    {
        /* Gupta et al suggest that specific enthalpy below 300 K
         * should be calculated assuming that the specific heat at
         * constant pressure is constant. That suggestion is used
         * here. Additionally, we apply the same idea to temperatures
         * above 30000 K. We will assume that specific heat is constant
         * above 30000 K.
         */
        if (T < T_REF) {
            number Cp = CpFromCurveFits(to!number(300.0), isp);
            number h = Cp*(T - T_REF) + _del_hf[isp];
            return h;
        }
        if (T <= T_REF && T < 300.0) {
            // For the short region between 298.15(=T_REF) to 300.0 K,
            // we just do a linear blend between the value at 298.15 K
            // and the value at 300.0 K. This is to ensure that the
            // reference point is correct at 298.15 and that the enthalpy
            // value matches up at 300.0 K correctly.
            number h_REF = _del_hf[isp];
            number h_300 = enthalpyFromCurveFits(to!number(300.0), isp);
            number w = T - T_REF;
            number h = (1.0 - w)*h_REF + w*h_300;
            return h;
        }
        if ( T > 30000.0) {
            number Cp = CpFromCurveFits(to!number(300000.0), isp);
            number h = Cp*(T - 30000.0) + enthalpyFromCurveFits(to!number(30000.0), isp);
            return h;
        }
        // For all other, determine coefficients and compute specific enthalpy.
        determineCoefficients(T, isp);
        number T2 = T*T;
        number T3 = T2*T;
        number T4 = T3*T;
        number h = _A[0] + _A[1]*T/2. + _A[2]*T2/3. + _A[3]*T3/4. + _A[4]*T4/5. + _A[5]/T;
        h *= (R_universal*T/_mol_masses[isp]);
        return h;
    }

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
        return CpFromCurveFits(Tve, isp) - _Cp_tr_rot[isp];
    }

    @nogc number vibSpecHeatConstV(in GasState Q, number Tve)
    {
        return Q.massf[Species.N2] * vibSpecHeatConstV(Tve, Species.N2);
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

static double[2] __mol_masses;
static double[2] del_hf;
static double[7][5][2] thermoCoeffs;
static double[string] A_11, B_11, C_11, D_11;
static double[string] A_22, B_22, C_22, D_22;

static this()
{
    __mol_masses[Species.N2] = 28.0134e-3;
    __mol_masses[Species.N] = 14.0067e-3;

    /**
     * See Table I in Gnoffo et al.
     */
    del_hf[Species.N2]  = 0.0;
    del_hf[Species.N]   = 112.951e3*4.184/__mol_masses[Species.N];

    thermoCoeffs[Species.N2] = 
        [
         [  0.3674826e+01, -0.1208150e-02,  0.2324010e-05, -0.6321755e-09, -0.2257725e-12, -0.1061160e+04,  0.23580e+01 ], // 300 -- 1000 K
         [  0.2896319e+01,  0.1515486e-02, -0.5723527e-06,  0.9980739e-10, -0.6522355e-14, -0.9058620e+03,  0.43661e+01 ], // 1000 -- 6000 K
         [     0.3727e+01,     0.4684e-03,    -0.1140e-06,     0.1154e-10,    -0.3293e-15,    -0.1043e+04,  0.46264e+01 ], // 6000 -- 15000 K
         [  0.9637690e+01, -0.2572840e-02,  0.3301980e-06, -0.1431490e-10,  0.2033260e-15, -0.1043000e+04, -0.37587e+02 ], // 15000 -- 25000 K
         [ -0.5168080e+01,  0.2333690e-02, -0.1295340e-06,  0.2787210e-11, -0.2135960e-16, -0.1043000e+04,  0.66217e+02 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[Species.N] = 
        [
         [  0.2503071e+01, -0.2180018e-04,  0.5420528e-07, -0.5647560e-10,  0.2099904e-13,  0.5609890e+05,  0.41676e+01 ], // 300 -- 1000 K
         [  0.2450268e+01,  0.1066145e-03, -0.7465337e-07,  0.1879652e-10, -0.1025983e-14,  0.5611600e+05,  0.42618e+01 ], // 1000 -- 6000 K
         [     0.2748e+01,    -0.3909e-03,     0.1338e-06,    -0.1191e-10,     0.3369e-15,     0.5609e+05,  0.28720e+01 ], // 6000 -- 15000 K
         [ -0.1227990e+01,  0.1926850e-02, -0.2437050e-06,  0.1219300e-10, -0.1991840e-15,  0.5609000e+05,  0.28469e+02 ], // 15000 -- 25000 K
         [  0.1552020e+02, -0.3885790e-02,  0.3228840e-06, -0.9605270e-11,  0.9547220e-16,  0.5609000e+05, -0.88120e+02 ]  // 25000 -- 30000 K
         ];

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
