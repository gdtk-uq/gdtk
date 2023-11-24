/**
 * Author: Brad Semple, Rowan G.
 *
 */

module gas.electronically_specific_gas;

import std.algorithm.iteration : sum;
import std.math;
import std.stdio;
import std.string;
import std.path : relativePath;
import std.conv : to;
import util.lua;
import util.lua_service;
import std.algorithm : canFind;
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.thermo.perf_gas_mix_eos;
import gas.electronic_species;


immutable double T_REF = 298.15; // K
static string[] molecularSpeciesNames = ["N2", "O2", "NO", "N2+", "O2+", "NO+"];

class ElectronicallySpecificGas : GasModel
{
public:
    int[] molecularSpecies;
    this(lua_State *L)
    {
        type_str = "ElectronicallySpecificGas";
        auto electronicCompFilename = getString(L,"electronic_species_filename");
        auto macroSpeciesFilename = getString(L,"macro_species_filename");

        //define some things for the macro species
        //Macrostate species will always appear in this order:
        //N, O, N2, O2, NO, N+, O+, N2+, O2+, NO+, e-
        //0, 1, 2,  3,  4,  5,  6,  7,   8,   9,   10
        auto L_macro = init_lua_State();
        doLuaFile(L_macro, macroSpeciesFilename);
        getArrayOfStrings(L_macro, "species", _macro_species_names);


        _n_macro_species = to!int(_macro_species_names.length);
        _macro_molef.length = _n_macro_species;
        macro_Q = GasState(_n_macro_species,1);
        _macro_particleMass.length = _n_macro_species;
        foreach (isp; 0 .. _n_macro_species) {
            _macro_particleMass[isp] = _macro_mol_masses[isp]/Avogadro_number;
            _macro_particleMass[isp] *= 1000.0; // kg -> g
        }

        _macro_R.length = _n_macro_species;
        foreach (isp; 0 .. _n_macro_species) {
            _macro_R[isp] = R_universal/_macro_mol_masses[isp];
        }
        // TODO 2019-04-15
        // RJG think about this. How do we find an effective electron
        // temperature in an electronically specific gas?
        _pgMixEOS = new PerfectGasMixEOS(_macro_R, false, -1, -1);

        _del_hf.length = _n_macro_species;
        foreach (isp; 0 .. _n_macro_species) {
            _del_hf[isp] = del_hf[isp];
        }

        _Cp_tr_rot.length = _n_macro_species;
        foreach (isp; 0 .. _n_macro_species) {
            _Cp_tr_rot[isp] = (5./2.)*_macro_R[isp];
            if (canFind(molecularSpeciesNames, _macro_species_names[isp])) {
                _Cp_tr_rot[isp] += _macro_R[isp];
                molecularSpecies ~= isp;
            }
        }
        // Setup storage of parameters for collision integrals.
        _A_11.length = _n_macro_species;
        _B_11.length = _n_macro_species;
        _C_11.length = _n_macro_species;
        _D_11.length = _n_macro_species;
        _Delta_11.length = _n_macro_species;
        _alpha.length = _n_macro_species;
        _A_22.length = _n_macro_species;
        _B_22.length = _n_macro_species;
        _C_22.length = _n_macro_species;
        _D_22.length = _n_macro_species;
        _Delta_22.length = _n_macro_species;
        _mu.length = _n_macro_species;
        foreach (isp; 0 .. _n_macro_species) {
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
            _Delta_11[isp].length = _n_macro_species;
            _Delta_22[isp].length = _n_macro_species;
            _alpha[isp].length = _n_macro_species;
            foreach (jsp; 0 .. isp+1) {
                string key = _macro_species_names[isp] ~ ":" ~ _macro_species_names[jsp];
                if (!(key in A_11)) {
                    // Just reverse the order, eg. N2:O2 --> O2:N2
                    key = _macro_species_names[jsp] ~ ":" ~ _macro_species_names[isp];
                }
                _A_11[isp][jsp] = A_11[key];
                _B_11[isp][jsp] = B_11[key];
                _C_11[isp][jsp] = C_11[key];
                _D_11[isp][jsp] = D_11[key];
                _A_22[isp][jsp] = A_22[key];
                _B_22[isp][jsp] = B_22[key];
                _C_22[isp][jsp] = C_22[key];
                _D_22[isp][jsp] = D_22[key];
                double M_isp = _macro_mol_masses[isp];
                double M_jsp = _macro_mol_masses[jsp];
                _mu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                _mu[isp][jsp] *= 1000.0; // convert kg/mole --> g/mole
                double M_ratio = M_isp/M_jsp;
                double numer = (1.0 - M_ratio)*(0.45 - 2.54*M_ratio);
                double denom = (1.0 + M_ratio)^^2;
                _alpha[isp][jsp] = 1.0 + numer/denom;
                _alpha[jsp][isp] = _alpha[isp][jsp];
            }
        }

        //define some things for the electronic species
        auto L_elec = init_lua_State();
        doLuaFile(L_elec, electronicCompFilename);
        _n_N_species = getInt(L_elec, "number_N_species");
        _n_O_species = getInt(L_elec, "number_O_species");
        _n_elec_species = _n_N_species + _n_O_species;

        lua_getglobal(L_elec, "electronic_species");
        foreach (isp; 0 .. _n_elec_species) {
            lua_rawgeti(L_elec, -1, isp);
            if (lua_isnil(L_elec, -1)) {
                string msg = format("There was an error when attempting to information about electronic-species %d.\n", isp);
                throw new Error(msg);
            }

            _electronicSpecies ~= createElectronicSpecies(L_elec);
            lua_pop(L_elec, 1);
            _elec_species_names ~= _electronicSpecies[$-1].name;
            _elec_mol_masses ~= _electronicSpecies[$-1].mol_mass;
            _elec_R ~= R_universal/_electronicSpecies[$-1].mol_mass;
            _elec_electronic_energy ~= _electronicSpecies[$-1].electronic_energy;

        }

        //define global properties for the general gas model
        _n_species = _n_macro_species + _n_elec_species - 2;
        _n_modes = 1;
        _species_names = _elec_species_names ~ _macro_species_names[2 .. $];
        _energy_mode_names.length = 1;
        _energy_mode_names[0] = "vibroelectronic";
        create_energy_mode_reverse_lookup();
        _mol_masses = _elec_mol_masses ~ _macro_mol_masses[2 .. $];
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "TwoTemperatureAir";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        Update_Macro_State(macro_Q, Q);

        _pgMixEOS.update_density(macro_Q);
        macro_Q.u = transRotEnergy(macro_Q);
        macro_Q.u_modes[0] = vibEnergy(macro_Q, macro_Q.T_modes[0]);
        Update_Electronic_State(Q, macro_Q);
    }

    override void update_thermo_from_rhou(ref GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        // We can compute T by direct inversion since the Cp in
        // in translation and rotation are fully excited,
        // and, as such, constant.
        number sumA = 0.0;
        number sumB = 0.0;
        foreach (isp; 0 .. _n_macro_species) {
            sumA += macro_Q.massf[isp]*(_Cp_tr_rot[isp]*T_REF - _del_hf[isp]);
            sumB += macro_Q.massf[isp]*(_Cp_tr_rot[isp] - _macro_R[isp]);
        }
        macro_Q.T = (macro_Q.u + sumA)/sumB;
        // Next, we can compute Q.T_modes by iteration.
        // We'll use a Newton method since the function
        // should vary smoothly at the polynomial breaks.
        macro_Q.T_modes[0] = vibTemperature(macro_Q);
        // Now we can compute pressure from the perfect gas
        // equation of state.
        _pgMixEOS.update_pressure(macro_Q);
        Update_Electronic_State(Q, macro_Q);
    }

    override void update_thermo_from_rhoT(ref GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        _pgMixEOS.update_pressure(macro_Q);
        macro_Q.u = transRotEnergy(macro_Q);
        macro_Q.u_modes[0] = vibEnergy(macro_Q, macro_Q.T_modes[0]);
        Update_Electronic_State(Q, macro_Q);
    }

    override void update_thermo_from_rhop(ref GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        // In this function, we assume that T_modes is set correctly
        // in addition to density and pressure.
        _pgMixEOS.update_temperature(macro_Q);
        macro_Q.u = transRotEnergy(macro_Q);
        macro_Q.u_modes[0] = vibEnergy(macro_Q, macro_Q.T_modes[0]);
        Update_Electronic_State(Q, macro_Q);
    }

    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        throw new GasModelException("update_thermo_from_ps not implemented in ElectronicallySpecificGas.");
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        throw new GasModelException("update_thermo_from_hs not implemented in ElectronicallySpecificGas.");
    }

    override void update_sound_speed(ref GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        // We compute the frozen sound speed based on an effective gamma
        number R = gas_constant(Q);
        macro_Q.a = sqrt(gamma(Q)*R*macro_Q.T);
        Update_Electronic_State(Q, macro_Q);
    }

    override void update_trans_coeffs(ref GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        massf2molef(macro_Q, _macro_molef);
        // Computation of transport coefficients via collision integrals.
        // Equations follow those in Gupta et al. (1990)
        double kB = Boltzmann_constant;
        number T = macro_Q.T;
        number mylogT = log(macro_Q.T);
        foreach (isp; 0 .. _n_macro_species) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_22[isp][jsp]*(mylogT)^^2 + _B_22[isp][jsp]*mylogT + _C_22[isp][jsp];
                number pi_Omega_22 = exp(_D_22[isp][jsp])*pow(T, expnt);
                _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_22;
                _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
            }
        }
        number sumA = 0.0;
        number sumB;
        foreach (isp; 0 .. _n_macro_species) {
            sumB = 0.0;
            foreach (jsp; 0 .. _n_macro_species) {
                sumB += _macro_molef[jsp]*_Delta_22[isp][jsp];
            }
            sumA += _macro_particleMass[isp]*_macro_molef[isp]/sumB;
        }
        macro_Q.mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)

        // k = k_tr + k_rot
        sumA = 0.0;
        foreach (isp; 0 .. _n_macro_species) {
            sumB = 0.0;
            foreach (jsp; 0 .. _n_macro_species) {
                sumB += _alpha[isp][jsp]*_macro_molef[jsp]*_Delta_22[isp][jsp];
            }
            sumA += _macro_molef[isp]/sumB;
        }
        double kB_erg = 1.38066e-16; // erg/K
        number k_tr = 2.3901e-8*(15./4.)*kB_erg*sumA;
        k_tr *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        foreach (isp; 0 .. _n_macro_species) {
            foreach (jsp; 0 .. isp+1) {
                number expnt = _A_11[isp][jsp]*(mylogT)^^2 + _B_11[isp][jsp]*mylogT + _C_11[isp][jsp];
                number pi_Omega_11 = exp(_D_11[isp][jsp])*pow(T, expnt);
                _Delta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*_mu[isp][jsp]/(to!double(PI)*_R_U_cal*T))*pi_Omega_11;
                _Delta_11[jsp][isp] = _Delta_11[isp][jsp];
            }
        }
        number k_rot = 0.0;
        number k_vib = 0.0;
        foreach (isp; molecularSpecies) {
            sumB = 0.0;
            foreach (jsp; 0 .. _n_macro_species) {
                sumB += _macro_molef[jsp]*_Delta_11[isp][jsp];
            }
            k_rot += _macro_molef[isp]/sumB;
            number Cp_vib = vibSpecHeatConstV(macro_Q.T_modes[0], isp);
            k_vib += (Cp_vib*_macro_mol_masses[isp]/R_universal)*_macro_molef[isp]/sumB;
        }
        k_rot *= 2.3901e-8*kB_erg;
        k_rot *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        macro_Q.k = k_tr + k_rot;

        k_vib *= 2.3901e-8*kB_erg;
        k_vib *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        macro_Q.k_modes[0] = k_vib;
        Update_Electronic_State(Q, macro_Q);
    }

    override number dudT_const_v(in GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        number Cv = 0.0;
        number Cv_tr_rot, Cv_vib;
        foreach (isp; 0 .. _n_macro_species) {
            Cv_tr_rot = transRotSpecHeatConstV(isp);
            Cv_vib = vibSpecHeatConstV(macro_Q.T_modes[0], isp);
            Cv += macro_Q.massf[isp] * (Cv_tr_rot + Cv_vib);
        }
        return Cv;
    }
    override number dhdT_const_p(in GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        // Using the fact that internal structure specific heats
        // are equal, that is, Cp_vib = Cv_vib
        number Cp = 0.0;
        number Cp_vib;
        foreach (isp; 0 .. _n_macro_species) {
            Cp_vib = vibSpecHeatConstV(macro_Q.T_modes[0], isp);
            Cp += macro_Q.massf[isp] * (_Cp_tr_rot[isp] + Cp_vib);
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
        Update_Macro_State(macro_Q, Q);
        return  mass_average(macro_Q, _macro_R);
    }
    override number internal_energy(in GasState Q)
    {
        return Q.u + Q.u_modes[0];
    }
    override number enthalpy(in GasState Q)
    {
        Update_Macro_State(macro_Q, Q);
        number e = transRotEnergy(macro_Q) + vibEnergy(macro_Q, macro_Q.T_modes[0]);
        number R = gas_constant(Q);
        number h = e + R*macro_Q.T;
        return h;
    }
    override number entropy(in GasState Q)
    {
        throw new GasModelException("entropy not implemented in ElectronicallySpecificGas.");
    }

    @nogc number vibEnergy(number Tve, int isp)
    {
        number h_at_Tve = enthalpyFromCurveFits(Tve, isp);
        number h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Tve - T_REF) - _del_hf[isp];
        debug {
            // writeln("Tve= ", Tve, " h_at_Tve= ", h_at_Tve, " h_ve= ", h_ve);
        }
        return h_ve;
    }

    @nogc number vibEnergy(in GasState macro_Q, number Tve)
    {
        //Update_Macro_State(macro_Q, Q);
        number e_ve = 0.0;
        foreach (isp; molecularSpecies) {
            e_ve += macro_Q.massf[isp] * vibEnergy(Tve, isp);
        }
        return e_ve;
    }


private:
    //storage for electronically specific species
    number _elecNumden0;
    number _elecNumden;
    int _n_N_species;
    int _n_O_species;
    int _n_elec_species;
    string[] _elec_species_names;
    double[] _elec_mol_masses;
    double[] _elec_electronic_energy;
    double[] _elec_R;
    ElectronicSpecies[] _electronicSpecies;
    GasState macro_Q;

    //storage for macro species
    int _n_macro_species;
    string[] _macro_species_names;


    //Storage from two_temperature_air.d file
    PerfectGasMixEOS _pgMixEOS;
    double _R_U_cal = 1.987; // cal/(mole.K)
    number[] _macro_molef; // will be getting mole-fractions from outside, so may be complex
    double[] _macro_particleMass;
    double[] _macro_R;
    double[] _del_hf;
    double[] _Cp_tr_rot;
    number[7] _A; // working storage of coefficients
    number[][] _A_11, _B_11, _C_11, _D_11, _Delta_11, _alpha;
    number[][] _A_22, _B_22, _C_22, _D_22, _Delta_22, _mu;

    @nogc
    void Update_Macro_State(ref GasState macro_state, in GasState Q)
    {
        macro_state.rho = Q.rho;
        macro_state.p = Q.p;
        macro_state.p_e = Q.p_e;
        macro_state.a = Q.a;
        macro_state.T = Q.T;
        macro_state.u = Q.u;
        macro_state.u_modes[] = Q.u_modes[];
        macro_state.T_modes[] = Q.T_modes[];
        macro_state.mu = Q.mu;
        macro_state.k = Q.k;
        macro_state.k_modes[] = Q.k_modes[];
        macro_state.sigma = Q.sigma;
        macro_state.quality = Q.quality;

        number N_massf_sum = 0.0;
        foreach (isp; 0 .. _n_N_species) {
            N_massf_sum += Q.massf[isp];
        }
        macro_state.massf[0] = N_massf_sum;

        number O_massf_sum = 0.0;
        foreach (isp; _n_N_species .. _n_N_species + _n_O_species) {
            O_massf_sum += Q.massf[isp];
        }
        macro_state.massf[1] = O_massf_sum;

        foreach ( isp; 2 .. _n_macro_species) {
            macro_state.massf[isp] = Q.massf[_n_elec_species+isp-2];
        }
    }

    @nogc
    void Update_Electronic_State(ref GasState Q, in GasState macro_state)
    {
        Q.rho = macro_state.rho;
        Q.p = macro_state.p;
        Q.p_e = macro_state.p_e;
        Q.a = macro_state.a;
        Q.T = macro_state.T;
        Q.u = macro_state.u;
        Q.u_modes[] = macro_state.u_modes[];
        Q.T_modes[] = macro_state.T_modes[];
        Q.mu = macro_state.mu;
        Q.k = macro_state.k;
        Q.k_modes[] = macro_state.k_modes[];
        Q.sigma = macro_state.sigma;
        Q.quality = macro_state.quality;

        number N_massf_sum = 0.0;
        foreach (isp; 0 .. _n_N_species) {
            N_massf_sum += Q.massf[isp];
        }
        if (N_massf_sum == 0.0){
            Q.massf[0] = macro_state.massf[0];
            foreach (isp; 1 .. _n_N_species) {
                Q.massf[isp] = 0.0;
            }
        } else {
            number N_massf_factor = macro_state.massf[0]/N_massf_sum;
            foreach (isp; 0 .. _n_N_species) {
                Q.massf[isp] *= N_massf_factor;
            }
        }

        number O_massf_sum = 0.0;
        foreach (isp; _n_N_species .. _n_N_species + _n_O_species) {
            O_massf_sum += Q.massf[isp];
        }
        if (O_massf_sum == 0.0) {
            Q.massf[_n_N_species] = macro_state.massf[1];
            foreach (isp; _n_N_species + 1 .. _n_N_species + _n_O_species) {
                Q.massf[isp] = 0.0;
            }
        } else {
            number O_massf_factor = macro_state.massf[1]/O_massf_sum;
            foreach (isp; _n_N_species .. _n_N_species + _n_O_species) {
                Q.massf[isp] *= O_massf_factor;
            }
        }
    }

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
        Cp *= (R_universal/_macro_mol_masses[isp]);
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
        h *= (R_universal*T/_macro_mol_masses[isp]);
        return h;
    }

    @nogc number transRotEnergy(in GasState macro_Q)
    {
        number e_tr_rot = 0.0;
        foreach (isp; 0 .. _n_macro_species) {
            number h_tr_rot = _Cp_tr_rot[isp]*(macro_Q.T - T_REF) + _del_hf[isp];
            e_tr_rot += macro_Q.massf[isp]*(h_tr_rot - _macro_R[isp]*macro_Q.T);
        }
        return e_tr_rot;
    }

    @nogc number vibSpecHeatConstV(number Tve, int isp)
    {
        return CpFromCurveFits(Tve, isp) - _Cp_tr_rot[isp];
    }

    @nogc number vibSpecHeatConstV(in GasState macro_Q, number Tve)
    {
        number Cv_vib = 0.0;
        foreach (isp; molecularSpecies) {
            Cv_vib += macro_Q.massf[isp] * vibSpecHeatConstV(Tve, isp);
        }
        return Cv_vib;
    }

    @nogc number transRotSpecHeatConstV(int isp)
    {
        return to!number(_Cp_tr_rot[isp] - _macro_R[isp]);
    }

    @nogc number transRotSpecHeatConstV(in GasState macro_Q)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. _n_macro_species) {
            Cv += macro_Q.massf[isp]*transRotSpecHeatConstV(isp);
        }
        return Cv;
    }

    @nogc number vibTemperature(in GasState macro_Q)
    {
        int MAX_ITERATIONS = 20;
        // We'll keep adjusting our temperature estimate
        // until it is less than TOL.
        double TOL = 1.0e-6;
        // Take the supplied T_modes[0] as the initial guess.
        number T_guess = macro_Q.T_modes[0];
        number f_guess = vibEnergy(macro_Q, T_guess) - macro_Q.u_modes[0];
        // Before iterating, check if the supplied guess is
        // good enough. Define good enough as 1/100th of a Joule.
        double E_TOL = 0.01;
        if (fabs(f_guess) < E_TOL) {
            // Given temperature is good enough.
            return macro_Q.T_modes[0];
        }

        // Begin iterating.
        int count = 0;
        number Cv, dT;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            Cv = vibSpecHeatConstV(macro_Q, T_guess);
            dT = -f_guess/Cv;
            T_guess += dT;
            if (fabs(dT) < TOL) {
                break;
            }
            f_guess = vibEnergy(macro_Q, T_guess) - macro_Q.u_modes[0];
            count++;
        }

        if (count == MAX_ITERATIONS) {
            string msg = "The 'vibTemperature' function failed to converge.\n";
            debug {
                msg ~= format("The final value for Tvib was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= macro_Q.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }

        return T_guess;
    }
}

double[11] _macro_mol_masses;
static double[11] del_hf;
static double[7][5][11] thermoCoeffs;
static double[string] A_11, B_11, C_11, D_11;
static double[string] A_22, B_22, C_22, D_22;


version(electronically_specific_gas_test) {
    int main()
    {
        import util.msg_service;

        auto L = init_lua_State();
        doLuaFile(L,"sample-data/electronic-and-macro-species.lua");

        auto gm = new ElectronicallySpecificGas(L);
        auto gd = GasState(gm.n_species,1);

        gd.p=666.0;
        gd.T = 293;
        gd.T_modes[0] = 293;
        gd.massf = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        gm.update_thermo_from_pT(gd);

        assert(isClose(-89784.2, gd.u, 1.0e-6), failedUnitTest());
        assert(isClose(0.00787412, gd.rho, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_rhou(gd);
        assert(isClose(666.0, gd.p, 1.0e-6), failedUnitTest());
        assert(isClose(293, gd.T, 1.0e-6), failedUnitTest());

        return 0;
    }
}
static this()
{
    _macro_mol_masses[0] = 14.0067e-3;
    _macro_mol_masses[1] = 0.01599940;
    _macro_mol_masses[2] = 28.0134e-3;
    _macro_mol_masses[3] = 0.03199880;
    _macro_mol_masses[4] = 0.03000610;
    _macro_mol_masses[5] = 14.0061514e-3;
    _macro_mol_masses[6] = 15.9988514e-3;
    _macro_mol_masses[7] = 28.0128514e-3;
    _macro_mol_masses[8] = 31.9982514e-3;
    _macro_mol_masses[9] = 30.0055514e-3;
    _macro_mol_masses[10] = 0.000548579903e-3;

    /**
     * See Table I in Gnoffo et al.
     */
    del_hf[0]   = 112.951e3*4.184/_macro_mol_masses[0];
    del_hf[1]   = 59.544e3*4.184/_macro_mol_masses[1];
    del_hf[2]  = 0.0;
    del_hf[3]  = 0.0;
    del_hf[4]  = 21.6009e3*4.184/_macro_mol_masses[4];
    del_hf[5]  = 449.709e3*4.184/_macro_mol_masses[5];
    del_hf[6]  = 374.867e3*4.184/_macro_mol_masses[6];
    del_hf[7] = 364.9392e3*4.184/_macro_mol_masses[7];
    del_hf[8] = 280.2099e3*4.184/_macro_mol_masses[8];
    del_hf[9] = 237.3239e3*4.184/_macro_mol_masses[9];
    del_hf[10]  =  0.0;

    thermoCoeffs[2] =
        [
         [  0.3674826e+01, -0.1208150e-02,  0.2324010e-05, -0.6321755e-09, -0.2257725e-12, -0.1061160e+04,  0.23580e+01 ], // 300 -- 1000 K
         [  0.2896319e+01,  0.1515486e-02, -0.5723527e-06,  0.9980739e-10, -0.6522355e-14, -0.9058620e+03,  0.43661e+01 ], // 1000 -- 6000 K
         [     0.3727e+01,     0.4684e-03,    -0.1140e-06,     0.1154e-10,    -0.3293e-15,    -0.1043e+04,  0.46264e+01 ], // 6000 -- 15000 K
         [  0.9637690e+01, -0.2572840e-02,  0.3301980e-06, -0.1431490e-10,  0.2033260e-15, -0.1043000e+04, -0.37587e+02 ], // 15000 -- 25000 K
         [ -0.5168080e+01,  0.2333690e-02, -0.1295340e-06,  0.2787210e-11, -0.2135960e-16, -0.1043000e+04,  0.66217e+02 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[3] =
        [
         [  0.3625598e+01, -0.1878218e-02,  0.7055454e-05, -0.6763513e-08,  0.2155599e-11, -0.1047520e+04,  0.43628e+01 ], // 300 -- 1000 K
         [  0.3621953e+01,  0.7361826e-03, -0.1965222e-06,  0.3620155e-10, -0.2894562e-14, -0.1201980e+04,  0.38353e+01 ], // 1000 -- 6000 K
         [     0.3721e+01,     0.4254e-03,    -0.2835e-07,     0.6050e-12,    -0.5186e-17,    -0.1044e+04,  0.23789e+01 ], // 6000 -- 15000 K
         [  0.3486660e+01,  0.5238420e-03, -0.3912340e-07,  0.1009350e-11, -0.8871830e-17, -0.1044000e+04,  0.48179e+01 ], // 15000 -- 25000 K
         [  0.3961980e+01,  0.3944550e-03, -0.2950580e-07,  0.7397450e-12, -0.6420930e-17, -0.1044000e+04,  0.13985e+01 ]
         ];

    thermoCoeffs[0] =
        [
         [  0.2503071e+01, -0.2180018e-04,  0.5420528e-07, -0.5647560e-10,  0.2099904e-13,  0.5609890e+05,  0.41676e+01 ], // 300 -- 1000 K
         [  0.2450268e+01,  0.1066145e-03, -0.7465337e-07,  0.1879652e-10, -0.1025983e-14,  0.5611600e+05,  0.42618e+01 ], // 1000 -- 6000 K
         [     0.2748e+01,    -0.3909e-03,     0.1338e-06,    -0.1191e-10,     0.3369e-15,     0.5609e+05,  0.28720e+01 ], // 6000 -- 15000 K
         [ -0.1227990e+01,  0.1926850e-02, -0.2437050e-06,  0.1219300e-10, -0.1991840e-15,  0.5609000e+05,  0.28469e+02 ], // 15000 -- 25000 K
         [  0.1552020e+02, -0.3885790e-02,  0.3228840e-06, -0.9605270e-11,  0.9547220e-16,  0.5609000e+05, -0.88120e+02 ]  // 25000 -- 30000 K
         ];

    thermoCoeffs[1] =
        [
         [  0.2946428e+01, -0.1638166e-02,  0.2421031e-05, -0.1602843e-08,  0.3890696e-12,  0.2914760e+05,  0.35027e+01 ], // 300 -- 1000 K
         [  0.2542059e+01, -0.2755061e-04, -0.3102803e-08,  0.4551067e-11, -0.4368051e-15,  0.2923080e+05,  0.49203e+01 ], // 1000 -- 6000 K
         [     0.2546e+01,    -0.5952e-04,     0.2701e-07,    -0.2798e-11,     0.9380e-16,    0.29150e+05,  0.50490e+01 ], // 6000 -- 15000 K
         [ -0.9787120e-02,  0.1244970e-02, -0.1615440e-06,  0.8037990e-11, -0.1262400e-15,  0.2915000e+05,  0.21711e+02 ], // 15000 -- 25000 K
         [  0.1642810e+02, -0.3931300e-02,  0.2983990e-06, -0.8161280e-11,  0.7500430e-16,  0.2915000e+05, -0.94358e+02 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[4] =
        [
         [  0.4045952e+01, -0.3418178e-02,  0.7981919e-05, -0.6113931e-08,  0.1591907e-11,  0.9745390e+04,  0.51497e+01 ], // 300 -- 1000 K
         [  0.3189000e+01,  0.1338228e-02, -0.5289932e-06,  0.9591933e-10, -0.6484793e-14,  0.9828330e+04,  0.66867e+01 ], // 1000 -- 6000 K
         [     0.3845e+01,     0.2521e-03,    -0.2658e-07,     0.2162e-11,    -0.6381e-16,  0.9764000e+04,  0.31541e+01 ], // 6000 -- 15000 K
         [  0.4330870e+01, -0.5808630e-04,  0.2805950e-07, -0.1569410e-11,  0.2410390e-16,  0.9764000e+04,  0.10735e+00 ], // 15000 -- 25000 K
         [  0.2350750e+01,  0.5864300e-03, -0.3131650e-07,  0.6049510e-12, -0.4055670e-17,  0.9764000e+04,  0.14026e+02 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[5] =
        [
         [     0.2727e+01,    -0.2820e-03,     0.1105e-06,    -0.1551e-10,     0.7847e-15,     0.2254e+06,  0.36450e+01 ], // 300 -- 1000 K
         [     0.2727e+01,    -0.2820e-03,     0.1105e-06,    -0.1551e-10,     0.7847e-15,     0.2254e+06,  0.36450e+01 ], // 1000 -- 6000 K , yes same as previous entry
         [     0.2499e+01,    -0.3725e-05,     0.1147e-07,    -0.1102e-11,     0.3078e-16,     0.2254e+06,  0.49500e+01 ], // 6000 -- 15000 K
         [  0.2385610e+01,  0.8349470e-04, -0.5881510e-08,  0.1884970e-12, -0.1611950e-17,     0.2254e+06,  0.56462e+01 ], // 15000 -- 25000 K
         [  0.2228570e+01,  0.1245820e-03, -0.8763570e-08,  0.2620400e-12, -0.2167420e-17,     0.2554e+06,  0.67811e+01 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[6] =
        [
         [  0.2498479e+01,  0.1141097e-04, -0.2976139e-07,  0.3224653e-10, -0.1237551e-13,  0.1879490e+06,  0.43864e+01 ], // 300 -- 1000 K
         [  0.2506048e+01, -0.1446424e-04,  0.1244604e-07, -0.4685847e-11,  0.6554887e-15,  0.1879470e+06,  0.43480e+01 ], // 1000 -- 6000 K
         [     0.2944e+01,    -0.4108e-03,     0.9156e-07,    -0.5848e-11,     0.1190e-15,  0.1879000e+06,  0.17500e+01 ], // 6000 -- 15000 K
         [  0.1278400e+01,  0.4086590e-03, -0.2173100e-07,  0.3325180e-12,  0.6316040e-18,  0.1879000e+06,  0.12761e+02 ], // 15000 -- 25000 K
         [  0.1288860e+01,  0.4334250e-03, -0.2675820e-07,  0.6215900e-12, -0.4513150e-17,  0.1879000e+06,  0.12604e+02 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[7] =
        [
         [  0.3397000e+01,  0.4525000e-03,  0.1272000e-06, -0.3879000e-10,  0.2459000e-14,  0.1826000e+06,  0.36535e+01 ], // 300 -- 1000 K
         [  0.3397390e+01,  0.4524870e-03,  0.1272300e-06, -0.3879340e-10,  0.2458950e-14,  0.1826000e+06,  0.42050e+01 ], // 1000 -- 6000 K
         [  0.3369950e+01,  0.8628820e-03, -0.1275510e-06,  0.8087120e-11, -0.1879660e-15,  0.1826000e+06,  0.40730e+01 ], // 6000 -- 15000 K
         [  0.4394250e+01,  0.1886760e-03, -0.7127180e-08, -0.1751090e-12,  0.6717580e-17,  0.1826000e+06, -0.23693e+01 ], // 15000 -- 25000 K
         [  0.3949290e+01,  0.3679480e-03, -0.2691020e-07,  0.6711050e-12, -0.5824370e-17,  0.1826000e+06,  0.65472e+00 ] // 25000 -- 30000 K
         ];

    thermoCoeffs[8] =
        [
         [  0.3243000e+01,  0.1174000e-02, -0.3900000e-06,  0.5437000e-10, -0.2392000e-14,  0.1400000e+06,  0.59250e+01 ], // 300 -- 1000 K
         [  0.3242980e+01,  0.1173910e-02, -0.3900420e-06,  0.5437260e-10, -0.2392320e-14,  0.1400000e+06,  0.59250e+01 ], // 1000 -- 6000 K
         [  0.5168650e+01, -0.8619690e-03,  0.2041410e-06, -0.1300410e-10,  0.2494210e-15,  0.1400000e+06, -0.52960e+01 ], // 6000 -- 15000 K
         [ -0.2801710e+00,  0.1667410e-02, -0.1210740e-06,  0.3211290e-11, -0.2834890e-16,  0.1400000e+06,  0.31013e+02 ], // 15000 -- 25000 K
         [  0.2044550e+01,  0.1031320e-02, -0.7404630e-07,  0.1925750e-11, -0.1746100e-16,  0.1400000e+06,  0.14310e+02 ] // 25000 -- 30000 K
         ];
    thermoCoeffs[9] =
        [
         [  0.3668506e+01, -0.1154458e-02,  0.2175561e-05, -0.4822747e-09, -0.2784791e-12,  0.1180340e+06,  0.37852e+01 ],
         [  0.2888549e+01,  0.1521712e-02, -0.5753124e-06,  0.1005108e-09, -0.6604429e-14,  0.1181920e+06,  0.51508e+01 ],
         [  0.2214170e+01,  0.1776060e-02, -0.4303860e-06,  0.4173770e-10, -0.1282890e-14,  0.1181920e+06,  0.83904e+01 ],
         [ -0.3324050e+01,  0.2441960e-02, -0.1905720e-06,  0.6858000e-11, -0.9911240e-16,  0.1181920e+06, -0.11079e+02 ],
         [ -0.4348760e+01,  0.2401210e-02, -0.1445990e-06,  0.3381320e-11, -0.2825510e-16,  0.1181920e+06,  0.65896e+02 ]
         ];

    thermoCoeffs[10] =
        [
         [  0.2500000e+01,            0.0,            0.0,            0.0,            0.0, -0.7453750e+03, -0.11734e+02],
         [  0.2500000e+01,            0.0,            0.0,            0.0,            0.0, -0.7453750e+03, -0.11734e+02],
         [  0.2508e+01   , -0.6332e-05   ,  0.1364e-08   , -0.1094e-12   ,  0.2934e-17   , -0.7450000e+03, -0.11734e+02],
         [  0.250010e+01 , -0.311281e-09 ,  0.357207e-13 , -0.1603670e-17,  0.250707e-22 , -0.7450000e+03, -0.11734e+02],
         [  0.250010e+01 ,  0.301577e-09 , -0.226204e-13 ,  0.667344e-18 , -0.689169e-23 , -0.7450000e+03, -0.11734e+02]
         ];


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

// Everything from here is old. Will delete when appropriate
/*
class ElectronicallySpecificGas : GasModel {
public:
    this(lua_State *L)
    {
        auto nElecSpecies = getInt(L, "number_electronic_species");
        _n_species = cast(uint) nElecSpecies;

        _n_modes = 1;
        _species_names.length = _n_species;
        _numden.length = _n_species;

        _s1 = getDouble(L, "s1");
        _T1 = getDouble(L, "T1");
        _p1 = getDouble(L, "p1");

        lua_getfield(L, "electronic_species");
        foreach (isp; 0 .. _n_species) {
            lua_rawgeti(L, -1, isp);
            if (lua_isnil(L, -1)) {
                string msg = format("There was an error when attempting to information about electronic-species %d.\n", isp);
                throw new Error(msg);
            }
            _electronicSpecies ~= createElectronicSpecies(L);
            lua_pop(L, 1);
            _species_names[isp] = _electronicSpecies[$-1].name;
            _mol_masses ~= _electronicSpecies[$-1].mol_mass;
            _R ~= R_universal/_electronicSpecies[$-1].mol_mass;
            _lowerLevel ~= _electronicSpecies[$-1].lowerLevel;
            _upperLevel ~= _electronicSpecies[$-1].upperLevel;
            _group_degeneracy ~= _electronicSpecies[$-1].group_degeneracy;
            _dof ~= _electronicSpecies[$-1].dof;
            _electronic_energy ~= _electronicSpecies[$-1].electronic_energy;
        }
        lua_pop(L, 1);
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "ElectronicSpeciesGas=()";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.rho = Q.p/(R_mix*Q.T);
        auto Cv_heavy = heavy_Cv(Q);
        Q.u = Cv_heavy*Q.T;
        auto Cv_electron = electron_Cv(Q);
        Q.u_modes[0]=Cv_electron*Q.T_modes[0];
    }

    override void update_thermo_from_rhou(ref GasState Q)
    {
        auto Cv_heavy = heavy_Cv(Q);
        Q.T = (Q.u)/Cv_heavy;
        auto Cv_electron = electron_Cv(Q);
        Q.T_modes[0] = Q.u_modes[0]/Cv_electron;
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
    }

    override void update_thermo_from_rhoT(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.p = Q.rho*R_mix*Q.T;
        auto Cv_heavy = heavy_Cv(Q);
        Q.u = Cv_heavy*Q.T;
        auto Cv_electron = electron_Cv(Q);
        Q.u_modes[0] = Cv_electron*Q.T_modes[0];
    }

    override void update_thermo_from_rhop(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        Q.T = Q.p/(Q.rho*R_mix);
        auto Cv_heavy = heavy_Cv(Q);
        Q.u = Cv_heavy*Q.T;
        auto Cv_electron = electron_Cv(Q);
        Q.u_modes[0] = Cv_electron*Q.T_modes[0];
    }

    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        Q.T = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
        update_thermo_from_pT(Q);
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        Q.T = h / _Cp;
        Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.T/_T1)));
        update_thermo_from_pT(Q);
    }

    override void update_sound_speed(ref GasState Q)
    {
        auto R_mix = gas_constant(Q);
        auto Cv_mix = Cv(Q);
        auto gamma_mix = (R_mix/Cv_mix) + 1;
        Q.a = sqrt(gamma_mix*R_mix*Q.T);
    }

    override void update_trans_coeffs(ref GasState Q)
    {
        Q.mu = 0.0;
        Q.k = 0.0;
    }

    override number dudT_const_v(in GasState Q) const
    {
        number Cv = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv += Q.massf[isp] * (_dof[isp]/2.0) * _R[isp];
        }
        return Cv;
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return gas_constant(Q) + dudT_const_v(Q);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        number R = gas_constant(Q);
        return R*Q.T;
    }
    override number gas_constant(in GasState Q) const
    {
        number R_mix = 0.0;
        foreach (isp; 0 .. _n_species) {
            R_mix += Q.massf[isp]*_R[isp];
        }
        return R_mix;
    }
    override number internal_energy(in GasState Q)
    {
        auto uNoneq = energyInNoneq(Q);
        return Q.u + Q.u_modes[0] + uNoneq;
    }
    override number enthalpy(in GasState Q) const
    {
        return Q.u + sum(Q.u_modes) + Q.p/Q.rho;
    }
    override number entropy(in GasState Q) const
    {
        return _s1 + _Cp * log(Q.T/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    ElectronicSpecies[] _electronicSpecies;
    double[] _R;
    int[] _lowerLevel;
    int[] _upperLevel;
    int[] _group_degeneracy;
    int[] _dof;
    number[] _numden; //  #/cm^3

    // Thermodynamic constants
    double _Rgas; // J/kg/K
    number _Cv1; // J/kg/K
    number _Cv2; //J/kg/K
    double _Cv;
    double _Cp;
    double _gamma;   // ratio of specific heats
    // Reference values for entropy
    double _s1;  // J/kg/K
    double _T1;  // K
    double _p1;  // Pa

    //physical cconstants
    double _pi = 3.14159265359; //honestly, this should probably be defined in physical constants
    double _me = 9.10938356e-28; //electron mass in g
    double _kb = 1.3807e-16; //Boltzmann constant in cm^2 g s^-1 K^-1
    double _e = 4.8032e-10; // electron charge, cm^(3/2) g s^-2 K^-1


    //create array of coefficients for appleton and Bray energy relaxation method
    //Rows in order N2, N, O2, O
    //each element: a,b,c in equation cross_section = a + b*Te + c*Te^2
    double[][] AB_coef = [[7.5e-20, 5.5e-24, -1e-28],
                            [5.0e-20, 0.0, 0.0],
                            [2.0e-20, 6.0e-24, 0.0],
                            [1.2e-20, 1.7e-24, -2e-29]];

    number _cs;
    number _Nsum;
    number _Osum;
    number _sumterm;

    @nogc number energyInNoneq(ref const(GasState) Q) const
    {
        number uNoneq = 0.0;
        foreach (isp; 0 .. _n_species) {
            uNoneq += Q.massf[isp] * _electronicSpecies[isp].electronic_energy;
        }
        return uNoneq;
    }

    @nogc number heavy_Cv(ref GasState Q)
    {
        _Cv1 = 0.0;
        foreach (isp; 0 .. _n_species-3){
            _Cv1 += Q.massf[isp] * (_dof[isp]/2.0) * _R[isp];
        }
        foreach (isp; _n_species-2 .. _n_species) {
            _Cv1 += Q.massf[isp] * (_dof[isp]/2.0) * _R[isp];
        }
        return _Cv1;
    }

    @nogc number electron_Cv(ref GasState Q)
    {
        number elec_cv = Q.massf[_n_species-3]*(_dof[_n_species-3]/2.0) * _R[_n_species-3];
        return elec_cv;
    }



    @nogc number electronMassf(ref GasState Q) {
        return Q.massf[_n_species-3];
    }

    @nogc number heavyMassf(ref GasState Q) {
        return 1.0 - electronMassf(Q);
    }
}


version(electronically_specific_gas_test) {
    int main()
    {
        import util.msg_service;

        string filename = "sample-data/electronic_composition.lua";

        auto L = init_lua_State();
        doLuaFile(L, relativePath(filename));
        auto gm = new ElectronicallySpecificGas(L);
        auto gd = GasState(gm.n_species,1);

        // gd.massf[] = 0.0;
        // gd.massf[0] = 0.037041674288877; //initialises massf of NI
        // gd.massf[9] = 0.010577876366622; //initialises massf of OI
        // gd.massf[19] = 0.74082290750449; //N2
        // gd.massf[20] = 0.21155752733244; //O2
        // gd.massf[18] = 1.0 - (gd.massf[0] + gd.massf[9] + gd.massf[19] + gd.massf[20]); //tiny massf for free electron
        // gd.massf = [0.0313603, 0.00492971, 0.000741705, 1.06916e-06, 4.90114e-07,
        //                 2.46998e-07, 9.58454e-08, 6.6456e-07, 6.41328e-06, 0.010005,
        //                 0.000565079, 8.59624e-06, 2.58411e-07, 9.00322e-08, 5.80925e-08,
        //                 3.67871e-08, 9.06483e-08, 4.16313e-07, 1.4773e-08, 0.740823, 0.211558];
        gd.massf[] = 0;
        gd.massf[gm.n_species-3] = 1e-8;
        gd.massf[gm.n_species-2] = 0.78;
        gd.massf[gm.n_species-1] = 1.0 - 0.78 - 1e-8;

        gd.p = 100000.0;
        gd.T = 300.0;
        gd.T_modes[0]=4000.0;
        gm.update_thermo_from_pT(gd);

        double massfsum=0.0;
        foreach(number eachmassf;gd.massf) {
            massfsum += eachmassf;
        }
        assert(isClose(massfsum, 1.0, 1e-2), failedUnitTest());

        return 0;
    }

}
*/
