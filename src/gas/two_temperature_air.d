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

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.thermo.perf_gas_mix_eos;


immutable double T_REF = 298.15; // K
static bool[string] molecularSpecies;
static double[string] del_hf;
static double[7][5][string] thermoCoeffs;

static this()
{
    molecularSpecies["N2"] = true;
    molecularSpecies["O2"] = true;
    molecularSpecies["NO"] = true;
    molecularSpecies["N2+"] = true;
    molecularSpecies["O2+"] = true;
    molecularSpecies["NO+"] = true;

    /**
     * See Table B1 in Gupta et al.
     */
    del_hf["N"]   = 112.973*4184;
    del_hf["O"]   = 59.553*4184;
    del_hf["N2"]  = 0.0;
    del_hf["O2"]  = 0.0;
    del_hf["NO"]  = 21.580*4184;
    del_hf["N+"]  = 449.840*4184;
    del_hf["O+"]  = 374.949*4184;
    del_hf["N2+"] = 360.779*4184;
    del_hf["O2+"] = 279.849*4184;
    del_hf["e-"]  =  0.0;

    thermoCoeffs["N2"] = 
	[
	 [  0.36748e+01, -0.12081e-02,  0.23240e-05, -0.63218e-09, -0.22577e-12, -0.10430e+04,  0.23580e+01 ], // 300 -- 1000 K
	 [  0.32125e+01,  0.10137e-02, -0.30467e-06,  0.41091e-10, -0.20170e-14, -0.10430e+04,  0.43661e+01 ], // 1000 -- 6000 K
	 [  0.31811e+01,  0.89745e-03, -0.20216e-06,  0.18266e-10, -0.50334e-15, -0.10430e+04,  0.46264e+01 ], // 6000 -- 15000 K
	 [  0.96377e+01, -0.25728e-02,  0.33020e-06, -0.14315e-10,  0.20333e-15, -0.10430e+04, -0.37587e+02 ], // 15000 -- 25000 K
	 [ -0.51681e+01,  0.23337e-02, -0.12953e-06,  0.27872e-11, -0.21360e-16, -0.10430e+04,  0.66217e+02 ] // 25000 -- 30000 K
	 ];

    thermoCoeffs["O2"] = 
	[
	 [  0.36146e+01, -0.18598e-02,  0.70814e-05, -0.68070e-08,  0.21628e-11, -0.10440e+04,  0.43628e+01 ], // 300 -- 1000 K
	 [  0.35949e+01,  0.75213e-03, -0.18732e-06,  0.27913e-10, -0.15774e-14, -0.10440e+04,  0.38353e+01 ], // 1000 -- 6000 K
	 [  0.38599e+01,  0.32510e-03, -0.92131e-08, -0.78684e-12,  0.29426e-16, -0.10440e+04,  0.23789e+01 ], // 6000 -- 15000 K
	 [  0.34867e+01,  0.52384e-03, -0.39123e-07,  0.10094e-11, -0.88718e-17, -0.10440e+04,  0.48179e+01 ], // 15000 -- 25000 K
	 [  0.39620e+01,  0.39446e-03, -0.29506e-07,  0.73975e-12, -0.64209e-17, -0.10440e+04,  0.13985e+01 ]
	 ];

    thermoCoeffs["N"] = 
	[
	 [  0.25031e+01, -0.21800e-04,  0.54205e-07, -0.56476e-10,  0.20999e-13,  0.56130e+05,  0.41676e+01 ], // 300 -- 1000 K
	 [  0.24820e+01,  0.69258e-04, -0.63065e-07,  0.18387e-10, -0.11747e-14,  0.56130e+05,  0.42618e+01 ], // 1000 -- 6000 K
	 [  0.27480e+01, -0.39090e-03,  0.13380e-06, -0.11910e-10,  0.33690e-15,  0.56130e+05,  0.28720e+01 ], // 6000 -- 15000 K
	 [ -0.12280e+01,  0.19268e-02, -0.24370e-06,  0.12193e-10, -0.19918e-15,  0.56130e+05,  0.28469e+02 ], // 15000 -- 25000 K
	 [  0.15520e+02, -0.38858e-02,  0.32288e-06, -0.96053e-11,  0.95472e-16,  0.56130e+05, -0.88120e+02 ]  // 25000 -- 30000 K
	 ];

    thermoCoeffs["O"] = 
	[
	 [  0.28236e+01, -0.89478e-03,  0.83060e-06, -0.16837e-09, -0.73205e-13,  0.29150e+05,  0.35027e+01 ], // 300 -- 1000 K
	 [  0.25421e+01, -0.27551e-04, -0.31028e-08,  0.45511e-11, -0.43681e-15,  0.29150e+05,  0.49203e+01 ], // 1000 -- 6000 K
	 [  0.25460e+01, -0.59520e-04,  0.27010e-07, -0.27980e-11,  0.93800e-16,  0.29150e+05,  0.50490e+01 ], // 6000 -- 15000 K
	 [ -0.97871e-02,  0.12450e-02, -0.16154e-06,  0.80380e-11, -0.12624e-15,  0.29150e+05,  0.21711e+02 ], // 15000 -- 25000 K
	 [  0.16428e+02, -0.39313e-02,  0.29840e-06, -0.81613e-11,  0.75004e-16,  0.29150e+05, -0.94358e+02 ], // 25000 -- 30000 K
	 ];
	
    thermoCoeffs["NO"] =
	[
	 [  0.35887e+01, -0.12479e-02,  0.39786e-05, -0.28651e-08,  0.63015e-12,  0.97640e+04,  0.51497e+01 ], // 300 -- 1000 K
	 [  0.32047e+01,  0.12705e-02, -0.46603e-06,  0.75007e-10, -0.42314e-14,  0.97640e+04,  0.66867e+01 ], // 1000 -- 6000 K
	 [  0.38543e+01,  0.23409e-03, -0.21354e-07,  0.16689e-11, -0.49070e-16,  0.97640e+04,  0.31541e+01 ], // 6000 -- 15000 K
	 [  0.43309e+01, -0.58086e-04,  0.28059e-07, -0.15694e-11,  0.24104e-16,  0.97640e+04,  0.10735e+00 ], // 15000 -- 25000 K
	 [  0.23507e+01,  0.58643e-03, -0.31316e-07,  0.60495e-12, -0.40557e-17,  0.97640e+04,  0.14026e+02 ], // 25000 -- 30000 K
	 ];
}

class TwoTemperatureAir : GasModel {
public:
    this(string model)
    {
	switch (model) {
	case "5-species":
	    _n_species = 5;
	    _n_modes = 1;
	    _species_names.length = 5;
	    _mol_masses.length = 5;
	    _species_names[0] = "N"; _mol_masses[0] = 14.0067e-3;
	    _species_names[1] = "O"; _mol_masses[1] = 0.01599940;
	    _species_names[2] = "N2"; _mol_masses[2] = 28.0134e-3;
	    _species_names[3] = "O2"; _mol_masses[3] = 0.03199880;
	    _species_names[4] = "NO"; _mol_masses[4] = 0.03000610;
	    create_species_reverse_lookup();
	    break;
	case "7-species":
	    _n_species = 7;
	    _n_modes = 1;
	    _species_names.length = 7;
	    _mol_masses.length = 7;
	    _species_names[0] = "N"; _mol_masses[0] = 14.0067e-3;
	    _species_names[1] = "O"; _mol_masses[1] = 0.01599940; 
	    _species_names[2] = "N2"; _mol_masses[2] = 28.0134e-3; 
	    _species_names[3] = "O2"; _mol_masses[3] = 0.03199880;
	    _species_names[4] = "NO"; _mol_masses[4] = 0.03000610;
	    _species_names[5] = "NO+"; _mol_masses[5] = 30.0055514e-3;
	    _species_names[6] = "e-"; _mol_masses[6] = 0.000548579903e-3;
	    create_species_reverse_lookup();
	    break;
	case "11-species":
	    _n_species = 11;
	    _n_modes = 1;
	    _species_names.length = 11;
	    _mol_masses.length = 11;
	    _species_names[0] = "N"; _mol_masses[0] = 14.0067e-3;
	    _species_names[1] = "O"; _mol_masses[1] = 0.01599940; 
	    _species_names[2] = "N2"; _mol_masses[2] = 28.0134e-3;
	    _species_names[3] = "O2"; _mol_masses[3] = 0.03199880;
	    _species_names[4] = "NO"; _mol_masses[4] = 0.03000610;
	    _species_names[5] = "N+"; _mol_masses[5] = 14.0061514e-3;
	    _species_names[6] = "O+"; _mol_masses[6] = 15.9988514e-3;
	    _species_names[7] = "N2+"; _mol_masses[7] = 28.0128514e-3;
	    _species_names[8] = "O2+"; _mol_masses[8] = 31.9982514e-3;
	    _species_names[9] = "NO+"; _mol_masses[9] = 30.0055514e-3;
	    _species_names[10] = "e-"; _mol_masses[10] = 0.000548579903e-3;
	    create_species_reverse_lookup();
	    break;
	default:
	    string errMsg = format("The model name '%s' is not a valid selection for the two-temperature air model.", model);
	    throw new Error(errMsg);
	}
	
	_R.length = _n_species;
	foreach (isp; 0 .. _n_species) {
	    _R[isp] = R_universal/_mol_masses[isp];
	}
	_pgMixEOS = new PerfectGasMixEOS(_R);

	_del_hf.length = _n_species;
	foreach (isp; 0 .. _n_species) {
	    _del_hf[isp] = del_hf[_species_names[isp]];
	}

	_Cp_tr_rot.length = _n_species;
	foreach (isp; 0 .. _n_species) {
	    _Cp_tr_rot[isp] = (5./2.)*_R[isp];
	    if (_species_names[isp] in molecularSpecies) {
		_Cp_tr_rot[isp] += _R[isp];
		_molecularSpecies ~= isp;
	    }
	}
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "TwoTemperatureAir";
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
	double sumA = 0.0;
	double sumB = 0.0;
	foreach (isp; 0 .. _n_species) {
	    sumA += Q.massf[isp]*(_Cp_tr_rot[isp] - _del_hf[isp]);
	    sumB += Q.massf[isp]*_Cp_tr_rot[isp];
	}
	Q.T = (Q.u - sumA)/sumB;
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

    override void update_thermo_from_ps(GasState Q, double s)
    {
	throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureAir.");
    }

    override void update_thermo_from_hs(GasState Q, double h, double s)
    {
	throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureAir.");
    }

    override void update_sound_speed(GasState Q)
    {
	
	
    }

    override void update_trans_coeffs(GasState Q)
    {

    }

    override double dudT_const_v(in GasState Q)
    {
	return 0.0;
    }
    override double dhdT_const_p(in GasState Q)
    {
	return 0.0;
    }
    override double dpdrho_const_T(in GasState Q)
    {
	return 0.0;
    }
    override double gas_constant(in GasState Q)
    {
	return  mass_average(Q, _R);
    }
    override double internal_energy(in GasState Q)
    {
	return Q.u + Q.u_modes[0]; // all internal energy
    }
    override double enthalpy(in GasState Q)
    {
	return 0.0;
    }
    override double entropy(in GasState Q)
    {
	throw new GasModelException("entropy not implemented in TwoTemperatureNitrogen.");
    }

private:
    PerfectGasMixEOS _pgMixEOS;
    double[] _R;
    double[] _del_hf;
    double[] _Cp_tr_rot;
    double[7] _A; // working storage of coefficients
    int[] _molecularSpecies;

    /**
     * In order to get smooth variations in thermodynamic properties
     * at the edges of the polynomial breaks, Gupta et al. suggest
     * taking a linear average of the coefficient values near the
     * the polynomial breaks. They describe this process in words
     * on p.14 of their report, and give a Fortran implementation in
     * Appendix C (pp 40 & 41).
     */
    void determineCoefficients(double T, int isp)
    {
	string spName = species_name(isp);
	if (T < 800.0) {
	    _A[] = thermoCoeffs[spName][0][];
	}
	if (T >= 800.0 && T <= 1200.0) {
	    double wB = (1./400.0)*(T - 800.0);
	    double wA = 1.0 - wB;
	    foreach (i; 0 .. 7) {
		_A[isp] = wA*thermoCoeffs[spName][0][i] + wB*thermoCoeffs[spName][1][i];
	    }
	}
	if (T > 1200.0 && T < 5500.0) {
	    _A[] = thermoCoeffs[spName][1][];
	}
	if (T >= 5500.0 && T <= 6500.0) {
	    double wB = (1./1000.0)*(T - 5500.0);
	    double wA = 1.0 - wB;
	    foreach (i; 0 .. 7) {
		_A[isp] = wA*thermoCoeffs[spName][1][i] + wB*thermoCoeffs[spName][2][i];
	    }
	}
	if (T > 6500.0 && T < 14500.0) {
	    _A[] = thermoCoeffs[spName][2][];
	}
	if (T >= 14500.0 && T <= 15500.0) {
	    double wB = (1./1000.0)*(T - 14500.0);
	    double wA = 1.0 - wB;
	    foreach (i; 0 .. 7) {
		_A[isp] = wA*thermoCoeffs[spName][2][i] + wB*thermoCoeffs[spName][3][i];
	    }
	}
	if (T > 15500.0 && T < 24500.0) {
	    _A[] = thermoCoeffs[spName][3][];
	}
	if (T >= 24500.0 && T <= 25500.0) {
	    double wB = (1./1000.0)*(T - 14500.0);
	    double wA = 1.0 - wB;
	    foreach (i; 0 .. 7) {
		_A[isp] = wA*thermoCoeffs[spName][3][i] + wB*thermoCoeffs[spName][4][i];
	    }
	}
	if ( T > 25500.0) {
	    _A[] = thermoCoeffs[spName][3][];
	}
    }

    double CpFromCurveFits(double T, int isp)
    {
	/* Assume that Cp is constant off the edges of the curve fits.
	 * For T < 300.0, this is a reasonable assumption to make.
	 * For T > 30000.0, that assumption might be questionable.
	 */
	if (T < 300.0) {
	    determineCoefficients(300.0, isp);
	}
	if (T > 30000.0) {
	    determineCoefficients(30000.0, isp);
	}
	// For all other cases, use supplied temperature
	determineCoefficients(T, isp);
	double T2 = T*T;
	double T3 = T2*T;
	double T4 = T3*T;
	double Cp = (_A[0] + _A[1]*T + _A[2]*T2 + _A[3]*T3 + _A[4]*T4);
	Cp *= (R_universal/_mol_masses[isp]);
	return Cp;

    }

    double enthalpyFromCurveFits(double T, int isp)
    {
	/* Gupta et al suggest that specific enthalpy below 300 K
	 * should be calculated assuming that the specific heat at
	 * constant pressure is constant. That suggestion is used
	 * here. Additionally, we apply the same idea to temperatures
	 * above 30000 K. We will assume that specific heat is constant
	 * above 30000 K.
	 */
	if (T < 300.0) {
	    double Cp = CpFromCurveFits(300.0, isp);
	    double h = Cp*(T - T_REF) + _del_hf[isp];
	}
	if ( T > 30000.0) {
	    double Cp = CpFromCurveFits(300000.0, isp);
	    double h = Cp*(T - 30000.0) + enthalpyFromCurveFits(30000.0, isp);
	}
	// For all other, determine coefficients and compute specific enthalpy.
	determineCoefficients(T, isp);
	double T2 = T*T;
	double T3 = T2*T;
	double T4 = T3*T;
	double h = _A[0] + _A[1]*T/2. + _A[2]*T2/3. + _A[3]*T3/4. + _A[4]*T4/5. + _A[5]/T;
	h *= (R_universal*T/_mol_masses[isp]);
	return h;
    }

    double vibEnergy(double Tve, int isp)
    {
	double h_at_Tve = enthalpyFromCurveFits(Tve, isp);
	double h_ve = h_at_Tve - _Cp_tr_rot[isp]*(Tve - T_REF) - _del_hf[isp];
	return h_ve;
    }

    double vibEnergy(in GasState Q, double Tve)
    {
	double e_ve = 0.0;
	foreach (isp; _molecularSpecies) {
	    e_ve += Q.massf[isp] * vibEnergy(Tve, isp);
	}
	return e_ve;
    }

    double transRotEnergy(in GasState Q)
    {
	double e_tr_rot = 0.0;
	foreach (isp; 0 .. _n_species) {
	    e_tr_rot += Q.massf[isp]*(_Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp]);
	}
	return e_tr_rot;
    }

    double vibSpecHeatConstV(double Tve, int isp)
    {
	return CpFromCurveFits(Tve, isp) - _Cp_tr_rot[isp];
    }

    double vibSpecHeatConstV(in GasState Q, double Tve)
    {
	double Cv_vib = 0.0;
	foreach (isp; _molecularSpecies) {
	    Cv_vib += Q.massf[isp] * vibSpecHeatConstV(Tve, isp);
	}
	return Cv_vib;
    }

    double vibTemperature(in GasState Q)
    {
	int MAX_ITERATIONS = 20;
	// We'll keep adjusting our temperature estimate
	// until it is less than TOL.
	double TOL = 1.0e-6;
	
	// Take the supplied T_modes[0] as the initial guess.
	double T_guess = Q.T_modes[0];
	double f_guess = vibEnergy(Q, T_guess) - Q.u_modes[0];
	// Before iterating, check if the supplied guess is
	// good enough. Define good enough as 1/100th of a Joule.
	double E_TOL = 0.01;
	if (fabs(f_guess) < E_TOL) {
	    // Given temperature is good enough.
	    return Q.T_modes[0];
	}

	// Begin iterating.
	int count = 0;
	double Cv, dT;
	foreach (iter; 0 .. MAX_ITERATIONS) {
	    Cv = vibSpecHeatConstV(Q, T_guess);
	    dT = -f_guess/Cv;
	    T_guess += dT;
	    if (fabs(dT) < TOL) {
		break;
	    }
	    count++;
	}
	
	if (count == (MAX_ITERATIONS-1)) {
	    string msg = "The 'vibTemperature' function failed to converge.\n";
	    msg ~= format("The final value for Tvib was: %12.6f\n", T_guess);
	    msg ~= "The supplied GasState was:\n";
	    msg ~= Q.toString() ~ "\n";
	    throw new GasModelException(msg);
	}

	return T_guess;
    }

}

