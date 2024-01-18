module gas.thermo.energy_modes;

import std.algorithm;
import std.math;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas;
import gas.thermo.thermo_model;
import gas.thermo.cea_thermo_curves;

immutable double T_REF = 298.15;

interface InternalEnergy {
    // NOTE: These might need an entire gas state rather than just T
    // depending on the model
    @nogc abstract number energy(number T);
    @nogc abstract number Cv(number T);
    @nogc abstract number entropy(number T);
}


// Vibrational energy modes

// A non-truncated harmonic oscillator. Care should be taken
// that the temperature is lower than the dissociation temperature
// as this model makes no effort to account for the molecules breaking apart.
// Ref: Vincenti & Kruger chapter IV, section 12
class HarmonicOscillator : InternalEnergy {
public:
    this(double theta_v, double R)
    {
        m_theta_v = theta_v;
        mR = R;
    }

    @nogc override number energy(number T)
    {
        // Ref: Vincenti & Kruger eq. 12.12
        return mR * m_theta_v / ( exp(m_theta_v / T) - 1 );
    }

    @nogc override number Cv(number T)
    {
        // Ref: Vincenti & Kruger eq. 12.13
        number exp_t = exp(m_theta_v / T);
        number tmp = exp_t / ((exp_t - 1) * (exp_t - 1));
        return mR * m_theta_v * m_theta_v / (T*T) * tmp;
    }

    @nogc override number entropy(number T)
    {
        return mR * (-log(1 - exp(-m_theta_v / T)) + (m_theta_v / T) / (exp(m_theta_v / T) - 1.0));
    }

private:
    double m_theta_v;
    double mR;
}

class TruncatedHarmonicOscillator : InternalEnergy {
    /*
    A truncated harmonic oscillator internal energy 
    */

    this (double theta_v, double theta_D, double R) {
        _theta_v = theta_v;
        _theta_D = theta_D;
        _R = R;
    }

    @nogc override number energy(number T) {
        // Equation 5.7 of Rowan's thesis
        number vib = _theta_v / (exp(_theta_v / T) - 1.);
        number diss = _theta_D / (exp(_theta_D / T) - 1.);
        return _R * (vib - diss);
    }

    @nogc override number Cv(number T){
        // Equation 5.9 of Rowan's thesis
        number exp_v = exp(_theta_v/T);
        number exp_D = exp(_theta_D/T);
        number tmp_vib = exp_v / ((exp_v - 1.) * (exp_v - 1.));
        number tmp_diss = exp_D / ((exp_D - 1.) * (exp_D - 1.));
        return _R * (_theta_v*_theta_v/(T*T)*tmp_vib - _theta_D*_theta_D/(T*T)*tmp_diss);
    }

    @nogc override number entropy(number T){
        throw new Error("TruncatedHarmonicOscillator.entropy is not implemented yet");
    }

private:
    double _theta_v; // characteristic vibration temperature (K)
    double _theta_D; // characteristic dissociation temperature (K)
    double _R;       // specific gas constant (J/kg/K)
}


class ZeroVibration : InternalEnergy
{
public:
    this () {}
    @nogc override number energy(number T) { return to!number(0.0); }
    @nogc override number Cv(number T) { return to!number(0.0); }
    @nogc override number entropy(number T) { return to!number(0.0); }
}

/* 
Compute the vibrational part of the cea energy curve
by subtracting out the known trans/rot and electronic
components.
 */    
class CEAThermoCurveVib : InternalEnergy {
public:
    this(lua_State* L, double R) {
        _electronic_component = create_electronic_energy_model(L, R); 
        lua_getfield(L, -1, "thermoCoeffs");
        _thermo_curve = new CEAThermoCurve(L, R);
        lua_pop(L, 1);

        // since this is a vibration mode, we'll assume that it is a molecule
        string molecule_type = getString(L, -1, "molecule_type");
        _Cp = (molecule_type == "linear") ? (7.0/2.)*R : (8.0/2.)*R;
    }

    @nogc override number energy(number T) {
        number h = _thermo_curve.eval_h(T);
        return h - _Cp*(T) - (_electronic_component.energy(T));
    }

    @nogc override number Cv(number T) { 
        return _thermo_curve.eval_Cp(T) - _Cp - _electronic_component.Cv(T);
    }

    @nogc override number entropy(number T) {
        throw new Error("CEAThermoCurveVib.entropy is not implemented yet");
    }

private:
    InternalEnergy _electronic_component;
    CEAThermoCurve _thermo_curve;
    double _Cp;
}


// Electronic Energy modes
class TwoLevelElectronic : InternalEnergy
{
public:
    this(double theta_e, double R, int g0, int g1)
    {
        m_g0 = g0;
        m_g1 = g1;
        m_theta_e = theta_e;
        mR = R;
    }

    @nogc override number energy(number T)
    {
        // Ref: Vincenti & Kruger Eq. 11.3
        double g1g0 = m_g1/m_g0;
        number numerator = g1g0 * exp(-m_theta_e / T);
        m_denom = 1 + g1g0 * exp(-m_theta_e/T);
        return mR * m_theta_e * numerator / m_denom;
    }

    @nogc override number Cv(number T)
    {
        // Ref: Vincenti & Kruger Eq. 11.4
        double g1g0 = m_g1/m_g0;
        return  energy(T) * m_theta_e / ((T*T) * m_denom);
    }

    // copy of Eilmer3's Two_level_electronic::s_eval_entropy_from_T
    @nogc override number entropy(number T)
    {
        number tmp = m_g1 * exp(-m_theta_e / T);
        number Q = m_g0 + tmp;
        return mR * (log(Q) + m_theta_e / T / (m_g0 / tmp + 1.0));
    }

private:
    // store the denominator, so that Cv doesn't need to re-compute it
    // after calling energy
    // We don't have to worry about this having an old value, because
    // it will be re-computed for each calculation
    number m_denom;
    number m_theta_e;

    // degeneracies
    int m_g0, m_g1;
    double mR;
}

class MultiLevelElectronic : InternalEnergy
{
public:
    this (double[] theta_e, double R, int[] g)
    {
        mR = R;
        m_theta_e = theta_e;
        m_g = g;
        m_n_levels = min(g.length, theta_e.length);
    }

    @nogc number energy(number T) {
        // Copy of Eilmer3's Multi_level_electronic::s_eval_energy_from_T
        // compute energy with no checking of temperature
        number numerator = 0.0;
        number denominator = 0.0;
        number tmp;
        foreach (i; 0 .. m_n_levels){
            tmp = m_g[i] * exp(-m_theta_e[i] / T);
            numerator += m_theta_e[i] * tmp;
            denominator += tmp;
        }
        return mR * numerator / denominator;
    }

    @nogc number Cv(number T) {
        number tmp;
        number u = 0.0;
        number v = 0.0;
        number u_dash_star = 0.0;
        number v_dash_star = 0.0;
        foreach( i; 0 .. m_n_levels ) {
            tmp = m_g[i] * exp(-m_theta_e[i] / T);
            u += m_theta_e[i] * tmp;
            v += tmp;
            u_dash_star += tmp * m_theta_e[i] * m_theta_e[i];
            v_dash_star += m_theta_e[i] * tmp;
        }
        return mR / (T * T) * (u_dash_star * v - u * v_dash_star) / (v * v);
    }


    // Ref: Bottin, B. "Aerothermodynamic model of an inductively-coupled plasma wind tunnel"
    // copy of Eilmer3's Multi_level_electronic::s_eval_entropy_from_T
    @nogc override number entropy(number T){
        number tmpA = 0.0;
        number tmpB = 0.0;
        number tmpC;
        foreach ( i; 0 .. m_n_levels ) {
            tmpC = m_g[i] * exp(-m_theta_e[i] / T);
            tmpA += m_theta_e[i] * tmpC;
            tmpB += tmpC;
        }
        return mR * log(tmpB) + mR/T*tmpA/tmpB;
    }

private:
    double[] m_theta_e;
    double mR;
    int[] m_g;
    ulong m_n_levels;


}

class ZeroEnergy : InternalEnergy
{
public:
    this () {}
    @nogc override number energy(number T) { return to!number(0.0); }
    @nogc override number Cv(number T) { return to!number(0.0); }
    @nogc override number entropy(number T) { return to!number(0.0); }
}

InternalEnergy create_vibrational_energy_model(lua_State * L, double R)
{
    lua_getfield(L, -1, "vib_data");
    string model = getString(L, -1, "model");
    double theta_v = getDouble(L, -1, "theta_v");
    double theta_D;
    if (model == "truncated-harmonic"){
        theta_D = getDouble(L, -1, "theta_D");
    }
    lua_pop(L, 1); // vib data
    switch (model)
    {
        case "harmonic":
            return new HarmonicOscillator(theta_v, R);
        case "truncated-harmonic":
            return new TruncatedHarmonicOscillator(theta_v, theta_D, R);
        case "from-cea-thermo-curve":
            return new CEAThermoCurveVib(L, R);
        case "frozen":
            return new ZeroEnergy();
        default:
            throw new Error("Vibrational internal energy model not recognised");
    }
}

InternalEnergy create_electronic_energy_model(lua_State * L, double R)
{
    lua_getfield(L, -1, "electronic_levels");
    string model = getString(L, -1, "model");
    double[] theta_e;
    getArrayOfDoubles(L, -1, "Te", theta_e);
    // convert cm^(-1) -> K
    // note the speed of light constant is in m/s, but we need it in cm/s
    // here, hence the multiply by 100
    theta_e[] = theta_e[] * speed_of_light * 100.0 * Plancks_constant / Boltzmann_constant;
    int[] g;
    getArrayOfInts(L, -1, "g", g);

    lua_pop(L, 1);
    switch (model){
        case "two-level":
            return new TwoLevelElectronic(theta_e[0], R, g[0], g[1]);
        case "multi-level":
            return new MultiLevelElectronic(theta_e, R, g);
        case "frozen":
            return new ZeroEnergy();
        default:
            throw new Error("Electronic internal energy model not recognised");
    }
}
