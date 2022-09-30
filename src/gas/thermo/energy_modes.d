module gas.thermo.energy_modes;

import std.algorithm;
import std.math;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas;
import gas.thermo.thermo_model;

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
class HarmonicOscillatorInternalEnergy : InternalEnergy {
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


class ZeroVibrationInternalEnergy : InternalEnergy
{
public:
    this () {}
    @nogc override number energy(number T) { return to!number(0.0); }
    @nogc override number Cv(number T) { return to!number(0.0); }
    @nogc override number entropy(number T) { return to!number(0.0); }
}


// Electronic Energy modes
class TwoLevelElectronicInternalEnergy : InternalEnergy
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

class MultiLevelElectronicInternalEnergy : InternalEnergy
{
public:
    this (double[] theta_e, double R, int[] g, double T_min, double T_low, double T_blend)
    {
        mR = R;
        m_theta_e = theta_e;
        m_g = g;
        m_n_levels = min(g.length, theta_e.length);
        m_min_temp = to!number(T_min);
        m_blend_temp = to!number(T_blend);
        m_low_temp = to!number(T_low);
        m_linear_energy_min = _energy(m_min_temp);
        m_linear_energy_max = _energy(m_low_temp);
        m_Cv_low = (m_linear_energy_max - m_linear_energy_min)/(m_low_temp - m_min_temp);
        m_T_blend_high = m_low_temp + 0.5*m_blend_temp;
        m_T_blend_low = m_low_temp - 0.5*m_blend_temp;
    }

    @nogc override number energy(number T) {
        if (T > m_T_blend_high) {
            return _energy(T);
        }
        number e_low = m_linear_energy_min + m_Cv_low*(T-m_min_temp);
        if (T < m_T_blend_low) {
            return e_low;
        }
        number e_high = _energy(T);
        number wB = (T - m_T_blend_low)/m_blend_temp;
        number wA = 1.0 - wB;
        return wA*e_low + wB*e_high;
    }

    @nogc override number Cv(number T){
        if ( T > m_T_blend_high ) {
            return _Cv(T);
        }
        if ( T < m_T_blend_low ) {
            return m_Cv_low;
        }
        number Cv_high = _Cv(T);
        number Cv_low = m_Cv_low;

        number wB = (T - m_T_blend_low)/m_blend_temp;
        number wA = 1.0 - wB;

        return wA*Cv_low + wB*Cv_high;
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
    number m_linear_energy_min, m_linear_energy_max, m_Cv_low;

    number m_blend_temp = 50;
    number m_low_temp = 400;
    number m_min_temp = 10;

    number m_T_blend_low, m_T_blend_high;

    @nogc number _energy(number T) {
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

    @nogc number _Cv(number T) {
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
}

class ZeroElectronicInternalEnergy : InternalEnergy
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
    lua_pop(L, 1);
    switch (model)
    {
        case "harmonic":
            return new HarmonicOscillatorInternalEnergy(theta_v, R);
        case "frozen":
            return new ZeroVibrationInternalEnergy();
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
    double T_low = getDoubleWithDefault(L, -1, "T_low", 0.0);
    double T_blend = getDoubleWithDefault(L, -1, "T_blend", 0.0);
    double T_min = getDoubleWithDefault(L, -1, "T_min", 0.1);

    lua_pop(L, 1);
    switch (model){
        case "two-level":
            return new TwoLevelElectronicInternalEnergy(theta_e[0], R, g[0], g[1]);
        case "multi-level":
            return new MultiLevelElectronicInternalEnergy(theta_e, R, g, T_min, T_low, T_blend);
        case "frozen":
            return new ZeroElectronicInternalEnergy();
        default:
            throw new Error("Electronic internal energy model not recognised");
    }
}
