module gas.thermo.three_temperature_gas;

import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.algorithm;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;
import nm.bracketing;

import gas;
import gas.thermo.thermo_model;
import gas.thermo.energy_modes;
import gas.thermo.cea_thermo_curves;

immutable double T_REF = 298.15; // K (temperature the heat of formation is given at)
immutable ulong VIB = 0; // The vibrational energy goes in index 0
immutable ulong EE = 1; // the electron/electronic energy goes in index 1

class ThreeTemperatureGasMixture : ThermodynamicModel {
public:
    this(lua_State *L, string[] speciesNames)
    {
        mNSpecies = to!int(speciesNames.length);
        mR.length = mNSpecies;
        mCpTR.length = mNSpecies;
        m_vibration_internal_energy.length = mNSpecies;
        m_electronic_internal_energy.length = mNSpecies;
        m_reference_electronic_energy.length = mNSpecies;
        m_reference_vibration_energy.length = mNSpecies;
        mHf.length = mNSpecies;
        mM.length = mNSpecies;
        ms.length = mNSpecies;
        m_tr_entropy_constant.length = mNSpecies;
        foreach (isp, spName; speciesNames) {
            if (spName == "e-") mElectronIdx = to!int(isp);
            lua_getglobal(L, "db");
            lua_getfield(L, -1, spName.toStringz);
            mM[isp] = getDouble(L, -1, "M");
            mR[isp] = R_universal/mM[isp];
            //mHf[isp] = getDouble(L, -1, "Hf") / mM[isp]; // J/mol -> J/kg
            lua_getfield(L, -1, "thermoCoeffs");
            mCurves ~= new CEAThermoCurve(L, mR[isp]);
            lua_pop(L, 1);
            mHf[isp] = mCurves[isp].eval_h(to!number(T_REF));
            number tmp_m = R_universal / mR[isp] / Avogadro_number;
            number tmp = pow(2.0 * 3.14159265 * tmp_m * Boltzmann_constant / Plancks_constant / Plancks_constant, 1.5) * Boltzmann_constant;
            m_tr_entropy_constant[isp] = mR[isp] * (log(tmp) + 2.5);
            string type = getString(L, -1, "type");
            switch (type) {
            case "electron":
                m_vibration_internal_energy[isp] = new ZeroVibration();
                m_electronic_internal_energy[isp] = new ZeroElectronic();
                mCpTR[isp] = 0.0;
                break;
            case "atom" :
                mCpTR[isp] = (5./2.)*mR[isp];
                m_vibration_internal_energy[isp] = new ZeroVibration();
                break;
            case "molecule":
                string molType = getString(L, -1, "molecule_type");
                mCpTR[isp] = (molType == "linear") ? (7./2.)*mR[isp] : (8./2.)*mR[isp];
                m_vibration_internal_energy[isp] = create_vibrational_energy_model(L, mR[isp]);
                break;
            default:
                string msg = "ThreeTemperatureGas: error trying to match particle type.\n";
                throw new GasModelException(msg);
            }
            if (type != "electron"){
                m_electronic_internal_energy[isp] = create_electronic_energy_model(L, mR[isp]);
            }
            m_reference_electronic_energy[isp] = m_electronic_internal_energy[isp].energy(to!number(T_REF));
            m_reference_vibration_energy[isp] = m_vibration_internal_energy[isp].energy(to!number(T_REF));
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
    }

    override void updateFromPU(ref GasState gs) {}
    
    @nogc
    override void updateFromPT(ref GasState gs)
    {
        // update the density, electron pressure, and energies. Assume that the pressure
        // and all temperatures are already set correctly

        // first, update the density
        number denom = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            number T = (isp == mElectronIdx) ? gs.T_modes[EE] : gs.T;
            denom += gs.massf[isp] * mR[isp] * T;
        }
        gs.rho = gs.p/denom;

        // now update energies
        gs.u = trans_rot_energy_mixture(gs);
        gs.u_modes[VIB] = vib_energy_mixture(gs);
        gs.u_modes[EE] = electron_electronic_energy_mixture(gs);

        // update electron pressure, if needed
        if (mElectronIdx != -1)
            gs.p_e = gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[EE];
    }

    @nogc
    override void updateFromRhoU(ref GasState gs)
    {
        // Given the assumption that the translation/rotation mode is fully excited
        // the temperature can be found by direct inversion
        number Cv = trans_rot_Cv_mixture(gs);
        number formation_energy = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            formation_energy += gs.massf[isp] * (mHf[isp] - mCpTR[isp] * T_REF);
        }
        gs.T = (gs.u - formation_energy) / Cv;

        // iterate to find vibrational temperature
        gs.T_modes[VIB] = vibration_temperature(gs);

        // iterate to find electron/electronic temperature
        gs.T_modes[EE] = electron_electronic_temperature(gs);

        // finally, we can calculate p
        update_pressure(gs);
    }

    @nogc
    override void updateFromRhoT(ref GasState gs)
    {
        update_pressure(gs);
        gs.u = trans_rot_energy_mixture(gs);
        gs.u_modes[VIB] = vib_energy_mixture(gs);
        gs.u_modes[EE] = electron_electronic_energy_mixture(gs);
    }

    @nogc
    override void updateFromRhoP(ref GasState gs)
    {
        // First, calculate the temperature
        number p_heavy = gs.p;
        if (mElectronIdx != -1)
            p_heavy -= gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[EE];
        number denom = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            denom += gs.massf[isp] * gs.rho * mR[isp];
        }
        gs.T = p_heavy / denom;
        gs.u = trans_rot_energy_mixture(gs);

        // We assume that u_modes is set correctly
        gs.T_modes[VIB] = vibration_temperature(gs);
        gs.T_modes[EE] = electron_electronic_temperature(gs);
        if (mElectronIdx != -1)
            gs.p_e = gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[EE];
    }

    @nogc override void updateFromPS(ref GasState gs, number s)
    {
        throw new GasModelException("ThreeTemperatureGas: updateFromPS not implemented");
    }

    @nogc override void updateFromHS(ref GasState gs, number h, number s)
    {
        throw new GasModelException("ThreeTemperatureGas: updatefromHS not implemented");
    }

    @nogc override void updateSoundSpeed(ref GasState gs)
    {
        // compute the frozen sound speed
        number gamma = dhdTConstP(gs) / dudTConstV(gs);
        gs.a = sqrt(gamma*gs.p/gs.rho);
    }

    // Methods related to computing thermo derivatives.
    @nogc override number dudTConstV(in GasState gs)
    {
        number Cv = 0.0;
        number Cv_tr, Cv_vib, Cv_ee;
        foreach (isp; 0 .. mNSpecies){
            Cv_tr = trans_rot_Cv_species(isp);
            Cv_vib = vib_Cv_species(gs.T_modes[VIB], isp);
            Cv_ee = electron_electronic_Cv_species(gs.T_modes[EE], isp);
            Cv += fmax(0.0, gs.massf[isp]) * (Cv_tr + Cv_vib + Cv_ee);
        }
        return Cv;
    }
    @nogc override number dhdTConstP(in GasState gs)
    {
        number Cp = 0.0;
        number Cp_tr, Cv_vib, Cv_ee;
        foreach (isp; 0 .. mNSpecies){
            Cp_tr = mCpTR[isp];
            Cv_vib = vib_Cv_species(gs.T_modes[VIB], isp);
            Cv_ee = electron_electronic_Cv_species(gs.T_modes[EE], isp);
            Cp += fmax(0.0, gs.massf[isp]) * (Cp_tr + Cv_vib + Cv_ee);
        }
        return Cp;
    }
    @nogc override number dpdrhoConstT(in GasState gs)
    {
        throw new GasModelException("dpdrhoConstT not implemented for three temperature gas");
    }

    // Methods related to single properties of the mixture.
    @nogc override number gasConstant(in GasState gs)
    {
        return mass_average(gs, mR);
    }

    @nogc override number gasConstant(in GasState gs, int isp)
    {
        return to!number(mR[isp]);
    }

    @nogc override number internalEnergy(in GasState gs)
    {
        number e_tr = trans_rot_energy_mixture(gs);
        number e_v = vib_energy_mixture(gs);
        number e_ee = electron_electronic_energy_mixture(gs);
        return e_tr + e_v + e_ee;
    }

    @nogc
    override number energyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        switch (imode) {
            case -1:
                return trans_rot_energy_species(gs.T, isp);
            case 0:
                return vib_energy_species(gs.T_modes[VIB], isp);
            case 1:
                return electron_electronic_energy_species(gs.T_modes[EE], isp);
            default:
                throw new GasModelException("Invalid energy mode for three temperature model");
        }
    }

    @nogc override number enthalpy(in GasState gs)
    {
       return internalEnergy(gs) + gs.p / gs.rho;
    }

    @nogc
    override void enthalpies(in GasState gs, number[] hs)
    {
        foreach(isp; 0 .. mNSpecies){
            hs[isp] = enthalpyPerSpecies(gs, isp);
        }
    }

    @nogc override number enthalpyPerSpecies(in GasState gs, int isp)
    {
        number h_tr = mCpTR[isp]*(gs.T - T_REF) + mHf[isp];
        number h_v = vib_energy_species(gs.T_modes[VIB], isp);
        number h_e = electron_electronic_energy_species(gs.T_modes[EE], isp);
        if (isp == mElectronIdx) {
            h_e += (gs.T_modes[EE] - T_REF) * mR[mElectronIdx];
        }
        return h_tr + h_v + h_e;
    }

    @nogc override number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        switch (imode) {
            case -1:
                return mCpTR[isp] * (gs.T - T_REF) + mHf[isp];
            case 0:
                return vib_energy_species(gs.T_modes[VIB], isp);
            case 1:
                number h_e = electron_electronic_energy_species(gs.T_modes[EE], isp);
                if (isp == mElectronIdx) {
                    h_e += (gs.T_modes[EE] - T_REF) * mR[mElectronIdx] + mHf[mElectronIdx];
                }
                return h_e;
            default:
                throw new GasModelException("Invalid energy mode for three temperature model");
        }
    }

    @nogc override number entropy(in GasState gs)
    {
        foreach (isp; 0 .. mNSpecies ) {
            ms[isp] = entropyPerSpecies(gs, isp) - mR[isp]*log(gs.p/P_atm);
        }
        return mass_average(gs, ms);
    }

    @nogc override number entropyPerSpecies(in GasState gs, int isp)
    {
        return mCurves[isp].eval_s(gs.T);
    }

    @nogc override number cpPerSpecies(in GasState gs, int isp)
    {
        return mCpTR[isp] + vib_Cv_species(gs.T, isp) + electron_electronic_Cv_species(gs.T, isp);
    }
    @nogc override void GibbsFreeEnergies(in GasState gs, number[] gibbs_energies)
    {
        // Note that this uses the transrotational temperature only. We typically only
        // call this routine when computing the backward reaction rates, in which case
        // we WANT to only use T, but still, there's some danger here. (NNG)
        number T = gs.T;
        number logT = log(T);
        number logp = log(gs.p/P_atm);

        foreach(isp; 0 .. mNSpecies){
            number h = enthalpyPerSpecies(gs, isp);
            number s = mCurves[isp].eval_s(T, logT) - mR[isp]*logp;
            gibbs_energies[isp] = h - T*s;
        }
    }

private:
    int mNSpecies;
    int mElectronIdx = -1;
    double[] mM; // molar masses
    double[] mCpTR;
    double[] mR;
    number[] mHf;
    InternalEnergy[] m_vibration_internal_energy;
    InternalEnergy[] m_electronic_internal_energy;
    number[] m_reference_vibration_energy;
    number[] m_reference_electronic_energy;
    number[] ms;
    number[] m_tr_entropy_constant;
    CEAThermoCurve[] mCurves;

    @nogc
    number energy_in_mode(in GasState gs, int mode)
    {
        switch (mode) {
            case -1:
                return trans_rot_energy_mixture(gs);
            case 0:
                return vib_energy_mixture(gs);
            case 1:
                return electron_electronic_energy_mixture(gs);
            default:
                throw new GasModelException("Invalid energy mode for three temperature model");
        }
    }

    @nogc
    number trans_rot_Cv_species(int isp)
    {
        // the contribution from free electrons is accounted for in the
        // electron/electronic Cv
        if (isp == mElectronIdx) return to!number(0.0);
        return to!number(mCpTR[isp] - mR[isp]);
    }

    @nogc
    number trans_rot_Cv_mixture(in GasState gs)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            Cv += fmax(0.0, gs.massf[isp]) * trans_rot_Cv_species(isp);
        }
        return Cv;
    }

    @nogc
    number trans_rot_energy_species(number T, int isp)
    {
        // free electron translational energy is accounted for in the
        // electron/electronic energy
        if (isp == mElectronIdx) return to!number(0.0);

        // Include the formation enthalpy in the translation/rotation
        // mode for consistency with the two temperature model
        return mCpTR[isp] * (T-T_REF) + mHf[isp] - mR[isp] * T;
    }

    @nogc
    number trans_rot_energy_mixture(in GasState gs)
    {
        number e_tr = 0.0;
        foreach (isp; 0 .. mNSpecies)
        {
            if (isp == mElectronIdx) continue;
            e_tr += fmax(0.0, gs.massf[isp]) * trans_rot_energy_species(gs.T, isp);
        }
        return e_tr;
    }

    @nogc
    number vib_Cv_species(number Tv, int isp)
    {
        if (isp == mElectronIdx) return to!number(0.0);

        // For internal modes, Cv = Cp, so we can evaluate Cp
        number cp_at_Tv = mCurves[isp].eval_Cp(Tv);
        number cp = cp_at_Tv - mCpTR[isp] - electron_electronic_Cv_species(Tv, isp);
        return cp;
    }

    @nogc
    number vib_Cv_mixture(in GasState gs, number Tv)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            Cv += fmax(0.0, gs.massf[isp]) * vib_Cv_species(Tv, isp);
        }
        return Cv;
    }

    @nogc
    number vib_Cv_mixture(in GasState gs)
    {
        return vib_Cv_mixture(gs, gs.T_modes[VIB]);
    }

    @nogc
    number vib_energy_species(number Tv, int isp)
    {
        //return m_vibration_internal_energy[isp].energy(Tv);
        if (isp == mElectronIdx) return to!number(0.0);
        number h_at_Tv = mCurves[isp].eval_h(Tv);
        number h_v = h_at_Tv - mCpTR[isp]*(Tv - T_REF) - mHf[isp] - electron_electronic_energy_species(Tv, isp);
        return h_v - m_reference_vibration_energy[isp];
    }

    @nogc
    number vib_energy_mixture(in GasState gs, number Tv)
    {
        number e_v = 0.0;
        foreach (isp; 0 .. mNSpecies)
        {
            e_v += fmax(0.0, gs.massf[isp]) * vib_energy_species(Tv, isp);
        }
        return e_v;
    }

    @nogc
    number vib_energy_mixture(in GasState gs)
    {
        return vib_energy_mixture(gs, gs.T_modes[VIB]);
    }

    @nogc
    number electron_electronic_energy_species(number Tee, int isp)
    {
        // free electrons have no electronic component, but they do contribute
        // to the energy via their translation degrees of freedom
        if (isp == mElectronIdx) return (3.0/2.0) * (Tee-T_REF) * mR[mElectronIdx];
        return m_electronic_internal_energy[isp].energy(Tee) - m_reference_electronic_energy[isp];
    }

    @nogc
    number electron_electronic_energy_mixture(in GasState gs, number Tee)
    {
        number e_ee = 0.0;
        foreach (isp; 0 .. mNSpecies){
            e_ee += fmax(0.0, gs.massf[isp]) * electron_electronic_energy_species(Tee, isp);
        }
        return e_ee;
    }

    @nogc
    number electron_electronic_energy_mixture(in GasState gs)
    {
        return electron_electronic_energy_mixture(gs, gs.T_modes[EE]);
    }

    @nogc number electron_electronic_Cv_species(number Tee, int isp)
    {
        // free electrons have no electronic component, but they do contribute
        // to Cv via their translation degrees of freedom
        if (isp == mElectronIdx) return to!number(3.0/2.0) * mR[isp];
        return m_electronic_internal_energy[isp].Cv(Tee);
    }

    @nogc number electron_electronic_Cv_mixture(in GasState gs, number Tee)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            Cv += fmax(0.0, gs.massf[isp]) * electron_electronic_Cv_species(Tee, isp);
        }
        return Cv;
    }

    @nogc number electron_electronic_Cv_mixture(in GasState gs)
    {
        return electron_electronic_Cv_mixture(gs, gs.T_modes[EE]);
    }

    @nogc
    void update_pressure(ref GasState gs)
    {
        gs.p = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            number T = (isp == mElectronIdx) ? gs.T_modes[EE] : gs.T;
            gs.p += gs.rho * fmax(0.0, gs.massf[isp]) * mR[isp] * T;
        }

        // set electron pressure whilst we're here
        if (mElectronIdx != 1) {
            gs.p_e = gs.rho * fmax(0.0, gs.massf[mElectronIdx]) * mR[mElectronIdx] * gs.T_modes[EE];
        }
    }

    @nogc
    number vibration_temperature(ref GasState gs)
    {
        const int MAX_ITERATIONS = 40;
        const double TOL = 1e-8;
        const double IMAGINARY_TOL = 1e-15;

        // Take the supplied T_modes[0] as the initial guess.
        number T_guess = gs.T_modes[VIB];
        number u0 = vib_energy_mixture(gs, T_guess);
        number f_guess =  u0 - gs.u_modes[VIB];

        // Begin iterating.
        int count = 0;
        number Cv, dT;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            Cv = vib_Cv_mixture(gs, T_guess);
            dT = -f_guess/Cv;
            T_guess += dT;
            version(complex_numbers) {
                if (fabs(dT) < TOL && fabs(dT.im) < IMAGINARY_TOL) {
                    break;
                }
            }
            else {
                if (fabs(dT) < TOL) {
                    break;
                }
            }
            f_guess = vib_energy_mixture(gs, T_guess) - gs.u_modes[VIB];
            count++;
        }

        if ((count == MAX_ITERATIONS)&&(fabs(dT)>1e-3)) {
            string msg = "The 'vibTemperature' function failed to converge.\n";
            debug {
                msg ~= format("The final value for Tvib was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= gs.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }
        return T_guess;
    }


    @nogc
    number electron_electronic_temperature(ref GasState gs)
    {
        const int MAX_ITERATIONS = 200;
        const double TOL = 1.0e-8;
        const double IMAGINARY_TOL = 1.0e-15;
        immutable number T_MIN = to!number(10.0);
        immutable number T_MAX = to!number(200_000.0);

        // we're aiming to find a temperature corresponding to this energy
        number target_u = gs.u_modes[EE];

        // guess that the new temperature is within +/- 5 K of gs.T_modes[imode]
        // this guess will be adjusted later by a bracketing algorithm
        number Ta = gs.T_modes[EE] - to!number(5.0);
        number Tb = gs.T_modes[EE] + to!number(5.0);

        number zero_func(number T) {
            number u = electron_electronic_energy_mixture(gs, T);
            return u - target_u;
        }

        // try to bracket the root
        if (bracket!(zero_func, number)(Ta, Tb, T_MIN, T_MAX) == -1){
            // If bracketing fails, there probably isn't any energy in this mode...
            // If our energy is already close, we'll just accept the temperature
            // that was already there. It is probably a pre-shock free stream
            // value, so has been explicitly set by the user.
            if (fabs(zero_func(Tb)) < 1e-1) {
                return gs.T_modes[EE];
            }

            // otherwise there's not much we can do...
            throw new GasModelException("Undefined electron/electronic temperature");
        }

        // the function evaluation at the bounds
        number fa, fb;
        fa = zero_func(Ta);
        fb = zero_func(Tb);
        number dT = fabs(Tb - Ta);
        number dT_old = dT;


        // use gs.T_modes[EE] as the initial guess
        number T = gs.T_modes[EE];

        // the function evaluation and derivative at the current guess
        number f = zero_func(T);
        number df = electron_electronic_Cv_mixture(gs, T);

        bool converged = false;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            bool newton_unstable = ((((T-Tb)*df-f) * ((T-Ta)*df-f)) > 0.0);
            bool newton_slow = (fabs(2.0*f) > fabs(dT_old*df));
            if (newton_unstable || newton_slow) {
                // use bisection
                dT_old = dT;
                dT = 0.5*(Tb-Ta);
                T = Ta + dT;
            }
            else {
                // use Newton
                dT_old = dT;
                dT = f/df;
                T -= dT;
            }

            // check for convergence
            version(complex_numbers) {
                if ((fabs(dT) < TOL) && fabs(dT.im) < IMAGINARY_TOL) {
                    converged = true;
                    break;
                }
            }
            else {
                if (fabs(dT) < TOL) {
                    converged = true;
                    break;
                }
            }

            // evaluate function for next iteration
            f = zero_func(T);
            df = electron_electronic_Cv_mixture(gs, T);

            // maintain the brackets on the root
            if (f < 0.0) {
                Ta = T;
            }
            else {
                Tb = T;
            }
        }

        if (!converged) {
            string msg = "ThreeTemperatureGas: Electron/electronic temperature failed to converge.\n";
            debug {
                msg ~= format("The final value was: %.16f\n", T);
                msg ~= "The supplied GasState was:\n";
                msg ~= gs.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }
        return T;
    }
}

version(three_temperature_gas_test)
{
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        doLuaFile(L, "sample-data/N2_3T.lua");
        string[] speciesNames;
        getArrayOfStrings(L, "species", speciesNames);
        auto tm = new ThreeTemperatureGasMixture(L, speciesNames);
        lua_close(L);
        auto gs = new GasState(1, 2);

        gs.T = 300.0; gs.T_modes[] = 300.0;
        gs.massf[0] = 1.0;
        gs.p = 1e5;
        tm.updateFromPT(gs);

        // known low temperature values to make sure things aren't severely broken
        assert(approxEqualNumbers(to!number(1.934130659e3), tm.enthalpy(gs), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(-8.7107292374465e4), tm.internalEnergy(gs), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(-8.711961231946e4), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(7.4247638645e2), tm.dudTConstV(gs), 1.0e-6), failedUnitTest());
        return 0;
    }
}
