/**
 * thermo/two_temperature_gas.d
 *
 * Author: Rowan G.
 * History: 2021-03-15 -- first refactor from other pieces
 **/

module gas.thermo.two_temperature_gas;

import std.stdio;
import std.math;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;
import nm.brent;
import nm.bracketing : bracket;

import gas;
import gas.thermo.thermo_model;
import gas.thermo.cea_thermo_curves;

immutable double T_REF = 298.15; // K (based on CEA reference temperature value)

/**
 * The TwoTemperatureGasMixture class provides a model for thermodynamic behaviour
 * of a two-temperature gas mixture.
 *
 * The assumption in this model is that the translational and rotational energy modes
 * are described by one temperature (T_tr) and the vibrational and electronic energy
 * modes are described my a second temperature (T_ve). The first of these appears in
 * in the GasState as T, and the second as T_modes[0].
 *
 * This model is useful for air mixtures, nitrogen mixtures and carbon dioxide flows.
 * Do not confuse this two temperature model with specialised two temperature models
 * for argon and hydrogen/helium mixtures. In those, the partitioning of energy is
 * focussed on separating the electron as having a distinct temperature from the
 * heavy particles.
 */
class TwoTemperatureGasMixture : ThermodynamicModel {
public:

    this(lua_State *L, string[] speciesNames)
    {
        mNSpecies = to!int(speciesNames.length);
        mR.length = mNSpecies;
        mDel_hf.length = mNSpecies;
        mCpTR.length = mNSpecies;
        ms.length = mNSpecies;
        foreach (isp, spName; speciesNames) {
            if (spName == "e-") mElectronIdx = to!int(isp);
            lua_getglobal(L, "db");
            lua_getfield(L, -1, spName.toStringz);
            double m = getDouble(L, -1, "M");
            mR[isp] = R_universal/m;
            lua_getfield(L, -1, "thermoCoeffs");
            mCurves ~= new CEAThermoCurve(L, mR[isp]);
            lua_pop(L, 1);
            mDel_hf[isp] = mCurves[isp].eval_h(to!number(T_REF));
            string type = getString(L, -1, "type");
            switch (type) {
            case "electron":
                mCpTR[isp] = 0.0;
                break;
            case "atom" :
                mCpTR[isp] = (5./2.)*mR[isp];
                break;
            case "molecule":
                string molType = getString(L, -1, "molecule_type");
                mCpTR[isp] = (molType == "linear") ? (7./2.)*mR[isp] : (8./2.)*mR[isp];
                break;
            default:
                string msg = "TwoTemperatureGas: error trying to match particle type.\n";
                throw new Error(msg);
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
    }

    override void updateFromPU(ref GasState gs) {}
    
    @nogc
    override void updateFromPT(ref GasState gs)
    {
        updateDensity(gs);
        gs.u = transRotEnergyMixture(gs);
        gs.u_modes[0] = vibElecEnergyMixture(gs, gs.T_modes[0]);
        if (mElectronIdx != -1) gs.p_e = gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[0];
    }

    @nogc
    override void updateFromRhoU(ref GasState gs)
    {
        // We can compute T by direct inversion since the Cp in
        // in translation and rotation are fully excited,
        // and, as such, constant.
        number sumA = 0.0;
        number sumB = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            if (isp == mElectronIdx) continue;
            sumA += gs.massf[isp]*(mCpTR[isp]*T_REF - mDel_hf[isp]);
            sumB += gs.massf[isp]*(mCpTR[isp] - mR[isp]);
        }
        gs.T = (gs.u + sumA)/sumB;
        // Next, we can compute T_modes[0] by iteration.
        // We'll use a Newton method since the function
        // should vary smoothly at the polynomial breaks
        // given that we blend the coefficients in those regions.
        gs.T_modes[0] = vibElecTemperature(gs);
        // Now we can compute pressure from the perfect gas
        // equation of state.
        updatePressure(gs);
    }

    @nogc
    override void updateFromRhoT(ref GasState gs)
    {
        updatePressure(gs);
        gs.u = transRotEnergyMixture(gs);
        gs.u_modes[0] = vibElecEnergyMixture(gs, gs.T_modes[0]);
    }

    @nogc
    override void updateFromRhoP(ref GasState gs)
    {
        // In this method, we assume that u_modes[0] is set correctly
        // in addition to density and pressure.
        gs.T_modes[0] = vibElecTemperature(gs);
        updateTemperatureFromRhoP(gs);
        gs.u = transRotEnergyMixture(gs);
        if (mElectronIdx != -1) gs.p_e = gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[0];
    }

    @nogc
    override void updateFromPS(ref GasState gs, number s)
    {
        double TOL = 1.0e-6;
        number delT = 100.0;
        number T1 = gs.T;
        // It's possible that T is 'nan' if it's never been set, so:
        if (isNaN(T1))
            T1 = 300.0; // Just set at some value to get started.
        number Tsave = T1;
        number T2 = T1 + delT;

        number zeroFun(number T) {
            gs.T = T;
            gs.T_modes[0] = T;
            number s_guess = entropy(gs);
            return s - s_guess;
        }

        number Tmin = 10.0;
        number Tmax = 200000.0;
        if (bracket!(zeroFun,number)(T1, T2, Tmin, Tmax) == -1) {
            string msg = "The 'bracket' function failed to find temperature values\n";
            debug {
                msg ~= "that bracketed the zero function in GasMixtureThermo.update_thermo_from_ps().\n";
                msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
            }
            throw new Exception(msg);
        }

        if (T1 < T_MIN)
            T1 = T_MIN;

        try {
            gs.T = nm.brent.solve!(zeroFun,number)(T1, T2, TOL);
            gs.T_modes[0] = gs.T;
        }
        catch (Exception e) {
            string msg = "There was a problem iterating to find temperature\n";
            debug {
                msg ~= "in function GasMixtureThermo.update_thermo_from_ps().\n";
                msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                msg ~= format("The target entropy value was: %12.6f\n", s);
                msg ~= format("The GasState is currently:\n");
                msg ~= gs.toString();
                msg ~= "The message from the brent.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }
        updateFromPT(gs);
    }

    @nogc
    override void updateFromHS(ref GasState gs, number h, number s)
    {
        /*
            TODO: This implementation enforces that T_modes[0]==T, since that
            is what we expect at nozzle inflow boundaries, which is what
            this routine is typically used for. Beware however, that
            this may not be true in general.
            @author: NNG (23/11/22), copied from therm_perf_gas_mix.d
        */

        // We do this in two stages.
        // First, from enthalpy we compute temperature.
        double TOL = 1.0e-6;
        number delT = 100.0;
        number T1 = gs.T;
        // It's possible that T is 'nan' if it's never been set, so:
        if (isNaN(T1))
            T1 = 300.0; // Just set at some value to get started.
        number Tsave = T1;
        number T2 = T1 + delT;

        number zeroFun(number T)
        {
            gs.T = T;
            gs.T_modes[0] = T;
            number h_guess = enthalpy(gs);
            return h - h_guess;
        }

        if (bracket!(zeroFun,number)(T1, T2) == -1) {
            string msg = "The 'bracket' function failed to find temperature values";
            debug {
                msg ~= "\nthat bracketed the zero function in GasMixtureThermo.update_thermo_from_hs().\n";
                msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
            }
            throw new Exception(msg);
        }

        if (T1 < T_MIN)
            T1 = T_MIN;

        try {
            gs.T = nm.brent.solve!(zeroFun,number)(T1, T2, TOL);
            gs.T_modes[0] = gs.T;
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find temperature";
            debug {
                msg ~= "\nin function GasMixtureThermo.update_thermo_from_hs().\n";
                msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                msg ~= format("The target enthalpy value was: %12.6f\n", h);
                msg ~= format("The GasState is currently:\n");
                msg ~= gs.toString();
                msg ~= "The message from the brent.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }

        // Second, we can iterate to find the pressure that gives
        // correct entropy.
        TOL = 1.0e-3;
        number delp = 1000.0;
        number p1 = gs.p;
        number psave = p1;
        number p2 = p1 + delp;

        number zeroFun2(number p)
        {
            gs.p = p;
            number s_guess = entropy(gs);
            return s - s_guess;
        }

        if ( bracket!(zeroFun2,number)(p1, p2) == -1 ) {
            string msg = "The 'bracket' function failed to find pressure values";
            debug {
                msg ~= "\nthat bracketed the zero function in GasMixtureThermo.update_thermo_from_hs().\n";
                msg ~= format("The final values are: p1 = %12.6f and p2 = %12.6f\n", p1, p2);
            }
            throw new Exception(msg);
        }

        if ( p1 < 0.0 )
            p1 = 1.0;

        try {
            gs.p = nm.brent.solve!(zeroFun2,number)(p1, p2, TOL);
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find pressure";
            debug {
                msg ~= "\nin function GasMixtureThermo.update_thermo_from_hs().\n";
                msg ~= format("The initial pressure guess was: %12.6f\n", psave);
                msg ~= format("The target entropy value was: %12.6f\n", s);
                msg ~= format("The GasState is currently:\n");
                msg ~= gs.toString();
                msg ~= "The message from the brent.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }
        updateFromPT(gs);
    }

    @nogc
    override void updateSoundSpeed(ref GasState gs)
    {
        // We compute the frozen sound speed
        number gamma = dhdTConstP(gs)/dudTConstV(gs);
        gs.a = sqrt(gamma*gs.p/gs.rho);
    }

    @nogc
    override number dudTConstV(in GasState gs)
    {
        number Cv = 0.0;
        number Cv_tr_rot, Cv_vib;
        foreach (isp; 0 .. mNSpecies) {
            Cv_tr_rot = transRotCvPerSpecies(isp);
            Cv_vib = vibElecCvPerSpecies(gs.T_modes[0], isp);
            Cv += gs.massf[isp] * (Cv_tr_rot + Cv_vib);
        }
        return Cv;
    }

    @nogc
    override number dhdTConstP(in GasState gs)
    {
        // Using the fact that internal structure specific heats
        // are equal, that is, Cp_vib = Cv_vib
        number Cp = 0.0;
        number Cp_vib;
        foreach (isp; 0 .. mNSpecies) {
            Cp_vib = vibElecCvPerSpecies(gs.T_modes[0], isp);
            Cp += gs.massf[isp] * (mCpTR[isp] + Cp_vib);
        }
        return Cp;
    }

    @nogc
    override number dpdrhoConstT(in GasState gs)
    {
        number sum = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            number T = (isp == mElectronIdx) ? gs.T_modes[0] : gs.T;
            sum += gs.massf[isp] * mR[isp] * T;
        }
        return sum;
    }

    @nogc
    override number gasConstant(in GasState gs)
    {
        return mass_average(gs, mR);
    }

    @nogc
    override number gasConstant(in GasState gs, int isp)
    {
        return to!number(mR[isp]);
    }

    @nogc
    override number internalEnergy(in GasState gs)
    {
        number Tve = gs.T_modes[0];
        number logTve = log(Tve);
        number u_tr = transRotEnergyMixture(gs);
        number u_ve = vibElecEnergyMixture(gs, Tve, logTve);
        return u_tr + u_ve;
    }

    @nogc
    override number energyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        if (imode == 0) {
            return vibElecEnergyPerSpecies(gs.T_modes[0], isp);
        }
        return transRotEnergyPerSpecies(gs.T, isp);
    }

    @nogc
    override number enthalpy(in GasState gs)
    {
        number logTve = log(gs.T_modes[0]);
        number u = transRotEnergyMixture(gs) + vibElecEnergyMixture(gs, gs.T_modes[0], logTve);
        return u + gs.p/gs.rho;
    }

    @nogc
    override void enthalpies(in GasState gs, number[] hs)
    {
        number Tve = gs.T_modes[0];
        number logTve = log(Tve);
        foreach(isp; 0 .. mNSpecies){
            number h_tr = mCpTR[isp]*(gs.T - T_REF) + mDel_hf[isp];
            number h_ve = vibElecEnergyPerSpecies(Tve, logTve, isp);
            hs[isp] = h_tr + h_ve;
        }
    }

    @nogc
    override number enthalpyPerSpecies(in GasState gs, int isp)
    {
        number h_tr = mCpTR[isp]*(gs.T - T_REF) + mDel_hf[isp];
        number h_ve = vibElecEnergyPerSpecies(gs.T_modes[0], isp);
        return h_tr + h_ve;
    }

    @nogc
    override number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        if (imode == 0) {
            return vibElecEnergyPerSpecies(gs.T_modes[0], isp);
        }
        return mCpTR[isp]*(gs.T - T_REF) + mDel_hf[isp];
    }

    @nogc
    override number entropy(in GasState gs)
    {
        number logT = log(gs.T);
        foreach ( isp; 0 .. mNSpecies ) {
            ms[isp] = mCurves[isp].eval_s(gs.T, logT) - mR[isp]*log(gs.p/P_atm);
        }
        return mass_average(gs, ms);
    }

    @nogc
    override number entropyPerSpecies(in GasState gs, int isp)
    {
        return mCurves[isp].eval_s(gs.T);
    }

    @nogc
    number vibElecEnergyPerSpecies(number Tve, number logTve, int isp)
    {
        // The electron possess energy only in translation.
        // We put this contribution in the electronic energy since
        // its translation is governed by the vibroelectronic temperature.
        if (isp == mElectronIdx) return (3./2.)* mR[isp] * Tve;
        // For heavy particles
        number h_at_Tve = mCurves[isp].eval_h(Tve, logTve);
        number h_ve = h_at_Tve - mCpTR[isp]*(Tve - T_REF) - mDel_hf[isp];
        return h_ve;
    }

    @nogc
    number vibElecEnergyPerSpecies(number Tve, int isp)
    {
        number logTve = log(Tve);
        return vibElecEnergyPerSpecies(Tve, logTve, isp);
    }

    @nogc
    number vibElecEnergyMixture(in GasState gs, number Tve, number logTve)
    {
        number e_ve = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            e_ve += gs.massf[isp] * vibElecEnergyPerSpecies(Tve, logTve, isp);
        }
        return e_ve;
    }

    @nogc
    number vibElecEnergyMixture(in GasState gs, number Tve)
    {
        number logTve = log(Tve);
        return vibElecEnergyMixture(gs, Tve, logTve);
    }
    @nogc
    override number cpPerSpecies(in GasState gs, int isp)
    {
        // Using the fact that internal structure specific heats
        // are equal, that is, Cp_vib = Cv_vib
        return mCpTR[isp] + vibElecCvPerSpecies(gs.T_modes[0], isp);
    }
    @nogc override void GibbsFreeEnergies(in GasState gs, number[] gibbs_energies)
    {
        number T = gs.T;
        number logT = log(T);
        number Tve = gs.T_modes[0];
        number logTve = log(Tve);
        number logp = log(gs.p/P_atm);

        foreach(isp; 0 .. mNSpecies){
            number h_tr = mCpTR[isp]*(gs.T - T_REF) + mDel_hf[isp];
            number h_ve = vibElecEnergyPerSpecies(Tve, logTve, isp);
            number h = h_tr + h_ve;
            number s = mCurves[isp].eval_s(T, logT) - mR[isp]*logp;
            gibbs_energies[isp] = h - T*s;
        }
    }


private:
    double[] mR;
    number[] ms;
    double[] mCpTR;
    number[] mDel_hf;
    CEAThermoCurve[] mCurves;
    int mNSpecies;
    int mElectronIdx = -1; // Set to this in case never set in neutrals-only simulations

    @nogc
    void updateDensity(ref GasState gs)
    {
        number denom = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            number T = (isp == mElectronIdx) ? gs.T_modes[0] : gs.T;
            denom += gs.massf[isp] * mR[isp] * T;
        }
        gs.rho = gs.p/denom;
    }

    @nogc
    void updatePressure(ref GasState gs)
    {
        gs.p = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            number T = (isp == mElectronIdx) ? gs.T_modes[0] : gs.T;
            gs.p += gs.rho * gs.massf[isp] * mR[isp] * T;
        }
        // Also set electron pressure, while we're computing pressures.
        if (mElectronIdx != -1) {
            gs.p_e = gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[0];
        }
    }

    @nogc
    void updateTemperatureFromRhoP(ref GasState gs)
    {
        // This assumes the T_modes[0] is known, and we're only trying to determine T
        number pHeavy = gs.p;
        if (mElectronIdx != -1) {
            pHeavy -= gs.rho * gs.massf[mElectronIdx] * mR[mElectronIdx] * gs.T_modes[0];
        }
        number denom = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            if (isp == mElectronIdx) continue;
            denom += gs.rho * gs.massf[isp] * mR[isp];
        }
        gs.T = pHeavy/denom;
    }

    @nogc
    number transRotEnergyPerSpecies(number T, int isp)
    {
        if (isp == mElectronIdx) return to!number(0.0);
        number h = mCpTR[isp] * (T - T_REF) + mDel_hf[isp];
        number e = h - mR[isp]*T;
        return e;
    }

    @nogc
    number transRotEnergyMixture(in GasState gs)
    {
        number e = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            if (isp == mElectronIdx) continue;
            e += gs.massf[isp] * transRotEnergyPerSpecies(gs.T, isp);
        }
        return e;
    }

    @nogc
    number vibElecCvPerSpecies(number Tve, int isp)
    {
        // electron as special case
        if (isp == mElectronIdx) return to!number((3./2.)*mR[isp]);

        // all other heavy species
        // Why Cp and not Cv?
        // We are computing an "internal" Cv here, the Cv_ve.
        // For internal energy storage Cv_int = Cp_int,
        // so we can make the calculation using Cp relations.
        return mCurves[isp].eval_Cp(Tve) - mCpTR[isp];
    }

    @nogc
    number vibElecCvMixture(in GasState gs, number Tve)
    {
        // Why pass in Tve, when we could dip into gs.T_modes[0]?
        // There are time when we'd like to use a value for Tve other
        // than T_modes[0]. For example, if we want to compute this
        // value when in thermal equilibrium, we can pass in Tve = gs.T.
        number Cv_ve = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            Cv_ve += gs.massf[isp] * vibElecCvPerSpecies(Tve, isp);
        }
        return Cv_ve;
    }

    @nogc
    number transRotCvPerSpecies(int isp)
    {
        // special case for electrons
        if (isp == mElectronIdx) return to!number(0.0);
        // all other heavy species
        return to!number(mCpTR[isp] - mR[isp]);
    }

    @nogc
    number transRotCvMixture(in GasState gs)
    {
        number Cv_tr = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            Cv_tr += gs.massf[isp] * transRotCvPerSpecies(isp);
        }
        return Cv_tr;
    }

    version(complex_numbers) {
        // For the complex numbers version of the code we need a Newton's method with a fixed number of iterations.
        // An explanation can be found in:
        //     E. J. Nielsen and W. L. Kleb
        //     Efficient Construction of Discrete Adjoint Operators on Unstructured Grids by Using Complex Variables
        //     AIAA Journal, vol. 44, no. 4, 2006.
        //
        // Changed Oct 2024 by NNG to match single species T from u.
        //   1.) We can now exit early if the real component is converged, this saves lots of time.
        //   2.) We always do at least one iteration, to ensure the complex component is set.
        //   3.) We throw an exception if the real part of dT isn't small at the end. This helps the
        //       physicality check make sure that the gas state has been set to a sensible value.
        //   TODO: We don't really need a different routine for complex and real now.
        @nogc
        number vibElecTemperature(ref GasState gs)
        {
            immutable int MAX_ITERATIONS = 10;
            immutable double TOL = 1e-10;

            // Take the supplied T_modes[0] as the initial guess.
            number T_guess = gs.T_modes[0];
            number logT_guess = log(T_guess);
            number u0 = vibElecEnergyMixture(gs, T_guess, logT_guess);
            number f_guess =  u0 - gs.u_modes[0];

            // Begin iterating.
            int count = 0;
            number Cv, dT;
            foreach (iter; 0 .. MAX_ITERATIONS) {
                Cv = vibElecCvMixture(gs, T_guess);
                dT = -f_guess/Cv;
                T_guess += dT;
                if (fabs(dT.re)/T_guess.re < TOL) {
                    return T_guess;
                }
                logT_guess = log(T_guess);
                f_guess = vibElecEnergyMixture(gs, T_guess, logT_guess) - gs.u_modes[0];
                count++;
            }
            if (fabs(dT)>1e-3) {
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
    } else {
        @nogc
        number vibElecTemperature(ref GasState gs)
        {
            int MAX_ITERATIONS = 40;

            // Take the supplied T_modes[0] as the initial guess.
            number T_guess = gs.T_modes[0];
            number u0 = vibElecEnergyMixture(gs, T_guess);
            number f_guess =  u0 - gs.u_modes[0];

            // Before iterating, check if the supplied guess is
            // good enough. Define good enough as 1/100th of a Joule.
            double E_TOL = 1e-5;
            if (fabs(f_guess) < E_TOL) {
                return gs.T_modes[0];
            }

            // We'll keep adjusting our temperature estimate
            // until it is less than TOL.
            double TOL = 1.0e-9;

            // Begin iterating.
            int count = 0;
            number Cv, dT;
            foreach (iter; 0 .. MAX_ITERATIONS) {
                Cv = vibElecCvMixture(gs, T_guess);
                dT = -f_guess/Cv;
                T_guess += dT;
                if (fabs(dT) < TOL) {
                    break;
                }
                f_guess = vibElecEnergyMixture(gs, T_guess) - gs.u_modes[0];
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
    } // version()
}

version(two_temperature_gas_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        doLuaFile(L, "sample-data/five-species-air.lua");
        string[] speciesNames;
        getArrayOfStrings(L, "species", speciesNames);
        auto tm = new TwoTemperatureGasMixture(L, speciesNames);
        lua_close(L);
        auto gs = GasState(5, 1);

        gs.p = 1.0e6;
        gs.T = 2000.0;
        gs.T_modes[0] = gs.T;
        gs.massf = [to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2)];
        tm.updateFromPT(gs);
        assert(approxEqualNumbers(to!number(11801825.6), gs.u + gs.u_modes[0], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.2840117), gs.rho, 1.0e-6), failedUnitTest());

        tm.updateFromRhoU(gs);
        assert(approxEqualNumbers(to!number(2000.0), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(2000.0), gs.T_modes[0], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1e6), gs.p, 1.0e-6), failedUnitTest());

        return 0;
    }
}
