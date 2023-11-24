/*
 * thermo/therm_perg_gas_mix.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2021-02-09
 */

module gas.thermo.therm_perf_gas_mix;

import std.math;
import std.string;
import std.conv;
import std.stdio;
import util.lua;
import util.lua_service;
import nm.nm_exception : NumericalMethodException;
import ntypes.complex;
import nm.number;
import nm.brent;
import nm.newton;
import nm.bracketing : bracket;

import gas;
import gas.thermo.thermo_model;
import gas.thermo.cea_thermo_curves;

class ThermPerfGasMixture : ThermodynamicModel {
public:

    this(lua_State *L, string[] speciesNames)
    {
        n_species = speciesNames.length;
        mR.length = n_species;
        mVals.length = n_species;

        foreach (isp, spName; speciesNames) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, spName.toStringz);
            double m = getDouble(L, -1, "M");
            mR[isp] = R_universal/m;
            lua_getfield(L, -1, "thermoCoeffs");
            mCurves ~= new CEAThermoCurve(L, mR[isp]);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
    }

    override void updateFromPU(ref GasState gs)
    {
        update_temperature_from_energy(gs);
        update_density(gs);
    }

    
    @nogc
    override void updateFromPT(ref GasState gs)
    {
        update_density(gs);
        update_energy(gs);
    }
    @nogc
    override void updateFromRhoU(ref GasState gs)
    {
        gs.T = update_temperature_from_energy(gs);
        update_pressure(gs);
    }
    @nogc
    override void updateFromRhoT(ref GasState gs)
    {
        update_energy(gs);
        update_pressure(gs);
    }
    @nogc
    override void updateFromRhoP(ref GasState gs)
    {
        update_temperature_from_rhop(gs);
        update_energy(gs);
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
        update_energy(gs);
        update_density(gs);
    }

    @nogc
    override void updateFromHS(ref GasState gs, number h, number s)
    {
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
        update_energy(gs);
        update_density(gs);
    }

    @nogc
    override void updateSoundSpeed(ref GasState gs)
    {
        // Reference:
        // Cengel and Boles (1998)
        // Thermodynamics: an Engineering Approach, 3rd edition
        // McGraw Hill
        // Equation 16-10 on p. 849

        // "frozen" sound speed
        gs.a = sqrt(gamma(gs)*dpdrhoConstT(gs));
    }

    @nogc
    override number dudTConstV(in GasState gs)
    {
        // Using Cv = Cp - R
        foreach (isp, ref Cv; mVals) {
            Cv = mCurves[isp].eval_Cp(gs.T) - mR[isp];
        }
        number Cv = mass_average(gs, mVals);
        return Cv;
    }

    @nogc
    override number dhdTConstP(in GasState gs)
    {
        foreach (isp, ref Cp; mVals) {
            Cp = mCurves[isp].eval_Cp(gs.T);
        }
        number Cp = mass_average(gs, mVals);
        return Cp;
    }

    @nogc
    override number dpdrhoConstT(in GasState gs)
    {
        number Rmix = gasConstant(gs);
        return Rmix*gs.T;
    }

    @nogc
    override number gasConstant(in GasState gs)
    {
        number Rmix = mass_average(gs, mR);
        return Rmix;
    }

    @nogc
    override number gasConstant(in GasState gs, int isp)
    {
        return to!number(mR[isp]);
    }

    @nogc
    override number internalEnergy(in GasState gs)
    {
        number logT = log(gs.T);
        foreach (isp, ref u; mVals) {
            u = mCurves[isp].eval_h(gs.T, logT) - mR[isp]*gs.T;
        }
        number u = mass_average(gs, mVals);
        return u;
    }

    @nogc
    override number energyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        /* This is a single-T gas. */
        return mCurves[isp].eval_h(gs.T) - mR[isp]*gs.T;
    }

    @nogc
    override number enthalpy(in GasState gs)
    {
        number logT = log(gs.T);
        foreach (isp, ref h; mVals) {
            h = mCurves[isp].eval_h(gs.T, logT);
        }
        number h = mass_average(gs, mVals);
        return h;
    }

    @nogc
    override void enthalpies(in GasState gs, number[] hs)
    {
        number logT = log(gs.T);
        foreach (isp; 0 .. n_species) {
            hs[isp] = mCurves[isp].eval_h(gs.T, logT);
        }
    }

    @nogc
    override number enthalpyPerSpecies(in GasState gs, int isp)
    {
        return mCurves[isp].eval_h(gs.T);
    }

    @nogc
    override number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        /* This is a single-T gas. */
        return mCurves[isp].eval_h(gs.T);
    }

    @nogc
    override number entropy(in GasState gs)
    {
        foreach (isp, ref s; mVals) {
            s = mCurves[isp].eval_s(gs.T) - mR[isp]*log(gs.p/P_atm);
        }
        number s = mass_average(gs, mVals);
        return s;
    }

    @nogc
    override number entropyPerSpecies(in GasState gs, int isp)
    {
        return mCurves[isp].eval_s(gs.T) - mR[isp]*log(gs.p/P_atm);
    }

    @nogc override number cpPerSpecies(in GasState gs, int isp)
    {
        return mCurves[isp].eval_Cp(gs.T);
    }
    @nogc override void GibbsFreeEnergies(in GasState gs, number[] gibbs_energies)
    {
        number T = gs.T;
        number logT = log(gs.T);
        number logp = log(gs.p/P_atm);
        foreach(isp; 0 .. n_species){
            number h = mCurves[isp].eval_h(T, logT);
            number s = mCurves[isp].eval_s(T, logT) - mR[isp]*logp;
            gibbs_energies[isp] = h - T*s;
        }
    }

private:
    size_t n_species;
    double[] mR;
    CEAThermoCurve[] mCurves;
    number[] mVals; // a private work array

    @nogc
    number gamma(ref GasState gs)
    {
        return dhdTConstP(gs)/dudTConstV(gs);
    }

    @nogc
    void update_pressure(ref GasState gs)
    {
        number Rmix = gasConstant(gs);
        gs.p = gs.rho*Rmix*gs.T;
    }

    @nogc
    void update_density(ref GasState gs)
    {
        number Rmix = gasConstant(gs);
        gs.rho = gs.p/(Rmix*gs.T);
    }

    @nogc
    void update_temperature_from_rhop(ref GasState gs)
    {
        number Rmix = gasConstant(gs);
        gs.T = gs.p/(Rmix*gs.rho);
    }

    @nogc
    void update_energy(ref GasState gs)
    {
        gs.u = internalEnergy(gs);
    }

    @nogc
    number internal_energy(size_t nsp, in number[] massf, in number T, in number logT)
    {
        number u = 0.0;
        foreach (isp; 0 .. nsp) {
            u += massf[isp]*(mCurves[isp].eval_h(T, logT) - mR[isp]*T);
        }
        return u;
    }

    @nogc
    number constant_volume_specific_heat(size_t nsp, number[] massf, number T)
    {
        number Cv = 0.0;
        foreach (isp; 0 .. nsp) {
            Cv += massf[isp]*(mCurves[isp].eval_Cp(T) - mR[isp]);
        }
        return Cv;
    }

    @nogc
    number update_temperature_from_energy(ref GasState gs)
    {
    /*
        Bespoke newton's method by NNG. This was turning into a major
        bottleneck for multispecies flow, so it was time for some serious
        corner-cutting...

        @author: Nick Gibbons
    */
        immutable double TOL = 1e-10;
        immutable number u_tgt = gs.u;
        immutable size_t nsp = n_species;

        number T_guess = gs.T;
        number logT_guess = log(T_guess);
        number u = internal_energy(nsp, gs.massf, T_guess, logT_guess);
        number f_guess =  u - u_tgt;

        // Do at least one iteration, in case we need to set the imaginary components of T
        number Cv, dT;
        immutable int MAX_ITERATIONS = 10;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            Cv = constant_volume_specific_heat(nsp, gs.massf, T_guess);
            dT = -f_guess/Cv;
            T_guess += dT;
            if (fabs(dT.re)/T_guess.re < TOL) {
                return T_guess;
            }
            logT_guess = log(T_guess);
            f_guess = internal_energy(nsp, gs.massf, T_guess, logT_guess) - u_tgt;
        }

        if (fabs(dT)>1e-3) {
            string msg = "update_temperature_from_energy function failed to converge.\n";
            debug {
                msg ~= format("The final value for T was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= gs.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }
        return T_guess;
    }
}

version(therm_perf_gas_mix_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        doLuaFile(L, "sample-data/therm-perf-5-species-air.lua");
        string[] speciesNames;
        getArrayOfStrings(L, "species", speciesNames);
        auto tm = new ThermPerfGasMixture(L, speciesNames);
        lua_close(L);
        auto gs = GasState(5, 0);

        gs.p = 1.0e6;
        gs.T = 2000.0;
        gs.massf = [to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2)];
        tm.updateFromPT(gs);
        assert(approxEqualNumbers(to!number(11801825.6), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.2840117), gs.rho, 1.0e-6), failedUnitTest());

        gs.rho = 2.0;
        gs.u = 14.0e6;
        tm.updateFromRhoU(gs);
        assert(approxEqualNumbers(to!number(3373757.4), gs.p, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(4331.944), gs.T, 1.0e-6), failedUnitTest());

        gs.T = 10000.0;
        gs.rho = 1.5;
        tm.updateFromRhoT(gs);
        assert(approxEqualNumbers(to!number(5841068.3), gs.p, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(20340105.9), gs.u, 1.0e-6), failedUnitTest());

        gs.rho = 10.0;
        gs.p = 5.0e6;
        tm.updateFromRhoP(gs);
        assert(approxEqualNumbers(to!number(11164648.5), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1284.012), gs.T, 1.0e-6), failedUnitTest());

        gs.p = 1.0e6;
        number s = 10000.0;
        tm.updateFromPS(gs, s);
        assert(approxEqualNumbers(to!number(2560.118), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(12313952.52), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.00309), gs.rho, 1.0e-6), failedUnitTest());

        s = 11000.0;
        number h = 17.0e6;
        tm.updateFromHS(gs, h, s);
        assert(approxEqualNumbers(to!number(5273.103), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(14946629.7), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(0.4603513), gs.rho, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(945271.84), gs.p, 1.0e-4), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = gd.u; // copy unperturbed value, but we don't really need it
            double ih = 1.0e-20;
            gs.T += complex(0.0,ih);
            tm.updateFromRhoT(gs);
            double myCv = gs.u.im/ih;
            assert(isClose(myCv, tm.dudTConstV(gs).re), failedUnitTest());
        }

        return 0;
    }
}

