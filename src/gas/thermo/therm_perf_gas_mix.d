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
import util.lua;
import util.lua_service;
import nm.nm_exception : NumericalMethodException;
import nm.complex;
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
        auto n_species = speciesNames.length;
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

    @nogc
    override void updateFromPT(GasState gs)
    {
        update_density(gs);
        update_energy(gs);
    }
    @nogc
    override void updateFromRhoU(GasState gs)
    {
        update_temperature_from_energy(gs);
        update_pressure(gs);
    }
    @nogc
    override void updateFromRhoT(GasState gs)
    {
        update_energy(gs);
        update_pressure(gs);
    }
    @nogc
    override void updateFromRhoP(GasState gs)
    {
        update_temperature_from_rhop(gs);
        update_energy(gs);
    }
    @nogc
    override void updateFromPS(GasState gs, number s)
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

        if (bracket!(zeroFun,number)(T1, T2) == -1) {
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
    override void updateFromHS(GasState gs, number h, number s)
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
    override void updateSoundSpeed(GasState gs)
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
    override number internalEnergy(in GasState gs)
    {
        foreach (isp, ref u; mVals) {
            u = mCurves[isp].eval_h(gs.T) - mR[isp]*gs.T;
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
        foreach (isp, ref h; mVals) {
            h = mCurves[isp].eval_h(gs.T);
        }
        number h = mass_average(gs, mVals);
        return h;
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

    
    
private:
    double[] mR;
    CEAThermoCurve[] mCurves;
    number[] mVals; // a private work array

    @nogc
    number gamma(GasState gs)
    {
        return dhdTConstP(gs)/dudTConstV(gs);
    }
    
    @nogc
    void update_pressure(GasState gs)
    {
        number Rmix = gasConstant(gs);
        gs.p = gs.rho*Rmix*gs.T;
    }

    @nogc
    void update_density(GasState gs)
    {
        number Rmix = gasConstant(gs);
        gs.rho = gs.p/(Rmix*gs.T);
    }

    @nogc
    void update_temperature_from_rhop(GasState gs)
    {
        number Rmix = gasConstant(gs);
        gs.T = gs.p/(Rmix*gs.rho);
    }
    
    @nogc
    void update_energy(GasState gs)
    {
        gs.u = internalEnergy(gs);
    }

    @nogc
    void update_temperature_from_energy(GasState gs)
    {
        number Tsave = gs.T; // Keep a copy for diagnostics purpose.
        double TOL = 1.0e-6;
        //
        // PJ 2020-01-01
        // Hack to cover over low energy problems with large ionization fractions.
        // Should not be needed if things are being solved well.
        //   number u_original = Q.u;
        //   Q.T = to!number(T_MIN+1.0);
        //   update_energy(Q);
        //   if (u_original > Q.u) { Q.u = u_original; }
        // End of hack.
        //
        // The "target" energy is the value we will iterate to find.
        // We search (via a numerical method) for the temperature
        // value that gives this target energy.
        number e_tgt = gs.u;
        // delT is the initial guess for a bracket size.
        // We set this quite agressivley at 100 K hoping to
        // keep the number of iterations required to a small
        // value. In general, we are taking small timesteps
        // so the value of temperature from the previous step
        // should not change too much. If it does, there should
        // be enough robustness in the bracketing and
        // the function-solving method to handle this.
        number delT = 1000.0;
        version(complex_numbers) {
            number T1 = nm.complex.fmax(gs.T - 0.5*delT, to!number(T_MIN));
        } else {
            double T1 = std.math.fmax(gs.T - 0.5*delT, T_MIN);
        }
        number T2 = T1 + delT;

        number zeroFn(number T)
        {
            gs.T = T;
            update_energy(gs);
            /*
            debug {
                writefln("T= %.3f  u= %.3e  Fn= %.6e ", T, Q.u, e_tgt - Q.u);
            }
            */
            return e_tgt - gs.u;
        }
        
        number dzdT(number T)
        {
            // We evaluate Cv. And return the negative.
            gs.T = T;
            return -1.0*dudTConstV(gs);
        }

        try {
            gs.T = nm.newton.solve!(zeroFn, dzdT)(Tsave, T1, T2, TOL);
        }
        catch (NumericalMethodException e) {
            // If we fail, we'll have a second attempt with an extended temperature range.
            try {
                gs.T = solve!(zeroFn, dzdT)(Tsave, to!number(T_MIN), to!number(T_MAX), TOL);
            }
            catch (NumericalMethodException e) {
                string msg = "There was a problem iterating to find temperature";
                debug {
                    msg ~= "\nin function ThermallyPerfectGasMix.update_temperature().\n";
                    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                    msg ~= format("The target energy value was: %12.6f\n", e_tgt);
                    msg ~= format("zeroFn(Tsave)=%g\n", zeroFn(Tsave));
                    msg ~= format("zeroFn(T_MIN)=%g\n", zeroFn(to!number(T_MIN)));
                    msg ~= format("zeroFn(T_MAX)=%g\n", zeroFn(to!number(T_MAX)));
                    msg ~= format("The GasState is currently:\n");
                    msg ~= gs.toString() ~ "\n";
                    msg ~= "The message from the solve function is:\n";
                    msg ~= e.msg;
                }
                // If we fail, leave temperature at value upon entry to method.
                gs.T = Tsave;
                update_energy(gs);
                throw new GasModelException(msg);
            }
        }
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
        getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", speciesNames);
        auto tm = new ThermPerfGasMixture(L, speciesNames);
        lua_close(L);
        auto gs = new GasState(5, 0);

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

