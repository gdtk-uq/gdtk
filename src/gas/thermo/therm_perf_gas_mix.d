/*
 * gas_mixture_thermo.d
 * 
 * Author: Rowan G. and Peter J.
 * Version: 2021-02-09
 */

module gas.thermo.gas_mixture_thermo;

import nm.brent : solve;
import nm.bracketing : bracket;

class ThermPerfGasMix : ThermoModel {
public:

    /* Bundle interacting parts here.

    
    @nogc
    override void update_thermo_from_pT(GasState gs)
    {
        mThermalEOS.update_density(gs);
        mCaloricEOS.update_energy(gs);
    }
    @nogc
    override void update_thermo_from_rhou(GasState gs)
    {
        mCaloricEOS.update_temperature(gs);
        mThermalEOS.update_pressure(gs);
    }
    @nogc
    override void update_thermo_from_rhoT(GasState gs)
    {
        mCaloricEOS.update_energy(gs);
        mThermalEOS.update_pressure(gs);
    }
    @nogc
    override void update_thermo_from_rhop(GasState gs)
    {
        mThermalEOS.update_temperature(gs);
        mCaloricEOS.update_energy(gs);
    }
    @nogc
    override void update_thermo_from_ps(GasState gs, number s)
    {
        double TOL = 1.0e-6;
        number delT = 100.0;
        number T1 = Q.T;
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
            gs.T = solve!(zeroFun,number)(T1, T2, TOL);
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
        mCaloricEOS.update_energy(Q);
        mThermalEOS.update_density(Q);
    }

    @nogc
    override void update_thermo_from_hs(GasState gs, number h, number s)
    {
        // We do this in two stages.
        // First, from enthalpy we compute temperature.
        double TOL = 1.0e-6;
        number delT = 100.0;
        number T1 = gs.T;
        // It's possible that T is 'nan' if it's never been set, so:
        if ( isNaN(T1) )
            T1 = 300.0; // Just set at some value to get started.
        number Tsave = T1;
        number T2 = T1 + delT;

        number zeroFun(number T)
        {
            gs.T = T;
            number h_guess = enthalpy(Q);
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
            gs.T = solve!(zeroFun,number)(T1, T2, TOL);
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
        number p1 = Q.p;
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
            gs.p = solve!(zeroFun2,number)(p1, p2, TOL);
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find pressure";
            debug {
                msg ~= "\nin function GasMixtureThermo.update_thermo_from_hs().\n";
                msg ~= format("The initial pressure guess was: %12.6f\n", psave);
                msg ~= format("The target entropy value was: %12.6f\n", s);
                msg ~= format("The GasState is currently:\n");
                msg ~= Q.toString();
                msg ~= "The message from the brent.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }
        mCaloricEOS.update_energy(Q);
        mThermalEOS.update_density(Q);
    }
     
    @nogc override void update_sound_speed(GasState Q)
    {
        // Reference:
        // Cengel and Boles (1998)
        // Thermodynamics: an Engineering Approach, 3rd edition
        // McGraw Hill
        // Equation 16-10 on p. 849

        // "frozen" sound speed
        gs.a = sqrt(gamma(gs)*dpdrho_const_T(gs));
    }

    
    

private:
    EVT_EOS mCaloricEOS;
    PVT_EOS mThermalEOS;
}

