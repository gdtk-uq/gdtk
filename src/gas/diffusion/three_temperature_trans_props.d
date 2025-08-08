/**
 * two_temperature_trans_props.d
 *
 * Author: Rowan G.
 * Date: 2021-03-16
 * History: 2021-03-16 -- code ported from two_temperature air model.
 *
 * References
 * ----------
 * Gnoffo, Gupta and Shin (1989)
 * Conservation equations and physical models for hypersonic air flows
 * in thermal and chemical nonequilibrium.
 * NASA TP-2867, February 1989, NASA Langley Research Center
 *
 * Gupta, Yos and Thompson (1989)
 * A review of reaction rate and thermodynamic and transport properties
 * for the 11-species air model for chemical and thermal nonequilibrium
 * calculations to 30,000 K.
 * NASA TM-101528, February 1989, NASA Langley Research Center
 **/

module gas.diffusion.three_temperature_trans_props;

import std.math;
import std.string;
import std.conv;
import std.algorithm;

import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas;
import gas.diffusion.transport_properties_model;
import gas.physical_constants;

class ThreeTemperatureTransProps : TransportPropertiesModel {
public:

    this(lua_State *L, string[] speciesNames)
    {
        mNSpecies = to!int(speciesNames.length);
        mElectronicMode = getInt(L, "electronic_mode");
        double[] molMasses;
        mMolMasses.length = mNSpecies;
        mParticleMass.length = mNSpecies;
        mMolef.length = mNSpecies;
        mCharge.length = mNSpecies;
        allocateArraysForCollisionIntegrals();

        foreach (isp, spName; speciesNames) {
            if (spName == "e-") mElectronIdx = to!int(isp);
            lua_getglobal(L, "db");
            lua_getfield(L, -1, spName.toStringz);
            double m = getDouble(L, -1, "M");
            mMolMasses[isp] = m;
            mParticleMass[isp] = 1000.0*m/Avogadro_number; // kg -> g
            mCharge[isp] = getInt(L, -1, "charge");
            string type = getString(L, -1, "type");
            if (type == "molecule") {
                mMolecularSpecies ~= to!int(isp);
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        // Fill in support data for collision integrals
        lua_getglobal(L, "db");
        lua_getfield(L, -1, "CIs");
        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                string key = speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                lua_getfield(L, -1, key.toStringz);
                if (lua_isnil(L, -1)) {
                    lua_pop(L, 1);
                    // Try reverse order, eg. N2:O2 --> O2:N2
                    key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                    lua_getfield(L, -1, key.toStringz);
                    if (lua_isnil(L, -1)) {
                        // There's still a problem, can't find entry for CI data.
                        string msg = format("Collision integral data for '%s' is missing.\n", key);
                        throw new Error(msg);
                    }
                }
                lua_getfield(L, -1, "pi_Omega_11");
                mA_11[isp][jsp] = getDouble(L, -1, "A");
                mB_11[isp][jsp] = getDouble(L, -1, "B");
                mC_11[isp][jsp] = getDouble(L, -1, "C");
                mD_11[isp][jsp] = getDouble(L, -1, "D");
                lua_pop(L, 1); // pop: pi_Omega_11
                lua_getfield(L, -1, "pi_Omega_22");
                mA_22[isp][jsp] = getDouble(L, -1, "A");
                mB_22[isp][jsp] = getDouble(L, -1, "B");
                mC_22[isp][jsp] = getDouble(L, -1, "C");
                mD_22[isp][jsp] = getDouble(L, -1, "D");
                lua_pop(L, 1); // pop: pi_Omega_22
                lua_pop(L, 1); // pop: collision pair.
            }
        }
        lua_pop(L, 1); // pop: cis
        lua_pop(L, 1); // pop: db

        // Compute alphas
        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. mNSpecies) {
                double M_isp = mMolMasses[isp];
                double M_jsp = mMolMasses[jsp];
                mMu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                mMu[isp][jsp] *= 1000.0; // convert kg/mole --> g/mole
                double M_ratio = M_isp/M_jsp;
                double numer = (1.0 - M_ratio)*(0.45 - 2.54*M_ratio);
                double denom = (1.0 + M_ratio)^^2;
                mAlpha[isp][jsp] = 1.0 + numer/denom;
                mIsCoulombCollision[isp][jsp] = (mCharge[isp]!=0)&&(mCharge[jsp]!=0);
            }
        }
    }

    void allocateArraysForCollisionIntegrals()
    {
        mMu.length = mNSpecies;
        mAlpha.length = mNSpecies;
        mA_11.length = mNSpecies;
        mB_11.length = mNSpecies;
        mC_11.length = mNSpecies;
        mD_11.length = mNSpecies;
        mDelta_11.length = mNSpecies;
        mA_22.length = mNSpecies;
        mB_22.length = mNSpecies;
        mC_22.length = mNSpecies;
        mD_22.length = mNSpecies;
        mDelta_22.length = mNSpecies;
        mIsCoulombCollision.length = mNSpecies;
        foreach (isp; 0 .. mNSpecies) {
            mMu[isp].length = mNSpecies;
            mAlpha[isp].length = mNSpecies;
            mA_11[isp].length = isp+1;
            mB_11[isp].length = isp+1;
            mC_11[isp].length = isp+1;
            mD_11[isp].length = isp+1;
            mDelta_11[isp].length = mNSpecies;
            mA_22[isp].length = isp+1;
            mB_22[isp].length = isp+1;
            mC_22[isp].length = isp+1;
            mD_22[isp].length = isp+1;
            mDelta_22[isp].length = mNSpecies;
            mIsCoulombCollision[isp].length = mNSpecies;
        }
    }

    @nogc
    override void updateTransProps(ref GasState gs, GasModel gm)
    {
        massf2molef(gs.massf, mMolMasses, mMolef);
        // Computation of transport coefficients via collision integrals.
        // Equations follow those in Gupta et al. (1990)
        computeDelta22(gs);
        computeDelta11(gs);

        // Compute mixture viscosity.
        number sumA = 0.0;
        number denom;
        foreach (isp; 0 .. mNSpecies) {
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_22[isp][jsp];
            }
            if (isp == mElectronIdx) continue;
            sumA += mParticleMass[isp]*mMolef[isp]/denom;
        }
        // Add term if electron present.
        if (mElectronIdx != -1) {
            // An additional term is required in the mixture viscosity.
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_22[mElectronIdx][jsp];
            }
            sumA += mParticleMass[mElectronIdx]*mMolef[mElectronIdx]/denom;

        }
        gs.mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)

        // Compute component thermal conductivities
        // k in transrotational = k_tr + k_rot
        // k in vibroelectronic = k_ve + k_E
        // 1. k_tr
        sumA = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                if (jsp != mElectronIdx) {
                    denom += mAlpha[isp][jsp]*mMolef[jsp]*mDelta_22[isp][jsp];
                }
                else {
                    denom += 3.54*mAlpha[isp][jsp]*mMolef[jsp]*mDelta_22[isp][jsp];
                }
            }
            if (isp == mElectronIdx) continue;
            sumA += mMolef[isp]/denom;
        }
        double kB_erg = 1.38066e-16; // erg/K
        number k_tr = 2.3901e-8*(15./4.)*kB_erg*sumA;
        k_tr *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        // 2. k_rot and k_vib at the same time
        // Assuming rotation is fully excited, eq (75) in Gnoffo
        // and vibration is partially excited eq (52a) in Gupta, Yos, Thomson
        number k_rot = 0.0;
        number k_vib = 0.0;

        // For k_vib, Cp(T_vib) needs to be evaluated
        foreach (isp; mMolecularSpecies) {
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_11[isp][jsp];
            }
            k_rot += mMolef[isp]/denom;

            // figure out how much vibrational excitation there is
            // The equation in the paper doesn't call for the min,
            // but it has been included to get the same result as
            // when the flow is fully excited
            number Cp = mMolMasses[isp]/R_universal*gm.Cp(gs, isp);
            k_vib += fmax((fmin(Cp.re, 9.0/2.0) - 7.0/2.0), 0.0) * mMolef[isp]/denom;
        }

        // electronic contribution
        number k_ee = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            // Free electrons don't contribute to the electronic thermal
            // conduction
            if (isp == mElectronIdx) continue;
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_11[isp][jsp];
            }
            number Cp = mMolMasses[isp]/R_universal*gm.Cp(gs, isp);
            k_ee += fmax((Cp.re - 9.0/2.0), 0.0) * mMolef[isp]/denom;
        }

        k_rot *= 2.3901e-8*kB_erg;
        k_rot *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        k_vib *= 2.3901e-8*kB_erg;
        k_vib *= (4.184/1.0e-2); // cal/cm.s.K) --> J/(m.s.K)
        k_ee *= 2.3901e-8*kB_erg;
        k_ee *= (4.184/1.0e-2); // cal/cm.s.K) --> J/(m.s.K)

        // Eq (76) in Gnoffo.
        gs.k = k_tr + k_rot;

        // 4. k_e
        number k_E = 0.0;
        if (mElectronIdx != -1) {
            // electron present.
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += 1.45*fmax(0.0, mMolef[jsp])*mDelta_22[mElectronIdx][jsp];
            }
            k_E = fmax(0.0, mMolef[mElectronIdx])/denom;
            k_E *= 2.3901e-8*(15./4.)*kB_erg;
            k_E *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        }
        gs.k_modes[0] = k_vib;
        gs.k_modes[1] = k_E;
        // gs.k_modes[mElectronicMode] += k_ee;
    }

    @nogc void binaryDiffusionCoefficients(ref const(GasState) gs, ref number[][] D)
    {
    /*
        Compute the binary diffusion coefficients using the Omega_11 collisions
        integrals and equations 42a from Gupta et al. 1990.

        Notes:
         - Prior to December 2021, this routine was returning diffusion coefficients
           that were 100 times too small. This was due to a unit conversion error
           related to Delta_11, which is in cm-s, not m-s. (NNG)

         - TODO: Think about "p" in this equation. Should it be total pressure, or
           bath pressure of binary interactors only? (RJG)
    */
        massf2molef(gs.massf, mMolMasses, mMolef);
        computeDelta11(gs);
        immutable double kB_erg = 1.38066e-16; // erg/K
        number p_cgs = gs.p*10.0; // CGS pressure in Baryes (?!)

        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                number T = (isp == mElectronIdx || jsp == mElectronIdx) ? gs.T_modes[$-1] : gs.T;
                number Dij = (kB_erg * T)/(p_cgs * mDelta_11[isp][jsp]); // cm^2/s
                Dij /= (100*100); // cm2/s -> m2/s
                D[isp][jsp] = Dij;
                D[jsp][isp] = Dij;
            }
        }
    }


private:
    int mNSpecies;
    int mElectronIdx = -1;
    int mElectronicMode;
    int[] mMolecularSpecies;
    double[] mMolMasses;
    double[] mParticleMass;
    int[] mCharge;
    // working array space
    number[] mMolef;
    number[][] mA_11, mB_11, mC_11, mD_11, mDelta_11, mAlpha;
    number[][] mA_22, mB_22, mC_22, mD_22, mDelta_22, mMu;
    bool[][] mIsCoulombCollision;

    @nogc
    number electronPressureCorrection(ref const(GasState) gs)
    {
    /*
        Compute a correction factor for the collision cross-sections of the charged species,
        based on the electron pressure. This is required as per Gupta, 1989, as explained on
        page 20.

        Notes:
         - Instead of computing p_em, the maximum allowed pressure, we just limit log(Lambda)
           to a minimum of one, which is equivalent, but neater.
         - This function becomes undefined for p_e==0.0, which can happen if no electrons are
           present. Analytically, this would be fine because later on the result will get
           multiplied by the electron mole fraction, which is also equal to zero. But numerically
           the divide be zero causes problems. To fix this we just limit p_e. It should be
           impossible for termTerm or logLambda to overflow as long as p_e is finite, so I chose
           a really small limit.
         - It would be nicer if we could solve this problem by switching to a series expansion
           as p_e gets small, but it messes up the code because you need the expression in terms
           of p_e*log(Lambda).

        @author: Nick Gibbons
    */
        number p_e = fmax(1e-32, gs.p_e);
        number p_e_atm = p_e/P_atm;
        number tempTerm = gs.T/1000.0/pow(p_e_atm, 0.25);
        number logLambda = 0.5*log(2.09e-2*pow(tempTerm, 4.0) + 1.52*pow(tempTerm, 8./3.));
        logLambda = fmax(logLambda, 1.0); // Equivalent to limiting by p_em, from equation 23g
        return logLambda;
    }

    @nogc
    void computeDelta11(ref const(GasState) gs)
    {
        number logLambda = 1.0;
        if (mElectronIdx != -1) logLambda = electronPressureCorrection(gs);

        double kB = Boltzmann_constant;
        number T_CI;
        number log_T_CI;

        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                 if (isp != mElectronIdx) {
                    // heavy-particle colliders: use transrotational temperature in calculation
                    T_CI = gs.T;
                    log_T_CI = log(T_CI);
                }
                else {
                    // collisions with electron: use electron temperature in calculation
                    // assume the electron temperature is the last modal temperature
                    T_CI = gs.T_modes[$-1];
                    log_T_CI = log(T_CI);
                }
                number expnt = mA_11[isp][jsp]*(log_T_CI)^^2 + mB_11[isp][jsp]*log_T_CI + mC_11[isp][jsp];
                number pi_Omega_11 = exp(mD_11[isp][jsp])*pow(T_CI, expnt);
                if (mIsCoulombCollision[isp][jsp])
                    pi_Omega_11 *= logLambda;

                mDelta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*mMu[isp][jsp]/(to!double(PI)*R_universal_cal*T_CI))*pi_Omega_11;
                mDelta_11[jsp][isp] = mDelta_11[isp][jsp];
            }
        }
    }

    @nogc
    void computeDelta22(ref const(GasState) gs)
    {
        number logLambda = 1.0;
        if (mElectronIdx != -1) logLambda = electronPressureCorrection(gs);

        double kB = Boltzmann_constant;
        number T_CI;
        number log_T_CI;

        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                if (isp != mElectronIdx) {
                    // heavy-particle colliders: use transrotational temperature in calculation
                    T_CI = gs.T;
                    log_T_CI = log(T_CI);
                }
                else {
                    // collisions with electron: use vibroelectronic temperature in calculation
                    T_CI = gs.T_modes[$-1];
                    log_T_CI = log(T_CI);
                }
                number expnt = mA_22[isp][jsp]*(log_T_CI)^^2 + mB_22[isp][jsp]*log_T_CI + mC_22[isp][jsp];
                number pi_Omega_22 = exp(mD_22[isp][jsp])*pow(T_CI, expnt);
                if (mIsCoulombCollision[isp][jsp])
                    pi_Omega_22 *= logLambda;

                mDelta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*mMu[isp][jsp]/(to!double(PI)*R_universal_cal*T_CI))*pi_Omega_22;
                mDelta_22[jsp][isp] = mDelta_22[isp][jsp];
            }
        }
    }
}


version(three_temperature_trans_props_test)
{
    int main()
    {
        import util.msg_service;
        import gas.composite_gas;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        GasModel gm = new CompositeGas("sample-data/N2-N.lua");
        doLuaFile(L, "sample-data/N2-N.lua");
        string[] speciesNames;
        getArrayOfStrings(L, "species", speciesNames);
        auto ttp = new ThreeTemperatureTransProps(L, speciesNames);
        lua_close(L);
        auto gs = GasState(2, 1);
        gs.p = 1.0e5;
        gs.T = 2000.0;
        gs.T_modes[0] = 3000.0;
        gs.massf[0] = 0.8;
        gs.massf[1] = 0.2;

        ttp.updateTransProps(gs, gm);

        import std.stdio;
        writefln("mu= %.6e  k= %.6e  k_v= %.6e\n", gs.mu, gs.k, gs.k_modes[0]);

        gs.T_modes[0] = 300.0;
        ttp.updateTransProps(gs, gm);
        writefln("mu = %.6e k = %.6e k_v = %.6e", gs.mu, gs.k, gs.k_modes[0]);

        //auto L = init_lua_State();
        //GasModel gm = new CompositeGas("gm-air11-2T.lua");
        //doLuaFile(L, "gm-air11-2T.lua");
        //string[] speciesNames;
        //getArrayOfStrings(L, "species", speciesNames);
        //auto ttp = new TwoTemperatureTransProps(L, speciesNames);
        //lua_close(L);
        //auto gs = GasState(11, 1);
        //gs.p = 0.1*101.35e3;
        //gs.T = 1000.0;
        //gs.T_modes[0] = 3000.0;

        //gs.massf[0] =     7.455871e-01; // N2
        //gs.massf[1] =     2.543794e-01; // O2
        //gs.massf[2] =     0.000000e+00; // N
        //gs.massf[3] =     1.310344e-10; // O
        //gs.massf[4] =     3.355965e-05; // NO
        //gs.massf[5] =     0.000000e+00; // N2+
        //gs.massf[6] =     0.000000e+00; // O2+
        //gs.massf[7] =     0.000000e+00; // N+
        //gs.massf[8] =     0.000000e+00; // O+
        //gs.massf[9] =     0.000000e+00; // NO+
        //gs.massf[10] =     0.000000e+00; // e-


        //gm.update_thermo_from_pT(gs);
        //ttp.updateTransProps(gs, gm);
        //import std.stdio;
        //writefln("mu= %.6e  k= %.6e  k_v= %.6e\n", gs.mu, gs.k, gs.k_modes[0]);
        return 0;
    }
}
