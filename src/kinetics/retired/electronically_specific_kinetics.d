/**
 * Authors: Brad Semple
 * Date: 2018-09-04
 *
 * This module provides an integrator for the master
 * equation for an electronic state-specific system
 *
 */

module kinetics.electronically_specific_kinetics;

import core.stdc.stdlib : exit;
import std.stdio;
import std.conv : to;
import std.string;
import std.math;
import std.algorithm;
import std.file;

import ntypes.complex;
import nm.number;
import nm.smla;
import util.lua;
import util.lua_service;
import gas;
import gas.physical_constants;
import gas.electronically_specific_gas;
import gas.electronic_species;
import gas.two_temperature_air;

import kinetics.thermochemical_reactor;
import kinetics.chemistry_update;
import kinetics.electronic_update;


immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
                                             // Yes, you read that right. Sometimes a step is succesful
                                             // but the timestep selection algorithm will suggest
                                             // a reduction. We limit that reduction to no more than 50%.
immutable double DT_REDUCTION_FACTOR = 10.0; // factor by which to reduce timestep
                                             // after a failed attempt
immutable double H_MIN = 1.0e-15; // Minimum allowable step size
enum ResultOfStep { success, failure };

final class ElectronicallySpecificKinetics : ThermochemicalReactor {
    this(string listOfFiles, GasModel gmodel)
    {

        auto L_filenames = init_lua_State();
        doLuaFile(L_filenames,listOfFiles);

        auto energyExchFile = getString(L_filenames, "energyExchFile");
        auto chemFile = getString(L_filenames,"chemFile");
        auto twotemperaturegmodel = getString(L_filenames,"macro_species_filename");
        auto ESK_N_Filename = getString(L_filenames,"ESK_N_Filename");
        auto ESK_O_Filename = getString(L_filenames,"ESK_O_Filename");

        _n_N_species = getInt(L_filenames, "number_N_species");
        _n_O_species = getInt(L_filenames, "number_O_species");
        _n_elec_species = _n_N_species + _n_O_species;

        auto L_TT = init_lua_State();
        doLuaFile(L_TT,twotemperaturegmodel);

        // Initialise much of the macro species data
        auto macro_gmodel = new TwoTemperatureAir(L_TT);
        _macroAirModel = cast(TwoTemperatureAir) macro_gmodel;
        if (_macroAirModel is null) {
            string errMsg = "Error in construction of TwoTemperatureAir.\n";
            errMsg ~= "The supplied gas model must be a TwoTemperature model.\n";
            throw new ThermochemicalReactorUpdateException(errMsg);
        }

        super(gmodel);
        _fullAirModel = cast(ElectronicallySpecificGas) gmodel;
        if (_fullAirModel is null) {
            string errMsg = "Error in construction of ElectronicallySpecificGas.\n";
            errMsg ~= "The supplied gas model must be an ElectronicallySpecific model.\n";
            throw new ThermochemicalReactorUpdateException(errMsg);
        }

        _n_macro_species = _macroAirModel.n_species;

        _macro_Q = GasState(macro_gmodel);
        _macro_Qinit = GasState(macro_gmodel);
        _macro_Q0 = GasState(macro_gmodel);
        _macro_chemUpdate = new ChemistryUpdate(chemFile, macro_gmodel);
        _macro_molef.length = _n_macro_species;
        _macro_numden.length = _macroAirModel.n_species;
        initModel(energyExchFile);

        //Create the electronic objects

        ES_N = new ElectronicUpdate(ESK_N_Filename, [0,_n_elec_species+3], _fullAirModel);

        ES_O = new ElectronicUpdate(ESK_O_Filename, [_n_N_species, _n_elec_species + 4], _fullAirModel);
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        // This section is the macro species updates
        //First, create test macro species gas state, _macro_Q from given state Q
        Update_Macro_State(_macro_Q, Q);

        number uTotal = _macro_Q.u + _macro_Q.u_modes[0];
        // 1. Perform chemistry update.

        _macro_chemUpdate(_macro_Q, tInterval, dtSuggest, params);

        debug {
            // writeln("--- 1 ---");
            // writefln("uTotal= %.12e u= %.12e uv= %.12e", uTotal, _macro_Q.u, _macro_Q.u_modes[0]);
            // writefln("T= %.12e  Tv= %.12e", _macro_Q.T, _macro_Q.T_modes[0]);
        }
        // Changing mass fractions does not change the temperature
        // of the temperatures associated with internal structure.
        // Now we can adjust the transrotational energy, given that
        // some of it was redistributed to other internal structure energy bins.
        _macro_Q.u_modes[0] = _macroAirModel.vibElecEnergy(_macro_Q, _macro_Q.T_modes[0]);
        _macro_Q.u = uTotal - _macro_Q.u_modes[0];
        try {
            _macroAirModel.update_thermo_from_rhou(_macro_Q);
            debug {
                // writeln("--- 2 ---");
                // writefln("uTotal= %.12e u= %.12e uv= %.12e", uTotal, _macro_Q.u, _macro_Q.u_modes[0]);
                // writefln("T= %.12e  Tv= %.12e", _macro_Q.T, _macro_Q.T_modes[0]);
            }
        }
        catch (GasModelException err) {
            string msg = "Call to update_thermo_from_rhou failed in electronically specific kinetics.";
            debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
        }

        // 2. Perform energy exchange update.
        // As prep work, compute mole fractions.
        // These values won't change during the energy update
        // since the composition does not change
        // during energy exhange.

        _macroAirModel.massf2molef(_macro_Q, _macro_molef);
        _macroAirModel.massf2numden(_macro_Q, _macro_numden);
        try {
            energyUpdate(_macro_Q, tInterval, dtSuggest);
            debug {
                // writeln("--- 3 ---");
                // writeln(_macro_Q);
            }
        }
        catch (GasModelException err) {
            string msg = "The energy update in the electronically specific kinetics module failed.\n";
            debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
        }

        Update_Electronic_State(Q,_macro_Q);
        //update the electronic distribution

        //Remember, gas model always initialised in this order:
        //N, O, N2, O2, NO, N+, O+, N2+, O2+, NO+, e-
        //input in electronic solver takes the form
        //[N1, N2, N3 .... Nn, N+, O1, O2, O3 ... On, O+, e-]

        ES_N.Update(Q, tInterval);
        ES_O.Update(Q, tInterval);

    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for electronically_specific_kinetics.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    ElectronicallySpecificGas _fullAirModel;
    TwoTemperatureAir _macroAirModel;
    ElectronicUpdate ES_N, ES_O;

    //macro definitions

    GasState _macro_Qinit, _macro_Q0, _macro_Q;
    ChemistryUpdate _macro_chemUpdate;
    number[] _macro_molef;
    int _n_macro_species;
    double _c2 = 1.0;
    number _uTotal;
    number _T_sh, _Tv_sh;
    double[] _A;
    number[] _macro_numden;
    double[] _particleMass;
    double[][] _mu;
    int[][] _reactionsByMolecule;
    // Numerics control
    int _maxSubcycles = 10000;
    int _maxAttempts = 3;
    double _energyAbsTolerance = 1.0e-9;
    double _energyRelTolerance = 1.0e-9;

    //definitions for electronics
    number initialEnergy;
    number[] _numden;
    number[] _numden_input;
    number[] _numden_output;
    int _n_N_species, _n_O_species, _n_elec_species;

    //macro functions
    void initModel(string energyExchFile)
    {
        // Read parameters from file.
        auto L = init_lua_State();
        doLuaFile(L, energyExchFile);

        lua_getglobal(L, "parameters");
        _T_sh = getDouble(L, -1, "post_shock_Ttr");
        _Tv_sh = getDouble(L, -1, "post_shock_Tve");
        lua_getfield(L, -1, "c2");
        if (!lua_isnil(L, -1)) {
            _c2 = luaL_checknumber(L, -1);
        }
        lua_pop(L, 1);
        lua_pop(L, 1);

        lua_getglobal(L, "numerics");
        if (!lua_isnil(L, -1)) {
            lua_getfield(L, -1, "maxSubcycles");
            if (!lua_isnil(L, -1)) {
                _maxSubcycles = to!int(luaL_checkinteger(L, -1));
            }
            lua_pop(L, 1);

            lua_getfield(L, -1, "maxAttempts");
            if (!lua_isnil(L, -1)) {
                _maxAttempts = to!int(luaL_checkinteger(L, -1));
            }
            lua_pop(L, 1);

            lua_getfield(L, -1, "absTolerance");
            if (!lua_isnil(L, -1)) {
                _energyAbsTolerance = to!double(luaL_checknumber(L, -1));
            }
            lua_pop(L, 1);

            lua_getfield(L, -1, "relTolerance");
            if (!lua_isnil(L, -1)) {
                _energyRelTolerance = to!double(luaL_checknumber(L, -1));
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);

        lua_close(L);

        int nSpecies = _macroAirModel.n_species;
        _A.length = nSpecies;
        _A[] = 0.0;
        foreach (isp; _macroAirModel.molecularSpecies) {
            _A[isp] = A_MW[_macroAirModel.species_name(isp)];
        }

        _mu.length = nSpecies;
        foreach (isp; 0 .. nSpecies) _mu[isp].length = nSpecies;
        foreach (isp; 0 .. nSpecies) {
            foreach (jsp; 0 .. nSpecies) {
                double M_isp = _macroAirModel.mol_masses[isp];
                double M_jsp = _macroAirModel.mol_masses[jsp];
                _mu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                _mu[isp][jsp] *= 1000.0; // convert kg/mole to g/mole
            }
        }

        _particleMass.length = nSpecies;
        foreach (isp; 0 .. nSpecies) {
            _particleMass[isp] = _macroAirModel.mol_masses[isp]/Avogadro_number;
        }

        _reactionsByMolecule.length = nSpecies;
        foreach (isp; _macroAirModel.molecularSpecies) {
            foreach (ir; 0 .. to!int(_macro_chemUpdate.rmech.n_reactions)) {
                if (_macro_chemUpdate.rmech.reactionHasParticipant(ir, isp)) {
                    _reactionsByMolecule[isp] ~= ir;
                }
            }
        }
    }

    @nogc
    void energyUpdate(ref GasState Q, double tInterval, ref double dtSuggest) //only pass macro Q
    {
        _uTotal = Q.u + Q.u_modes[0];
        // We borrow the algorithm from ChemistryUpdate.opCall()
        // Take a copy of what's passed in, in case something goes wrong
        _macro_Qinit.copy_values_from(Q);
        // 1. Sort out the time step for possible subcycling.
        double t = 0.0;
        double h;
        double dtSave;
        if ( dtSuggest > tInterval )
            h = tInterval;
        else if ( dtSuggest <= 0.0 )
            h = 0.01*tInterval;
        else
            h = dtSuggest;

        // 2. Now do the interesting stuff, increment energy change
        int cycle = 0;
        int attempt = 0;
        for ( ; cycle < _maxSubcycles; ++cycle) {
            /* Copy the last good timestep suggestion before moving on.
             * If we're near the end of the interval, we want
             * the penultimate step. In other words, we don't
             * want to store the fractional step of (tInterval - t)
             * that is taken as the last step.
             */
            dtSave = dtSuggest;
            h = fmin(h, tInterval - t);
            attempt= 0;
            for ( ; attempt < _maxAttempts; ++attempt) {
                ResultOfStep result = RKFStep(Q, h, dtSuggest);
                if (result == ResultOfStep.success) {
                    t += h;
                    /* We can now make some decision about how to
                     * increase the timestep. We will take some
                     * guidance from the ODE step, but also check
                     * that it's sensible. For example, we won't
                     * increase by more than 10% (or INCREASE_PERCENT).
                     * We'll also balk if the ODE step wants to reduce
                     * the stepsize on a successful step. This is because
                     * if the step was successful there shouldn't be any
                     * need (stability wise or accuracy related) to warrant
                     * a reduction. Thus if the step is successful and
                     * the dtSuggest comes back lower, let's just set
                     * h as the original value for the successful step.
                     */
                    double hMax = h*(1.0 + DT_INCREASE_PERCENT/100.0);
                    if (dtSuggest > h) {
                        /* It's possible that dtSuggest is less than h.
                         * When that occurs, we've made the decision to ignore
                         * the new timestep suggestion supplied by the ODE step.
                         * Our reasoning is that we'd like push for an aggressive
                         * timestep size. If that fails, then a later part of
                         * the timestepping algorithm will catch that and reduce
                         * the timestep.
                         *
                         * It's also possible that the suggested timestep
                         * is much larger than what we just used. For that
                         * case, we cap it at hMax. That's the following line.
                         */
                        h = fmin(dtSuggest, hMax);
                    }
                    break;
                }
                else { // in the case of failure...
                    /* We now need to make some decision about
                     * what timestep to attempt next.
                     */
                    h /= DT_REDUCTION_FACTOR;
                    if (h < H_MIN) {
                        string errMsg = "Hit the minimum allowable timestep in 2-T air energy exchange update.";
                        debug { errMsg ~= format("\ndt= %.4e", H_MIN); }
                        Q.copy_values_from(_macro_Qinit);
                        throw new ThermochemicalReactorUpdateException(errMsg);
                    }
                }
            } // end attempts at single subcycle.
            if (attempt == _maxAttempts) {
                string errMsg = "Hit maximum number of step attempts within a subcycle for 2-T air energy exchange update.";
                // We did poorly. Let's put the original GasState back in place,
                // and let the outside world know via an Exception.
                Q.copy_values_from(_macro_Qinit);
                throw new ThermochemicalReactorUpdateException(errMsg);
            }
            /* Otherwise, we've done well. */
            if (t >= tInterval) {
                // We've done enough work.
                // If we've only take one cycle, then we would like to use
                // dtSuggest rather than using the dtSuggest from the
                // penultimate step.
                if (cycle == 0) {
                    dtSave = dtSuggest;
                }
                break;
            }
        }
        if (cycle == _maxSubcycles) {
            // If we make it out here, then we didn't complete within
            // the maximum number of subscyles... we are taking too long.
            // Let's return the gas state to how it was before we failed
            // and throw an Exception.
            string errMsg = "Hit maximum number of subcycles while attempting 2-T air energy exchange update.";
            Q.copy_values_from(_macro_Qinit);
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        // At this point, it appears that everything has gone well.
        // Update energies and leave.
        try {
            Q.u = _uTotal - Q.u_modes[0];
            _macroAirModel.update_thermo_from_rhou(Q);
        }
        catch (GasModelException err) {
            string msg = "Call to update_thermo_from_rhou failed in two-temperature air kinetics.";
            debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
        }
        dtSuggest = dtSave;
    }

    @nogc
    number evalRelaxationTime(ref GasState Q, int isp)
    {
        number totalND = 0.0;
        number sum = 0.0;

        foreach (csp; 0 .. _macroAirModel.n_species) {
            // We exclude electrons in the vibrational relaxation time.
            if (_macroAirModel.species_name(csp) == "e-") {
                continue;
            }
            // 1. Compute Millikan-White value for each interaction
            if (_macro_molef[csp] >= SMALL_MOLE_FRACTION) {
                number nd = _macro_numden[csp];
                sum += nd * exp(_A[isp] * (pow(Q.T, -1./3) - 0.015*pow(_mu[isp][csp], 0.25)) - 18.42);
                totalND += nd;
            }
        }
        number tauMW = sum/totalND; // <-- p*tauMW (in atm.s)
        tauMW *= P_atm/Q.p;
        // 2. Compute Park value for high-temperature correction
        number nd = _macro_numden[isp];
        double kB = Boltzmann_constant;
        number cBar_s = sqrt(8*kB*Q.T/(to!double(PI)*_particleMass[isp]));
        double sigma_s = 1e-20; // Gnoffo p. 17 gives 1.0e-16 cm^2,
                                // converted to m^2 this is: 1.0e-20 m^2
        number tauP = 1.0/(sigma_s*cBar_s*nd);
        // 3. Combine M-W and Park values.
        number tauV = tauMW + tauP;
        return tauV;
        /* RJG, 2018-12-22
           Remove this correction for present. It is causing issues.
           Need to do some deeper research on this correction.
        // 4. Correct Landau-Teller at high temperatures by
        //    modifying tau.
        number s = 3.5*exp(-5000.0/_T_sh);
        number tauInv = (1.0/tauV)*(pow(fabs((_T_sh - Q.T_modes[0])/(_T_sh - _Tv_sh)), s-1.0));
        return 1.0/tauInv;
           END: RJG, 2018-12-22
        */
    }

    @nogc
    number evalRate(ref GasState Q)
    {
        number rate = 0.0;
        foreach (isp; _macroAirModel.molecularSpecies) {
            // Vibrational energy exchange via collisions.
            number tau = evalRelaxationTime(Q, isp);
            number evStar = _macroAirModel.vibElecEnergy(Q.T, isp);
            number ev = _macroAirModel.vibElecEnergy(Q.T_modes[0], isp);
            rate += Q.massf[isp] * (evStar - ev)/tau;
            // Vibrational energy change due to chemical reactions.
            number chemRate = 0.0;
            foreach (ir; _reactionsByMolecule[isp]) {
                chemRate += _macro_chemUpdate.rmech.rate(ir, isp);
            }
            chemRate *= _macroAirModel.mol_masses[isp]; // convert mol/m^3/s --> kg/m^3/s
            number Ds = _c2*ev;
            rate += Q.massf[isp]*chemRate*Ds;
        }
        return rate;
    }

    @nogc
    ResultOfStep RKFStep(ref GasState Q, double h, ref double hSuggest)
    {
        _macro_Q0.copy_values_from(Q);
        immutable double a21=1./5., a31=3./40., a32=9./40.,
            a41=3./10., a42=-9./10., a43 = 6./5.,
            a51=-11./54., a52=5./2., a53=-70./27., a54=35./27.,
            a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096.;
        immutable double b51=37./378., b53=250./621., b54=125./594., b56=512./1771.,
            b41=2825./27648., b43=18575./48384., b44=13525./55296., b45=277./14336., b46=1.0/4.0;

        number k1 = evalRate(Q);
        Q.u_modes[0] = _macro_Q0.u_modes[0] + h*(a21*k1);
        Q.u = _uTotal - Q.u_modes[0];
        _macroAirModel.update_thermo_from_rhou(Q);


        number k2 = evalRate(Q);
        Q.u_modes[0] = _macro_Q0.u_modes[0] + h*(a31*k1 + a32*k2);

        Q.u = _uTotal - Q.u_modes[0];
        _macroAirModel.update_thermo_from_rhou(Q);


        number k3 = evalRate(Q);
        Q.u_modes[0] = _macro_Q0.u_modes[0] + h*(a41*k1 + a42*k2 + a43*k3);

        Q.u = _uTotal - Q.u_modes[0];
        _macroAirModel.update_thermo_from_rhou(Q);


        number k4 = evalRate(Q);
        Q.u_modes[0] = _macro_Q0.u_modes[0] + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4);

        Q.u = _uTotal - Q.u_modes[0];
        _macroAirModel.update_thermo_from_rhou(Q);


        number k5 = evalRate(Q);
        Q.u_modes[0] = _macro_Q0.u_modes[0] + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5);

        Q.u = _uTotal - Q.u_modes[0];
        _macroAirModel.update_thermo_from_rhou(Q);


        number k6 = evalRate(Q);
        Q.u_modes[0] = _macro_Q0.u_modes[0] + h*(b51*k1 + b53*k3 + b54*k4 + b56*k6);

        Q.u = _uTotal - Q.u_modes[0];
        _macroAirModel.update_thermo_from_rhou(Q);


        // Compute error estimate.
        number errEst = Q.u_modes[0] - (_macro_Q0.u_modes[0] + h*(b41*k1 + b43*k3 + b44*k4 + b45*k5 + b46*k6));

        // And use error estimate as a means to suggest a new timestep.
        double atol = _energyAbsTolerance;
        double rtol = _energyRelTolerance;
        double sk = atol + rtol*fmax(fabs(_macro_Q0.u_modes[0].re), fabs(Q.u_modes[0].re));
        double err = fabs(errEst.re/sk);
        // Now use error as an estimate for new step size
        double scale = 0.0;
        const double maxscale = 10.0;
        const double minscale = 0.2;
        const double safe = 0.9;
        const double alpha = 0.2;
        if ( err <= 1.0 ) { // A succesful step...
            if ( err == 0.0 )
                scale = maxscale;
            else {
                scale = safe * pow(err, -alpha);
                if ( scale < minscale )
                    scale = minscale;
                if ( scale > maxscale )
                    scale = maxscale;
            }
            hSuggest = scale*h;
            return ResultOfStep.success;
        }
        // else, failed step
        scale = fmax(safe*pow(err, -alpha), minscale);
        hSuggest = scale*h;
        return ResultOfStep.failure;
    }

    //electronic functions

    @nogc
    void update_input_from_numden(number[] elec_numden, number[] full_numden)
    {
        foreach(int i;0 .. _n_N_species) { //number density in state solver in #/cm^3
            elec_numden[i] = full_numden[i] / 1e6;
        }
        elec_numden[_n_N_species] = full_numden[_n_elec_species + 3] / 1e6;

        foreach(int i;_n_N_species .. _n_elec_species) {
            elec_numden[i+1] = full_numden[i] / 1e6;
        }
        elec_numden[_n_elec_species+1] = full_numden[_n_elec_species + 4] / 1e6;
        elec_numden[$-1] = full_numden[$-1] / 1e6;
    }

    @nogc update_numden_from_output(number[] full_numden, number[] elec_numden)
    { // Remember - N1, N2, ..., O1, O2, ..., N2, O2, NO, N+, O+, N2+, O2+, NO+, e
        foreach (int i;0 .. _n_N_species) { // Do N
            full_numden[i] = elec_numden[i] * 1e6;
        }
        full_numden[_n_elec_species + 3] = elec_numden[_n_N_species] * 1e6; //Do N+

        foreach (int i; _n_N_species .. _n_elec_species) {
            full_numden[i] = elec_numden[i+1] * 1e6; //Do O
        }
        full_numden[_n_elec_species + 4] = elec_numden[_n_elec_species+1] * 1e6; //Do O+
        full_numden[$-1] = elec_numden[$-1] * 1e6;
    }


    @nogc
    number energyInNoneq(ref GasState Q)
    {
        return to!number(0);
    }

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
    void Update_Electronic_State(ref GasState Q, ref GasState macro_state)
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
        foreach (isp; _n_N_species .. _n_elec_species) {
            O_massf_sum += Q.massf[isp];
        }
        if (O_massf_sum == 0.0) {
            Q.massf[_n_N_species] = macro_state.massf[1];
            foreach (isp; _n_N_species + 1 .. _n_elec_species) {
                Q.massf[isp] = 0.0;
            }
        } else {
            number O_massf_factor = macro_state.massf[1]/O_massf_sum;
            foreach (isp; _n_N_species .. _n_elec_species) {
                Q.massf[isp] *= O_massf_factor;
            }
        }

        foreach (isp; _n_elec_species .. Q.massf.length) {
            Q.massf[isp] = macro_state.massf[isp - _n_elec_species + 2];
        }
    }
}


static double[string] A_MW; // A parameter in Millikan-White expression
static this()
{
    // Gnoffo et al. Table 1
    A_MW["N2"] = 220.0;
    A_MW["O2"] = 129.0;
    A_MW["NO"] = 168.0;
    A_MW["N2+"] = 220.0;
    A_MW["O2+"] = 129.0;
    A_MW["NO+"] = 168.0;
}

version(electronically_specific_kinetics_test) {
    int main()
    {
        //writeln("---------Start Test---------");
        import util.msg_service;

        auto L = init_lua_State();
        string filename = "../gas/sample-data/electronic-and-macro-species.lua";
        doLuaFile(L, filename);
        auto gm = new ElectronicallySpecificGas(L);
        auto gd = GasState(gm.n_species,1);
        lua_close(L);
        //writeln("----------Initialised Gas model and state-----------");
        ElectronicallySpecificKinetics esk = new ElectronicallySpecificKinetics("sample-input/ES_files.lua",gm);

        //writeln("--------Initialised kinetics-------------");

        gd.massf[] = 0;
        gd.massf[16] = 0.78;
        gd.massf[16 + 1] = 0.22;

        gd.p = 21478.7;
        gd.T = 25000.9;
        gd.T_modes[0] = 20000.0;
        gm.update_thermo_from_pT(gd);
        // writeln(gd);
        double _dt = 1e-9;
        double dtSuggest = 1e-10;
        number[maxParams] params;

        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        esk(gd, _dt, dtSuggest, params);
        // writeln("The final gas state after the kinetics steps is: ");
        // writeln(gd);

        double massfsum=0.0;
        foreach(number eachmassf;gd.massf) {
            massfsum += eachmassf;
        }
        assert(isClose(massfsum, 1.0, 1e-2), failedUnitTest());



        return 0;
    }
}

//Everything from here is old, will delete when appropriate
/*
final class ElectronicallySpecificKinetics : ThermochemicalReactor {
public:

    this(string ESK_N_Filename, string ESK_O_Filename, GasModel gmodel)
    {
        super(gmodel);
        _numden.length = gmodel.n_species;
        _numden_input.length = gmodel.n_species - 2;
        _numden_output.length = gmodel.n_species - 2;

        full_grouped_data.length = 2;

        foreach(int i;0 .. gmodel.n_species) {
            if (gmodel.species_name(i).length > 2) {
                if (to!(char[])(gmodel.species_name(i))[0..3] == "NI ") {
                    NInum += 1;
                } else if (to!(char[])(gmodel.species_name(i))[0..3] == "OI ") {
                    OInum += 1;
                }
            }
            if (gmodel.species_name(i) == "e") {
                eind = i;
            } else if (gmodel.species_name(i) == "NII") {
                NIIind = i;
            } else if (gmodel.species_name(i) == "OII") {
                OIIind = i;
            } else if (gmodel.species_name(i) == "N2") {
                N2ind = i;
            } else if (gmodel.species_name(i) == "O2") {
                O2ind=i;
            }
        }


        //Need to construct the energy data table and the reaction rate parameters from file.
        PopulateRateFits(ESK_N_Filename,ESK_O_Filename);
        kinetics.electronic_state_solver.Init(full_rate_fit, [NInum,OInum]);
    }

    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {

        //  Process for kinetics update is as follows
        //  1. Calculate initial energy in higher electronic states
        //  2. Pass mass fraction to number density vector and convert to #/cm^3 for the solver
        //  3. Solve for the updated electronic states over a time of tInterval (use dtSuggest in solver)
        //  4. Convert number density vector to mol/cm^3 for the chemistry dissociation
        //  5. Call molecule update with state Q --> updates numden vector
        //  6. Convert number density vector to #/m^3 and pass back to Q.massf
        //  8. Calculate final energy in higher electronic states --> update u_modes[0] --> update thermo

        // 1.
        initialEnergy = energyInNoneq(Q);

        // 2.
        foreach(int i; 0 .. _gmodel.n_species){ //give massf values for the required species
            _numden[i] = Q.massf[i];
        }
        _gmodel.massf2numden(Q, _numden); //Convert from mass fraction to number density

        foreach(int i;0 .. _gmodel.n_species - 2) { //number density in state solver in #/cm^3
            _numden_input[i] = _numden[i] / 1e6;
        }

        // 3.
        Electronic_Solve(_numden_input, _numden_output, Q.T_modes[0], tInterval, dtSuggest);
        // 4.
        foreach (int i; 0 .. _gmodel.n_species - 2) {//convert back to number density in #/m^3
            _numden[i] = _numden_output[i] * 1e6;
        }

        foreach(int i; 0 .. _gmodel.n_species){ //convert numden to mol/cm^3
            _numden[i] = _numden[i] / (Avogadro_number*1e6);
        }

        // 5.
        Molecule_Update(Q, tInterval);
        // 6.
        foreach(int i; 0 .. _gmodel.n_species){ //convert mol/cm^3 to numden
            _numden[i] = _numden[i] * Avogadro_number*1e6;
        }
        _gmodel.numden2massf(_numden,Q);
        // 8.
        finalEnergy = energyInNoneq(Q);
        Q.u -= finalEnergy-initialEnergy;
        _gmodel.update_thermo_from_rhou(Q);
    }

private:
    number[] _numden; //Total set of species including static N2 and O2
    number[] _numden_input; //Input for the electronic state solver, exclusive of N2 and O2
    number[] _numden_output; //Output for the state solver
    number initialEnergy;
    number finalEnergy;
    number N_sum;
    number O_sum;
    number N2_step_change;
    number O2_step_change;
    double[8][47][47][2] full_rate_fit;
    double[][][] full_grouped_data;

    //physical cconstants
    double _pi = 3.14159265359;
    double _me = 9.10938356e-28; //electron mass in g
    double _kb = 1.3807e-16; //Boltzmann constant in cm^2 g s^-1 K^-1
    double _e = 4.8032e-10; // electron charge, cm^(3/2) g s^-2 K^-1

    //define number of states
    int NInum;
    int OInum;
    int NIIind;
    int OIIind;
    int eind;
    int N2ind;
    int O2ind;


    //Arhenius coefficients for N2 and O2 dissociation
    //In order: N2, N, O2, O
    //within each group: scale, A, n, E
    double[][] ArrN2_fr = [[2.5, 1.920e+17, -0.50, 113100.00],
                            [1.0, 4.150e+22, -1.50, 113100.00],
                            [1.0, 1.920e+17, -0.50, 113100.00],
                            [1.0, 1.920e+17, -0.50, 113100.00]];

    double[][] ArrN2_br = [[2.5, 1.090e+16, -0.50, 0.0],
                            [1.0, 2.320e+21, -1.50, 0.0],
                            [1.0, 1.090e+16, -0.50, 0.0],
                            [1.0, 1.090e+16, -0.50, 0.0]];

    double[][] ArrO2_fr = [[2.0, 3.610e+18, -1.00, 59400.00],
                            [1.0, 3.610e+18, -1.00, 59400.00],
                            [9.0, 3.610e+18, -1.00, 59400.00],
                            [25.0, 3.610e+18, -1.00, 59400.00]];

    double[][] ArrO2_br = [[2.0, 3.010e+15, -0.50, 0.0],
                            [1.0, 3.010e+15, -0.50, 0.0],
                            [9.0, 3.010e+15, -0.50, 0.0],
                            [25.0, 3.010e+15, -0.50, 0.0]];

    @nogc
    number N2_fr(ref GasState Q) //must only be called when in conc mol/cm^3, not massf
    {
        number rate_coef(int i){
            return ArrN2_fr[i][0]*ArrN2_fr[i][1]*(Q.T^^ArrN2_fr[i][2])*exp(-ArrN2_fr[i][3]/Q.T);
        }

        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. NInum) {
            N_sum += _numden[i];
        }
        foreach(int i; NInum+1 .. NInum+1+OInum) {
            O_sum += _numden[i];
        }

        return rate_coef(0)*_numden[N2ind]*_numden[N2ind] + rate_coef(1)*_numden[N2ind]*N_sum + 
                rate_coef(2)*_numden[N2ind]*_numden[O2ind] + rate_coef(3)*_numden[N2ind]*O_sum;
    }
    @nogc
    number N2_br(ref GasState Q) //must only be called when in conc mol/cm^3, not massff
    {
        number rate_coef(int i){
            return ArrN2_br[i][0]*ArrN2_br[i][1]*(Q.T^^ArrN2_br[i][2])*exp(-ArrN2_br[i][3]/Q.T);
        }

        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. NInum) {
            N_sum += _numden[i];
        }
        foreach(int i; NInum+1 .. NInum+1+OInum) {
            O_sum += _numden[i];
        }

        return rate_coef(0)*N_sum*N_sum*_numden[N2ind] + rate_coef(1)*N_sum*N_sum*N_sum + 
                rate_coef(2)*N_sum*N_sum*_numden[O2ind] + rate_coef(3)*N_sum*N_sum*O_sum;
    }
    @nogc
    number O2_fr(ref GasState Q) //must only be called when in conc mol/cm^3, not massf
    {
        number rate_coef(int i){
            return ArrO2_fr[i][0]*ArrO2_fr[i][1]*(Q.T^^ArrO2_fr[i][2])*exp(-ArrO2_fr[i][3]/Q.T);
        }

        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. NInum) {
            N_sum += _numden[i];
        }
        foreach(int i; NInum+1 .. NInum+1+OInum) {
            O_sum += _numden[i];
        }

        return rate_coef(0)*_numden[O2ind]*_numden[N2ind] + rate_coef(1)*_numden[O2ind]*N_sum +
                rate_coef(2)*_numden[O2ind]*_numden[O2ind] + rate_coef(3)*_numden[O2ind]*O_sum;
    }
    @nogc
    number O2_br(ref GasState Q) //must only be called when in conc mol/cm^3, not massf
    {
        number rate_coef(int i){
            return ArrO2_br[i][0]*ArrO2_br[i][1]*(Q.T^^ArrO2_br[i][2])*exp(-ArrO2_br[i][3]/Q.T);
        }

        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. NInum) {
            N_sum += _numden[i];
        }
        foreach(int i; NInum+1 .. NInum+1+OInum) {
            O_sum += _numden[i];
        }

        return rate_coef(0)*O_sum*O_sum*_numden[N2ind] + rate_coef(1)*O_sum*O_sum*N_sum +
                rate_coef(2)*O_sum*O_sum*_numden[O2ind] + rate_coef(3)*O_sum*O_sum*O_sum;
    }
    @nogc
    number N2_rate(ref GasState Q)
    {
        return N2_fr(Q) - N2_br(Q);
    }
    @nogc
    number O2_rate(ref GasState Q)
    {
        return O2_fr(Q) - O2_br(Q);
    }
    @nogc
    void Molecule_Update(ref GasState Q, double tInterval)
    {
        double _dt = tInterval/100.0;
        foreach(int i; 0 .. 100){
            N2_step_change = N2_rate(Q);
            O2_step_change = O2_rate(Q);
            _numden[N2ind] += -_dt*N2_step_change;
            _numden[O2ind] += -_dt*O2_step_change;
            _numden[0] += _dt*2*N2_step_change;
            _numden[NInum+1] += _dt*2*O2_step_change;
        }
    }

    @nogc number energyInNoneq(ref const(GasState) Q)
    {
        number uNoneq = 0.0;
        foreach (int isp; 0 .. _gmodel.n_species) {
            uNoneq += Q.massf[isp] * _gmodel.electronic_energy[isp];
        }
        return uNoneq;
    }

    @nogc
    double[][] Import_2D(string filename) {
        debug{
            double[][] output_data;
            if (exists(filename)) { //check for existance of file
                File file=File(filename, "r");
                while (!file.eof()) {
                    output_data~=to!(double[])(split(strip(file.readln()),","));
                }
                if (output_data[$-1].length==0) { //accounts for the sometimes blank line at the end of csv files
                    output_data = output_data[0..$-1];
                }
            } else { //if no filename exists
                writeln("no such filename: ",filename);
            }
        return output_data;
        } else {
            throw new Error("Not implemented for nondebug build.");
        }
    }

    @nogc
    void PopulateRateFits(string Nfilename, string Ofilename)
    {
        debug{
            double[][] N_rate_fit = Import_2D(Nfilename);
            double[][] O_rate_fit = Import_2D(Ofilename);
            foreach (double[] row;N_rate_fit){
                full_rate_fit[0][to!int(row[0])][to!int(row[1])] = row[2 .. $];
            }
            foreach (double[] row;O_rate_fit){
                full_rate_fit[1][to!int(row[0])][to!int(row[1])] = row[2 .. $];
            }
        }
    }

}

version(electronically_specific_kinetics_test) {
    int main()
    {
        return 0;
    }
}
*/
/*
version(electronically_specific_kinetics_test)
{
    int main()
    {
        import util.msg_service;

        auto L = init_lua_State();
        string filename = "../gas/sample-data/electronic_composition.lua";
        doLuaFile(L, filename);
        auto gm = new ElectronicallySpecificGas(L);
        auto gd = GasState(gm.n_species,1);

        number initial_uNoneq=0.0;

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

        gd.p = 64502.1;
        gd.T = 54826.4;
        gd.T_modes[0] = 15000.0;
        gm.update_thermo_from_pT(gd);

        lua_close(L);

        double _dt = 1e-9;
        double _duration = 1e-7;
        double _t = 0.0;

        ElectronicallySpecificKinetics esk = new ElectronicallySpecificKinetics("sample-input/ESK-N.txt","sample-input/ESK-O.txt",gm);

        double dtSuggest = 1e-9;
        number[maxParams] params;
        while (_t < _duration+_dt) {
            esk(gd, _dt, dtSuggest, params);
            _t+=_dt;
        }
        double massfsum=0.0;
        foreach(number eachmassf;gd.massf) {
            massfsum += eachmassf;
        }
        assert(isClose(massfsum, 1.0, 1e-2), failedUnitTest());

        return 0;

    }
}

//arrhenius mol/cm^3


*/
