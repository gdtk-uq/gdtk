/**
 * Two-temperature air kinetics module.
 *
 * This module implements the two-temperature version
 * of the high-temperature air model described in
 * the report by Gupta et al.
 *
 * Author: Rowan G.
 * Major Edits by N. Gibbons (21/03/19)
 */

module kinetics.two_temperature_air_kinetics;

import std.stdio;
import std.conv;
import std.string;
import std.math;
import std.algorithm;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas;
import gas.two_temperature_air;
import kinetics.thermochemical_reactor;
import kinetics.chemistry_update;
import kinetics.reaction_mechanism;


immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
                                             // Yes, you read that right. Sometimes a step is succesful
                                             // but the timestep selection algorithm will suggest
                                             // a reduction. We limit that reduction to no more than 50%.
immutable double DT_REDUCTION_FACTOR = 10.0; // factor by which to reduce timestep
                                             // after a failed attempt
immutable double H_MIN = 1.0e-15; // Minimum allowable step size

final class TwoTemperatureAirKinetics : ThermochemicalReactor {
    ReactionMechanism rmech;
    ChemODEStep cstep;
    bool tightTempCoupling;
    int maxSubcycles;
    int maxAttempts;

    this(string chemFile, string energyExchFile, GasModel gmodel)
    {
        super(gmodel);
        _airModel = cast(TwoTemperatureAir) gmodel;
        if (_airModel is null) {
            string errMsg = "Error in construction of TwoTemperatureAirKinetics.\n";
            errMsg ~= "The supplied gas model must be a TwoTemperatureAir model.\n";
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        _Qinit = GasState(gmodel);
        _Q0 = GasState(gmodel);
        chemupdate_constructor(chemFile, gmodel);
        _molef.length = _airModel.n_species;
        _numden.length = _airModel.n_species;
        initModel(energyExchFile);
    }

    void chemupdate_constructor(string fname, GasModel gmodel)
    {
        // Allocate memory
        //_Qinit = GasState(gmodel.n_species, gmodel.n_modes);
        _conc0.length = gmodel.n_species;
        _concOut.length = gmodel.n_species;
        _conc_for_source_terms.length = gmodel.n_species + gmodel.n_modes;
        _rates_for_source_terms.length = gmodel.n_species + gmodel.n_modes;

        // Configure other parameters via Lua state.
        auto L = init_lua_State();
        doLuaFile(L, fname);

        // Check on species order before proceeding.
        lua_getglobal(L, "species");
        if (lua_isnil(L, -1)) {
            string errMsg = format("There is no species listing in your chemistry input file: '%s'\n", fname);
            errMsg ~= "Please be sure to use an up-to-date version of prep-chem to generate your chemistry file.\n";
            throw new Error(errMsg);
        }
        foreach (isp; 0 .. gmodel.n_species) {
            lua_rawgeti(L, -1, isp);
            auto sp = to!string(luaL_checkstring(L, -1));
            if (sp != gmodel.species_name(isp)) {
                string errMsg = "Species order is incompatible between gas model and chemistry input.\n";
                errMsg ~= format("Chemistry input file is: '%s'.\n", fname);
                errMsg ~= format("In gas model: index %d ==> species '%s'\n", isp, gmodel.species_name(isp));
                errMsg ~= format("In chemistry: index %d ==> species '%s'\n", isp, sp);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);

        lua_getglobal(L, "config");
        lua_getfield(L, -1, "tempLimits");
        lua_getfield(L, -1, "lower");
        double T_lower_limit = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, -1, "upper");
        double T_upper_limit = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_pop(L, 1);
        lua_pop(L, 1);

        lua_getglobal(L, "reaction");
        rmech = createReactionMechanism(L, gmodel, T_lower_limit, T_upper_limit);
        lua_pop(L, 1);

        lua_getglobal(L, "config");
        lua_getfield(L, -1, "odeStep");
        lua_getfield(L, -1, "method");
        string ode_method = to!string(lua_tostring(L, -1));
        lua_pop(L, 1);
        switch (ode_method) {
        case "rkf":
             lua_getfield(L, -1, "errTol");
             double err_tol = lua_tonumber(L, -1);
             lua_pop(L, 1);
             cstep = new RKFStep(gmodel, rmech, err_tol);
             break;
        case "alpha-qss":
             lua_getfield(L, -1, "eps1");
             double eps1 = lua_tonumber(L, -1);
             lua_pop(L, 1);
             lua_getfield(L, -1, "eps2");
             double eps2 = lua_tonumber(L, -1);
             lua_pop(L, 1);
             lua_getfield(L, -1, "delta");
             double delta = lua_tonumber(L, -1);
             lua_pop(L, 1);
             lua_getfield(L, -1, "maxIters");
             int max_iters = to!int(lua_tointeger(L, -1));
             lua_pop(L, 1);
             cstep = new AlphaQssStep(gmodel, rmech, eps1, eps2, delta, max_iters);
             break;
        default:
             string errMsg = format("ERROR: The odeStep '%s' is unknown.\n", ode_method);
             throw new Error(errMsg);
        }
        lua_pop(L, 1); // pops 'odeStep'

        lua_getfield(L, -1, "tightTempCoupling");
        tightTempCoupling = lua_toboolean(L, -1);
        lua_pop(L, 1);

        lua_getfield(L, -1, "maxSubcycles");
        maxSubcycles = to!int(luaL_checkinteger(L, -1));
        lua_pop(L, 1);

        lua_getfield(L, -1, "maxAttempts");
        maxAttempts = to!int(luaL_checkinteger(L, -1));
        lua_pop(L, 1);

        lua_pop(L, 1); // pops 'config'
        lua_close(L);
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        _Qinit.copy_values_from(Q);
        _gmodel.massf2conc(Q, _conc0);
        number uTotal = _gmodel.internal_energy(Q);
        _uTotal = Q.u + Q.u_modes[0];

        // 0. Evaluate the rate constants.
        //    It helps to have these computed before doing other setup work.
        rmech.eval_rate_constants(_gmodel, Q);
        // 1. Sort out the time step for possible subcycling.
        double t = 0.0;
        double h;
        double dtSave;
        if ( dtSuggest > tInterval )
            h = tInterval;
        else if ( dtSuggest <= 0.0 )
            h = rmech.estimateStepSize(_conc0);
        else
            h = dtSuggest;

        // 2. Now do the interesting stuff, increment species change
        int cycle = 0;
        int attempt = 0;
        for ( ; cycle < maxSubcycles; ++cycle ) {
            /* Copy the last good timestep suggestion before moving on.
             * If we're near the end of the interval, we want
             * the penultimate step. In other words, we don't
             * want to store the fractional step of (tInterval - t)
             * that is taken as the last step.
             */
            dtSave = dtSuggest;
            h = min(h, tInterval - t);
            attempt= 0;
            for ( ; attempt < maxAttempts; ++attempt ) {
                ResultOfStep result1 = cstep(_conc0, h, _concOut, dtSuggest);
                ResultOfStep result2 = EnergyStep(Q, h, dtSuggest);
                ResultOfStep result = ResultOfStep.failure;
                if (result1==ResultOfStep.success && result2==ResultOfStep.success){
                    result = ResultOfStep.success;
                }

                bool passesMassFractionTest = true;
                if ( result  == ResultOfStep.success ) {
                    /* We succesfully took a step of size h according to the ODE method.
                     * However, we need to test that that mass fractions have remained ok.
                     */
                    _gmodel.conc2massf(_concOut, Q);
                    auto massfTotal = sum(Q.massf);
                    if ( fabs(massfTotal - 1.0) > ALLOWABLE_MASSF_ERROR ) {
                        passesMassFractionTest = false;
                    }
                }

                if ( result == ResultOfStep.success && passesMassFractionTest ) {
                    t += h;
                    // Copy succesful values in concOut to conc0, ready for next step
                    _conc0[] = _concOut[];
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
                        h = min(dtSuggest, hMax);
                    }
                    break;
                }
                else { // in the case of failure...
                    /* We now need to make some decision about
                     * what timestep to attempt next. We follow
                     * David Mott's suggestion in his thesis (on p. 51)
                     * and reduce the timestep by a factor of 2 or 3.
                     * (The actual value is set as DT_REDUCTION_FACTOR).
                     * In fact, for the types of problems we deal with
                     * I have found that a reduction by a factor of 10
                     * is more effective.
                     */
                    h /= DT_REDUCTION_FACTOR;
                    if ( h < H_MIN ) {
                        string errMsg = "Hit the minimum allowable timestep in chemistry update.";
                        debug { errMsg ~= format("\ndt= %.4e", H_MIN); }
                        Q.copy_values_from(_Qinit);
                        throw new ThermochemicalReactorUpdateException(errMsg);
                    }
                }
            } // end attempts at single subcycle.
            if ( attempt == maxAttempts ) {
                string errMsg = "Hit maximum number of step attempts within a subcycle for chemistry update.";
                // We did poorly. Let's put the original GasState back in place,
                // and let the outside world know via an Exception.
                Q.copy_values_from(_Qinit);
                throw new ThermochemicalReactorUpdateException(errMsg);
            }
            /* Otherwise, we've done well.
             * If tight temperature coupling is requested, we can reevaluate
             * the temperature at this point. With regards to tight coupling,
             * we follow Oran and Boris's advice on pp. 140-141.
             * To paraphrase: solving a separate differential equation for
             * temperature is computationally expensive, however, it usually
             * suffices to reevaluate the temperature assuming that total internal
             * energy of the system has not changed but has been redistributed
             * among chemical states. Since the chemistry has not altered much,
             * the iteration for temperature should converge in one or two
             * iterations.
             *
             * My own additional argument for avoiding a temperature differential
             * equation is that it does not play nicely with the special
             * ODE methods for chemistry that exploit structure in the species
             * differential equations.
             */
            if ( tightTempCoupling ) {
                _gmodel.conc2massf(_conc0, Q);
                Q.u = _uTotal - Q.u_modes[0];
                _gmodel.update_thermo_from_rhou(Q);
                rmech.eval_rate_constants(_gmodel, Q);
            }

            if ( t >= tInterval ) { // We've done enough work.
                // If we've only take one cycle, then we would like to use
                // dtSuggest rather than using the dtSuggest from the
                // penultimate step.
                if (cycle == 0) {
                    dtSave = dtSuggest;
                }
                break;
            }
        }
        if ( cycle == maxSubcycles ) {
            // If we make it out here, then we didn't complete within
            // the maximum number of subscyles... we are taking too long.
            // Let's return the gas state to how it was before we failed
            // and throw an Exception.
            string errMsg = "Hit maximum number of subcycles while attempting chemistry update.";
            Q.copy_values_from(_Qinit);
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        // At this point, it appears that everything has gone well.
        // We'll tweak the mass fractions in case they are a little off.
        _gmodel.conc2massf(_concOut, Q);
        auto massfTotal = sum(Q.massf);
        foreach (ref mf; Q.massf) mf /= massfTotal;
        /*
         * RJG 2018-11-16
         * Let's remove this piece of complication.
         * Let's have the chemistry module *only* update mass fractions.
         * The caller can decide on the thermodynamic constraint.
         * This allows re-use as a chemistry updater in multi-temperature models.
         *
        if (_gmodel.n_modes >= 1) {
            // Changing mass fractions does not change internal energies
            // but it might have been redistributed.
            try {
                _gmodel.update_thermo_from_rhoT(Q);
            }
            catch (GasModelException err) {
                string msg = "The call to update_thermo_from_rhoT in the chemistry update failed.";
                debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
            }
            auto uOther = sum(Q.u_modes);
            Q.u = uTotal - uOther;
            try {
                _gmodel.update_thermo_from_rhou(Q);
            }
            catch (GasModelException err) {
                string msg = "The call to update_thermo_from_rhou in the chemistry update failed.";
                debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
            }
        }
        */
        dtSuggest = dtSave;
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for two_temperature_air_kinetics.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    // Stuff from chemistry update
    //GasState _Qinit;
    number[] _conc0;
    number[] _concOut;
    number[] _conc_for_source_terms;
    number[] _rates_for_source_terms;

    TwoTemperatureAir _airModel;
    GasState _Qinit, _Q0;
    double _c2 = 1.0;
    number _uTotal;
    number _T_sh, _Tv_sh;
    double[] _A;
    double[] _as, _bs, _cs;
    number[] _molef;
    number[] _numden;
    double[] _particleMass;
    double[][] _mu;
    int[][] _reactionsBySpecies;
    int[] _neutralidxs;
    int[] _ionidxs;
    int _electronidx;
    bool do_ET_exchange;
    // Numerics control
    int _maxSubcycles = 10000;
    int _maxAttempts = 3;
    double _energyAbsTolerance = 1.0e-09;
    double _energyRelTolerance = 1.0e-09;

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
                _energyAbsTolerance = to!int(luaL_checknumber(L, -1));
            }
            lua_pop(L, 1);

            lua_getfield(L, -1, "relTolerance");
            if (!lua_isnil(L, -1)) {
                _energyRelTolerance = to!int(luaL_checknumber(L, -1));
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);

        lua_close(L);

        int nSpecies = _airModel.n_species;
        _A.length = nSpecies;
        _A[] = 0.0;
        foreach (isp; _airModel.molecularSpecies) {
            _A[isp] = A_MW[_airModel.species_name(isp)];
        }
        foreach(isp; 0 .. _airModel.n_species){
            if (_airModel.charge[isp]>0.0) _ionidxs ~= isp;
            if (_airModel.charge[isp]<0.0) _electronidx = isp;
            if (_airModel.charge[isp]==0.0){
                string name = _airModel.species_name(isp);
                _as ~= as[name];
                _bs ~= bs[name];
                _cs ~= cs[name];
                _neutralidxs ~= isp;
            }
        }

        _mu.length = nSpecies;
        foreach (isp; 0 .. nSpecies) _mu[isp].length = nSpecies;
        foreach (isp; 0 .. nSpecies) {
            foreach (jsp; 0 .. nSpecies) {
                double M_isp = _airModel.mol_masses[isp];
                double M_jsp = _airModel.mol_masses[jsp];
                _mu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                _mu[isp][jsp] *= 1000.0; // convert kg/mole to g/mole
            }
        }
        if (_airModel.is_plasma) do_ET_exchange = true;

        _particleMass.length = nSpecies;
        foreach (isp; 0 .. nSpecies) {
            _particleMass[isp] = _airModel.mol_masses[isp]/Avogadro_number;
        }

        _reactionsBySpecies.length = nSpecies;
        foreach (isp; 0 .. nSpecies) {
            foreach (ir; 0 .. to!int(rmech.n_reactions)) {
                if (rmech.reactionHasParticipant(ir, isp)) {
                    _reactionsBySpecies[isp] ~= ir;
                }
            }
        }
    }

    @nogc
    number evalRelaxationTime(ref GasState Q, int isp)
    {

        number totalND = 0.0;
        number sum = 0.0;
        foreach (csp; 0 .. _airModel.n_species) {
            // We exclude electrons in the vibrational relaxation time.
            if (_airModel.species_name(csp) == "e-") {
                continue;
            }
            // 1. Compute Millikan-White value for each interaction
            if (_molef[csp] >= SMALL_MOLE_FRACTION) {
                number nd = _numden[csp];
                sum += nd * exp(_A[isp] * (pow(Q.T, -1./3) - 0.015*pow(_mu[isp][csp], 0.25)) - 18.42);
                totalND += nd;
            }
        }
        number tauMW = sum/totalND; // <-- p*tauMW (in atm.s)
        tauMW *= P_atm/Q.p;

        // 2. Compute Park value for high-temperature correction
        number nd = _numden[isp];
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
        _gmodel.massf2molef(Q, _molef);
        _gmodel.massf2numden(Q, _numden);
        number rate = 0.0;
        foreach (isp; _airModel.molecularSpecies) {
            // Vibrational energy exchange via collisions.
            number tau = evalRelaxationTime(Q, isp);
            number evStar = _airModel.vibElecEnergy(Q.T, isp);
            number ev = _airModel.vibElecEnergy(Q.T_modes[0], isp);
            rate += Q.massf[isp] * (evStar - ev)/tau;
        }
        rate += EnergyReactiveSourceTerm(Q);
        if(do_ET_exchange)
           rate += ElectronEnergyExchangeRate(Q);
        return rate;
    }

    @nogc
    number ElectronEnergyExchangeRate(ref GasState Q)
    {
        /*
        Appleton-Bray electron-translational energy exchange expression. Taken from Gnoffo, 1989
        equation (16) term 7, and (64)/(65) for the collision cross section data.

        The expression differs slightly from the paper, because we are using SI units. Their expression
        uses "esu" for the electric charge, which is defined so that the factor 1/4/pi/epsilon_0 is equal
        to one. See notes 21/03/24 for the derivation of the expression with this term inserted.

        @author: Nick Gibbons
        */
        number Te = Q.T_modes[0];
        number Tfit = fmin(30e3, Te);
        number me = _particleMass[_electronidx];
        number ne = _numden[_electronidx];
        number nues_on_Ms = 0.0;
        double pi = to!double(PI);

        // Neutral energy exchange from tabulated collision cross sections eqn. (65)
        foreach(isp; _neutralidxs){
            number sigma = _as[isp] + _bs[isp]*Tfit + _cs[isp]*Tfit*Tfit;
            number nues = _numden[isp]*sigma*sqrt(8*Boltzmann_constant*Te/pi/me);
            nues_on_Ms += nues/_airModel.mol_masses[isp];
        }

        // Ion energy exchange from analytical Coulomb collision cross section eqn. (64)
        // Note that Le is a weird inverse lengthscale related to the Debeye length.
        number Le = 4.0*PI*vacuum_permittivity*Boltzmann_constant*Te/electron_volt_energy/electron_volt_energy;
        number Lambda = log(Le*Le*Le/pi/fmax(ne,1.0));
        number A = 8.0/3.0*sqrt(pi/me/8.0)*electron_volt_energy/sqrt(4.0*pi*vacuum_permittivity);
        number B = sqrt(1.0/Le/Le/Le);

        foreach(isp; _ionidxs){
            number nues = A*_numden[isp]*B*Lambda;
            nues_on_Ms += nues/_airModel.mol_masses[isp];
        }

        // Unlike Gnoffo's equation 16, we are working directly on u_ve in J/kg,
        // rather than rho u_ve in J/m3. So we use massf[e-] instead of rho[e-]
        return Q.massf[_electronidx]*3.0*R_universal*(Q.T-Te)*nues_on_Ms;
    }

    @nogc
    number EnergyReactiveSourceTerm(ref GasState Q)
    {
        /*
        Since the different species have different amounts of vibronic energy, chemical reactions
        actually contribute to the ev source terms, albeit not by much.

        See derivation from 21/04/02
        @author: Nick Gibbons
        */

        number rate = 0.0;
        foreach (isp; 0 .. _airModel.n_species) {
            // Vibrational energy per kilogram of species isp
            number ev = _airModel.vibElecEnergy(Q.T_modes[0], isp);

            // Vibrational energy change due to chemical reactions.
            number chemRate = 0.0;
            foreach (ir; _reactionsBySpecies[isp]) {
                chemRate += rmech.rate(ir, isp);
            }
            chemRate *= _airModel.mol_masses[isp]; // convert mol/m^3/s --> kg/m^3/s
            rate += chemRate*ev;
        }
        rate /= Q.rho; // convert J/m3/s to J/kg(of mixture)/s

        return rate;
    }

    @nogc
    ResultOfStep EnergyStep(ref GasState Q, double h, ref double hSuggest)
    {
        _Q0.copy_values_from(Q);
        immutable double a21=1./5., a31=3./40., a32=9./40.,
            a41=3./10., a42=-9./10., a43 = 6./5.,
            a51=-11./54., a52=5./2., a53=-70./27., a54=35./27.,
            a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096.;
        immutable double b51=37./378., b53=250./621., b54=125./594., b56=512./1771.,
            b41=2825./27648., b43=18575./48384., b44=13525./55296., b45=277./14336., b46=1.0/4.0;

        number k1 = evalRate(Q);
        Q.u_modes[0] = _Q0.u_modes[0] + h*(a21*k1);

        Q.u = _uTotal - Q.u_modes[0];
        _airModel.update_thermo_from_rhou(Q);


        number k2 = evalRate(Q);
        Q.u_modes[0] = _Q0.u_modes[0] + h*(a31*k1 + a32*k2);

        Q.u = _uTotal - Q.u_modes[0];
        _airModel.update_thermo_from_rhou(Q);


        number k3 = evalRate(Q);
        Q.u_modes[0] = _Q0.u_modes[0] + h*(a41*k1 + a42*k2 + a43*k3);

        Q.u = _uTotal - Q.u_modes[0];
        _airModel.update_thermo_from_rhou(Q);


        number k4 = evalRate(Q);
        Q.u_modes[0] = _Q0.u_modes[0] + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4);

        Q.u = _uTotal - Q.u_modes[0];
        _airModel.update_thermo_from_rhou(Q);


        number k5 = evalRate(Q);
        Q.u_modes[0] = _Q0.u_modes[0] + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5);

        Q.u = _uTotal - Q.u_modes[0];
        _airModel.update_thermo_from_rhou(Q);


        number k6 = evalRate(Q);
        Q.u_modes[0] = _Q0.u_modes[0] + h*(b51*k1 + b53*k3 + b54*k4 + b56*k6);

        Q.u = _uTotal - Q.u_modes[0];
        _airModel.update_thermo_from_rhou(Q);


        // Compute error estimate.
        number errEst = Q.u_modes[0] - (_Q0.u_modes[0] + h*(b41*k1 + b43*k3 + b44*k4 + b45*k5 + b46*k6));

        // And use error estimate as a means to suggest a new timestep.
        double atol = _energyAbsTolerance;
        double rtol = _energyRelTolerance;
        double sk = atol + rtol*fmax(fabs(_Q0.u_modes[0].re), fabs(Q.u_modes[0].re));
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

}

static double[string] A_MW; // A parameter in Millikan-White expression
static double[string] as; // a parameter in Electron-Neutral Collision Cross Section Curve Fit
static double[string] bs; // b parameter in Electron-Neutral Collision Cross Section Curve Fit
static double[string] cs; // c parameter in Electron-Neutral Collision Cross Section Curve Fit

static this()
{
    // Gnoffo et al. Table 1
    A_MW["N2"] = 220.0;
    A_MW["O2"] = 129.0;
    A_MW["NO"] = 168.0;
    A_MW["N2+"] = 220.0;
    A_MW["O2+"] = 129.0;
    A_MW["NO+"] = 168.0;

    // Gnoffo et al. Table 5 (These appear to be in m2)
    as["N"]  = 5.0e-20; bs["N"]  = 0.0e+00; cs["N"]  =  0.0e+00;
    as["O"]  = 1.2e-20; bs["O"]  = 1.7e-24; cs["O"]  = -2.0e-29;
    as["N2"] = 7.5e-20; bs["N2"] = 5.5e-24; cs["N2"] = -1.0e-28;
    as["O2"] = 2.0e-20; bs["O2"] = 6.0e-24; cs["O2"] =  0.0e-00;
    as["NO"] = 1.0e-19; bs["NO"] = 0.0e+00; cs["NO"] =  0.0e+00;
}
