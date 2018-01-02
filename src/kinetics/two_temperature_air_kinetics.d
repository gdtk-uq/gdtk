/**
 * Two-temperature air kinetics module.
 *
 * This module implements the two-temperature version
 * of the high-temperature air model described in
 * the report by Gupta et al. 
 *
 * Author: Rowan G.
 */

module kinetics.two_temperature_air_kinetics;

import std.stdio;
import std.conv;
import std.string;
import std.math;
import util.lua;
import util.lua_service;
import gas;
import gas.two_temperature_air;
import kinetics.thermochemical_reactor;
import kinetics.chemistry_update;


immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
                                             // Yes, you read that right. Sometimes a step is succesful
                                             // but the timestep selection algorithm will suggest
                                             // a reduction. We limit that reduction to no more than 50%.
immutable double DT_REDUCTION_FACTOR = 10.0; // factor by which to reduce timestep
                                             // after a failed attempt
immutable double H_MIN = 1.0e-15; // Minimum allowable step size
enum ResultOfStep { success, failure };

final class TwoTemperatureAirKinetics : ThermochemicalReactor {
    this(string chemFile, string energyExchFile, GasModel gmodel)
    {
	super(gmodel);
	_airModel = cast(TwoTemperatureAir) gmodel;
	if (_airModel is null) {
	    string errMsg = "Error in construction of TwoTemperatureAirKinetics.\n";
	    errMsg ~= "The supplied gas model must be a TwoTemperatureAir model.\n";
	    throw new ThermochemicalReactorUpdateException(errMsg);
	}
	_Qinit = new GasState(gmodel);
	_Q0 = new GasState(gmodel);
	_chemUpdate = new ChemistryUpdate(chemFile, gmodel);
	_molef.length = _airModel.n_species;
	initModel(energyExchFile);
    }

    override void opCall(GasState Q, double tInterval,
			 ref double dtChemSuggest, ref double dtThermSuggest,
			 ref double[] params)
    {
	double dummyDouble;
	// 1. Perform chemistry update.
	_chemUpdate(Q, tInterval, dtChemSuggest, dummyDouble, params);

	// 2. Perform energy exchange update.
	// As prep work, compute mole fractions.
	// These values won't change during the energy update
	// since the composition does not change
	// during energy exhange.
	massf2molef(Q.massf, _airModel.mol_masses, _molef);
	energyUpdate(Q, tInterval, dtThermSuggest);
    }

private:
    TwoTemperatureAir _airModel;
    GasState _Qinit, _Q0;
    ChemistryUpdate _chemUpdate;
    double _eTotal;
    double _c2 = 1.0;
    double _T_sh, _Tv_sh;
    double[] _A;
    double[] _molef;
    double[] _particleMass;
    double[][] _mu;
    int[][] _reactionsByMolecule;
    // Numerics control
    int _maxSubcycles = 1000;
    int _maxAttempts = 3;
    double _energyAbsTolerance = 1.0e-2;
    double _energyRelTolerance = 1.0e-3;

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

	_particleMass.length = nSpecies;
	foreach (isp; 0 .. nSpecies) {
	    _particleMass[isp] = _airModel.mol_masses[isp]/Avogadro_number;
	}

	_reactionsByMolecule.length = nSpecies;
	foreach (isp; _airModel.molecularSpecies) {
	    foreach (ir; 0 .. to!int(_chemUpdate.rmech.n_reactions)) {
		if (_chemUpdate.rmech.reactionHasParticipant(ir, isp)) {
		    _reactionsByMolecule[isp] ~= ir;
		}
	    }
	}
    }

    void energyUpdate(GasState Q, double tInterval, ref double dtSuggest)
    {
	_eTotal = Q.u + Q.u_modes[0];
	// We borrow the algorithm from ChemistryUpdate.opCall()
	// Take a copy of what's passed in, in case something goes wrong
	_Qinit.copy_values_from(Q);

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
			string errMsg = format("Hit the minimum allowable timestep in 2-T air energy exchange update: dt= %.4e", H_MIN);
			Q.copy_values_from(_Qinit);
			throw new ThermochemicalReactorUpdateException(errMsg);
		    }
		}
	    } // end attempts at single subcycle.
	    if (attempt == _maxAttempts) {
		string errMsg = "Hit maximum number of step attempts within a subcycle for 2-T air energy exchange update.";
		// We did poorly. Let's put the original GasState back in place,
		// and let the outside world know via an Exception.
		Q.copy_values_from(_Qinit);
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
	    Q.copy_values_from(_Qinit);
	    throw new ThermochemicalReactorUpdateException(errMsg);
	}
	// At this point, it appears that everything has gone well.
	// Update energies and leave.
	dtSuggest = dtSave;
    }

    double evalRelaxationTime(GasState Q, int isp)
    {
	double tauInv = 0.0;
	foreach (csp; 0 .. _airModel.n_species) {
	    // We exclude electrons in the vibrational relaxation time.
	    if (_airModel.species_name(csp) == "e-") {
		continue;
	    }
	    // 1. Compute Millikan-White value for each interaction
	    double pTau = exp(_A[isp] * (pow(Q.T, -1./3) - 0.015*pow(_mu[isp][csp], 0.25)) - 18.42);
	    double pCollider = _molef[csp]*Q.p/P_atm;
	    double tau = pTau/pCollider;
	    tauInv += 1.0/tau;
	}
	double tauMW = 1.0/tauInv;
	// 2. Compute Park value for high-temperature correction
	double k = Boltzmann_constant;
	double nd = Q.p/(k*Q.T); // 1/m^3
	nd *= 1.0e-6; // convert 1/m^3 --> 1/cm^3
	double c = sqrt(8*k*Q.T/(PI*_particleMass[isp]));
	c *= 100.0; // convert m/s --> cm/s
	double sigmaDash = 3.0e-17; // cm^2
	double sigma = sigmaDash*((50000.0/Q.T)^^2); // cm^2
	double tauP = 1.0/(nd*c*sigma);
	// 3. Combine M-W and Park values.
	double tauV = tauMW + tauP;
	// 4. Correct Landau-Teller at high temperatures by
	//    modifying tau.
	double s = 3.5*exp(-5000.0/_T_sh);
	tauInv = (1.0/tauV)*(pow(fabs((_T_sh - Q.T_modes[0])/(_T_sh - _Tv_sh)), s-1.0));
	return 1.0/tauInv;
    }

    double evalRate(GasState Q)
    {
	double rate = 0.0;
	foreach (isp; _airModel.molecularSpecies) {
	    // Vibrational energy exchange via collisions.
	    double tau = evalRelaxationTime(Q, isp);
	    double evStar = _airModel.vibEnergy(Q.T, isp);
	    double ev = _airModel.vibEnergy(Q.T_modes[0], isp);
	    rate += Q.massf[isp] * (evStar - ev)/tau;
	    // Vibrational energy change due to chemical reactions.
	    double chemRate = 0.0;
	    foreach (ir; _reactionsByMolecule[isp]) {
		chemRate += _chemUpdate.rmech.rate(ir, isp);
	    }
	    chemRate *= _airModel.mol_masses[isp]; // convert mol/m^3/s --> kg/m^3/s
	    double Ds = _c2*ev;
	    rate += Q.massf[isp]*chemRate*Ds;
	}
	return rate;
    }

    ResultOfStep RKFStep(GasState Q, double h, ref double hSuggest)
    {
	_Q0.copy_values_from(Q);
	immutable double a21=1./5., a31=3./40., a32=9./40.,
	    a41=3./10., a42=-9./10., a43 = 6./5.,
	    a51=-11./54., a52=5./2., a53=-70./27., a54=35./27.,
	    a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096.;
	immutable double b51=37./378., b53=250./621., b54=125./594., b56=512./1771.,
	    b41=2825./27648., b43=18575./48384., b44=13525./55296., b45=277./14336., b46=1.0/4.0;

	double k1 = evalRate(Q);
	Q.u_modes[0] = _Q0.u_modes[0] + h*(a21*k1);
	Q.u = _eTotal - Q.u_modes[0];
	_airModel.update_thermo_from_rhou(Q);

	double k2 = evalRate(Q);
	Q.u_modes[0] = _Q0.u_modes[0] + h*(a31*k1 + a32*k2);
	Q.u = _eTotal - Q.u_modes[0];
	_airModel.update_thermo_from_rhou(Q);

	double k3 = evalRate(Q);
	Q.u_modes[0] = _Q0.u_modes[0] + h*(a41*k1 + a42*k2 + a43*k3);
	Q.u = _eTotal - Q.u_modes[0];
	_airModel.update_thermo_from_rhou(Q);

	double k4 = evalRate(Q);
	Q.u_modes[0] = _Q0.u_modes[0] + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4);
	Q.u = _eTotal - Q.u_modes[0];
	_airModel.update_thermo_from_rhou(Q);

	double k5 = evalRate(Q);
	Q.u_modes[0] = _Q0.u_modes[0] + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5);
	Q.u = _eTotal - Q.u_modes[0];
	_airModel.update_thermo_from_rhou(Q);

	double k6 = evalRate(Q);
	Q.u_modes[0] = _Q0.u_modes[0] + h*(b51*k1 + b53*k3 + b54*k4 + b56*k6);
	Q.u = _eTotal - Q.u_modes[0];
	_airModel.update_thermo_from_rhou(Q);

	// Compute error estimate.
	double errEst = Q.u_modes[0] - (_Q0.u_modes[0] + h*(b41*k1 + b43*k3 + b44*k4 + b45*k5 + b46*k6));

	// And use error estimate as a means to suggest a new timestep.
	double atol = _energyAbsTolerance;
	double rtol = _energyRelTolerance;
	double sk = atol + rtol*fmax(fabs(_Q0.u_modes[0]), fabs(Q.u_modes[0]));
	double err = errEst/sk;
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

