/**
 * This module provides routines to update a gas phase
 * chemistry system.
 *
 * Author: Rowan G. and Kyle D.
 */

module kinetics.chemistry_update;

import std.stdio;
import std.string;
import std.math;
import std.algorithm;
import util.lua;
import util.lua_service;
import gas;
import kinetics.reaction_mechanism;

immutable int MAX_SUBCYCLES = 10000; // maximum number of subcycles to perform over tInterval
immutable int MAX_STEP_ATTEMPTS = 3; // maximum number of attempts on any one step
immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
immutable double DT_REDUCTION_FACTOR = 2.0; // factor by which to reduce timestep
                                            // after a failed attempt
immutable double H_MIN = 1.0e-15; // Minimum allowable step size

static bool working_memory_allocated = false;
static GasState Qinit;
static double[] conc0;
static double[] concOut;
enum ResultOfStep { success, failure };

final class ReactionUpdateScheme {
    ReactionMechanism rmech;
    ChemODEStep cstep;
    bool tightTempCoupling;
    int maxSubcycles;
    int maxAttempts;

    this(string fname, GasModel gmodel)
    {
	auto L = init_lua_State(fname);
	lua_getglobal(L, "reaction");
	rmech = createReactionMechanism(L, gmodel);
	lua_close(L);
	// Just hard code selection of RKF for present.
	cstep = new RKFStep(gmodel, rmech, 1.0e-3);
	// cstep = new AlphaQssStep(gmodel, rmech);
	tightTempCoupling = false;
	maxSubcycles = MAX_SUBCYCLES;
	maxAttempts = MAX_STEP_ATTEMPTS;
    }

    void update_state(GasState Q, double tInterval, ref double dtSuggest, GasModel gmodel)
    {
	if ( update_chemistry(Q, tInterval, dtSuggest, gmodel, rmech, cstep,
			      tightTempCoupling, maxSubcycles, maxAttempts) != 0 ) {
	    string errMsg = "There was a problem with the chemistry update.";
	    throw new Exception(errMsg);
	}
    }

}

int update_chemistry(GasState Q, double tInterval, ref double dtSuggest,
		     GasModel gmodel, ReactionMechanism rmech, ChemODEStep cstep,
		     bool tightTempCoupling, int maxSubcycles=MAX_SUBCYCLES, int maxAttempts=MAX_STEP_ATTEMPTS)
{
    // 0. On entry take a copy of the GasState in case we bugger it up.
    if ( !working_memory_allocated ) {
	Qinit = new GasState(gmodel.n_species, gmodel.n_modes);
	conc0.length = gmodel.n_species;
	concOut.length = gmodel.n_species;
	working_memory_allocated = true;
    }
    Qinit.copy_values_from(Q);
    gmodel.massf2conc(Q, conc0);

    // 0. Evaluate the rate constants. 
    //    It helps to have these computed before doing other setup work.
    rmech.eval_rate_constants(Q);
    // 1. Sort out the time step for possible subcycling.
    double t = 0.0;
    double h;
    if ( dtSuggest > tInterval )
	h = tInterval;
    else if ( dtSuggest <= 0.0 )
        h = rmech.estimateStepSize(conc0);
    else
	h = dtSuggest;
    double dtSave = h;
    // 2. Now do the interesting stuff, increment species change
    //writeln("THIS IS THE FIRST dt = ", h);
    int cycle = 0;
    int attempt = 0;
    for ( ; cycle < maxSubcycles; ++cycle ) {
	/* Copy the last good timestep before moving on.
	 * If we're near the end of the interval, we want
	 * the penultimate step. In other words, we don't
	 * want to store the fractional step of (tInterval - t)
	 * that is taken as the last step.
	 */
	dtSave = h;
	h = min(h, tInterval - t);
	attempt= 0;
      	for ( ; attempt < maxAttempts; ++attempt ) {
	    if ( cstep(conc0, h, concOut, dtSuggest) == ResultOfStep.success ) {
		/* We succesfully took a step of size h
		 * so increment the total time.
		 */       
		t += h;
		 foreach ( i; 0..conc0.length ) conc0[i] = concOut[i];
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
		 */
		h /= DT_REDUCTION_FACTOR;
		if ( h < H_MIN ) {
		    Q.copy_values_from(Qinit);
		    return -1;
		}
	    }
	} // end attempts at single step.
	if ( attempt == MAX_STEP_ATTEMPTS ) {
	    writeln("WARNING: hit max step attempts within a subcycle.");
	    // We did poorly. Let the outside world know by returning -1.
	    // But put the original GasState back in place.
	    Q.copy_values_from(Qinit);
	    return -1;
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
	    gmodel.conc2massf(conc0, Q);
	    gmodel.update_thermo_from_rhoe(Q);
	    rmech.eval_rate_constants(Q);
	}

	if ( t >= tInterval ) { // We've done enough work.
	    break;
	}
    }
    if ( cycle == maxSubcycles ) {
	// If we make it out here, then we didn't complete within
	// the maximum number of subscyles... we are taking too long.
	// Let's return the gas state to how it was before we failed.
	writeln("WARNING: hit max number of subcycles.");
	Q.copy_values_from(Qinit);
	return -1;
    }
    // else all is well, so update GasState Q and leave.
    gmodel.conc2massf(concOut, Q);
    gmodel.update_thermo_from_rhoe(Q);
    dtSuggest = dtSave;
    return 0;
}

/++
 + ChemODEStep is an abstract class that provides the interface
 + for a chemistry ODE update step. This class provides just one
 + public service so we make that happen when the object is
 + called like a function with opCall().
 +
 + You might rightly ask "Why use a class when a function would do?"
 + It turns out that there is a lot of extra state that needs to be
 + stored and a lot of working array space. This is possible as
 + module level static data. However, it just plain looks neater
 + in a class.
 +/

class ChemODEStep
{
public:
    this(ReactionMechanism rmech)
    {
	_rmech = rmech;
    }
    abstract ResultOfStep opCall(double[] y0, double h, double[] yOut, ref double hSuggest);
private:
    ReactionMechanism _rmech;
}

/++
 + The Runge-Kutta-Fehlberg method, specialised to work
 + with a chemical kinetic ODE.
 +
 + References:
 + Fehlberg, E. (1969)
 + Low-order classical Runge-Kutta formulas with stepsize control
 + and their application to some heat transfer problems.
 + NASA Technical Report, R-315
 +
 + Cash, J. R. and Karp, A. H. (1990)
 + A Variable Order Runge-Kutta Method for Initial Value
 + Problems with Rapidly Varying Right-Hand Sides
 + ACM Transactions on Mathematical Software, 16:3, pp. 201--222
 +
 + Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007)
 + Numerical Recipes: The Art of Scientific Computing, Third Edition
 + Cambridge University Press, New York, USA
 +/

class RKFStep : ChemODEStep {
public:
    this(in GasModel gmodel, ReactionMechanism rmech, double tolerance)
    {
	super(rmech);
	_tolerance = tolerance;
	// Allocate working arrays
	_ndim = gmodel.n_species;
	_yTmp.length = _ndim;
	_yErr.length = _ndim;
	_k1.length = _ndim;
	_k2.length = _ndim;
	_k3.length = _ndim;
	_k4.length = _ndim;
	_k5.length = _ndim;
	_k6.length = _ndim;
    }

    override ResultOfStep opCall(double[] y0, double h, double[] yOut, ref double hSuggest)
    {
	// 0. Set up the constants associated with the update formula
	// We place them in here for visibility reasons... no one else
	// needs to be using the coefficients a21, a31, etc.
	immutable double a21=1./5., a31=3./40., a32=9./40.,
	    a41=3./10., a42=-9./10., a43 = 6./5.,
	    a51=-11./54., a52=5./2., a53=-70./27., a54=35./27.,
	    a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096.;
	immutable double b51=37./378., b53=250./621., b54=125./594., b56=512./1771.,
	    b41=2825./27648., b43=18575./48384., b44=13525./55296., b45=277./14336., b46=1.0/4.0;
	
	// 1. Apply the formula to evaluate intermediate points
	_rmech.eval_rates(y0, _k1);
	_yTmp[] = y0[] + h*(a21*_k1[]);

	_rmech.eval_rates(_yTmp, _k2);
	_yTmp[] = y0[] + h*(a31*_k1[] + a32*_k2[]);

	_rmech.eval_rates(_yTmp, _k3);
	_yTmp[] = y0[] + h*(a41*_k1[] + a42*_k2[] + a43*_k3[]);

	_rmech.eval_rates(_yTmp, _k4);
	_yTmp[] = y0[] + h*(a51*_k1[] + a52*_k2[] + a53*_k3[] + a54*_k4[]);

	_rmech.eval_rates(_yTmp, _k5);
	_yTmp[] = y0[] + h*(a61*_k1[] + a62*_k2[] + a63*_k3[] + a64*_k4[] + a65*_k5[]);

	_rmech.eval_rates(_yTmp, _k6);
	
	// 2. Compute new value and error esimate
	yOut[] = y0[] + h*(b51*_k1[] + b53*_k3[] + b54*_k4[] + b56*_k6[]);
	_yErr[] = yOut[] - (y0[] + h*(b41*_k1[] + b43*_k3[] + b44*_k4[] + b45*_k5[] + b46*_k6[]));

	// 3. Lastly, use error estimate as a means to suggest
	//    a new timestep.
	//    Compute error using tol as atol and rtol as suggested in
	//    Press et al. (2007)
	double err = 0.0;
	double sk = 0.0;
	double atol = _tolerance;
	double rtol = _tolerance;

	foreach ( i; 0.._ndim ) {
	    sk = atol + rtol*max(fabs(y0[i]), fabs(yOut[i]));
	    err += (_yErr[i]/sk)*(_yErr[i]/sk);
	}
	err = sqrt(err/_ndim);
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
		// We are NOT using the PI control version as
		// given by Press et al. as I do not want to store
		// the old error from previous integration steps.
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
	scale = max(safe*pow(err, -alpha), minscale);
	hSuggest = scale*h;
	return ResultOfStep.failure;
    }
private:
    int _ndim;
    double _tolerance;
    double[] _yTmp, _yErr;
    double[] _k1, _k2, _k3, _k4, _k5, _k6;
}

/++
 + The Alpha-QSS method, specialised to work
 + with a chemical kinetic ODE.
 + 
 + author: K.Damm 2015
 + adapted from R.Gollan's eilmer3 implementation
 +
 + References:
 + Mott, D. (1999)
 + Integrates a system of chemical kinetics ODEs which contains both
 + stiff and non-stiff equations. It is proven to be A-stable and has
 + second order accuracy. The method is single-step self starting,
 + thus has little start-up costs.
 +
 + S.R. Qureshi, R. Prosser (2007)
 + Delta term added into Mott algorithm to prevent excessively small
 + concentrations forcing the solver to take miniscule steps.
 + Default delta = 10.0e-10.
 +/

class AlphaQssStep : ChemODEStep {
public:
  this(in GasModel gmodel, ReactionMechanism rmech)
    {
	super(rmech);
	// Allocate working arrays
	_ndim = gmodel.n_species;
	_yp.length = _ndim;
	_yp0.length = _ndim;
	_yc.length = _ndim;
	_q0.length = _ndim;
	_L0.length = _ndim;
	_p0.length = _ndim;
	_qtilda.length = _ndim;
	_pbar.length = _ndim;
	_qp.length = _ndim;
        _Lp.length = _ndim;
	_pp.length = _ndim;
	_alpha.length = _ndim;
	_alphabar.length = _ndim;
    }

    override ResultOfStep opCall(double[] y0, double h, double[] yOut, ref double hSuggest)
    {
	// 0. Alpha-QSS Constant parameters (currently hard-coded)
        static immutable double _EPS1 = 0.001;
        static immutable double _EPS2 = _EPS1/2.0;
        static immutable int _MAX_ITER = 10;
        static immutable double _ZERO_EPS = 1.0e-50;
	static immutable double _delta = 1.0e-10;
	
	// 1. Predictor Step
	_rmech.eval_split_rates(y0, _q0, _L0);
	p_on_y(_L0, y0, _p0);
        alpha_compute(_p0, _alpha, h);
	update_conc(_yp, y0, _q0, _p0, _alpha, h);
	// put aside the first predictor values for
	// later use in the convergence test
	foreach ( i; 0.._ndim )
	  _yp0[i] = _yp[i];

	// 2. Corrector Step
	foreach ( corr; 0.._MAX_ITER) {
	    _rmech.eval_split_rates(_yp, _qp, _Lp);
	    p_on_y(_Lp, _yp, _pp);
	    alpha_compute(_pp, _alphabar, h);
	    foreach ( i; 0.._ndim ) {
		_qtilda[i] = _alphabar[i]*_qp[i] + (1-_alphabar[i])*_q0[i];
		_pbar[i] = 0.5*(_p0[i]+_pp[i]);
	    }

	    // actual corrector step
	    update_conc(_yc, y0, _qtilda, _pbar, _alphabar, h);

	    // test if converged or not
	    bool converged = test_converged(_yc, _yp0, h);

	    if ( converged ) {
		// the solver has converged, let's set the current
		// corrector to the yOut vector and suggest a new
		// time step.
		foreach ( i; 0.._ndim ) {
		    yOut[i] = _yc[i];
		}
		hSuggest = step_suggest(h, _yc, _yp0);
		return ResultOfStep.success;
	    }
	    // if we haven't converged yet make the corrector the new predictor and try again
	    foreach ( i; 0.._ndim ) {
		_yp[i] = _yc[i];
	    }
	}
	// if we have failed the convergence test too many times let's start again with a smaller h
	// This suggested h will be passed through generic step size  check on outer loop
        hSuggest = step_suggest(h, _yc, _yp0);
	return ResultOfStep.failure;
    }

private:
    int _ndim;
    double[] _yTmp, _yp, _yp0,  _yc, _q0, _L0, _p0, _pp, _qtilda, _pbar, _qp, _Lp, _alpha, _alphabar;
    immutable double _EPS1 = 0.001;
    immutable double _EPS2 = _EPS1/2.0;
    immutable int _MAX_ITER = 10;
    immutable double _ZERO_EPS = 1.0e-50;
    immutable double _delta = 1.0e-10;
    
    // Private functions.
    void p_on_y(double[] L, double[] y, double[] p_y) {
	foreach( i; 0.._ndim) 
	    p_y[i] = L[i] / (y[i] + _ZERO_EPS);
	return;
    }

    void alpha_compute(double[] p, double[] alpha, double h) {
        foreach ( i; 0.._ndim ) {
	    double r = 1.0/(p[i]*h+_ZERO_EPS);  // ZERO_EPS prevents division by 0 
	    alpha[i] = (180.0*r*r*r+60.0*r*r+11.0*r+1.0)/(360.0*r*r*r+60.0*r*r+12.0*r+1.0);
	}
    }

    void update_conc(double[] yTmp, double[] y0, double[] q, double[] p, double[] alpha, double h) {
        foreach ( i; 0.._ndim )
	    yTmp[i] = y0[i] + (h*(q[i]-p[i]*y0[i]))/(1.0+alpha[i]*h*p[i]);
    }

    bool test_converged(in double[] yc, in double[] yp, double h) {
	int flag = 0;
	double test = 0.0;
	foreach ( i; 0.._ndim ) {
	    if ( yc[i] < _ZERO_EPS )
		continue;
	    test = fabs(yc[i] - yp[i]);
	    // +delta from Qureshi and Prosser (2007)
	    if ( test >= (_EPS1 * (yc[i] + _delta)) )
	      ++flag;
	}
	if ( flag == 0 )
	    return true;
	else
	    return false;
    }

    double step_suggest(double h, double[] yc, double[] yp)
    {
	double test = 0.0;
	double sigma = 0.0;
	double h_new = 0.0;

	foreach (i; 0.._ndim) {
	    if ( yc[i] < _ZERO_EPS )
		continue;
	    test = fabs( yc[i] - yp[i]) / (_EPS2*(yc[i]+ _delta));
	    if ( test > sigma )
		sigma = test;
	}

	if ( sigma <= 0.0 ) {
	    h_new = h;
	}
	else {
	    double x = sigma;
	    x = x - 0.5*(x*x - sigma)/x;
	    x = x - 0.5*(x*x - sigma)/x;
	    x = x - 0.5*(x*x - sigma)/x;
	    h_new = h * ( (1.0 / x) + 0.005 );
	}
    return h_new;
    }
}

version(chemistry_update_test) {
    import std.stdio;
    import util.msg_service;
    import kinetics.rate_constant;
    import kinetics.reaction;
    import gas.therm_perf_gas;

    double numericalEstimate(double dt, double tMax, double[] conc0, ChemODEStep step)
    {
	double t = 0.0;
	double[] conc1;
	conc1.length = conc0.length;
	double dtDummy;
	int count = 0;
	while ( (tMax - t) > 1.0e-9 ) {
	    dt = min(dt, tMax - t);
	    step(conc0, dt, conc1, dtDummy);
	    t += dt;
	    conc0 = conc1.dup;
	    count++;
	}
	return conc1[2];
    }

    int main() {

	auto gmodel = new ThermallyPerfectGas("sample-input/H2-I2-HI.lua");
	auto rmech = createReactionMechanism("sample-input/H2-I2-inp.lua", gmodel);

	auto gd = new GasState(3, 1);
	gd.T[0] = 700.0;
	double c0 = 4.54;
	double[] conc0 = [c0, c0, 0.0];
	gd.p = 2.0*c0*R_universal*gd.T[0];
	double[] molef = [0.5, 0.5, 0.0];
	gmodel.molef2massf(molef, gd);
	gmodel.update_thermo_from_pT(gd);
	rmech.eval_rate_constants(gd);
	double tInterval = 60000.0;
	double analyticalVal = 7.1420197868416215;

	/* 1. Test RKFStep
	 * We'll check that the order of accuracy for the method
	 * is as expected. We'll choose a timestep and run the problem
	 * once and compare the numerical result to an analytical
	 * calculation. Then we'll run the problem with a reduced
	 * timestep h*pow(2.0, -1/5). For this 5th-order accurate
	 * method, we should expect to see the error drop by a factor
	 * of 2.0 if we reduce the timestep by 2.0**(-1/5).
	 */
	auto rkfStep = new RKFStep(gmodel, rmech, 1.0e-3);
	double dt = 287.175;
	double numVal0 = numericalEstimate(dt, tInterval, conc0, rkfStep);
	double err0 = analyticalVal - numVal0;
	// Reduce timestep and repeat test
	dt *= 2.0^^(-1./5.);
	conc0 = [c0, c0, 0.0];
	double numVal1 = numericalEstimate(dt, tInterval, conc0, rkfStep);
	double err1 = analyticalVal - numVal1;
	assert(approxEqual(1.96809, err0/err1, 1.0e-4), failedUnitTest());

	/* 2. alpha-QSS step
	 * [TODO] : Kyle, fix me please.
	 */

	/* 3. Test the complete update algorithm as used
	 *    by the flow solver. This might be stretching a
	 *    bit what a *unit* is. In this test, we'll exercise
	 *    the chemistry update routine by solving the
	 *    complete hydrogen-iodine system.
	 */
	double dtSuggest = 200.0;
	int result = update_chemistry(gd, tInterval, dtSuggest,
				      gmodel, rmech, rkfStep,
				      false, 50, 2);
	double[] conc;
	conc.length = 3;
	gmodel.massf2conc(gd, conc);
	assert(approxEqual(7.14201983840, conc[2], 1.0e-9), failedUnitTest());

	return 0;
    }
}


