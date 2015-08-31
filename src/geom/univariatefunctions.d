/** univariatefunctions.d
 * Univariate functions ported from the C++ code -- for use as clustering functions.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-02-20 first code
 */

module univariatefunctions;

import std.conv;
import std.math;

// Base class for functions of one variable.
class UnivariateFunction {
public:
    double opCall(double t) const
    {
	return t; // expect this to be overridden
    }
    override string toString() const
    {
	return "UnivariateFunction()";
    }
    double[] distribute_parameter_values(int n, double t0=0.0, double t1=1.0) const
    // Returns an array of parameter values, distributed from t0 though t1.
    // The subclass determines the form of distribution.
    {
	double[] tv;
	double dt = 1.0 / (n-1);
	foreach (i; 0 .. n) {
	    double t = opCall(dt*i); // first mapping is done by the specific function
	    tv ~= (1.0-t)*t0 + t*t1; // map to specified range and save
	}
	return tv;
    } // end distribute_parameter_values()
} // end class

class LinearFunction : UnivariateFunction {
public:
    double t0;
    double t1;
    this(double t0=0.0, double t1=1.0)
    {
	this.t0 = t0;
	this.t1 = t1;
    }
    this(const LinearFunction other)
    {
	t0 = other.t0;
	t1 = other.t1;
    }
    LinearFunction dup() const
    {
	return new LinearFunction(this);
    }
    override double opCall(double t) const
    {
	return (1.0-t) * t0 + t * t1;
    }
    override string toString() const
    {
	return "LinearFunction(t0=" ~ to!string(t0) ~ ", t1=" ~ to!string(t1) ~ ")";
    }
} // end class LinearFunction

class RobertsFunction : UnivariateFunction {
public:
    // Stretching parameters for original Robert's transform.
    bool end0;
    bool end1;
    double alpha = 0.5;
    double beta;
    bool reverse = false;
    bool cluster = true;
    this(bool end0, bool end1, double beta)
    {
	this.end0 = end0;
	this.end1 = end1;
	this.beta = beta;
	if (!end0 && !end1) cluster = false;
	if (beta <= 1.0) cluster = false;
	if (end0 && end1) alpha = 0.5;
	if (end0 && !end1) {
	    reverse = true;
	    alpha   = 0.0;
	}
	if (!end0 && end1) {
	    reverse = 0;
	    alpha   = 0.0;
	}
    }
    this(const RobertsFunction other)
    {
	this.end0 = end0;
	this.end1 = end1;
	beta = other.beta;
	alpha = other.alpha;
	reverse = other.reverse;
	cluster = other.cluster;
    }
    RobertsFunction dup() const
    {
	return new RobertsFunction(this);
    }
    override double opCall(double t) const
    {
	double tbar;
	if (reverse) t = 1.0 - t;
	if (cluster) {
	    tbar = roberts_original(t, alpha, beta);
	} else { 
	    tbar = t;
	}
	if (reverse) tbar = 1.0 - tbar;
	return tbar;
    }
    override string toString() const
    {
	return "RobertsFunction(end0=" ~ to!string(end0) ~ 
	    ", end1=" ~ to!string(end1) ~ ", beta=" ~ to!string(beta) ~ ")";
    }
} // end class RobertsFunction()

double roberts_original(double eta, double alpha, double beta)
// eta   : unstretched coordinate, 0 <= eta <= 1.0
// alpha : location of stretching
//         alpha = 0.5 : clustering of nodes at both extremes of eta
//         alpha = 0.0 : nodes will be clustered near eta=1.0
// beta  : stretching factor (more stretching as beta --> 1.0)
// Returns stretched coordinate, 0 <= roberts <= 1.0
{
    double lambda = (beta + 1.0) / (beta - 1.0);
    lambda = pow ( lambda, ((eta - alpha)/(1.0 - alpha)) );
    double etabar = (beta + 2.0 * alpha) * lambda - beta + 2.0 * alpha;
    etabar = etabar / ((2.0 * alpha + 1.0) * (1.0 + lambda));
    return etabar;
} // end roberts_original()


unittest {
    auto cf = new RobertsFunction(false, true, 1.1);
    assert(approxEqual(cf(0.1), 0.166167), "RobertsFunction");
    assert(approxEqual(cf(0.9), 0.96657), "RobertsFunction");
    auto cf2 = new LinearFunction(1.0, 0.0);
    assert(approxEqual(cf2(0.1), 0.9), "LinearFunction");
    assert(approxEqual(cf2(0.9), 0.1), "LinearFunction");
}
