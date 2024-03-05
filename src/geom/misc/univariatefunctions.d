/** univariatefunctions.d
 * Univariate functions ported from the C++ code -- for use as clustering functions.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-02-20 first code
 */

module geom.misc.univariatefunctions;

import std.conv;
import std.math;
import std.stdio;

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
    double[] distribute_parameter_values(size_t n, double t0=0.0, double t1=1.0) const
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

class QuadraticFunction : UnivariateFunction {
/*
   A cluster function based on a quadratic. This function is designed for
   linearly increasing or decreasing cell sizes. It's only parameter is
   "ratio" which is the ratio of the largest cell sizes at the end to the
   smallest cell sizes at the beginning.

   @author: Nick Gibbons
*/
    this(double ratio, bool reverse)
    {
        this.ratio = ratio;
        this.reverse = reverse;
        this.a = (ratio-1.0)/(ratio+1.0);
        this.b = 2/(ratio+1.0);
    }
    this(const QuadraticFunction other)
    {
        ratio = other.ratio;
        reverse = other.reverse;
    }
    QuadraticFunction dup() const
    {
        return new QuadraticFunction (this);
    }
    override double opCall(double t) const
    {
        if (reverse) t = 1.0 - t;
        double y = a*t*t + b*t;
        if (reverse) y = 1.0 - y;
        return y;
    }
    override string toString() const
    {
        return "QuadraticFunction(ratio=" ~ to!string(ratio) ~ ", t1=" ~ to!string(reverse) ~ ")";
    }
private:
    double a,b,ratio;
    bool reverse;
} // end class QuadraticFunction

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

class GeometricFunction : UnivariateFunction {
/*
   A cluster function based on geometric progression, starting with a fixed
   size and growing by a constant multiple until a threshold is reached.
   After this the cells are equally sized

   Notes: This function mimics how clustering is done in GridPro's "clu" utility
   @author: Nick Gibbons
*/
public:
    this(double a, double r, int N, bool reverse)
    {
        this.a = a;
        this.r = r;
        this.N = N;
        this.reverse = reverse;
        if (a<=0.0) throw new Error("Problematic parameters in GeometricFunction a= "~to!string(a));
        if (r<=1.0) throw new Error("Problematic parameters in GeometricFunction r= "~to!string(r));
        double ns = solve_for_nswitch(a, r, N);
        if ((ns<1.0) || (ns>N)) throw new Error("Problematic parameters in GeometricFunction ns= "~to!string(ns));
        this.ts = ns/(N-1);
        this.ns = ns;
    }

    this(const GeometricFunction other)
    {
        a = other.a;
        r = other.r;
        N = other.N;
        reverse = other.reverse;
        ts = other.ts;
        ns = other.ns;
    }

    GeometricFunction dup() const
    {
        return new GeometricFunction(this);
    }

    override double opCall(double t) const
    {
        if (reverse) t = 1.0 - t;
        double n = t*(N-1);
        double t1;
        if (t<ts){
            t1 = a*(1-pow(r,n))/(1-r);
        } else {
            t1 = a*(1-pow(r,ns))/(1-r) + (n-ns)*a*pow(r,(ns-1));
        }
        if (reverse) t1 = 1.0 - t1;
        return t1;
    }

private:
    double a,r,ts,ns;
    int N;
    bool reverse;

    double lambertW(double z, double wguess=1.0, double tol=1e-15){
        /*
            Lamberts W function: solve for w given z=w*exp(w)
            See: en.wikipedia.org/wiki/Lambert_W_function#Numerical_evaluation

            Notes: Adjustments were performed to the tol parameter on March 2024.
            This required a concerningly tight tolerance to function properly.
            In the future, perhaps revisit this to avoid possible numerical issues.
        */
        double w, expw, zt, wexpw;
        w = wguess;
        foreach(n; 0 .. 100){
            if (n>=99) throw new Error("lambertW function eval timed out");
            expw = exp(w);
            zt = w*expw;
            if (fabs(zt-z)<tol) break;
            wexpw = w*expw;
            w -= (wexpw - z)/(expw + wexpw);
        }
        return w;
    }

    double solve_for_nswitch(double a, double r, double N){
        /*
           Solve for the switching point where the geometric progression
           switches over to a linear function. (See derivation 20/07/2020)

           Effectively we're solving for ns given:
           1.0 = a*(1-r**ns)/(1-r) + (N-1-ns)*a*r**(ns-1)

        */
        double A = 1/a - 1/(1-r);
        double B = -1.0/(1-r);
        double M = N-1;

        double powarg = -M-B*r+1;
        double arg = -A*pow(r,powarg)*log(r);
        double wguess = (N/4 - M - B*r)*log(r);
        double w = lambertW(arg, wguess=wguess);
        double ns = w/log(r) + B*r + M;
        return ns;
    }

}

class GaussianFunction : UnivariateFunction {
/*
   A cluster function based on the Gaussian Normal Function. Designed for clustering
   in the middle of a block somewhere:
   https://en.wikipedia.org/wiki/Gaussian_function

   Inputs:
     m - The location of peak clustering, (0<m<1)
     s - The width of the peak
     ratio - smallest cell size/largest cell size

   @author: Nick Gibbons
*/
public:
    this(double m, double s, double ratio)
    {
        this.m = m;
        this.s = s;
        this.ratio = ratio;
        if ((m<=0.0) || (m>=1.0)) throw new Error("Problematic parameters in GaussianFunction m= "~to!string(m));
        if ((ratio<=0.0) || (ratio>=1.0)) throw new Error("Problematic parameters in GaussianFunction ratio= "~to!string(ratio));
        // See derivation 20/12/11 (NNG)
        this.a = 1.0/(1.0+sqrt(pi/2.0)*(1.0-ratio)*s*(erf((m-1.0)/sqrt(2.0)/s) - erf(m/sqrt(2.0)/s)));
        this.c = -sqrt(pi/2.0)*a*(1-ratio)*s*erf(m/sqrt(2.0)/s);
    }

    this(const GaussianFunction other)
    {
        m = other.m;
        s = other.s;
        ratio = other.ratio;
        a = other.a;
        c = other.c;
    }

    GaussianFunction dup() const
    {
        return new GaussianFunction(this);
    }

    override double opCall(double x) const
    {
        return a*x + sqrt(pi/2.0)*a*(1.0-ratio)*s*erf((m-x)/sqrt(2.0)/s) + c;
    }

private:
    double m,s,ratio,a,c;
    immutable double pi = 3.1415926535;
    immutable double p  = 0.3275911;
    immutable double a1 = 0.254829592;
    immutable double a2 =-0.284496736;
    immutable double a3 = 1.421413741;
    immutable double a4 =-1.453152027;
    immutable double a5 = 1.061405429;

    const double erf(double x){
        /*
            Approximate the gaussian error function using curve fitted polynomial
            "Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables", NIST
            See: https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions
            Note: The approximation here is only valid for x>=0. For x less than zero we exploit the
            fact that the error function is odd: f(x) = -f(-x)
        */

        double xabs = fabs(x);
        double t = 1.0/(1.0+p*xabs);
        double f = 1.0 - (a1*t + a2*t*t + a3*t*t*t + a4*t*t*t*t + a5*t*t*t*t*t)*exp(-(xabs*xabs));
        f = copysign(f, x);
        return f;
    }
}

class GaussGeomHybridFunction : UnivariateFunction {
/*
   A hybrid function that includes both a geometric growth beginning
   and a smooth gaussian cluster in the middle.

   Inputs:
      - A : nondimensional starting size for the geometric progression (1.0 is the full width)
      - R : Geometric progression growth factor, should be larger than 1
      - N : Number of node points along the boundary.
      - m : location of Gaussian cluster minimum
      - s : Gaussian cluster width, equivalent to a standard deviation
      - ratio: cell size at the minimum divided by baseline largest cell size
      - reverse : whether to put the Geometric progression at the beginning or end of the boundary

   @author: Nick Gibbons
*/
public:
    this(double A, double R, int N, double m, double s, double ratio, bool reverse)
    {
        this.A = A;
        this.R = R;
        this.N = N;
        this.m = m;
        this.s = s;
        this.ratio = ratio;
        this.reverse = reverse;
        if (A<=0.0) throw new Error("Problematic parameters in GaussGeomHybridFunction A= "~to!string(A));
        if (R<=1.0) throw new Error("Problematic parameters in GaussGeomHybridFunction R= "~to!string(R));
        if ((m<=0.0) || (m>=1.0)) throw new Error("Problematic parameters in GaussGeomHybridFunction m= "~to!string(m));
        if ((ratio<=0.0) || (ratio>=1.0)) throw new Error("Problematic parameters in GaussGeomHybridFunction ratio= "~to!string(ratio));

        double xs = solve_for_xswitch();
        if ((xs<0.0) || (xs>1.0)) throw new Error("Problematic parameters in GaussGeomHybridFunction xs= "~to!string(xs));
        this.xs = xs;

        this.a = dGeometricdx(xs)/dGaussiandx(1.0,xs);
        this.c = 1.0 - Gaussian(this.a, 0.0, 1.0);
    }

    this(const GaussGeomHybridFunction other)
    {
        A = other.A;
        R = other.R;
        N = other.N;
        m = other.m;
        s = other.s;
        ratio = other.ratio;
        reverse = other.reverse;
        xs = other.xs;
    }

    GaussGeomHybridFunction dup() const
    {
        return new GaussGeomHybridFunction(this);
    }

    override double opCall(double x) const
    {
        if (reverse) x = 1.0 - x;
        double y;
        if (x<xs){
            y = Geometric(x);
        } else {
            y = Gaussian(a, c, x);
        }
        if (reverse) y = 1.0 - y;
        return y;
    }

private:
    double A,R,m,s,ratio,xs;
    int N;
    double a,c;
    bool reverse;

    immutable double pi = 3.1415926535;
    immutable double p  = 0.3275911;
    immutable double a1 = 0.254829592;
    immutable double a2 =-0.284496736;
    immutable double a3 = 1.421413741;
    immutable double a4 =-1.453152027;
    immutable double a5 = 1.061405429;

    const double erf(double x){
        /*
            Approximate the gaussian error function using curve fitted polynomial
            "Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables", NIST
            See: https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions
            Note: The approximation here is only valid for x>=0. For x less than zero we exploit the
            fact that the error function is odd: f(x) = -f(-x)
        */

        double xabs = fabs(x);
        double t = 1.0/(1.0+p*xabs);
        double f = 1.0 - (a1*t + a2*t*t + a3*t*t*t + a4*t*t*t*t + a5*t*t*t*t*t)*exp(-(xabs*xabs));
        f = copysign(f, x);
        return f;
    }

    const double Gaussian(double a_, double c_, double x){
        /*
           Note that during setup, we don't yet know what a and c are. During the process
           of figuring them out we need to call this function with dummy values (like a=1)
           so they are labelled a_ instead of a to avoid name clashes.

        */
        return a_*x + sqrt(pi/2.0)*a_*(1.0-ratio)*s*erf((m-x)/sqrt(2.0)/s) + c_;
    }

    const double dGaussiandx(double a_, double x){
        return a_ - a_*(1.0-ratio)*exp(-(m-x)*(m-x)/2.0/s/s);
    }

    const double d2Gaussiandx2(double a_, double x){
        return -a_*(1.0-ratio)*exp(-(m-x)*(m-x)/2.0/s/s)*(m-x)/s/s;
    }

    const double Geometric(double x){
        double n = x*(N-1.0);
        return A*(1.0-pow(R,n))/(1.0-R);
    }

    const double dGeometricdx(double x){
        double n = x*(N-1.0);
        return A*(N-1.0)*log(R)*pow(R,n)/(R-1.0);
    }

    const double d2Geometricdx2(double x){
        double dEdx = dGeometricdx(x);
        return (N-1.0)*log(R)*dEdx;
    }

    const double SwitchEquationError(double xs){
        /*
           This function returns zero when the following three constraints are met:
            1.) 1.0 = Gaussian(1.0)
            2.) Gaussian(xswitch) = Geometric(xswitch)
            3.) dGaussiandx(xswitch) = dGeometricdx(xswitch)

           It is solved numerically to get a smooth crossover point between the two functions.
        */
        double E = Geometric(xs);
        double dEdx = dGeometricdx(xs);
        double G = Gaussian(1.0, 0.0, xs);
        double G1= Gaussian(1.0, 0.0, 1.0);
        double dGdx = dGaussiandx(1.0, xs);

        double error = dEdx*(G-G1)/dGdx + 1.0 - E;
        return error;
    }

    const double SwitchEquationDerivative(double xs){
        double E = Geometric(xs);
        double dEdx = dGeometricdx(xs);
        double d2Edx2 = d2Geometricdx2(xs);
        double G = Gaussian(1.0, 0.0, xs);
        double G1= Gaussian(1.0, 0.0, 1.0);
        double dGdx = dGaussiandx(1.0, xs);
        double d2Gdx2 = d2Gaussiandx2(1.0, xs);

        double de = (G-G1)*(dGdx*d2Edx2 - d2Gdx2*dEdx)/dGdx/dGdx;
        return de;
    }

    const double solve_for_xswitch(double guess=0.5, double tol=1e-8){
        /*
            One variable Newton's method to find the switching point between the two functions
        */
        double xs = guess;
        double error = 1e32;
        double f,dfdx,dxs;
        int iterations=0;

        while (error>tol){
            f = SwitchEquationError(xs);
            error = sqrt(f*f);
            dfdx = SwitchEquationDerivative(xs);
            dxs = -f/dfdx;
            xs += dxs;
            iterations += 1;
            if (iterations>99) throw new Exception("Newton solver for xswitch timed out!");
        }
        return xs;
    }
}

version(univariatefunctions_test) {
    import util.msg_service;
    int main() {
        auto cf = new RobertsFunction(false, true, 1.1);
        assert(isClose(cf(0.1), 0.166167, 1.0e-4), failedUnitTest());
        assert(isClose(cf(0.9), 0.96657, 1.0e-4), failedUnitTest());
        auto cf2 = new LinearFunction(1.0, 0.0);
        assert(isClose(cf2(0.1), 0.9, 1.0e-4), failedUnitTest());
        assert(isClose(cf2(0.9), 0.1, 1.0e-4), failedUnitTest());
        auto cf3 = new GeometricFunction(0.005, 1.1, 40, false);
        assert(isClose(cf3(0.1), 0.0225106, 1.0e-4), failedUnitTest());
        assert(isClose(cf3(0.9), 0.85256413, 1.0e-4), failedUnitTest());
        auto cf4 = new GaussianFunction(0.5, 0.1, 0.2);
        assert(isClose(cf4(0.1), 0.1250750, 1.0e-4), failedUnitTest());
        assert(isClose(cf4(0.9), 0.8749249, 1.0e-4), failedUnitTest());
        auto cf5 = new GaussGeomHybridFunction(0.01, 1.2, 40, 0.8, 0.1, 0.2, false);
        assert(isClose(cf5(0.1), 0.0518068, 1.0e-4), failedUnitTest());
        assert(isClose(cf5(0.9), 0.8984721, 1.0e-4), failedUnitTest());
        return 0;
    }
} // end univariatefunctions_test
