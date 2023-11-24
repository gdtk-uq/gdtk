/** nelmin.d
 *
 * Nelder-Mead simplex minimization of a nonlinear (multivariate) function.
 *
 * Author: Peter J.
 * Version: 2020-07-25 made more like the Python 2020 version.
 *
 * This D version has been adapted from the C++ version 2014-06-14.
 * The C++ version had been adpated from the Python version in Jan 2006.
 * The Python code had been adapted from the C-coded nelmin.c which was
 * adapted from the Fortran-coded nelmin.f which was, in turn, adapted
 * from the papers
 *
 *     J.A. Nelder and R. Mead (1965)
 *     A simplex method for function minimization.
 *     Computer Journal, Volume 7, pp 308-313.
 *
 *     R. O'Neill (1971)
 *     Algorithm AS47. Function minimization using a simplex algorithm.
 *     Applied Statistics, Volume 20, pp 338-345.
 *
 * For the 2020 update, we make some of the stepping closer to
 * the description given in the paper:
 *
 *    Donghoon Lee and Matthew Wiswall (2007)
 *    A parallel implementation of the simplec function minimization routine.
 *    Computational Economics 30:171-187
 *
 * and some examples are in
 *
 *    D.M. Olsson and L.S. Nelson (1975)
 *    The Nelder-Mead Simplex procedure for function minimization.
 *    Technometrics, Volume 17 No. 1, pp 45-51.
 *
 * For a fairly recent and popular incarnation of this minimizer,
 * see the amoeba function in the famous "Numerical Recipes" text.
 * The programming interface is via the minimize() function; see below.
 */

module nm.nelmin;
import std.conv;
import std.math;
import std.stdio;
import std.format;
import std.algorithm;
import ntypes.complex;

//-----------------------------------------------------------------------
// The public face of the minimizer...

/**
 * Locate a minimum of the objective function, f.
 *
 * Input:
 *     f       : user-supplied multivariate function double f(double x[])
 *     x       : vector of N coordinates
 *     dx      : vector of N increments to apply to x when forming
 *               the initial simplex.  Their magnitudes determine the size
 *               and shape of the initial simplex.
 *     tol     : the terminating limit for the standard-deviation
 *               of the simplex function values.
 *     P       : number of points to replace in parallel, each step.
 *     maxfe   : maximum number of function evaluations that we will allow
 *     n_check : number of steps between convergence checks
 *     delta   : magnitude of the perturbations for checking a local minimum
 *               and for the scale reduction when restarting
 *     Kreflect, Kextend, Kcontract: coefficients for locating the new vertex
 *
 * Output:
 *     Returns a flag to indicate if convergence was achieved.
 *     On return x contains the coordinates for the best x location,
 *     corresponding to min(f(x)),
 *     f_min      : the function value at that point,
 *     n_fe      : the number of function evaluations and
 *     n_restart : the number of restarts (with scale reduction).
 */
bool minimize(alias f, T)
    (ref T[] x, out T f_min, out int n_fe, out int n_restart, T[] dx,
     in double tol=1.0e-6, in int P=1, in int max_fe=400,
     in int n_check=20, in double delta=0.001,
     in double Kreflect=1.0, in double Kextend=2.0, in double Kcontract=0.5)
    if (is(typeof(f([0.0,0.0])) == double) ||
        is(typeof(f([Complex!double(0.0),Complex!double(0.0)])) == Complex!double))
{
    bool converged = false;
    size_t N = x.length;
    // Might have to check for NaNs as well.
    if (dx.length != x.length) {
        dx = x.dup;
        foreach (ref elem; dx) elem = 0.1;
    }
    auto smplx = new NMMinimizer!(f,T)(N, dx, P, Kreflect, Kextend, Kcontract);
    smplx.build_initial_simplex(x);
    while (!converged && (smplx.nfe < max_fe)) {
        smplx.take_steps(n_check);
        smplx.vertices.sort!(less);
        // Pick out the current best vertex.
        x[] = smplx.vertices[0].x[];
        f_min = smplx.vertices[0].fx;
        // Check the scatter of vertex values to see if we are
        // close enough to call it quits.
        auto mean_std_dev = smplx.f_statistics();
        if (mean_std_dev[1] < tol) {
            // All of the points are close together but we need to test more carefully.
            converged = smplx.test_for_minimum(delta);
            if (!converged) smplx.rescale(delta);
        }
    } // end while !converged
    n_fe = smplx.nfe;
    n_restart = smplx.nrestarts;
    return converged;
} // end minimize()

//-----------------------------------------------------------------------
// Keep the data tidy and conveniently accessible...

struct Vertex(T) {
    // Stores the coordinates as and array and the associated function value.
    T[] x;
    T fx;
}

// Following function used when sorting vertices into ascending order.
bool less(T)(Vertex!T a, Vertex!T b) { return a.fx < b.fx; }


class NMMinimizer(alias f, T)
    if (is(typeof(f([0.0,0.0])) == double) ||
        is(typeof(f([Complex!double(0.0),Complex!double(0.0)])) == Complex!double))
{
    //  Stores the (nonlinear) simplex and some stepping parameters.

    this(size_t N, T[] dx, int P, double Kreflect, double Kextend, double Kcontract)
    // In an N-dimensional problem, each vertex has an array of N coordinates
    // and the simplex consists of N+1 vertices.
    // and the set of step-sizes dx.
    // f is a user-specified objective function f(x).
    {
        assert(dx.length == N, "Incorrect length for dx.");
        this.dx = dx.dup;
        this.N = N;
        this.P = P;
        this.Kreflect = Kreflect;
        this.Kextend = Kextend;
        this.Kcontract = Kcontract;
        this.nrestarts = 0;
        this.nfe = 0;
    }

    void build_initial_simplex(T[] x)
    // Set up the vertices about the user-specified point, x,
    {
        assert(x.length == N, "Incorrect length for x.");
        foreach (i; 0 .. N+1) {
            auto x_new = x.dup;
            if (i >= 1) x_new[i-1] += dx[i-1];
            auto f_new = f(x_new);
            vertices ~= Vertex!T(x_new, f_new);
        }
        vertices.sort!(less);
        this.nfe = to!int(N+1);
        return;
    }

    void take_steps(int nsteps)
    // Take a few steps, updating the simplex.
    // Returns a sorted simplex, with the best point at index 0.
    {
        foreach (istep; 0 .. nsteps) {
            vertices.sort!(less);
            bool any_success = false;
            foreach (i; 0 .. P) {
                bool success = replace_vertex(N-i);
                if (success) any_success = true;
            }
            if (!any_success) {
                contract_about_zero_point();
            }
        }
        return;
    }

    void rescale(double ratio)
    // Rebuild the simplex about the lowest point for a restart.
    {
        vertices.sort!(less);
        foreach (ref elem; dx) { elem *= ratio; }
        auto vtx = vertices[0];
        foreach (i; 0 .. N) {
            auto x_new = vtx.x.dup;
            x_new[i] += dx[i];
            vertices[i+1].x[] = x_new[];
            vertices[i+1].fx = f(x_new); nfe += 1;
        }
        vertices.sort!(less);
        nrestarts += 1;
        return;
    }

    T[2] f_statistics()
    // Returns mean and standard deviation of the vertex fn values.
    {
        T mean = 0.0;
        foreach (vtx; vertices) { mean += vtx.fx; }
        mean /= vertices.length;
        T ss = 0.0;
        foreach (vtx; vertices) { ss += (vtx.fx - mean)^^2; }
        T std_dev = sqrt(ss/(vertices.length-1));
        return [mean, std_dev];
    }

    T[] centroid(size_t imax)
    // Returns the centroid of the subset of vertices up to and including imax.
    {
        T[] xmid = vertices[0].x.dup;
        imax = min(imax, N);
        foreach (i; 1 .. imax+1) { xmid[] += vertices[i].x[]; }
        foreach (ref elem; xmid) { elem /= (imax+1); }
        return xmid;
    }

    void contract_about_zero_point()
    // Contract the simplex about the vertex[0].
    {
        T[] x_con = vertices[0].x.dup;
        foreach (i; 1 .. N+1) {
            foreach (j, ref elem; vertices[i].x) { elem = 0.5*elem + 0.5*x_con[j]; }
            vertices[i].fx = f(vertices[i].x); nfe += 1;
        }
        vertices.sort!(less);
        return;
    }

    bool test_for_minimum(double delta)
    // Look around vertex 0 to see if it is a local minimum.
    // This is expensive, so we don't want to do it often.
    {
        bool is_minimum = true;
        auto f_min = vertices[0].fx;
        foreach (j; 0 .. N) {
            auto x_new = vertices[j].x.dup;
            x_new[j] += dx[j] * delta;
            auto f_new = f(x_new); nfe += 1;
            if (f_new < f_min) { is_minimum = false; break; }
            x_new[j] -= dx[j] * delta * 2.0;
            f_new = f(x_new); nfe += 1;
            if (f_new < f_min) { is_minimum = false; break; }
        }
        return is_minimum;
    }

    bool replace_vertex(size_t i)
    // Try to replace the worst point, i, in the simplex.
    // Returns true is there was a successful replacement.
    {
        auto f_min = vertices[0].fx;
        if (i <= (N-P)) {
            throw new Exception(format("i=%d seems not to be in the high points", i));
        }
        auto x_high = vertices[i].x.dup;
        auto f_high = vertices[i].fx;
        // Centroid of simplex excluding point(s) that we are replacing.
        auto x_mid = centroid(N-P);
        //
        // First, try moving away from worst point by
        // reflection through centroid.
        auto x_refl = x_mid.dup;
        foreach (j, ref e; x_refl) { e += Kreflect*(e - x_high[j]); }
        auto f_refl = f(x_refl); nfe += 1;
        if (f_refl < f_min) {
            // The reflection through the centroid is good,
            // try to extend in the same direction.
            auto x_ext = x_mid.dup;
            foreach (j, ref e; x_ext) { e += Kextend*(x_refl[j] - x_mid[j]); }
            auto f_ext = f(x_ext); nfe += 1;
            if (f_ext < f_refl) {
                // Keep the extension because it's best.
                vertices[i].x = x_ext[]; vertices[i].fx = f_ext;
                return true;
            } else {
                // Settle for the original reflection.
                vertices[i].x = x_refl[]; vertices[i].fx = f_refl;
                return true;
            }
        } else {
            // The reflection is not going in the right direction, it seems.
            // See how many vertices are worse than the reflected point.
            size_t count = 0;
            foreach (j; 0 .. N+1) { if (vertices[j].fx > f_refl) count += 1; }
            if (count <= 1) {
                // Not too many points are higher than the original reflection.
                // Try a contraction on the reflection-side of the centroid.
                auto x_con = x_mid.dup;
                foreach (j, ref e; x_con) { e = (1.0-Kcontract)*e + Kcontract*x_high[j]; }
                auto f_con = f(x_con); nfe += 1;
                if (f_con < f_high) {
                    // At least we haven't gone uphill; accept.
                    vertices[i].x = x_con[]; vertices[i].fx = f_con;
                    return true;
                }
            } else {
                // Retain the original reflection because there are many
                // original vertices with higher values of the objective function.
                vertices[i].x = x_refl[]; vertices[i].fx = f_refl;
                return true;
            }
        }
        // If we arrive here, we have not replaced the highest point.
        return false;
    } // end replace_vertex()

public:
    // In Python style, we have the following attributes of the simplex public
    // because we want to be able to access them from the minimize function.
    T[] dx;
    size_t N;
    size_t P;
    double Kreflect;
    double Kextend;
    double Kcontract;
    int nfe;
    int nrestarts;
    Vertex!T[] vertices;
} // end class NMMinimizer

//-----------------------------------------------------------------------

version(nelmin_test) {
    import std.math;
    import util.msg_service;
    import nm.number;
    int main() {
        // Test objective function 1.
        //   x is expected to be a list of coordinates.
        //   Returns a single float value.
        number testFunction1(number[] x)
        {
            number sum = 0.0;
            foreach (elem; x) sum += (elem - 1.0) * (elem - 1.0);
            return sum;
        }
        number[] x = [to!number(0.0), to!number(0.0), to!number(0.0), to!number(0.0)];
        number[] dx;
        number fx;
        int nfe, nres;
        bool conv_flag = minimize!(testFunction1,number)(x, fx, nfe, nres, dx);
        // writeln("x = ", x, " fx = ", fx);
        // writeln("convergence-flag = ", conv_flag);
        // writeln("number-of-fn-evaluations = ", nfe);
        // writeln("number-of-restarts = ", nres);
        assert(conv_flag, failedUnitTest());
        assert(approxEqualNumbers(x, [to!number(1.0), to!number(1.0),
                                      to!number(1.0), to!number(1.0)]),
               failedUnitTest());
        assert(approxEqualNumbers(fx, to!number(0.0)), failedUnitTest());
        assert((nfe > 425) && (nfe < 435), failedUnitTest()); // 428
        assert(nres < 6, failedUnitTest()); // 5

        // Test objective function 2.
        // Example 3.3 from Olsson and Nelson.
        number testFunction2(number[] x)
        {
            number x1 = x[0]; number x2 = x[1];   // rename to match the paper
            if ( (x1 * x1 + x2 * x2) > 1.0 ) {
                return to!number(1.0e38);
            } else {
                number yp = 53.69 + 7.26 * x1 - 10.33 * x2 + 7.22 * x1 * x1
                    + 6.43 * x2 * x2 + 11.36 * x1 * x2;
                number ys = 82.17 - 1.01 * x1 - 8.61 * x2 + 1.40 * x1 * x1
                    - 8.76 * x2 * x2 - 7.20 * x1 * x2;
                return -yp + fabs(ys - 87.8);
            }
        }
        number[] x2 = [to!number(0.0), to!number(0.0)];
        number[] dx2 = [to!number(0.5), to!number(0.5)];
        conv_flag = minimize!(testFunction2,number)(x2, fx, nfe, nres, dx2, 1.0e-4);
        // writeln("x = ", x2, " fx = ", fx);
        // writeln("convergence-flag = ", conv_flag);
        // writeln("number-of-fn-evaluations = ", nfe);
        // writeln("number-of-restarts =", nres);
        assert(conv_flag, failedUnitTest());
        assert(approxEqualNumbers(x2, [to!number(0.811295), to!number(-0.584636)]),
               failedUnitTest());
        assert(approxEqualNumbers(fx, to!number(-67.1077)), failedUnitTest());
        assert((nfe > 80) && (nfe < 90), failedUnitTest()); // 86
        assert(nres == 0, failedUnitTest());

        // Test objective function 3.
        // Example 3.5 from Olsson and Nelson; least-squares.
        double testFunction3(double[] z)
        {
            double[] x = [0.25, 0.50, 1.00, 1.70, 2.00, 4.00];
            double[] y = [0.25, 0.40, 0.60, 0.58, 0.54, 0.27];
            double a1 = z[0]; double a2 = z[1];
            double alpha1 = z[2]; double alpha2 = z[3];
            double sum_residuals = 0.0;
            foreach (i; 0 .. 6) {
                double t = x[i];
                double eta = a1 * exp(alpha1 * t) + a2 * exp(alpha2 * t);
                double r = y[i] - eta;
                sum_residuals += r * r;
            }
            return sum_residuals;
        }
        double[] x3 = [1.0, 1.0, -0.5, -2.5];
        double[] dx3 = [0.1, 0.1, 0.1, 0.1];
        double fx3;
        conv_flag = minimize!(testFunction3,double)(x3, fx3, nfe, nres, dx3, 1.0e-9, 2, 800);
        // writeln("x = ", x3, " fx = ", fx3);
        // writeln("convergence-flag = ", conv_flag);
        // writeln("number-of-fn-evaluations = ", nfe);
        // writeln("number-of-restarts = ", nres);
        assert(conv_flag, failedUnitTest());
        assert(approxEqualNumbers(x3, [1.80105, -1.84177, -0.463388, -1.20508]), failedUnitTest());
        assert(approxEqualNumbers(fx3, 0.000908953), failedUnitTest());
        assert((nfe > 500) && (nfe < 505), failedUnitTest()); // 503
        assert(nres == 0, failedUnitTest());

        return 0;
    }
} // end version(nelmin_test)
