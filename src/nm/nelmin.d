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
import nm.complex;

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
bool minimize(alias f, T)(ref T[] x,
                          out T f_min,
                          out int n_fe,
                          out int n_restart,
                          T[] dx,
                          in double tol=1.0e-6,
                          in int P=1,
                          in int max_fe=300,
                          in int n_check=20,
                          in double delta=0.001,
                          in double Kreflect=1.0,
                          in double Kextend=2.0,
                          in double Kcontract=0.5)
    if (is(typeof(f([0.0,0.0])) == double) ||
        is(typeof(f([Complex!double(0.0),Complex!double(0.0)])) == Complex!double))
{
    bool converged = false;
    n_fe = 0;
    n_restart = 0;
    size_t N = x.length;
    // Might have to check for NaNs as well.
    if ( dx.length != x.length ) {
        dx = x.dup;
        foreach (ref elem; dx) elem = 0.1;
    }
    // In an N-dimensional problem, each vertex is an array of N coordinates
    // and the simplex consists of N+1 vertices.
    // Set up the vertices about the user-specified vertex, x,
    // and the set of step-sizes dx.
    // f is a user-specified objective function f(x).
    T[][] vertex_list = new T[][N+1];
    T[] f_list = new T[N+1];
    foreach (i; 0 .. N+1) {
        auto p = x.dup;
        if ( i > 0 ) p[i-1] += dx[i-1];
        vertex_list[i] = p;
        f_list[i] = f(p);
        n_fe += 1;
    }

    // Utility functions for dealing with the simplex.

    // Returns the index of the lowest vertex, excluding the one specified.
    int lowest(int exclude=-1)
    {
        int indx;
        T lowest_f_value;
        if ( exclude == 0 ) indx = 1; else indx = 0;
        lowest_f_value = f_list[indx];
        for ( int i = 0; i <= N;  ++i ) {
            if ( i == exclude ) continue;
            if ( f_list[i] < lowest_f_value ) {
                lowest_f_value = f_list[i];
                indx = i;
            }
        }
        return indx;
    }

    // Returns the index of the highest vertex, excluding the one specified.
    int highest(int exclude=-1)
    {
        int indx;
        T highest_f_value;
        if ( exclude == 0 ) indx = 1; else indx = 0;
        highest_f_value = f_list[indx];
        for ( int i = 0; i <= N;  ++i ) {
            if ( i == exclude ) continue;
            if ( f_list[i] > highest_f_value ) {
                highest_f_value = f_list[i];
                indx = i;
            }
        }
        return indx;
    }

    // Returns the standard deviation of the vertex fn values.
    T std_dev()
    {
        int i;
        T sum = 0.0;
        for ( i = 0; i <= N; ++i ) {
            sum += f_list[i];
        }
        T mean = sum / (N + 1);
        sum = 0.0;
        for ( i = 0; i <= N; ++i ) {
            T diff = f_list[i] - mean;
            sum += diff * diff;
        }
        return sqrt(sum / N);
    }

    // Pick out the current minimum and rebuild the simplex about that point.
    void rescale(double ratio)
    {
        foreach ( idx; 0..dx.length) dx[idx] *= ratio;
        T[] p = vertex_list[lowest()].dup; // save to use below
        foreach (i; 0 .. N+1) {
            vertex_list[i][] = p.dup;
            if ( i >= 1 ) vertex_list[i][i-1] += dx[i-1];
            f_list[i] = f(vertex_list[i]);
            n_fe += 1;
        }
        return;
    }

    // Returns the centroid of all vertices excluding the one specified.
    T[] centroid(int exclude=-1)
    {
        T[] xmid = new T[N];
        foreach(ref elem; xmid) elem = 0.0;
        foreach (i; 0 .. N+1) {
            if (i == exclude ) continue;
            xmid[] += vertex_list[i][];
        }
        foreach( idx; 0..xmid.length) xmid[idx] /= N;
        return xmid;
    }

    // Contract the simplex about the vertex i_con.
    void contract_about_one_point(int i_con)
    {
        T[] p_con = vertex_list[i_con].dup;
        foreach (i; 0 .. N+1) {
            if ( i == i_con ) continue;
            T[] p = vertex_list[i].dup;
            foreach( idx; 0..p.length) p[idx] = 0.5 * (p[idx] + p_con[idx]);
            f_list[i] = f(p);
            n_fe += 1;
        }
        return;
    }

    // Perturb the minimum vertex and check that it is a local minimum.
    bool test_for_minimum(int i_min, double delta)
    {
        bool is_minimum = true;  // Assume it is true and test for failure.
        T f_min = f_list[i_min];
        T f_p;
        for ( int j = 0; j < N; ++j ) {
            // Check either side of the minimum, perturbing one coordinate at a time.
            T[] p = vertex_list[i_min].dup;
            p[j] += dx[j] * delta;
            f_p = f(p);
            n_fe += 1;
            if ( f_p < f_min ) { is_minimum = false; break; }
            p[j] -= dx[j] * delta * 2;
            f_p = f(p);
            n_fe += 1;
            if ( f_p < f_min ) { is_minimum = false; break; }
        }
        return is_minimum;
    }

    void replace_vertex(int i, T[] p, T fp)
    {
        vertex_list[i] = p.dup;
        f_list[i] = fp;
    }

    // The core of the minimizer...
    // The following (nested) functions require access to the simplex.
    void take_a_step()
    {
        // Try to move away from the worst point in the simplex.
        // The new point will be inserted into the simplex (in place).
        int i_low = lowest();
        int i_high = highest();
        T[] x_high = vertex_list[i_high].dup;
        T f_high = f_list[i_high];
        // Centroid of simplex excluding worst point.
        T[] x_mid = centroid(i_high);
        T f_mid = f(x_mid);
        n_fe += 1;

        // First, try moving away from worst point by
        // reflection through centroid
        T[] x_refl = x_mid.dup;
        foreach( idx; 0..x_refl.length) x_refl[idx] = (1.0+Kreflect) * x_mid[idx] - Kreflect * x_high[idx];
        T f_refl = f(x_refl);
        n_fe += 1;

        if ( f_refl < f_mid ) {
            // The reflection through the centroid is good,
            // try to extend in the same direction.
            T[] x_ext = x_mid.dup;
            foreach( idx; 0..x_ext.length) x_ext[idx] = Kextend * x_refl[idx] + (1.0-Kextend) * x_mid[idx];
            T f_ext = f(x_ext);
            n_fe += 1;
            if ( f_ext < f_refl ) {
                // Keep the extension because it's best.
                replace_vertex(i_high, x_ext, f_ext);
                return;
            } else {
                // Settle for the original reflection.
                replace_vertex(i_high, x_refl, f_refl);
                return;
            }
        } else {
            // The reflection is not going in the right direction, it seems.
            // See how many vertices are better than the reflected point.
            int count = 0;
            for ( int i = 0; i <= N; ++i ) {
                if ( f_list[i] > f_refl )  count += 1;
            }
            if ( count <= 1 ) {
                // Not too many points are higher than the original reflection.
                // Try a contraction on the reflection-side of the centroid.
                T[] x_con = x_mid.dup;
                foreach( idx; 0..x_con.length) x_con[idx] = (1.0-Kcontract) * x_mid[idx] + Kcontract * x_high[idx];
                T f_con = f(x_con);
                n_fe += 1;
                if ( f_con < f_high ) {
                    // At least we haven't gone uphill; accept.
                    replace_vertex(i_high, x_con, f_con);
                } else {
                    // We have not been successful in taking a single step.
                    // Contract the simplex about the current lowest point.
                    contract_about_one_point(i_low);
                }
                return;
            } else {
                // Retain the original reflection because there are many
                // vertices with higher values of the objective function.
                replace_vertex(i_high, x_refl, f_refl);
                return;
            }
        }
    } // end take_a_step()

    // Now, do some real work with all of these tools...

    // TODO make a better loop termination condition.
    while ( !converged && n_fe < (max_fe+N+2) ) {
        // Take some steps and then check for convergence.
        for ( int i = 0; i < n_check; ++i ) {
            take_a_step();
            // Pick out the current best vertex.
            int i_best = lowest();
            x[] = vertex_list[i_best][];
            f_min = f_list[i_best];
            // Check the scatter of vertex values to see if we are
            // close enough to call it quits.
            if ( std_dev() < tol ) {
                // All of the points are close together but we need to
                // test more carefully to see if we are at a true minimum.
                converged = test_for_minimum(i_best, delta);
                if ( !converged ) {
                    // The function evaluations are all very close together
                    // but we are not at a true minimum; rescale the simplex.
                    rescale(delta);
                    n_restart += 1;
                }
            }
        } // end for i
    } // end while !converged

    return converged;
} // end minimize()


version(nelmin_test) {
    import std.math;
    import util.msg_service;
    import nm.number;
    int main() {
        /** Test objective function 1.
         *  x is expected to be a list of coordinates.
         *  Returns a single float value.
         */
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
        assert(conv_flag, failedUnitTest());
        assert(approxEqualNumbers(x, [to!number(1.0), to!number(1.0), to!number(1.0), to!number(1.0)]), failedUnitTest());
        assert(approxEqualNumbers(fx, to!number(0.0)), failedUnitTest());
        assert((nfe > 355) && (nfe < 365), failedUnitTest()); // 361
        assert(nres < 6, failedUnitTest()); // 5
        // writeln("x = ", x, " fx = ", fx);
        // writeln("convergence-flag = ", conv_flag);
        // writeln("number-of-fn-evaluations = ", nfe);
        // writeln("number-of-restarts = ", nres);

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
        assert(conv_flag, failedUnitTest());
        assert(approxEqualNumbers(x2, [to!number(0.811295), to!number(-0.584636)]), failedUnitTest());
        assert(approxEqualNumbers(fx, to!number(-67.1077)), failedUnitTest());
        assert((nfe > 125) && (nfe < 135), failedUnitTest()); // 131
        assert(nres == 0, failedUnitTest());
        // writeln("x = ", x2, " fx = ", fx);
        // writeln("convergence-flag = ", conv_flag);
        // writeln("number-of-fn-evaluations = ", nfe);
        // writeln("number-of-restarts =", nres);

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
        conv_flag = minimize!(testFunction3,double)(x3, fx3, nfe, nres, dx3, 1.0e-9, 1, 800);
        assert(conv_flag, failedUnitTest());
        assert(approxEqualNumbers(x3, [1.80102, -1.84173, -0.463377, -1.20505]), failedUnitTest());
        assert(approxEqualNumbers(fx3, 0.000908953), failedUnitTest());
        assert((nfe > 715) && (nfe < 725), failedUnitTest()); // 722
        assert(nres == 0, failedUnitTest());
        // writeln("x = ", x3, " fx = ", fx3);
        // writeln("convergence-flag = ", conv_flag);
        // writeln("number-of-fn-evaluations = ", nfe);
        // writeln("number-of-restarts = ", nres);

        return 0;
    }
} // end version(nelmin_test)
