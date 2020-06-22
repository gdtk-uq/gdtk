#! /bin/env python3
"""
nelmin.py: Nelder-Mead simplex minimization of a nonlinear (multivariate) function.

This code has been adapted from the C-coded nelmin.c which was
adapted from the Fortran-coded nelmin.f which was, in turn, adapted
from the papers:

    J.A. Nelder and R. Mead (1965)
    A simplex method for function minimization.
    Computer Journal, Volume 7, pp 308-313.

    R. O'Neill (1971)
    Algorithm AS47. Function minimization using a simplex algorithm.
    Applied Statistics, Volume 20, pp 338-345.

and some examples are in:

   D.M. Olsson and L.S. Nelson (1975)
   The Nelder-Mead Simplex procedure for function minimization.
   Technometrics, Volume 17 No. 1, pp 45-51.

For a fairly recent and popular incarnation of this minimizer,
see the amoeba function in the famous "Numerical Recipes" text.
The programming interface is via the minimize() function; see below.

.. Author: PA Jacobs, School of Engineering, The University of Queensland

.. Version:
   07-Jan-04 Python2 flavour for the cfcfd project.
   2020-06-22 Port to Python3 for the DGD project.

Example transcript::

    $ python ~/e3bin/cfpylib/nm/nelmin.py
    Begin nelmin self-test...
    ---------------------------------------------------
    test 1: simple quadratic with zero at (1,1,...)
    x= [1.0000000651818255, 1.000000021516589, 0.9999999925111813]
    fx= 4.76771637998e-15
    convergence-flag= 0
    number-of-fn-evaluations= 300
    number-of-restarts= 3
    ---------------------------------------------------
    test 2: Example 3.3 in Olsson and Nelson f(0.811,-0.585)=-67.1
    x= [0.8112948625421268, -0.5846355980866065]
    fx= -67.1077410608
    convergence-flag= 1
    number-of-fn-evaluations= 127
    number-of-restarts= 0
    ---------------------------------------------------
    test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares
    f(1.801, -1.842, -0.463, -1.205)=0.0009
    x= [1.8010249070374442, -1.8417283432511073, -0.46337704853342615, -1.205053720578973]
    fx= 0.000908952916125
    convergence-flag= 1
    number-of-fn-evaluations= 618
    number-of-restarts= 0
    ---------------------------------------------------
    Done.
"""

from copy import copy
from math import sqrt

#-----------------------------------------------------------------------
# The public face of the minimizer...

def minimize(f, x, dx=None, tol=1.0e-6,
             maxfe=300, n_check=20, delta=0.001,
             Kreflect=1.0, Kextend=2.0, Kcontract=0.5):
    """
    Locate a minimum of the objective function, f.

    :param f: user-specified function f(x)
    :param x: list of N coordinates
    :param dx: list of N increments to apply to x when forming
               the initial simplex.  Their magnitudes determine the size
               and shape of the initial simplex.
    :param tol: the terminating limit for the standard-deviation
                of the simplex function values.
    :param maxfe: maximum number of function evaluations that we will allow
    :param n_check: number of steps between convergence checks
    :param delta: magnitude of the perturbations for checking a local minimum
                  and for the scale reduction when restarting
    :param Kreflect:
    :param Kextend:
    :param Kcontract: coefficients for locating the new vertex

    :Returns: a tuple consisting of::

        [0] a list of coordinates for the best x location, corresponding to min(f(x)),
        [1] the function value at that point,
        [2] a flag to indicate if convergence was achieved
        [3] the number of function evaluations and
        [4] the number of restarts (with scale reduction)
    """
    converged = 0
    N = len(x)
    if dx == None:
        dx = [0.1] * N
    smplx = NMSimplex(x, dx, f)

    while (not converged) and (smplx.nfe < maxfe):
        # Take some steps and then check for convergence.
        for i in range(n_check):
            take_a_step(smplx, Kreflect, Kextend, Kcontract)
        # Pick out the current best vertex.
        i_best = smplx.lowest()
        x_best = copy(smplx.get_vertex(i_best))
        f_best = smplx.f_list[i_best]
        # Check the scatter of vertex values to see if we are
        # close enough to call it quits.
        mean, stddev = smplx.f_statistics()
        # Append some statistics to a log file as we go along the optimisation
        # loop, so that the user knows that the optimiser is still running.
        fo = open('nelmin.log','a')
        fo.write('x = [')
        for i in range(len(x_best)):
            fo.write('%f, ' % x_best[i])
        fo.write(']\n')
        fo.write('Number of function evaluations = %d \n' % smplx.nfe)
        fo.write('Target standard deviation  = %f \n' % tol)
        fo.write('Current standard deviation = %f \n' % stddev)
        fo.write('---------------------------------------------------\n')
        fo.close()
        #
        if stddev < tol:
            # All of the points are close together but we need to
            # test more carefully.
            converged = smplx.test_for_minimum(i_best, delta)
            if not converged:
                # The function evaluations are all very close together
                # but we are not at a true minimum; rescale the simplex.
                smplx.rescale(delta)
    return x_best, f_best, converged, smplx.nfe, smplx.nrestarts


#-----------------------------------------------------------------------
# Use a class to keep the data tidy and conveniently accessible...

class NMSimplex:
    """
    Stores the (nonlinear) simplex as a list of lists.

    In an N-dimensional problem, each vertex is a list of N coordinates
    and the simplex consists of N+1 vertices.
    """
    def __init__(self, x, dx, f):
        """
        Initialize the simplex.

        Set up the vertices about the user-specified vertex, x,
        and the set of step-sizes dx.
        f is a user-specified objective function f(x).
        """
        self.N = len(x)
        self.vertex_list = []
        self.f_list = []
        self.dx = copy(dx)
        self.f = f
        self.nfe = 0
        self.nrestarts = 0
        for i in range(self.N + 1):
            p = copy(x)
            if i >= 1: p[i-1] += dx[i-1]
            self.vertex_list.append(p)
            self.f_list.append(f(p))
            self.nfe += 1

    def rescale(self, ratio):
        """
        Pick out the current minimum and rebuild the simplex about that point.
        """
        i_min = self.lowest()
        for i in range(self.N):
            self.dx[i] *= ratio
        x = self.get_vertex(i_min)
        self.vertex_list = []
        self.f_list = []
        for i in range(self.N + 1):
            p = copy(x)
            if i >= 1: p[i-1] += self.dx[i-1]
            self.vertex_list.append(p)
            self.f_list.append(self.f(p))
            self.nfe += 1
        self.nrestarts += 1
        return

    def get_vertex(self, i):
        return copy(self.vertex_list[i])

    def replace_vertex(self, i, x, fvalue):
        self.vertex_list[i] = copy(x)
        self.f_list[i] = fvalue
        return

    def lowest(self, exclude=-1):
        """
        Returns the index of the lowest vertex, excluding the one specified.
        """
        if exclude == 0:
            indx = 1
        else:
            indx = 0
        lowest_f_value = self.f_list[indx]
        for i in range(self.N + 1):
            if i == exclude: continue
            if self.f_list[i] < lowest_f_value:
                lowest_f_value = self.f_list[i]
                indx = i
        return indx

    def highest(self, exclude=-1):
        """
        Returns the index of the highest vertex, excluding the one specified.
        """
        if exclude == 0:
            indx = 1
        else:
            indx = 0
        highest_f_value = self.f_list[indx]
        for i in range(self.N + 1):
            if i == exclude: continue
            if self.f_list[i] > highest_f_value:
                highest_f_value = self.f_list[i]
                indx = i
        return indx

    def f_statistics(self):
        """
        Returns mean and standard deviation of the vertex fn values.
        """
        sum = 0.0
        for i in range(self.N + 1):
            sum += self.f_list[i]
        mean = sum / (self.N + 1)
        sum = 0.0
        for i in range(self.N +1):
            diff = self.f_list[i] - mean
            sum += diff * diff
        std_dev = sqrt(sum / self.N)
        return mean, std_dev

    def centroid(self, exclude=-1):
        """
        Returns the centroid of all vertices excluding the one specified.
        """
        xmid = [0.0]*self.N
        for i in range(self.N + 1):
            if i == exclude: continue
            for j in range(self.N):
                xmid[j] += self.vertex_list[i][j]
        for j in range(self.N):
            xmid[j] /= self.N
        return xmid

    def contract_about_one_point(self, i_con):
        """
        Contract the simplex about the vertex i_con.
        """
        p_con = self.vertex_list[i_con]
        for i in range(self.N + 1):
            if i == i_con: continue
            p = self.vertex_list[i]
            for j in range(self.N):
                p[j] = 0.5 * (p[j] + p_con[j])
            self.f_list[i] = self.f(p)
            self.nfe += 1
        return

    def test_for_minimum(self, i_min, delta):
        """
        Perturb the minimum vertex and check that it is a local minimum.
        """
        is_minimum = 1  # Assume it is true and test for failure.
        f_min = self.f_list[i_min]
        for j in range(self.N):
            # Check either side of the minimum, perturbing one
            # coordinate at a time.
            p = self.get_vertex(i_min)
            p[j] += self.dx[j] * delta
            f_p = self.f(p)
            self.nfe += 1
            if f_p < f_min:
                is_minimum = 0
                break
            p[j] -= self.dx[j] * delta * 2
            f_p = self.f(p)
            self.nfe += 1
            if f_p < f_min:
                is_minimum = 0
                break
        return is_minimum

#--------------------------------------------------------------------
# The core of the minimizer...

def take_a_step(smplx, Kreflect, Kextend, Kcontract):
    """
    Try to move away from the worst point in the simplex.

    The new point will be inserted into the simplex (in place).
    """
    i_low = smplx.lowest()
    i_high = smplx.highest()
    x_high = smplx.vertex_list[i_high]
    f_high = smplx.f_list[i_high]
    # Centroid of simplex excluding worst point.
    x_mid = smplx.centroid(i_high)
    f_mid = smplx.f(x_mid)
    smplx.nfe += 1

    # First, try moving away from worst point by
    # reflection through centroid
    x_refl = create_new_point(1.0+Kreflect, x_mid, -Kreflect, x_high)
    f_refl = smplx.f(x_refl)
    smplx.nfe += 1
    if f_refl < f_mid:
        # The reflection through the centroid is good,
        # try to extend in the same direction.
        x_ext = create_new_point(Kextend, x_refl, 1.0-Kextend, x_mid)
        f_ext = smplx.f(x_ext)
        smplx.nfe += 1
        if f_ext < f_refl:
            # Keep the extension because it's best.
            smplx.replace_vertex(i_high, x_ext, f_ext)
        else:
            # Settle for the original reflection.
            smplx.replace_vertex(i_high, x_refl, f_refl)
    else:
        # The reflection is not going in the right direction, it seems.
        # See how many vertices are better than the reflected point.
        count = 0
        for i in range(smplx.N+1):
            if smplx.f_list[i] > f_refl: count += 1
        if count <= 1:
            # Not too many points are higher than the original reflection.
            # Try a contraction on the reflection-side of the centroid.
            x_con = create_new_point(1.0-Kcontract, x_mid, Kcontract, x_high)
            f_con = smplx.f(x_con)
            smplx.nfe += 1
            if f_con < f_high:
                # At least we haven't gone uphill; accept.
                smplx.replace_vertex(i_high, x_con, f_con)
            else:
                # We have not been successful in taking a single step.
                # Contract the simplex about the current lowest point.
                smplx.contract_about_one_point(i_low)
        else:
            # Retain the original reflection because there are many
            # vertices with higher values of the objective function.
            smplx.replace_vertex(i_high, x_refl, f_refl)
    return

def create_new_point(c1, p1, c2, p2):
    """
    Create a new N-dimensional point as a weighting of points p1 and p2.
    """
    p_new = []
    for j in range(len(p1)):
        p_new.append(c1 * p1[j] + c2 * p2[j])
    return p_new

#--------------------------------------------------------------------

def test_fun_1(x):
    """
    Test objective function 1.

    x is expected to be a list of coordinates.
    Returns a single float value.
    """
    n = len(x)
    sum = 0.0
    for i in range(n):
        sum += (x[i] - 1.0) * (x[i] - 1.0)
    return sum

def test_fun_2(x):
    """
    Test objective function 2.

    Example 3.3 from Olsson and Nelson.
    """
    x1, x2 = x   # rename to match the paper
    if (x1 * x1 + x2 * x2) > 1.0:
        return 1.0e38
    else:
        yp = 53.69 + 7.26 * x1 - 10.33 * x2 + 7.22 * x1 * x1 \
             + 6.43 * x2 * x2 + 11.36 * x1 * x2
        ys = 82.17 - 1.01 * x1 - 8.61 * x2 + 1.40 * x1 * x1 \
             - 8.76 * x2 * x2 - 7.20 * x1 * x2
        return -yp + abs(ys - 87.8)

from math import exp

def test_fun_3(z):
    """
    Test objective function 3.

    Example 3.5 from Olsson and Nelson; least-squares.
    """
    x = [0.25, 0.50, 1.00, 1.70, 2.00, 4.00]
    y = [0.25, 0.40, 0.60, 0.58, 0.54, 0.27]
    a1, a2, alpha1, alpha2 = z
    sum_residuals = 0.0
    for i in range(len(x)):
        t = x[i]
        eta = a1 * exp(alpha1 * t) + a2 * exp(alpha2 * t)
        r = y[i] - eta
        sum_residuals += r * r
    return sum_residuals

#--------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin nelmin self-test...")

    print("---------------------------------------------------")
    print("test 1: simple quadratic with zero at (1,1,...)")
    x, fx, conv_flag, nfe, nres = minimize(test_fun_1, [0.0, 0.0, 0.0])
    print("x=", x)
    print("fx=", fx)
    print("convergence-flag=", conv_flag)
    print("number-of-fn-evaluations=", nfe)
    print("number-of-restarts=", nres)

    print("---------------------------------------------------")
    print("test 2: Example 3.3 in Olsson and Nelson f(0.811,-0.585)=-67.1")
    x, fx, conv_flag, nfe, nres = minimize(test_fun_2,
                                           [0.0, 0.0], [0.5, 0.5],
                                           1.0e-4)
    print("x=", x)
    print("fx=", fx)
    print("convergence-flag=", conv_flag)
    print("number-of-fn-evaluations=", nfe)
    print("number-of-restarts=", nres)

    print("---------------------------------------------------")
    print("test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares")
    print("f(1.801, -1.842, -0.463, -1.205)=0.0009")
    x, fx, conv_flag, nfe, nres = minimize(test_fun_3,
                                           [1.0, 1.0, -0.5, -2.5],
                                           [0.1, 0.1, 0.1, 0.1],
                                           1.0e-9, 800)
    print("x=", x)
    print("fx=", fx)
    print("convergence-flag=", conv_flag)
    print("number-of-fn-evaluations=", nfe)
    print("number-of-restarts=", nres)

    print("---------------------------------------------------")
    print("Done.")
