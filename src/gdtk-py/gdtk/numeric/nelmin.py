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

For the 2020 update, we make some of the stepping closer to
the description given in the paper:

    Donghoon Lee and Matthew Wiswall (2007)
    A parallel implementation of the simplec function minimization routine.
    Computational Economics 30:171-187

Some examples are in:

   D.M. Olsson and L.S. Nelson (1975)
   The Nelder-Mead Simplex procedure for function minimization.
   Technometrics, Volume 17 No. 1, pp 45-51.

For a fairly recent and popular incarnation of this minimizer,
see the amoeba function in the famous "Numerical Recipes" text.
The programming interface is via the minimize() function; see below.

Author:
   PA Jacobs, School of Engineering, The University of Queensland

Version:
   2004-01-07 Python2 flavour for the cfcfd project.
   2020-06-22 Port to Python3 for the DGD project.
   2020-06-25 Concurrent evaluation of the candidate points.
   2021-06-07 Dan Smith added option to read the initial simplex.

Example transcript::

$ python3 nelmin.py
Begin nelmin self-test...
test 1: simple quadratic with zero at (1,1,...)
  x= [0.99999759305115321, 0.99999586570882093, 0.99999901474880748]
  fx= 2.38564862168e-11
  convergence-flag= True
  number-of-fn-evaluations= 202
  number-of-restarts= 2
------------------------------------------------------------
test 2: Example 3.3 in Olsson and Nelson f(0.811,-0.585)=-67.1
  x= [0.81129486254212679, -0.58463559808660648]
  fx= -67.1077410608
  convergence-flag= True
  number-of-fn-evaluations= 86
  number-of-restarts= 0
------------------------------------------------------------
test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares, P=1, n_workers=1
  f(1.801, -1.842, -0.463, -1.205)=0.0009
  calculation time (with 0.1 sec sleeps)= 54.83996343612671
  x= [1.8007773284762463, -1.8414976281967217, -0.46335920559222765, -1.2051840585794604]
  fx= 0.000908952987659
  convergence-flag= True
  number-of-fn-evaluations= 547
  number-of-restarts= 0
------------------------------------------------------------
test 4: Example 3.5 in Olsson and Nelson, nonlinear least-squares, P=2, n_workers=2
  f(1.801, -1.842, -0.463, -1.205)=0.0009
  calculation time (with 0.1 sec sleeps)= 26.823304653167725
  x= [1.8011390988049902, -1.841849174052093, -0.46339847883903756, -1.205042800777345]
  fx= 0.000908952812173
  convergence-flag= True
  number-of-fn-evaluations= 490
  number-of-restarts= 0
------------------------------------------------------------
Done.
"""

import math
import numpy
from collections import namedtuple
import concurrent.futures
import json

#-----------------------------------------------------------------------
# The public face of the simple minimizer...

def minimize(f, x, dx=None, options={}):
    """
    Locate a minimum of the objective function, f.

    f: user-specified function f(x)
    x: list of N coordinates
    dx: optional list of N increments to apply to x when forming the initial simplex.
        These increments determine the size and shape of the initial simplex.
    options: a dictionary with entries
        tol: (default 1.0e-6) the terminating limit for the standard-deviation
            of the simplex function values.
        P: (default 1) number of points to replace in parallel, each step.
        n_workers: (default 1) number of concurrent threads in pool
            We are using thread pool to do the concurrent evaluation of the objective
            function because typical use will involve running the Eilmer flow simulation
            code as a subprocess.  This will look like slow IO and we will get
            the benefit of running on multiple cores.
        maxfe: (default 300) maximum number of function evaluations that we will allow
            Note that this is a soft limit that will likely be exceeded by about 3N
            because N+1 evaluations are needed to initialize the simplex and then
            2N are required in the convergence check after taking a number of steps.
        n_check: (default 20) number of steps between convergence checks
        delta: (default 0.001) magnitude of the perturbations for checking a local minimum
            and for the scale reduction when restarting
        Kreflect: (default 1.0)
        Kextend: (default 2.0)
        Kcontract: (default 0.5) coefficients for locating the new vertex
        initial_simplex_fname: Filename directing to json file with desired initial simplex.
            Useful for restarting an optimisation calculation which was cut short for some reason.
        print_messages: (default False) If True, print messages as we go.

    Returns a namedtuple consisting of:
        x, a list of coordinates for the best x location, corresponding to min(f(x)),
        fun, the function value at that point,
        success, a flag to indicate if convergence was achieved
        nfe, the number of function evaluations and
        nrestarts, the number of restarts (with scale reduction)
        vertices, the full simplex of points and function values, useful for a restart
    """
    # Check input parameters.
    assert callable(f), "A function was expected for f."
    N = len(x)
    if dx is None: dx = numpy.array([0.1]*N)
    # Default values for options.
    opts = {'tol':1.0e-6, 'maxfe':300, 'n_check':20, 'delta':0.001,
            'Kreflect':1.0, 'Kextend':2.0, 'Kcontract':0.5,
            'P':1, 'n_workers':1,
            'initial_simplex_fname':None,
            'print_messages':False}
    # Overwrite default values.
    for k in options.keys():
        if k in opts.keys():
            opts[k] = options[k]
        else:
            raise RuntimeError(f'option {k} is not available')
    # Check that we have a consistent value for the number of steps between checks.
    opts['n_check'] = min(int(opts['maxfe']/opts['P']), opts['n_check'])
    # Start working.
    global workers
    if opts['n_workers'] > 1:
        workers = concurrent.futures.ThreadPoolExecutor(max_workers=opts['n_workers'])
    else:
        workers = None
    # Get to work.
    converged = False
    smplx = NelderMeadMinimizer(f, N, dx, opts['P'], opts['Kreflect'], opts['Kextend'], opts['Kcontract'])
    if opts['initial_simplex_fname'] is None:
        if opts['print_messages']: print('Build the initial simplex.')
        smplx.build_initial_simplex(x)
        # Save the initial version of the simplex in case something happens.
        smplx.dump_simplex('initial_simplex.json')
    else:
        if opts['print_messages']: print('Load the initial simplex from file.')
        smplx.load_simplex(opts['initial_simplex_fname'])
    while (not converged) and (smplx.nfe < opts['maxfe']):
        if opts['print_messages']: print('Take some steps with the simplex.')
        smplx.take_steps(opts['n_check'])
        # Save the latest version of the simplex in case something happens.
        smplx.dump_simplex('latest_simplex.json')
        x_best = list(smplx.vertices[0].x.copy())
        f_best = smplx.vertices[0].f
        # Check the scatter of vertex values to see if we are
        # close enough to call it quits.
        mean, stddev = smplx.f_statistics()
        if stddev < opts['tol']:
            # All of the points are close together but we need to test more carefully.
            converged = smplx.test_for_minimum(opts['delta'])
            if not converged: smplx.rescale(opts['delta'])
        # print("smplx.nfe=", smplx.nfe) # for debug
    #
    if workers:
        workers.shutdown()
        workers = None
    Result = namedtuple('Result', ['x', 'fun', 'success', 'nfe', 'nrestarts', 'vertices'])
    return Result(x_best, f_best, converged, smplx.nfe, smplx.nrestarts, smplx.vertices)


#-----------------------------------------------------------------------
# Use classes to keep the data tidy and conveniently accessible...

class Vertex:
    """
    Stores the coordinates as an array together with the associated function value.
    """
    __slots__ = 'x', 'f'

    def __init__(self, x, f):
        self.x = x.copy()
        self.f = f

    def __str__(self):
        return "Vertex(x=%s, f=%s)" % (self.x, self.f)


class NelderMeadMinimizer:
    """
    Stores the (nonlinear) simplex as a list of lists.

    In an N-dimensional problem, each vertex is a list of N coordinates
    and the simplex consists of N+1 vertices.
    """
    __slots__ = 'dx', 'f', 'N', 'P', 'vertices', 'nfe', 'nrestarts', \
                'Kreflect', 'Kextend', 'Kcontract'

    def __init__(self, f, N, dx=None, P=1, Kreflect=1.0, Kextend=2.0, Kcontract=0.5):
        """
        Initialize the minimizer.

        f is a user-specified objective function f(x).
        P is the number of points to be replaced in parallel, each step.
        """
        self.N = N
        self.P = P
        self.Kreflect = Kreflect
        self.Kextend = Kextend
        self.Kcontract = Kcontract
        self.dx = numpy.array(dx) if dx is not None else numpy.array([0.1]*N)
        assert callable(f), "A function was expected for f."
        self.f = f
        self.nfe = 0
        self.nrestarts = 0
        return

    def build_initial_simplex(self, x):
        """
        Set up the vertices about the user-specified vertex, x, using step-sizes dx.
        """
        if self.N != len(x):
            raise RuntimeError(f'N={self.N} does not match length of data {len(x)}')
        x = numpy.array(x) # since it may be a list or array
        xs = []
        for i in range(self.N + 1):
            x_new = x.copy()
            if i >= 1: x_new[i-1] += self.dx[i-1]
            xs.append(x_new)
        fs = self.evaluate_candidate_points(xs)
        self.vertices = [Vertex(item[0], item[1]) for item in zip(xs, fs)]
        self.vertices.sort(key=lambda v: v.f)
        return

    def set_simplex(self, vertex_list):
        """
        Sets the full simplex from the supplied list of Vertex objects.
        """
        self.vertices = [Vertex(v.x, v.f) for v in vertex_list]
        return

    def load_simplex(self, filename):
        """
        Load full simplex from a JSON file.
        """
        lines = open(filename,'r').readlines()
        dicts = [json.loads(line.strip()) for line in lines]
        self.vertices = [Vertex(d['x'], d['f']) for d in dicts]
        return

    def dump_simplex(self, filename):
        """
        Dump full simplex to a JSON file.
        """
        with open(filename,'w') as fp:
            for v in self.vertices:
                fp.write(json.dumps({'x':list(v.x), 'f':v.f})+'\n')
        return

    def take_steps(self, nsteps):
        """
        Take some steps, updating the simplex.
        On return, the best point is vertex[0].
        """
        global workers
        for istep in range(nsteps):
            self.vertices.sort(key=lambda v: v.f)
            if workers:
                # Concurrent replacements
                my_futures = [workers.submit(self.replace_vertex, self.N-i) for i in range(self.P)]
                success = [fut.result() for fut in my_futures]
            else:
                # Serial replacements
                success = [self.replace_vertex(self.N-i) for i in range(self.P)]
            if not any(success):
                # Contract the simplex about the current lowest point.
                self.contract_about_zero_point()
        self.vertices.sort(key=lambda v: v.f)
        return

    def evaluate_candidate_points(self, xs):
        """
        Evaluate the objective function for a list of candidate points.
        """
        global workers
        if workers:
            my_futures = [workers.submit(self.f, x) for x in xs]
            fs = [fut.result() for fut in my_futures]
        else:
            fs = [self.f(x) for x in xs]
        self.nfe += len(fs)
        return fs

    def rescale(self, ratio):
        """
        Rebuild the simplex about the lowest point for a restart.
        """
        self.vertices.sort(key=lambda v: v.f)
        self.dx *= ratio
        vtx = self.vertices[0]
        xs = []
        for i in range(self.N):
            x_new = vtx.x.copy()
            x_new[i] += self.dx[i]
            xs.append(x_new)
        fs = self.evaluate_candidate_points(xs)
        self.vertices = [vtx,] + [Vertex(item[0], item[1]) for item in zip(xs, fs)]
        self.nrestarts += 1
        self.vertices.sort(key=lambda v: v.f)
        return

    def f_statistics(self):
        """
        Returns mean and standard deviation of the vertex fn values.
        """
        f_list = [v.f for v in self.vertices]
        mean = sum(f_list)/(len(f_list))
        ss = sum([(v-mean)**2 for v in f_list])
        std_dev = math.sqrt(ss/(len(f_list)-1))
        return mean, std_dev

    def centroid(self, imax):
        """
        Returns the centroid of the subset of vertices up to and including imax.
        """
        xmid = numpy.array([0.0]*self.N)
        imax = min(imax, self.N)
        for i in range(imax+1): xmid += self.vertices[i].x
        xmid /= (imax+1)
        return xmid

    def contract_about_zero_point(self):
        """
        Contract the simplex about the vertex[0].
        """
        x_con = self.vertices[0].x
        xs = [0.5*self.vertices[i].x + 0.5*x_con for i in range(1, self.N+1)]
        fs = self.evaluate_candidate_points(xs)
        self.vertices = [self.vertices[0],] + [Vertex(item[0], item[1]) for item in zip(xs, fs)]
        self.vertices.sort(key=lambda v: v.f)
        return

    def test_for_minimum(self, delta):
        """
        Look around vertex 0 to see if it is a local minimum.

        This is expensive, so we don't want to do it often.
        """
        f_min = self.vertices[0].f
        xs = []
        for j in range(self.N):
            # Check either side of the candidate minimum,
            # perturbing one coordinate at a time.
            x_new = self.vertices[0].x.copy()
            x_new[j] += self.dx[j] * delta
            xs.append(x_new)
            x_new[j] -= self.dx[j] * delta * 2
            xs.append(x_new)
        fs = self.evaluate_candidate_points(xs)
        is_minimum = all([f >= f_min for f in fs])
        return is_minimum

    def replace_vertex(self, i):
        """
        Try to replace the worst point, i, in the simplex.

        Returns True is there was a successful replacement.
        """
        f_min = self.vertices[0].f
        assert i > (self.N-self.P), ("i=%d, seems not to be in the high points" % i)
        x_high = numpy.array(self.vertices[i].x.copy())
        f_high = self.vertices[i].f
        # Centroid of simplex excluding point(s) that we are replacing.
        x_mid = self.centroid(self.N-self.P)
        #
        # First, try moving away from worst point by
        # reflection through centroid.
        x_refl = x_mid + self.Kreflect*(x_mid - x_high)
        f_refl = self.f(x_refl); self.nfe += 1
        if f_refl < f_min:
            # The reflection through the centroid is good,
            # try to extend in the same direction.
            x_ext = x_mid + self.Kextend*(x_refl - x_mid)
            f_ext = self.f(x_ext); self.nfe += 1
            if f_ext < f_refl:
                # Keep the extension because it's best.
                self.vertices[i] = Vertex(x_ext, f_ext)
                return True
            else:
                # Settle for the original reflection.
                self.vertices[i] = Vertex(x_refl, f_refl)
                return True
        else:
            # The reflection is not going in the right direction, it seems.
            # See how many vertices are worse than the reflected point.
            count = 0
            for j in range(self.N+1):
                if self.vertices[j].f > f_refl: count += 1
            if count <= 1:
                # Not too many points are higher than the original reflection.
                # Try a contraction on the reflection-side of the centroid.
                x_con = (1.0-self.Kcontract)*x_mid + self.Kcontract*x_high
                f_con = self.f(x_con); self.nfe += 1
                if f_con < f_high:
                    # At least we haven't gone uphill; accept.
                    self.vertices[i] = Vertex(x_con, f_con)
                    return True
            else:
                # Retain the original reflection because there are many
                # original vertices with higher values of the objective function.
                self.vertices[i] = Vertex(x_refl, f_refl)
                return True
        #
        # If we arrive here, we have not replaced the highest point.
        return False

#--------------------------------------------------------------------
# Self test follows.

import time

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
    time.sleep(0.1) # To make this calculation seem expensive.
    return sum_residuals


if __name__ == '__main__':
    print("Begin nelmin self-test...")
    def pretty_print(r):
        print("  x=", r.x)
        print("  fx=", r.fun)
        print("  convergence-flag=", r.success)
        print("  number-of-fn-evaluations=", r.nfe)
        print("  number-of-restarts=", r.nrestarts)
        print(60*"-")
        return
    #
    print("test 1: simple quadratic with zero at (1,1,...)")
    result = minimize(test_fun_1, [0.0, 0.0, 0.0])
    pretty_print(result)
    #
    print("test 2: Example 3.3 in Olsson and Nelson f(0.811,-0.585)=-67.1")
    result = minimize(test_fun_2, [0.0, 0.0], [0.5, 0.5], options={'tol':1.0e-4})
    pretty_print(result)
    #
    print("test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares, P=1, n_workers=1")
    print("  f(1.801, -1.842, -0.463, -1.205)=0.0009")
    start_time = time.time()
    result = minimize(test_fun_3, [1.0, 1.0, -0.5, -2.5], [0.1, 0.1, 0.1, 0.1],
                      options={'tol':1.0e-9, 'maxfe':800, 'print_messages':True})
    print("  calculation time (with 0.1 sec sleeps)=", time.time()-start_time)
    pretty_print(result)
    #
    print("test 4: Example 3.5 in Olsson and Nelson, nonlinear least-squares, P=2, n_workers=2")
    print("  f(1.801, -1.842, -0.463, -1.205)=0.0009")
    start_time = time.time()
    result = minimize(test_fun_3, [1.0, 1.0, -0.5, -2.5], [0.1, 0.1, 0.1, 0.1],
                      options={'tol':1.0e-9, 'P':2, 'maxfe':800, 'n_workers':2})
    print("  calculation time (with 0.1 sec sleeps)=", time.time()-start_time)
    pretty_print(result)
    #
    print("Done.")
