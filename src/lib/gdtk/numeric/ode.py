"""
ode.py: Integrate a set of first-order ODEs.

Author: Peter J.

Version:
  2003-10-05 implementation with lists for storage
  2005-02-21 use Numeric arrays for storage and manipulation of the state data.
  2020-07-10 Changed to accumulating the sample data and adapted to Python3.

Running this module as a Python script gives me the following transcript::

Start sample integrations...
(1) Constant derivatives:
t1= 10.0
y1= [ 10. -20.  30.]
err1= [5.55111512e-17 1.11022302e-16 0.00000000e+00]
(2) Second-order linear ODE:
t1= 6.283185307179586
y1= [-5.64489861e-15  1.00000000e+00]
err1= [6.16715441e-11 6.16715459e-11]
Done.
"""

import math
import numpy as np

def ode_integrate(t0, tlast, nstep, f, n, y0):
    """
    Steps the set of ODEs until independent variable, t, reaches tlast.

    This function coordinates the work of integrating a system
    of first-order differential equations of the form:
      y' = f(t, y)
      y(t=t0) = y0
    The actual work is done by rkf45_step, a more specialised
    stepping function, that appears below.

    t0: is the starting value of the independent variable
    tlast: the desired finishing value for x
    nstep: number of steps to take to arrive at tlast
    f: a callable function that returns the derivative of y wrt t
      The signature of this function is f(t, y, n) where
      t is a float value, y is an array of float values
      and n is an integer specifying the number of equations.
    n: the number of dependent variables (in y)
    y0: an array of starting values for the dependent variables
      It is assumed that the y-elements are indexed 0..n-1

    Lists of t, y, and error estimates for y values are returned in a tuple.
    """
    assert callable(f)
    assert n <= len(y0)
    assert t0 <= tlast
    assert nstep >= 1
    ts = np.linspace(t0, tlast, nstep+1)
    ys = []
    ys.append(y0.copy())  # Set up a new list so we don't mangle y0 itself.
    err_sums = []
    err_sums.append(np.array([0.0]*n))
    for i in range(nstep):
        t, y, err = rkf45_step(ts[i], ts[i+1]-ts[i], f, n, ys[-1])
        ys.append(y)
        err_sums.append(err_sums[-1]+err)
    # Return the accumulated data: times, the set of y-values and an error estimate.
    return ts, ys, err_sums


def rkf45_step(t0, h, f, n, y0):
    """
    Single-step the set of ODEs by the Runge-Kutta-Fehlberg method.

    t0: is the starting value of the independent variable
    h: the requested step size
    f: a callable function that returns the derivative of y wrt t
      The signature of this function is f(t, y, n) where
      t is a float value, y is a list (or array) or float values
      and n is an integer specifying the number of equations.
    n: the number of dependent variables (in y)
    y0: an array of starting values for the dependent variables
      It is assumed that the y-elements are indexed 0..n-1

    Final values of t, y, and error estimates for y values are returned in a tuple.
    """
    # Build up the sample point information.
    k1 = f(t0, y0.copy(), n)
    k2 = f(t0+h/4.0, y0+0.25*h*k1, n)
    k3 = f(t0+3.0*h/8.0, y0+3.0*h*k1/32.0+9.0*h*k2/32.0, n)
    k4 = f(t0+12.0*h/13.0, y0+1932.0*h*k1/2197.0-7200.0*h*k2/2197.0+7296.0*h*k3/2197.0, n)
    k5 = f(t0+h, y0+439.0*h*k1/216.0-8.0*h*k2+3680.0*h*k3/513.0-845.0*h*k4/4104.0, n)
    k6 = f(t0+h/2.0, y0-8.0*h*k1/27.0+2.0*h*k2-3544.0*h*k3/2565.0+1859.0*h*k4/4104.0-11.0*h*k5/40.0, n)
    # Now, do the integration as a weighting of the sampled data.
    t1 = t0+h
    y1 = y0+16.0*h*k1/135.0+6656.0*h*k3/12825.0+28561.0*h*k4/56430.0-9.0*h*k5/50.0+2.0*h*k6/55.0
    err = abs(h*k1/360.0-128.0*h*k3/4275.0-2197.0*h*k4/75240.0+h*k5/50.0+2.0*h*k6/55.0)
    return t1, y1, err

#----------------------------------------------------------------------

if __name__ == "__main__":
    def f_sample_1(t, y, n):
        "Constant derivatives"
        assert n == 3
        return np.array([1.0, -2.0, 3.0])

    def f_sample_2(t, y, n):
        "Second-order linear ODE with sine solution"
        assert n == 2
        return np.array([y[1], -y[0]])

    print("Start sample integrations...")
    print("(1) Constant derivatives:")
    ts, ys, errs = ode_integrate(0.0, 10.0, 8, f_sample_1, 3, np.array([0.0]*3))
    print("t1=", ts[-1])
    print("y1=", ys[-1])
    print("err1=", errs[-1])
    assert all(np.isclose(ys[-1], np.array([10.0, -20.0, 30.0]))), "Constant derivative test"

    print("(2) Second-order linear ODE:")
    ts, ys, errs = ode_integrate(0.0, 2.0*math.pi, 600, f_sample_2, 2, np.array([0.0, 1.0]))
    print("t1=", ts[-1])
    print("y1=", ys[-1])
    print("err1=", errs[-1])
    assert all(np.isclose(ys[-1], np.array([0.0, 1.0]))), "Linear second-order ODE test"
    print("Done.")

