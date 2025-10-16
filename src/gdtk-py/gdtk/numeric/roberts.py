#! /bin/env python3
"""
roberts.py: Node distribution and coordinate stretching functions.

These functions should behave the same as the C code functions.

.. Author PA Jacobs

.. Version:
   2005-03-22: Python2 port from C code.
   2020-04-06: Python3 port for l1d4
   2022-11-08: Allow arrays of ordinate values to be transformed.
"""

import numpy as np

def roberts(eta, alpha, beta):
    """
    Computes the stretched coordinate in the range [0.0..1.0]
    using the boundary-layer-like transformation devised by Roberts.

    :param eta: unstretched coordinate, 0 <= eta <= 1.0
    :param beta: stretching factor (more stretching as beta --> 1.0)
    :param alpha: location of stretching:

        | alpha = 0.5: clustering of nodes at both extremes of eta
        | alpha = 0.0: nodes will be clustered near eta=1.0

    Works for both scalars and arrays.
    """
    lmbda = (beta + 1.0) / (beta - 1.0)
    lmbda = np.power(lmbda, ((eta - alpha)/(1.0 - alpha)))
    etabar = (beta + 2.0 * alpha) * lmbda - beta + 2.0 * alpha
    etabar = etabar / ((2.0 * alpha + 1.0) * (1.0 + lmbda))
    return etabar

def roberts_1(eta, end1, end2, beta):
    """
    Compute a stretched ordinate.

    eta: unstretched ordinate, scalar or array
    end1: (bool) cluster flag for end 1:
        False: points are not clustered to end 1
        True: points ARE clustered to end 1
    end2: (bool) cluster flag for end 2:
        False: points are not clustered to end 2
        True: points ARE clustered to end 2
    beta: grid stretching parameter:
        1 < beta < +inf : points are clustered
        The closer to 1, the more the clustering.
        beta < 1 for no clustering.

    Returns the stretched ordinate.
    """
    # Decide on stretching parameters for Robert's transform.
    alpha = 0.5
    reverse = False
    cluster = True
    if ((not end1) and (not end2)) or (beta < 1.0):
        cluster = False
    if end1 and end2:
        alpha = 0.5
    if end1 and (not end2):
        reverse = True
        alpha = 0.0
    if (not end1) and end2:
        reverse = False
        alpha = 0.0
    if cluster:
        if reverse: eta = 1.0 - eta
        etabar = roberts(eta, alpha, beta)
        if reverse: etabar = 1.0 - etabar
    else:
        etabar = eta
    return etabar

def distribute_points_1(t1, t2, n, end1, end2, beta):
    """
    Generate a set of n+1 points nonuniformly distributed from t1 to t2.

    t1: parameter value  1
    t2: parameter value  2
    n: number of intervals (the number of points is n+1)
    end1: (bool) cluster flag for end 1:
        False: points are not clustered to end 1
        True: points ARE clustered to end 1
    end2: (bool) cluster flag for end 2:
        False: points are not clustered to end 2
        True: points ARE clustered to end 2
    beta: grid stretching parameter:
        1 < beta < +inf : points are clustered
        The closer to 1, the more the clustering.
        beta < 1 for no clustering.

    Returns the array of n+1 distributed values.
    """
    # Compute the grid points as an array.
    etabar = roberts_1(np.linspace(0.0, 1.0, n+1), end1, end2, beta)
    # Compute the parameter value within the given end-points.
    x = (1.0 - etabar) * t1 + etabar * t2
    return x

#------------------------------------------------------------------

if __name__ == "__main__":
    print("Begin roberts demo...")
    a = 0.5
    b = 1.1
    print("Scalar use...")
    print("eta  roberts")
    for eta in np.arange(0.0, 1.0, 0.1):
        print(eta, roberts(eta, a, b))
    print("Vector use...")
    eta = np.arange(0.0, 1.0, 0.1)
    print("eta=", eta)
    print("roberts=", roberts(eta, a, b))
    print("Distribute points...")
    print("x=", distribute_points_1(0.0, 1.0, 5, True, False, 1.1))
    print("x=", distribute_points_1(0.0, 1.0, 5, False, True, 1.1))
    print("x=", distribute_points_1(0.0, 1.0, 5, True, True, 1.1))
    print("x=", distribute_points_1(0.0, 1.0, 5, False, False, 1.1))
    print("Done.")
