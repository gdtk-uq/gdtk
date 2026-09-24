# nelmin_test.py
# Derived from the self-test code in nelmin.py

from gdtk.numeric.nelmin import minimize, NelderMeadMinimizer
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
    print("Begin nelmin test...")
    def pretty_print(r):
        print("  x=", r.x)
        print("  fx=", r.fun)
        print("  convergence-flag=", r.success)
        print("  number-of-fn-evaluations=", r.nfe)
        print("  number-of-restarts=", r.nrestarts)
        print("  vertices=", [str(v) for v in r.vertices])
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
                      options={'tol':1.0e-9, 'maxfe':800})
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
    print("test 1 repeated with more-direct control of stepping")
    print("  simple quadratic with zero at (1,1,...)")
    nmm1 = NelderMeadMinimizer(test_fun_1, 3)
    nmm1.build_initial_simplex([0.0, 0.0, 0.0])
    nmm1.dump_simplex('test_1_simplex.json')
    print("  Start a new minimizer.")
    nmm2 = NelderMeadMinimizer(test_fun_1, 3)
    nmm2.load_simplex('test_1_simplex.json')
    print("  Initial simplex, loaded from JSON file.")
    for v in nmm2.vertices: print("    x=", v.x, "f=", v.f)
    nmm2.take_steps(20)
    print("  Simplex after 20 steps.")
    for v in nmm2.vertices: print("    x=", v.x, "f=", v.f)
    nmm2.take_steps(20)
    print("  Simplex after a further 20 steps.")
    for v in nmm2.vertices: print("    x=", v.x, "f=", v.f)
    #
    print("Done.")
