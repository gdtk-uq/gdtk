# least_squares.py
"""
Least-squares model fitting with user-specified basis.
The model function is of the form:
p(x) = alpha[0]*phi[0](x) + alpha[1]*phi[1](x) + ... + alpha[m-1]*phi[m-1](x)

PeterJ, mech2700 demo code, 09-October-2013
        mech3750 adaption, 2014-07-27
        2022-11-26 Port to the GDTK project
"""

import numpy as np

def eval_model(x, phi, alpha):
    """
    Evaluates the model function p(x).
    Input:
        x: point at which to evaluate model (may be an array if
           the user-supplied functions are suitable.)
        phi: sequence of functions
        alpha: sequence of coefficients in the polynomial
    """
    m = len(alpha); assert len(phi) == m
    p = 0.0
    for j in range(m):
        p += alpha[j] * phi[j](x)
    return p

def fit_model(phi, xdata, ydata):
    """
    Determine the coefficients to fit a data set in least-squares sense.
    Input:
        phi : sequence of basis functions defining model
        xdata : sequence of x-coordinates
        ydata : sequence of y-coordinates
    Returns:
        the single-dimensioned array of coefficients.
    """
    m = len(phi); N = len(xdata); assert len(ydata) == N
    assert m < N
    vdm = np.zeros((m,m), float)
    rhs = np.zeros((m,), float)
    for j in range(m): # for each constraint equation
        for k in range(m):
            for i in range(N):
                vdm[j,k] += phi[j](xdata[i]) * phi[k](xdata[i])
        for i in range(N):
            rhs[j] += ydata[i] * phi[j](xdata[i])
    return np.linalg.solve(vdm,rhs), np.linalg.cond(vdm)


if __name__ == '__main__':
    print("Start general least-squares demo...")
    T_degreeC  = [ 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0]
    reading_mV = [0.01, 0.12, 0.24, 0.38, 0.51,  0.67,  0.84,  1.01,  1.15,  1.31]
    if 1:
        m = 3
        model_basis = [lambda x, k=j: x**k for j in range(m)]
    else:
        def x0(x): return 1.0
        def x1(x): return x
        def x2(x): return x*x
        model_basis = [x0, x1, x2]
    poly_degree = len(model_basis) - 1
    alpha, cond = fit_model(model_basis, reading_mV, T_degreeC)
    print('alpha=', alpha, '\ncondition number=', cond)
    modelx = np.linspace(0.0, 1.32, 100)
    modely = eval_model(modelx, model_basis, alpha)
    import pylab as plt
    plt.plot(reading_mV,T_degreeC,'o', modelx, modely, '-')
    plt.title("Thermocouple data and model polynomial of degree %d"%(poly_degree,))
    plt.xlabel("Thermocouple voltage, mV")
    plt.ylabel("Temperature, degreeC")
    plt.show()
    print("Done.")
