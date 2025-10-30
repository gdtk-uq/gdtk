# polyomial.py
"""
Simple implementation of polynomial interpolation using supplied basis.
The model equation is
    p(x) = a[0]*phi[0](x) + a[1]*phi[1](x) + ... + a[n-1]*phi[n-1](x)

PeterJ, MECH3750, October 2011
2015-08-10 Python3 version
2022-11-26 Ported to the GDTK project.
"""

def eval_polynomial(x, phi, a):
    """
    Simple evaluation of polynomial p(x).
    """
    n = len(a)
    assert len(phi) == n
    p = 0.0
    for j in range(n):
        p += a[j] * phi[j](x)
    return p

from numpy import array, zeros, linalg

def fit_interpolating_polynomial(phi, xi, yi):
    """
    Determine the polynomial coefficients to interpolate a data set.

    phi : sequence of basis functions
    xi : sequence of x-coordinates
    yi : sequence of y-coordinates

    Returns the single-dimensioned array of coefficients.
    """
    n = len(xi); assert n == len(yi)
    vdm = zeros((n,n), float)
    for i in range(n):
        for j in range(n):
            vdm[i,j] = phi[j](xi[i])
    return linalg.solve(vdm,yi), linalg.cond(vdm)

#----------------------------------------------------

if __name__ == '__main__':
    print("Start monomial demo...")
    print("Fit the coefficients.")
    xi = [3.2, 2.7, 1.0, 4.8, 5.6]
    yi = [22.0, 17.8, 14.2, 5.6, 51.7]
    n = len(xi)
    monomial_basis = [lambda x, k=j: x**k for j in range(n)]
    alpha, cond = fit_interpolating_polynomial(monomial_basis, xi, yi)
    print('alpha=', alpha, '\ncondition number=', cond)
    x = 3.0
    print('p(',x,')=', eval_polynomial(x, monomial_basis, alpha))
    #
    from pylab import linspace, plot, title, xlabel, ylabel, show
    x_plot = linspace(min(xi), max(xi), 100)
    y_plot = [eval_polynomial(x, monomial_basis, alpha) for x in x_plot]
    plot(x_plot, y_plot, '-')
    plot(xi, yi, 'o')
    title("Interpolatory polynomial with order n-1=%d" % (n-1,))
    xlabel('x'); ylabel('y')
    show()
    print("Done.")
