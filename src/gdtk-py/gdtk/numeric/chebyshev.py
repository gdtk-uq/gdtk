# chebyshev.py
"""
Simple implementation of polynomial interpolation using Chebyshev basis.
PeterJ, MECH3750, October 2011
2015-08-10 Python3 version
2022-11-26 Port to GDTK project
"""

from gdtk.numeric.polynomial import fit_interpolating_polynomial, eval_polynomial

def make_chebyshev_basis(n):
    """
    Returns a list of functions T[0]..T[n-1].
    After T[0] and T[1], each function is built on top of previous 2.

    """
    flist = [lambda x: 1.0, lambda x: x] # T[0] and T[1]
    for j in range(2,n):
        def T(x, j=j, Tlist=flist):
            return 2 * x * Tlist[j-1](x) - Tlist[j-2](x)
        flist.append(T)
    return flist

if __name__ == '__main__':
    print("Begin Chebyshev interpolation demo...")
    xi = [3.2, 2.7, 1.0, 4.8, 5.6]
    yi = [22.0, 17.8, 14.2, 5.6, 51.7]
    chebyshev_basis = make_chebyshev_basis(len(xi))
    x = 3.0
    alpha, cond = fit_interpolating_polynomial(chebyshev_basis, xi, yi)
    print('alpha=', alpha, '\ncondition number=', cond)
    print('p(',x,')=', eval_polynomial(x, chebyshev_basis, alpha))
    print("Done")
