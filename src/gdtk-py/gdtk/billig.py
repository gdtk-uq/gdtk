"""
billig.py: Fred Billig's correlations for hypersonic shock-wave shapes.

These are a direct implementation of equations 5.36, 5.37 and 5.38
from J.D. Anderson's text Hypersonic and High Temperature Gas Dynamics

.. Author: PJ

.. Version: 19-June-2005
            2022-11-26 Port to the GDTK project and Python3
"""

from math import exp, sqrt, pow, tan, asin
from gdtk.ideal_gas_flow import beta_obl, beta_cone2

def delta_over_R(M, axi):
    """
    Calculates the normalised stand-off distance.
    """
    if axi:
        # Spherical nose
        d_R = 0.143 * exp(3.24/(M*M))
    else:
        # Cylindrical nose
        d_R = 0.386 * exp(4.67/(M*M))
    return d_R


def Rc_over_R(M, axi):
    """
    Calculates the normalised radius of curvature of the shock.
    """
    if axi:
        # Spherical nose
        Rc_R = 1.143 * exp(0.54/pow(M-1, 1.2))
    else:
        # Cylindrical nose
        Rc_R = 1.386 * exp(1.8/pow(M-1, 0.75))
    return Rc_R


def x_from_y(y, M, theta=0.0, axi=False, R_nose=1.0):
    """
    Determine the x-coordinate of a point on the shock wave.

    :param y: y-coordinate of the point on the shock wave
    :param M: free-stream Mach number
    :param theta: angle (in radians wrt free-stream direction)
                  of the downstream surface
    :param axi: (bool) axisymmetric flag:

                | == False : cylinder-wedge
                | == True  : sphere-cone

    :param R_nose: radius of the forebody (either cylinder or sphere)

    It is assumed that, for the ideal gas, gamma=1.4.
    That's the only value relevant to the data used for
    Billig's correlations.
    """
    Rc = R_nose * Rc_over_R(M, axi)
    d = R_nose * delta_over_R(M, axi)
    if theta == 0.0:
        beta = asin(1.0/M)
    else:
        beta = beta_cone2(M, theta) if axi else beta_obl(M, theta)
    tan_beta = tan(beta)
    cot_beta = 1.0/tan_beta
    x = R_nose + d - Rc * cot_beta**2 * (sqrt(1 + (y*tan_beta/Rc)**2) - 1)
    return x

def y_from_x(x, M, theta=0.0, axi=False, R_nose=1.0):
    """
    Determine the y-coordinate of a point on the shock wave.

    :param x: x-coordinate of the point on the shock wave
    :param M: free-stream Mach number
    :param theta: angle (in radians wrt free-stream direction)
                  of the downstream surface
    :param axi: (bool) axisymmetric flag:

                | == False : cylinder-wedge
                | == True  : sphere-cone

    :param R_nose: radius of the forebody (either cylinder or sphere)

    It is assumed that, for the ideal gas, gamma=1.4.
    That's the only value relevant to the data used for
    Billig's correlations.
    """
    Rc = R_nose * Rc_over_R(M, axi)
    d = R_nose * delta_over_R(M, axi)
    if axi:
        # Use shock angle on a cone
        beta = beta_cone2(M, theta)
    else:
        # Use shock angle on a wedge
        beta = beta_obl(M, theta)
    tan_beta = tan(beta)
    cot_beta = 1.0/tan_beta
    tmpA = (x - R_nose - d)/(-Rc * cot_beta**2) + 1
    y = sqrt( ((tmpA**2 - 1) * Rc**2) / (tan_beta**2) )
    return y

#------------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin demo of Billig's correlations.")
    print("Compare with Fig 5.31 in Anderson's text.")
    M_inf = 4.0
    for y in [0.0, 0.5, 1.0, 2.0]:
        print("x=", x_from_y(y, M_inf, 0.0, 1), "y=", y)
    print("Done.")
