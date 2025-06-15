"""
A module to generate M-flow.

M-flow is an internal, axisymmetric, conical flow field described by
Moelder in 1967. It turns free stream flow inwardly through
a conical shock.

References:

Moelder, S. (1967)
Internal, Axisymmetric, Conical Flow.
AIAA Journal, 5(7), pp. 1252--1255

.. Author: RJG
.. Version: 2025-06-15
"""

from gdtk.ideal_gas_flow import beta_obl, M2_obl, theta_obl
import numpy as np
from scipy.integrate import solve_ivp
from math import pi, radians, cos, sin, tan, degrees, atan2

def cot(x): return 1/tan(x)

class MFlow():

    def __init__(self, M1, delta, gamma=1.4):
        """Initialise an M-flow calculation."""
        self._M1 = M1
        self._delta = delta
        self._gamma = gamma
        # M2 is flow behind weak shock with deflection delta
        self._beta = beta_obl(M1, delta, gamma)
        self._M2 = M2_obl(M1, self._beta, delta, gamma)
        self._theta2 = pi - self._beta
        self._u2 = -self._M2*cos(self._beta - delta)
        self._v2 = -self._M2*sin(self._beta - delta)

    def xy_from_rtheta(self, r, theta):
        return r*np.cos(theta), r*np.sin(theta)

    def dydx(self, theta, u, v):
        numer = (u/v)*np.sin(theta) + np.cos(theta)
        denom = (u/v)*np.cos(theta) - np.sin(theta)
        return numer/denom

    def generate_contour(self, r2=1.0, max_dtheta=2.0e-4*pi, rtol=1.0e-3, atol=1.0e-3):
        """Generate the M-flow contour by integration."""

        def fODE(theta, Y):
            u, v, r = Y
            tmp = (u + v*cot(-theta))/(v*v - 1)
            gm1_2 = (self._gamma - 1)/2
            du_dtheta = (v + gm1_2*u*v*tmp)
            dv_dtheta = -u + (1 + gm1_2*v*v)*tmp
            dr_dtheta = r*u/v
            return np.array([-du_dtheta, -dv_dtheta, -dr_dtheta])

        def event(theta, Y):
            u, v, r = Y
            return v + 1
        event.terminal = True

        theta = self._theta2
        Y = np.array([self._u2, self._v2, r2])
        
        self._ode_soln = solve_ivp(fODE, (-theta, -theta+pi/4), Y, method='DOP853', dense_output=True, events=event, max_step=max_dtheta, rtol=rtol, atol=atol)

        self._thetas = -self._ode_soln.t.copy()
        self._us = self._ode_soln.y[0,:].copy()
        self._vs = self._ode_soln.y[1,:].copy()
        self._rs = self._ode_soln.y[2,:].copy()
        self._Ms = np.sqrt(self._us**2 + self._vs**2)
        self._xs, self._ys = self.xy_from_rtheta(self._rs, self._thetas)
        self._dydxs = self.dydx(self._thetas, self._us, self._vs)
        return

    def write_contour(self, filename, x_at_zero=True, scale=1.0):
        x_start = self._xs[0] if x_at_zero else 0.0
        with open(filename, "w") as f:
            f.write(f"# x\ty\n")
            for i in range(len(self._xs)):
                f.write(f"{scale*(self._xs[i] - x_start)}\t{scale*self._ys[i]}\n")

        return

    def write_mach_distribution(self, filename, x_at_zero=True, scale=1.0):
        x_start = self._xs[0] if x_at_zero else 0.0
        with open(filename, "w") as f:
            f.write(f"x\tMach\n")
            for i in range(len(self._xs)):
                f.write(f"{scale*(self._xs[i] - x_start)}\t{self._Ms[i]}\n")

        return

if __name__ == '__main__':

    M1 = 5.0
    theta_s = radians(140)
    beta = pi - theta_s
    delta = theta_obl(M1, beta)
    mFlow = MFlow(M1, delta)
    mFlow.generate_contour()
    mFlow.write_contour('m-flow-contour.dat')
    mFlow.write_mach_distribution('mach-distribution.dat')
