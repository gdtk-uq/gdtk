"""
A module to generate a Busemann diffuser.

The original idea of an axially-symmetric coverging diffuser for supersonic flow
is due to Adolph Busemann (1944).
The diffuser has a single terminating conical shock.
This leads to a highly efficient compression field because the shock to turn
the flow occurs where Mach number is lowest.

The implementation here borrows from the description in Moelder (2019) and
an earlier work of his: Moelder (2003).

References:

Busemann, A. (1944)
Die achsensymmetrische kegelige ueberschallstroemung.
Luftfahrtforschung. 19(4):137--144

Moelder, S. (2019)
The Busemann Air Intake for Hypersonic Speeds.
Chapter 5 in: Hypersonic Vehicles --- Past, Present and Future Developments.
Ed. Pezzella and Viviani, IntechOpen.

Moelder, S. (2003)
A Benchmark for Internal Flow CFD Codes
Computational Fluid Dynamics Journal, 12(2):47, pp.408--414

.. Author: RJG (and helpful discussions with Reece Otto)

.. Version: 2024-07-03
            2025-02-15  Use scpiy ODE integrator
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from math import sin, cos, tan, pi, sqrt, asin, fabs, atan2, degrees
from gdtk.ideal_gas_flow import (
    theta_obl,
    M2_obl,
    p02_p01_obl,
    p2_p1_obl,
    T2_T1_obl,
    p0_p,
    T0_T,
)
from gdtk.numeric.spline import CubicSpline
from collections import namedtuple

# Properties are:
#   M1 : Mach number at entrance to diffuser, free-stream Mach number
#   M2 : Mach number at end of isentropic compression, BEFORE terminating shock
#   M3 : Mach number BEHIND terminating shock, Mach number of flow entering throat
#   Pi : total pressure recovery, ratio of pt3/pt2
#   PR : pressure ratio
#   TR : temperature ratio
#   CR : contraction ratio
BDProperties = namedtuple("BDProperties", ["M1", "M2", "M3", "Pi", "PR", "TR", "CR"])


def cot(x):
    return 1 / tan(x)


class BusemannDiffuser:
    def __init__(self, M2, theta_23, gamma=1.4):
        """Initialise a Busemann diffuser integration from M2 and theta_23."""
        self._M2 = M2
        self._theta_23 = theta_23
        self._gamma = 1.4
        self._delta = theta_obl(M2, theta_23, gamma)
        self._theta2 = theta_23 - self._delta
        self._u2 = M2 * cos(theta_23)
        self._v2 = -M2 * sin(theta_23)
        self._M3 = M2_obl(M2, theta_23, self._delta, gamma)
        self._Pi = p02_p01_obl(M2, theta_23, gamma)
        self._M1 = None  # we can't set this until we've integrated the flow field
        self._ode_soln = None
        return

    @property
    def theta2(self):
        return self._theta2

    @property
    def theta23(self):
        return self._theta_23

    @property
    def M2(self):
        return self._M2

    @property
    def start(self):
        assert self._M1, "Busemann contour not yet computed"
        return self._xs[-1], self._ys[-1]

    @property
    def end(self):
        assert self._M1, "Busemann contour not yet computed"
        return self._xs[0], self._ys[0]

    @property
    def properties(self):
        # compute ratios in isentopic portion
        p2_p1 = p0_p(self._M1, self._gamma) / p0_p(self._M2, self._gamma)
        T2_T1 = T0_T(self._M1, self._gamma) / T0_T(self._M2, self._gamma)
        # compute ratios across terminating shock
        p3_p2 = p2_p1_obl(self._M2, self._theta_23, self._gamma)
        T3_T2 = T2_T1_obl(self._M2, self._theta_23, self._gamma)
        # compute ratios from capture to throat
        PR = p3_p2 * p2_p1
        TR = T3_T2 * T2_T1
        CR = self._ys[-1] ** 2 / self._ys[0] ** 2
        return BDProperties(self._M1, self._M2, self._M3, self._Pi, PR, TR, CR)

    def xy_from_rtheta(self, r, theta):
        return r * np.cos(theta), r * np.sin(theta)

    def dydx(self, theta, u, v):
        numer = (u / v) * np.sin(theta) + np.cos(theta)
        denom = (u / v) * np.cos(theta) - np.sin(theta)
        return numer / denom

    def ode_soln(self, theta):
        return self._ode_soln.sol(theta)

    @property
    def last_characteristic(self):
        return float(self._ode_soln.t_events[0][0])

    def theta_from_x(self, x):
        def fZERO(theta):
            u, v, r = self.ode_soln(theta)
            x_test, y = self.xy_from_rtheta(r, theta)
            return x - x_test

        theta0 = self._thetas[-1]
        theta1 = self._thetas[0]
        sol = root_scalar(fZERO, bracket=[theta0, theta1], method="brentq")
        return sol.root

    def generate_contour(self, r2=1.0, max_dtheta=0.01 / pi, rtol=1.0e-3, atol=1.0e-6):
        """Generate the Busemann diffuser contour by integration."""

        def fODE(theta, Y):
            u, v, r = Y
            tmp = (u + v * cot(theta)) / (v * v - 1)
            gm1_2 = (self._gamma - 1) / 2
            du_dtheta = v + gm1_2 * u * v * tmp
            dv_dtheta = -u + (1 + gm1_2 * v * v) * tmp
            dr_dtheta = r * u / v
            return np.array([du_dtheta, dv_dtheta, dr_dtheta])

        def last_characteristic(theta, Y):
            u, v, r = Y
            M = sqrt(u * u + v * v)
            mu = asin(1 / M)
            vx = u * cos(theta) - r * v * sin(theta)
            vy = u * sin(theta) + r * v * cos(theta)
            delta = atan2(vy, vx)
            return theta - (pi - (fabs(mu) + fabs(delta)))

        def singularity(theta, Y):
            u, v, r = Y
            return u * sin(theta) + v * cos(theta)

        singularity.terminal = True

        # Set initial conditions
        theta = self._theta2
        Y = np.array([self._u2, self._v2, r2])
        self._ode_soln = solve_ivp(
            fODE,
            (theta, pi),
            Y,
            method="DOP853",
            dense_output=True,
            events=[last_characteristic, singularity],
            max_step=max_dtheta,
            rtol=rtol,
            atol=atol,
        )

        self._thetas = self._ode_soln.t.copy()
        self._us = self._ode_soln.y[0, :].copy()
        self._vs = self._ode_soln.y[1, :].copy()
        self._rs = self._ode_soln.y[2, :].copy()
        self._Ms = np.sqrt(self._us**2 + self._vs**2)
        self._xs, self._ys = self.xy_from_rtheta(self._rs, self._thetas)
        self._dydxs = self.dydx(self._thetas, self._us, self._vs)

        # Set M1 now that we've completed integration
        self._M1 = self._Ms[-1]

        return

    def write_contour(self, filename, number_points, scale=1.0):
        assert self._M1, "Busemann contour not yet computed."
        idxs = self._collect_sample_indices(number_points)
        with open(filename, "w") as f:
            f.write("x\ty\n")
            for i in idxs:
                f.write(f"{scale * self._xs[i]}\t{scale * self._ys[i]}\n")

        return

    def write_wall_properties(self, filename, number_points=-1, scale=1.0):
        assert self._M1, "Busemann contour not yet computed."
        if number_points > 0:
            idxs = self._collect_sample_indices(number_points)
        else:
            idxs = range(len(self._xs))
        with open(filename, "w") as f:
            f.write("x\ty\tdydx\ttheta\tr\tu\tM\n")
            for i in idxs:
                f.write(
                    f"{scale * self._xs[i]}\t"
                    f"{scale * self._ys[i]}\t"
                    f"{self._dydxs[i]}\t"
                    f"{self._thetas[i]}\t"
                    f"{scale * self._rs[i]}\t"
                    f"{self._us[i]}\t"
                    f"{self._Ms[i]}\n"
                )
        return

    def contour_as_spline(self, number_points, scale=1.0):
        assert self._M1, "Busemann contour not yet computed."
        idxs = self._collect_sample_indices(number_points)
        xs = []
        ys = []
        for i in idxs:
            xs.append(scale * self._xs[i])
            ys.append(scale * self._ys[i])
        return CubicSpline(xs, ys)

    def _collect_sample_indices(self, number_points):
        total_points = len(self._xs)
        idxs = np.round(np.linspace(0, total_points - 1, number_points)).astype(int)
        return idxs[::-1]  # reversed


if __name__ == "__main__":
    M2 = 3.0
    k = 1.4
    theta_23 = asin(k / M2)
    bd = BusemannDiffuser(M2, theta_23)
    bd.generate_contour()

    print(bd.properties)
    print(f"last centred characteristic at {degrees(bd.last_characteristic)} deg")
