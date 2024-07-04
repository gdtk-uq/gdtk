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
"""
import numpy as np
from math import sin, cos, tan, pi, sqrt, asin, ceil
from gdtk.ideal_gas_flow import theta_obl, M2_obl, p02_p01_obl
from gdtk.numeric.ode import rkf45_step
from gdtk.numeric.spline import CubicSpline
from collections import namedtuple

# Properties are:
#   M1 : Mach number at entrance to diffuser, free-stream Mach number
#   M2 : Mach number at end of isentropic compression, BEFORE terminating shock
#   M3 : Mach number BEHIND terminating shock, Mach number of flow entering throat
#   Pi : total pressure recovery, ratio of pt3/pt2
BDProperties = namedtuple('BDProperties', ['M1', 'M2', 'M3', 'Pi'])

def cot(x): return 1/tan(x)


class BusemannDiffuser():

    def __init__(self, M2, theta_23, gamma=1.4):
        """Initialise a Busemann diffuser integration from M2 and theta_23."""
        self._M2 = M2
        self._theta_23 = theta_23
        self._gamma = 1.4
        self._delta = theta_obl(M2, theta_23, gamma)
        self._theta2 = theta_23 - self._delta
        self._u2 = M2*cos(theta_23)
        self._v2 = -M2*sin(theta_23)
        self._M3 = M2_obl(M2, theta_23, self._delta, gamma)
        self._Pi = p02_p01_obl(M2, theta_23, gamma)
        self._M1 = None # we can't set this until we've integrated the flow field
        return

    def properties(self):
        return BDProperties(self._M1, self._M2, self._M3, self._Pi)
         
    
    def _xy_from_rtheta(self, r, theta):
        return r*cos(theta), r*sin(theta)

    def generate_contour(self, r2=1.0, dtheta=pi/180.0):
        """Generate the Busemann diffuser contour by integration."""

        def fODE(theta, Y, n):
            u, v, r = Y
            tmp = (u + v*cot(theta))/(v*v - 1)
            gm1_2 = (self._gamma - 1)/2
            du_dtheta = v + gm1_2*u*v*tmp
            dv_dtheta = -u + (1 + gm1_2*v*v)*tmp
            dr_dtheta = r*u/v
            return np.array([du_dtheta, dv_dtheta, dr_dtheta])

        # Initialise storage for the integration values and derived values
        self._us = [self._u2]
        self._vs = [self._v2]
        self._Ms = [self._M2]
        self._rs = [r2]
        self._thetas = [self._theta2]
        x, y = self._xy_from_rtheta(r2, self._theta2)
        self._xs = [x]
        self._ys = [y]

        # Set initial conditions
        n = 3
        theta = self._theta2
        Y = np.array([self._u2, self._v2, r2])
        while (self._vs[-1] < -1.0):
            theta, Y, err = rkf45_step(theta, dtheta, fODE, n, Y)
            self._thetas.append(theta)
            self._us.append(Y[0])
            self._vs.append(Y[1])
            self._rs.append(Y[2])
            M = sqrt(Y[0]**2 + Y[1]**2)
            self._Ms.append(M)
            x, y = self._xy_from_rtheta(self._rs[-1], self._thetas[-1])
            self._xs.append(x)
            self._ys.append(y)

        # Remove final point.
        self._thetas.pop()
        self._us.pop()
        self._vs.pop()
        self._rs.pop()
        self._Ms.pop()
        self._xs.pop()
        self._ys.pop()

        # Set M1 now that we've completed integration
        self._M1 = self._Ms[-1]
            
        return

    def write_contour(self, filename, number_points, scale=1.0):
        assert self._xs, "Busemann contour not yet computed."
        idxs = self._collect_sample_indices(number_points)
        with open(filename, "w") as f:
            f.write("x\ty\n")
            for i in idxs:
                f.write(f"{scale*self._xs[i]}\t{scale*self._ys[i]}\n")

        return

    def contour_as_spline(self, number_points, scale=1.0):
        assert self._xs, "Busemann contour not yet computed."
        idxs = self._collect_sample_indices(number_points)
        xs = []
        ys = []
        for i in idxs:
            xs.append(scale*self._xs[i])
            ys.append(scale*self._ys[i])
        return CubicSpline(xs, ys)

    def _collect_sample_indices(self, number_points):
        total_points = len(self._xs)
        n_total_spaces = total_points - 1
        n_sample_spaces = (number_points - 2) + 1 # remove end points (-2); count spaces (+1)
        step = ceil(n_total_spaces/n_sample_spaces)
        idxs = [0] + list(range(step, total_points, step)) + [total_points-1]
        return idxs[::-1] # reversed
    
    
                
            
        
    

        
            
            
        
        


