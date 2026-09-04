"""
A module to generate a Prandtl-Meyer contour.

Two-dimensional flow is isentropically compressed, and redirected with
the freestream through a terminating planar shock.

The implementation of Prandtl-Meyer theory here is based on the
description in Moradian and Timofeev (2012), and Moradian and
Timofeev (2018).

References:

N. Moradian and E. Timofeev. (2018)
Starting Characteristics of Prandtl-Meyer Scramjet Intakes
with Overboard Spillage. 
Journal of Propulsion and Power, 34(1):189-197

N. Moradian and E. Timofeev. (2012)
Limiting Contractions for Starting Prandtl-Meyer-Type Scramjet
Inlets with Overboard Spillage.
28th International Symposium on Shock Waves, pages 307-312, Berlin, Heidelberg.

Author: Brock Duffy
Date: August 2024

Update: RJG, 2025-02-10
        Some clean up to integrate with gdtk python library
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from gdtk.numeric.spline import CubicSpline
from gdtk.ideal_gas_flow import PM1, theta_obl


class PrandtlMeyerInlet():
    """Prandtl-Meyer contour."""

    def __init__(self, M1, M2, gamma=1.4):
        """
        M1: initial Mach number of the flow.
        M2: final Mach number of the flow after isentropic compression.
        gamma: ratio of specific heats for gas
        """
        # Chosen parameters.
        self.M1 = M1
        self.M2 = M2
        self.gamma = gamma

        self.focal = [0,0]

        Ai = 1 # [m^2], inlet capture area (assumes width into the page, W = 1m).
        Si = self.M1 * Ai
        alphai = np.arcsin(Ai/Si)
        Xi = Si * np.cos(alphai)
        
        self.Si = Si
        self.Ai = Ai

    def generate_contour(self, number_points=200):
        """Generate a Prandtl-Meyer contour and populate
        the x and y coordinate arrays."""
        self.x = []
        self.y = []
        Ms = np.linspace(self.M1, self.M2, number_points-1) # because we append a point later on
        
        for i in range(len(Ms)):
            delta_nu = PM1(self.M1, self.gamma) - PM1(Ms[i], self.gamma)
            alpha = np.arcsin(1/Ms[i])
            beta = np.pi/2 - alpha
            arg = beta + (np.pi/2 - delta_nu)
            S = self.Si * self.f(Ms[i]) / self.f(self.M1)
            
            self.x.append(S*np.cos(arg))
            self.y.append(S*np.sin(arg))
        
        # Append straight section.
        As = self.Ai*self.f(self.M2)*self.M1 / (self.f(self.M1)*self.M2)
        sigma2 = self.compute_beta(delta_nu, self.M2) # Oblique shock angle at cowl.
        L = As / np.sin(sigma2)
        self.x.extend(L*np.cos(sigma2-delta_nu))
        self.y.extend(L*np.sin(sigma2-delta_nu))

        # Shift contour to the positive xy-plane.
        self.x = list(self.x - self.x[0])

        # Generate gradients.
        self._compute_dydx()
    
    def _compute_dydx(self):
        """Compute gradient dy/dx."""
        self.dydx = []
        for i in range(len(self.y)):
            if i == 0: # One sided difference.
                y_curr = self.y[i]
                y_next = self.y[i+1]
                dx = self.x[i+1] - self.x[i]
                self.dydx.append((y_next - y_curr) / dx)
            elif i == len(self.y)-1: # One sided difference.
                y_prev = self.y[i-1]
                y_curr = self.y[i]
                dx = self.x[i] - self.x[i-1]
                self.dydx.append((y_curr - y_prev) / dx)
            else: # Central difference.
                y_prev = self.y[i-1]
                y_next = self.y[i+1]
                dx = self.x[i+1] - self.x[i-1]
                self.dydx.append((y_next - y_prev) / dx)
   
    def f(self, mach):
        """Used to calculate the radial distance, S,
        from the focal point to the contour."""
        k = (self.gamma - 1) / (self.gamma + 1)
        pwr = 1 / (2*k)
        return (k*(mach**2 - 1) + 1)**pwr
    
    def compute_beta(self, theta, mach):
        """Return beta, the oblique shock angle which turns a M=mach
        flow by an angle of theta."""
        betaguess = theta
        def f_zero(beta):
            return theta_obl(mach, beta, self.gamma) - theta
        return fsolve(f_zero, betaguess)

    def write_contour(self, filename, capture_area=1.0, aspect_ratio=1):
        """Write contour data to a file."""
        assert self.x, "Prandtl-Meyer contour not yet computed."
        y_scale = (capture_area/(4*aspect_ratio))**0.5
        x_scale = y_scale
        with open(filename, 'w') as f:
            f.write("x\ty\n")
            for i in range(len(self.x)):
                f.write(f"{x_scale*self.x[i]:.6e}\t{y_scale*self.y[i]:.6e}\n")

    def contour_as_spline(self, capture_area=1.0, aspect_ratio=1):
        assert self.x, "Prandtl-Meyer contour not yet computed."
        y_scale = (capture_area/(4*aspect_ratio))**0.5
        x_scale = y_scale
        xs = []
        ys = []
        for x, y in zip(self.x, self.y):
            xs.append(x_scale*x)
            ys.append(y_scale*y)
        return CubicSpline(xs, ys)

if __name__ == "__main__":

    def plot_contour(inlet):
        """Plot the Prandtl-Meyer contour."""
        assert inlet.x, "Prandtl-Meyer contour not yet computed."
        fig, ax = plt.subplots(figsize=(16,6))
        ax.plot(inlet.x, inlet.y, 'k')
        #ax.plot(self.x, self.y, 'k') # Basic plotting.
        ax.set_title('Prandtl-Meyer Contour')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid()
        plt.axis('equal')

        fig.tight_layout()
        ### fig.savefig('plot_prandtl_meyer_contour.png', dpi=1000)

        plt.show()
    # The below is used for testing purposes.
    M1 = 8.33
    M2 = 6.20

    pm = PrandtlMeyerInlet(M1, M2)
    pm.generate_contour(number_points=1000)
    plot_contour(pm)
