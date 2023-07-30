#! /usr/bin/env python3
# self_similar_bl.py
"""
Set up a boundary-layer flow profile.

The calculation is based on the self-similar solution for a flat plate,
as described in section 6.5 of J.D. Anderson's text
Hypersonic and High Temperature Gas Dynamics.

Usage:
Edit the parameters describing the flow conditions and the streamwise position
and then run this program to generate the InflowBoundaryLayerBC.data file
for use by the chicken flow solver.

Author: PJ
06-Mar-2011 Original code.
2022-11-26 Update for the GDTK project.
2023-07-29 Adjust for DJM's crossing-shock experiment
"""

import sys, os
import math, numpy
import gdtk.numeric.ode as ode
import gdtk.sutherland as sutherland
import gdtk.numeric.nelmin as nelmin

print("self_similar_bl: begin...")

print("Gas constants for ideal air")
R = 287.1  # J/kg.K
gamma = 1.4
C_p = gamma/(gamma-1.0) * R
Pr = 0.71  # to match Schetz' boundary-layer code
print("R=", R, "gamma=", gamma, "C_p=", C_p, "Pr=", Pr)

print("Conditions just outside of boundary layer to match T4 Mach 6 nozzle exit condition.")
p_e = 3531.5   # Pa
u_e = 2158.2   # m/s
T_e = 340.2    # degrees K
h_e = C_p * T_e
T_wall = 300.0
h_wall = C_p * T_wall
print("p_e=", p_e, "u_e=", u_e, "T_e=", T_e, "h_e=", h_e)
print("Condition at wall")
print("T_wall=", T_wall, "h_wall=", h_wall)

rho_e = p_e / (R * T_e) # ideal equation of state
mu_e = sutherland.mu(T_e,'Air')
k_e = mu_e * C_p / Pr
print("rho_e=", rho_e, "mu_e=", mu_e, "k_e=", k_e)

# Choose a position along the plate
x = 0.2  # metres
xi = rho_e * u_e * mu_e * x
print("x=", x, "xi=", xi)

def C(g):
    """
    Ratio of density.viscosity product at points in the boundary layer.
    """
    T = g * h_e / C_p
    T = max(T, 100.0) # to avoid difficult values
    rho = p_e / (R * T)
    mu = sutherland.mu(T,'Air')
    return rho*mu/(rho_e*mu_e)

def Cd(g, gd):
    """
    Finite-difference estimate of dC/deta.
    """
    deltag = 0.01 * g  # something not too big, not too small
    C0 = C(g)
    C1 = C(g+deltag)
    dCdg = (C1 - C0)/deltag
    return dCdg * gd

def odes(eta, z, n):
    """
    Functions defining the differential equations from Anderson's text.

    Elements of state vector z:
    f    stream function
    fd   normalized velocity df/deta = u/u_e
    fdd  d(fd)/deta
    g    normalized enthalpy h/h_e
    gd   dg/deta
    y    spatial coordinate through the boundary layer
    """
    assert len(z) == 6
    f, fd, fdd, g, gd, y = z
    fddd = 1.0/C(g) * (-f*fdd - Cd(g,gd)*fdd)
    gdd = Pr/C(g) * (-gd*(Cd(g,gd)/Pr+f) - C(g)*u_e*u_e/h_e*fdd*fdd)
    yd = math.sqrt(2*xi)/u_e * h_e/p_e * (gamma-1.0)/gamma * g
    dzdeta = numpy.array([fd, fdd, fddd, gd, gdd, yd])
    return dzdeta

def integrate_through_bl(fdd, gd):
    """
    Start at wall and integrate the differential equations through the BL.

    Input: fdd, gd: guesses for these elements
    Returns: eta and state vector values through the boundary layer (and beyond)
    """
    # Starting values at the wall.
    f = 0.0
    fd = 0.0
    g = h_wall/h_e
    y = 0.0
    z0 = numpy.array([f, fd, fdd, g, gd, y])
    # Integrate through the boundary layer to a large value of eta.
    eta = numpy.linspace(0.0, 5.0, 500)
    deta = eta[1] - eta[0]
    z_list = [z0]
    for i in range(1,len(eta)):
        eta1, z1, err = ode.rkf45_step(eta[i-1], deta, odes, 6, z_list[i-1])
        z_list.append(z1)
    return eta, z_list

def objective(params):
    """
    Evaluate the guessed values for fdd and gd by returning a measure
    of the error at the outer edge of the boundary layer.
    """
    assert len(params) == 2
    fdd, gd = params
    eta, z = integrate_through_bl(fdd, gd)
    f, fd, fdd, g, gd, y = z[-1]
    penalty_value = abs(fd - 1.0) + abs(g - 1.0)
    return penalty_value

if 1:
    # Determine the best initial values of fdd and gd.
    r = nelmin.minimize(objective, [0.5, 1.0])
    print("params=", r.x, "obj=", r.fun, "nfe=", r.nfe, "nrestart=", r.nrestarts)
    # Integrate this particular boundary layer.
    fdd, gd = r.x
else:
    # I happened to prepare this solution a little earlier...
    # with obj= 1.77e-08 nfe= 155 nrestart= 0
    fdd, gd = 0.4108115969934917, 2.1259086459379164
eta, z = integrate_through_bl(fdd, gd)

# Write the data for plotting (with gnuplot, maybe) and curve fitting.
fp = open('profile.data', 'w')
fp.write("# eta f fd fdd g gd y p T rho u\n")
ys = []; Ts = []; rhos = []; us = []
for i in range(len(eta)):
    f, fd, fdd, g, gd, y = z[i]
    h = g * h_e; T = h / C_p; rho = p_e / (R * T); u = fd * u_e
    ys.append(y); Ts.append(T); rhos.append(rho); us.append(u)
    fp.write("%f %f %f %f %f %f %f %f %f %f %f\n" %
             (eta[i], f, fd, fdd, g, gd, y, p_e, T, rho, u))
fp.close()

# Put a spline through a number of the sampled points.
from gdtk.numeric.spline import CubicSpline
y_sample = []; T_sample = []; rho_sample = []; u_sample = []
N = len(ys)
M = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
for i in M:
    j = int((N-1)*i)
    y_sample.append(ys[j]);
    T_sample.append(Ts[j])
    u_sample.append(us[j])
sT = CubicSpline(y_sample, T_sample)
su = CubicSpline(y_sample, u_sample)
#
# Write the spline data for use as a Chicken InflowBoundaryLayerBC
fp = open('InflowBoundaryLayerBC.data', 'w')
fp.write('%g\n' % p_e)
fp.write('%d\n' % len(y_sample))
for i in range(len(y_sample)):
    fp.write('%g %g %g\n' % (y_sample[i], T_sample[i], u_sample[i]))
fp.close()

from pylab import linspace, plot, title, xlabel, ylabel, show, subplots
y_plot = linspace(min(ys), max(ys), 100)
T_plot = [sT(y) for y in y_plot]
rho_plot = [p_e/(sT(y)*R) for y in y_plot]
u_plot = [su(y) for y in y_plot]
fig, axs = subplots(3, sharex=True)
axs[0].plot(ys, Ts, '-')
axs[0].plot(y_plot, T_plot, '-r')
axs[1].plot(ys, rhos, '-')
axs[1].plot(y_plot, rho_plot, '-r')
axs[2].plot(ys, us, '-')
axs[2].plot(y_plot, u_plot, '-r')
fig.suptitle("Spline approximation of boundary-layer properties.")
axs[0].set(ylabel ='T, degree K')
axs[1].set(ylabel ='rho, kg/m**3')
axs[2].set(xlabel='y, m', ylabel ='velx, m/s')
show()

print("Done.")
