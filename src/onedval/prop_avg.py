# Author: Rowan J. Gollan
# Place: UQ, Brisbane, Queensland, Australia
# Date: 26-Jun-2012
#
# This module provides methods to give
# one-dimensionalised (averaged) flow properties
# based on a variety of techniques.
#

import sys

from math import sqrt, pow
from cfpylib.geom.minimal_geometry import Vector, dot
from gaspy import Gas_data, set_massf, PC_P_atm
from cfpylib.nm.zero_solvers import secant
from scipy.optimize import minimize, brute
from numpy import median
from copy import copy

DEBUG = False
N_RETRIES = 3

def area(cells):
    A = 0.0
    for c in cells:
        A = A + c.area()
    return A

def avg_pos(cells, var_map):
    x = 0.0; y = 0.0; z = 0.0
    xlabel = var_map['x']
    ylabel = var_map['y']
    zlabel = var_map['z']
    A = area(cells)
    for c in cells:
        dA = c.area()
        x += c.get(xlabel)*dA
        y += c.get(ylabel)*dA
        z += c.get(zlabel)*dA
    x /= A
    y /= A
    z /= A
    return Vector(x, y, z)

def compute_fluxes(cells, var_map, species, gmodel, special_fns):
    f_mass = 0.0
    f_mom = Vector(0.0, 0.0, 0.0)
    f_energy = 0.0
    f_sp = [0.0,]*len(species)
    nsp = gmodel.get_number_of_species()
    N = Vector(0.0, 0.0, 0.0)
    A = 0.0
    rholabel = var_map['rho']
    plabel = var_map['p']
    ulabel = var_map['u']
    vlabel = var_map['v']
    wlabel = var_map['w']
    Tlabel = var_map['T']
    Q = Gas_data(gmodel)
    special_fluxes = copy(special_fns)
    for k in special_fluxes:
        special_fluxes[k] = 0.0
    for c in cells:
        dA = c.area()
        A += dA
        n = c.normal()
        rho = c.get(rholabel)
        p = c.get(plabel)
        vel = Vector(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        T = c.get(Tlabel)
        # Add increments
        u_n = dot(vel, n)
        f_mass = f_mass + rho*u_n*dA
        f_mom = f_mom + (rho*u_n*vel + p*n)*dA
        if gmodel.get_number_of_species() > 1:
            massfs = {}
            for isp, sp in enumerate(species):
                splabel = var_map.get(sp, sp)
                massf = c.get(splabel)
                f_sp[isp] = f_sp[isp] + rho*massf*u_n*dA
                massfs[sp] = massf
            set_massf(Q, gmodel, massfs)
        else:
            Q.massf[0] = 1.0
            f_sp[0] = f_sp[0] + rho*u_n*dA
        Q.rho = c.get(rholabel)
        Q.T[0] = c.get(Tlabel)
        gmodel.eval_thermo_state_rhoT(Q)
        h = gmodel.mixture_enthalpy(Q)
        h0 = h + 0.5*abs(vel)*abs(vel)
        f_energy = f_energy + rho*u_n*h0*dA
    
    # Process any special fns.
    # Special fns may like to know total flux and area
    fluxes = {'mass':abs(f_mass), 'mom':f_mom, 'energy':abs(f_energy), 'species':map(abs, f_sp)}
    
    for c in cells:
	dA = c.area()
        rho = c.get(rholabel)
        vel = Vector(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        n = c.normal()
        u_n = dot(vel, n)
        for l,f in special_fns.iteritems():
            special_fluxes[l] += f(c, rho, u_n, dA, A, fluxes)
    output = {'mass': f_mass, 'mom': f_mom, 'energy': f_energy, 'species':f_sp}
    for l,f in special_fluxes.iteritems():
        output[l] = special_fluxes[l]
    return output

def area_weighted_avg(cells, props, var_map):
    phis = dict.fromkeys(props, 0.0)
    area = 0.0
    for c in cells:
        dA = c.area()
        area = area + dA
        for p in props:
            label = p
            if p in var_map:
                label = var_map[p]
            phis[p] = phis[p] + c.get(label)*dA
    
    for p in props:
        phis[p] = phis[p]/area
    
    return phis

def mass_flux_weighted_avg(cells, props, var_map):
    phis = dict.fromkeys(props, 0.0)
    f_mass = 0.0
    rholabel = var_map['rho']
    ulabel = var_map['u']
    vlabel = var_map['v']
    wlabel = var_map['w']
    Tlabel = var_map['T']
    for c in cells:
        dA = c.area()
        rho = c.get(rholabel)
        vel = Vector(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        n = c.normal()
        w = rho*dot(vel, n)
        f_mass = f_mass + w*dA
        for p in props:
            label = p
            if p in var_map:
                label = var_map[p]
            phis[p] = phis[p] + c.get(label)*w*dA
    
    for p in props:
        phis[p] = phis[p]/f_mass

    return phis


def stream_thrust_avg(cells, props, var_map, species, gmodel):
    flag = 'success'
    f_mass = 0.0
    f_mom = Vector(0.0, 0.0, 0.0)
    f_energy = 0.0
    f_sp = [0.0,]*len(species)
    nsp = gmodel.get_number_of_species()
    N = Vector(0.0, 0.0, 0.0)
    A = 0.0
    rholabel = var_map['rho']
    plabel = var_map['p']
    ulabel = var_map['u']
    vlabel = var_map['v']
    wlabel = var_map['w']
    Tlabel = var_map['T']
    Q = Gas_data(gmodel)
    for c in cells:
        dA = c.area() 
        n = c.normal()
        rho = c.get(rholabel)
        p = c.get(plabel)
        vel = Vector(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        T = c.get(Tlabel)
        # Add increments
        u_n = dot(vel, n)
        f_mass = f_mass + rho*u_n*dA
        f_mom = f_mom + (rho*u_n*vel + p*n)*dA
        if gmodel.get_number_of_species() > 1:
            massfs = {}
            for isp, sp in enumerate(species):
                splabel = var_map.get(sp, sp)
                massf = c.get(splabel)
                f_sp[isp] = f_sp[isp] + rho*massf*u_n*dA
                massfs[sp] = massf
            set_massf(Q, gmodel, massfs)
        else:
            Q.massf[0] = 1.0
            f_sp[0] = f_sp[0] + rho*u_n*dA
        Q.rho = c.get(rholabel)
        Q.T[0] = c.get(Tlabel)
        gmodel.eval_thermo_state_rhoT(Q)
        h = gmodel.mixture_enthalpy(Q)
        h0 = h + 0.5*abs(vel)*abs(vel)
        f_energy = f_energy + rho*u_n*h0*dA
        A = A + dA
        N = N + n*dA

    N = N / A
    f_mom_s = dot(f_mom, N)
    Q = Gas_data(gmodel)
    if nsp > 1:
        massf = [ f_isp/f_mass for f_isp in f_sp ]
        set_massf(Q, gmodel, massf)
    else:
        Q.massf[0] = 1.0
    
    LARGE_PENALTY = 1.0e6

    def f_to_minimize(x):
        rho, T, u = x
        if rho < 0.0 or T < 0.0:
            # Give a big penalty
            return LARGE_PENALTY
        # Use equation of state to compute other thermo quantities
        Q.rho = rho
        Q.T[0] = T
        flag = gmodel.eval_thermo_state_rhoT(Q)
        if flag != 0:
            # If there are problems, then these are NOT good values
            # so return a large error
            return LARGE_PENALTY
        h = gmodel.mixture_enthalpy(Q)
        p = Q.p
        # Compute errors
        fmass_err = abs(f_mass - rho*u*A)/(abs(f_mass))
        fmom_err = abs(f_mom_s - (rho*u*u + p)*A)/(abs(f_mom_s))
        fe_err = abs(f_energy - (rho*u*A*(h + 0.5*u*u)))/(abs(f_energy))
        # Total error is the sum
        if DEBUG:
            print "DEBUG: fmass_err= ", fmass_err
            print "DEBUG: fmom_err= ", fmom_err
            print "DEBUG: fe_err= ", fe_err
            print "DEBUG: total error= ", fmass_err + fmom_err + fe_err
        return fmass_err + fmom_err + fe_err

    # Find bounds for answer.
    rhos = [c.get(rholabel) for c in cells]
    rho_min = min(rhos)
    rho_max = max(rhos)
    rho_mid = median(rhos)
    Ts = [c.get(Tlabel) for c in cells]
    T_min = min(Ts)
    T_max = max(Ts)
    T_mid = median(Ts)
    us = [sqrt(c.get(ulabel)**2 + c.get(vlabel)**2 + c.get(wlabel)**2) for c in cells]
    u_min = min(us)
    u_max = max(us)
    u_mid = median(us)
    # ------------ This is a complicated way to get an initial guess ----------- #
#    aw_props = area_weighted_avg(cells, ['rho', 'p', 'T', 'u', 'v', 'w', 'M'], var_map)
    # Use the Nelder-Mead minimiser with area-weighted averages as starting guess
#    u = sqrt(aw_props['u']**2 + aw_props['v']**2 + aw_props['w']**2)
#    rho = aw_props['rho']
#    f_mass_guess = rho*u*A
    # Scale guesses for rho and u such that f_mass is correct initially
#    rho = rho*sqrt(f_mass/f_mass_guess)
#    u = u*sqrt(f_mass/f_mass_guess)
#    p = f_mom_s/A - rho*u*u
#    Q.rho = rho
#    Q.p = p
#    gmodel.eval_thermo_state_rhop(Q)
#    T = Q.T[0]
    # ------------- Just use median values for initial guess -------------- #
    guess = [rho_mid, T_mid, u_mid]
    result = minimize(f_to_minimize, guess, method='Nelder-Mead', options={'ftol':1.0e-6})
    rho, T, u = result.x
    if DEBUG:
        print "DEBUG: First pass of optimiser... result= ", result.success
        print "DEBUG: rho= ", rho, " T= ", T, " u= ", u
        print "DEBUG: f_mass= ", f_mass
        print "DEBUG: computed from props= ", rho*u*A
        print "DEBUG: f_mom_s= ", f_mom_s
        print "DEBUG: computed from props= ", A*(rho*u*u + p)
        
    if rho < rho_min or rho > rho_max or T < T_min or T > T_max or u < u_min or u > u_max:
        print "Something went seriously wrong with the optimiser because at least one"
        print "of the returned values is outside the range of data."
        print "rho_min= ", rho_min, " rho_max= ", rho_max, " computed rho= ", rho
        print "T_min= ", T_min, " T_max= ", T_max, " computed T= ", T
        print "u_min= ", u_min, " u_max= ", u_max, " computed u= ", u
        result.success = False
        # We need a new initial guess
        aw_props = area_weighted_avg(cells, ['rho', 'p', 'T', 'u', 'v', 'w', 'M'], var_map)
        u = sqrt(aw_props['u']**2 + aw_props['v']**2 + aw_props['w']**2)
        rho = aw_props['rho']
        Q.rho = rho
        Q.p = p
        gmodel.eval_thermo_state_rhop(Q)
        T = Q.T[0]

    if not result.success:
        for i in range(N_RETRIES):
            print "The optimiser did not converge. Attempting another go."
            print "Retry attempt: %d/%d" % (i+1, N_RETRIES)
            print "Retrying with the last (unconverged) values as the initial guess."
            guess = [rho, T, u]
            print "rho= ", rho, " T= ", T, " u= ", u
            result = minimize(f_to_minimize, guess, method='Nelder-Mead', options={'ftol':1.0e-6})
            rho, T, u = result.x
            if rho < rho_min or rho > rho_max or T < T_min or T > T_max or u < u_min or u > u_max:
                print "Something went seriously wrong with the optimiser because at least one"
                print "of the returned values is outside the range of data."
                print "rho_min= ", rho_min, " rho_max= ", rho_max, " computed rho= ", rho
                print "T_min= ", T_min, " T_max= ", T_max, " computed T= ", T
                flag = 'failed'
                result.success = False
            if result.success:
                break
    
    if not result.success:
        # We'll set this as a failed state and the user can decide whether
        # to keep it or not by using the 'skip_bad_slices' option.
        flag = 'failed'
        # Now use the brute force search
        print "The Nelder-Mead optimiser has run into troubles."
        print "Now we'll use a brute force approach in the region of the"
        print "median values."
        delta = 0.15
        rho_l = (1.0 - delta)*rho_mid; rho_h = (1.0+delta)*rho_mid
        T_l = (1.0 - delta)*T_mid; T_h = (1.0+delta)*T_mid
        u_l = (1.0-delta)*u_mid; u_h = (1.0+delta)*u_mid
        x0 = brute(f_to_minimize, ((rho_l, rho_h), (T_l, T_h), (u_l, u_h)))
        flag = 'success'
        rho, T, u = x0
       
    Q.rho = rho; Q.T[0] = T
    gmodel.eval_thermo_state_rhoT(Q)
    h = gmodel.mixture_enthalpy(Q)
    M = u/Q.a
    p = Q.p
    gamma = gmodel.gamma(Q)
    T0 = T*(1.0 + (gamma - 1.0) * 0.5 * M**2)
    p0 = p * pow(T0/T, gamma/(gamma - 1.0))
    hr = 0.0
    # temporary dummy gas object
    Qd = Gas_data(gmodel)
    Qd.p = PC_P_atm
    Qd.T[0] = 298.15 # K <-- DO NOT ADJUST THIS VALUE
                     # CEA curve fits are constructed
                     # such that h(298.15) == h_f
    for isp, f in enumerate(f_sp):
        hr += f*gmodel.enthalpy(Qd, isp) # <-- is h_f because evaluated at T=298.15 K
    
    # Create a dictionary of all possible requests
    vals = {'rho':rho,
            'p': p,
            'p0': p0,
            'T': T,
            'T0': T0,
            'a': Q.a,
            'h': h,
            'h0': h + 0.5*u*u,
            's': gmodel.mixture_entropy(Q),
            'R': gmodel.R(Q),
            'Cp': gmodel.Cp(Q),
            'Cv': gmodel.Cv(Q),
            'gamma': gamma,
            'u': u,
            'M': M,
            'hr': hr}

    phis = dict.fromkeys(props, 0.0)
    
    for k in phis:
        if not k in vals:
            print "Do not know what to do for flux-conserved average of:", k
            print "Skipping this request."
            continue
        phis[k] = vals[k]
    
    return phis, flag
    
    
    

    
    
        

