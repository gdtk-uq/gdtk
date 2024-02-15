"""
Example to verify source term implementation in Puffin. 
The simulations look at 2-D channels (rectangular of circular), where
hear, mass, and momentum are added to impact the flow state. 

In addittion to analytical equations are solved to provide reference 
results for comparison. 

Caution: 
Excess source term values can lead to chocking and flow velocities 
less than Mach 1. In these cases the solutions loose validity. 

Authors: Ingo Jahn
Created: 2024-02-14
Last Modified: 2024-02-14
"""




import numpy as np
from scipy.optimize import root
from gdtk.geom.xpath import XPath

gas_file = "ideal-air-gas-model.lua"
init_gas_model(gas_file)
gas1 = GasState(config.gmodel)
gas1.p = 1.0e3
gas1.T = 300.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 6.0
V1 = M1 * gas1.a
gamma = gas1.gamma

print(f"################################")
print(f"### Gas Model and Flow State ###")
print(f"################################")
print(f"Gas Model: {gas_file}")
print(f"Mach No: {M1} [-]")
print(f"P1: {gas1.p} [Pa]")
print(f"T1: {gas1.T} [K]")
print(f"V1: {V1} [m/s]")
print(f"rho1: {gas1.rho} [kg/m3]")
print(f"gamma: {gamma} [-]")
print()


length = 1.2
height = 0.1
delta = 0.01

# create multiple parallel channel to simualte effect of different source terms
path_0 = {}
path_1 = {}
for i in range(4):
    h_low = i*(height+delta)
    h_upp = h_low + height
    path_0[i] = XPath().moveto(0, h_low).lineto(length, h_low)
    path_1[i] = XPath().moveto(0, h_upp).lineto(length, h_upp)

# set boundaries at top and bottom of each channel as slip walls
def bc_path_0(x): return 0   # slip wall
def bc_path_1(x): return 0   # slip wall

# define source term distributions. Source terms are active between x=0.1 and x=1.1.
source_mass = 5
source_xmom = 5e2
source_tE = 5e6
# Source term for mass
def rho_func(x): 
    if 0.1 < x and x < 1.1: 
        return source_mass 
    else:
        return 0
# Source term for x-momentum
def xmo_func(x): 
    if 0.1 < x and x < 1.1: 
        return source_xmom
    else:
        return 0
# Source term for thermal energy
def tE_func(x): 
    if 0.1 < x and x < 1.1: 
        return source_tE  
    else:
        return 0

# configure Puffin settings
config.max_step_relax = 40
config.max_x = length
config.dx = config.max_x/1000
config.add_user_supplied_source_terms = True

# Reference Case - No source terms
st0 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=path_0[0], y1=path_1[0], bc0=bc_path_0, bc1=bc_path_1,
                 ncells=20)

# Add mass
st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=path_0[1], y1=path_1[1], bc0=bc_path_0, bc1=bc_path_1,
                 add_rho=rho_func,
                 ncells=20)

# Add x-momentum 
st2 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=path_0[2], y1=path_1[2], bc0=bc_path_0, bc1=bc_path_1,
                 add_xmo=xmo_func,
                 ncells=20)

# Add thermal energy
st3 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=path_0[3], y1=path_1[3], bc0=bc_path_0, bc1=bc_path_1,
                 add_tE=tE_func,
                 ncells=20)


# Solving Analytical Solution
# Here we simulate the effect of the source terms by solving the conservation 
# equations across the region of mass, momentum, and energy addition.
# Calculations below assume source terms act over a region with length 1m. 

# Upstream Conditions
V1 = V1
rho1 = gas1.rho
T1 = gas1.T
x1 = (rho1, T1, V1)
gas2 = GasState(config.gmodel)

def conservation_fun(x2, x1, gas1, gas2, d_rho, d_xmom, d_tE):
    (rho2, T2, V2) = x2
    (rho1, T1, V1) = x1
    # use Equation of state to decode pressure and enthalpy
    gas1.rho = rho1
    gas1.T = T1
    gas1.update_thermo_from_rhoT()
    P1 = gas1.p
    h1 = gas1.enthalpy
    gas2.rho = rho2
    gas2.T = T2
    gas2.update_thermo_from_rhoT()
    P2 = gas2.p
    h2 = gas2.enthalpy
    # Define errors (Note - these are normalised with respect to input fluxes)
    error_mass = (V1*rho1 + d_rho - V2*rho2) / (V1*rho1)
    error_xmom = (P1 + V1*rho1*V1 + d_xmom - (P2 + V2*rho2*V2)) / (P1 + V1*rho1*V1)
    error_tE = (V1*rho1* (h1 + 0.5*V1*V1) + d_tE - V2*rho2* (h2+0.5*V2*V2)) / (V1*rho1* (h1 + 0.5*V1*V1))
    return [error_mass, error_xmom, error_tE] 

print(f'###########################################')
print(f'### Analytical Results (ideal gas only) ###')
print(f'###########################################')
print(f"Inlet:")
print(f"M1: {M1} [-]; V1: {V1} [m/s]; P1: {gas1.p} [Pa];  T1: {gas1.T} [Pa];  rho1: {gas1.rho} [kg/m3]")
print()
print(f"st0 - No Source terms")
x0 = (gas1.rho, gas1.T, V1)
args = (x1, gas1, gas2, 0, 0, 0)
sol = root(conservation_fun, x0, args=args, method='lm')
print(f"Solver Succes: {sol.success}")
rho2 = sol.x[0]
T2 = sol.x[1]
V2 = sol.x[2]
gas2.rho = rho2
gas2.T = T2
gas2.update_thermo_from_rhoT()
print(f"M2: {V2 / gas2.a} [-]; V2: {V2} [m/s]; P2: {gas2.p} [Pa];  T2: {T2} [Pa];  rho2: {rho2} [kg/m3]")
print()

print(f"st1 - Mass Source: {source_mass} [kg/m3]")
x0 = (gas1.rho, gas1.T, V1)
args = (x1, gas1, gas2, source_mass, 0, 0)
sol = root(conservation_fun, x0, args=args, method='lm')
print(f"Solver Succes: {sol.success}")
rho2 = sol.x[0]
T2 = sol.x[1]
V2 = sol.x[2]
gas2.rho = rho2
gas2.T = T2
gas2.update_thermo_from_rhoT()
print(f"M2: {V2 / gas2.a} [-]; V2: {V2} [m/s]; P2: {gas2.p} [Pa];  T2: {T2} [Pa];  rho2: {rho2} [kg/m3]")
print()

print(f"st2 - X-mom Source: {source_xmom} [Pa/m3]")
x0 = (gas1.rho, gas1.T, V1)
args = (x1, gas1, gas2, 0, source_xmom, 0)
sol = root(conservation_fun, x0, args=args, method='lm')
print(f"Solver Succes: {sol.success}")
rho2 = sol.x[0]
T2 = sol.x[1]
V2 = sol.x[2]
gas2.rho = rho2
gas2.T = T2
gas2.update_thermo_from_rhoT()
print(f"M2: {V2 / gas2.a} [-]; V2: {V2} [m/s]; P2: {gas2.p} [Pa];  T2: {T2} [Pa];  rho2: {rho2} [kg/m3]")
print()

print(f"st3 - tE Source: {source_tE} [J/m3]")
x0 = (gas1.rho, gas1.T, V1)
args = (x1, gas1, gas2, 0, 0, source_tE)
sol = root(conservation_fun, x0, args=args, method='lm')
print(f"Solver Succes: {sol.success}")
rho2 = sol.x[0]
T2 = sol.x[1]
V2 = sol.x[2]
gas2.rho = rho2
gas2.T = T2
gas2.update_thermo_from_rhoT()
print(f"M2: {V2 / gas2.a} [-]; V2: {V2} [m/s]; P2: {gas2.p} [Pa];  T2: {T2} [Pa];  rho2: {rho2} [kg/m3]")
print()
print(F"DONE.")
print()