# nenzfr2.py
# Compute the state-to-state processes and the nonequilibrium nozzle flow
# for a particular shot of the the T4 shock tunnel.
# Part A: state-to-state calculations of the shock-tube conditions,
#   to the sonic flow condition at the nozzle throat.
# Part B: steady-state space-marching along the supersonic part of the nozzle.
#
# Run with the commands:
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air5species-gas-model.lua .
# $ cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .
# $ cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
# $ prep-gas air-5sp-1T.inp air-5sp-1T.lua
# $ prep-chem air-5sp-1T.lua GuptaEtAl-air-reactions.lua air-5sp-1T-reactions.lua
# $ python3 nenzfr2.py > wilson-10572.txt
#
# The output file (wilson-10572.txt) will be in a format ready for gnuplot.
#
# Author: Peter J.
#
# Versions:
# 2019-11-30, state-to-state calculation of shock tube conditions
# 2019-12-27, nozzle flow calculation adapted from supersonic-diffuser example
#

#------------------------------------------------------------------

print("# Nenzfr2 proof-of-concept script, with Wilson's T4 shot 10572.")
print("# Part A. State-to-state calculation of the shock tube processes.")
import math
from eilmer.gas import GasModel, GasState, GasFlow
debug = False

gmodel = GasModel('cea-air5species-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 80.0e3 # Pa
state1.T = 300.0 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("# Initial test gas:")
print("#   state1: %s" % state1)

print("# Normal shock, given shock speed")
vs = 2890.0
print("#   vs=%g" % vs)
state2 = GasState(gmodel)
flow = GasFlow(gmodel)
v2, vg = flow.normal_shock(state1, vs, state2)
print("#   v2=%g vg=%g" % (v2, vg))
print("#   state2: %s" % state2)

print("# Reflected shock")
state5 = GasState(gmodel)
vr_b = flow.reflected_shock(state2, vg, state5)
print("#   vr_b=%g" % vr_b)
print("#   state5: %s" % state5)

print("# Expand from stagnation (with ratio of pressure to match observation).")
state5s = GasState(gmodel)
v5s = flow.expand_from_stagnation(state5, 37.47e6/state5.p, state5s)
print("#   v5s=%g Mach=%g" % (v5s, v5s/state5s.a))
print("#   state5s: %s" % state5s)
print("#   (h5s-h1)=%g" % (state5s.enthalpy - state1.enthalpy))

print("# Expand to throat condition (Mach 1.0001) assuming equilibrium.")
state6 = GasState(gmodel)
v6 = flow.expand_to_mach(state5s, 1.0001, state6)
print("#   v6=%g Mach=%g" % (v6, v6/state6.a))
print("#   state6: %s" % state6)
print("#   ceaSavedData=%s" % state6.ceaSavedData)

print("# Expand (with equilibrium thermochem) to something like Mach 6 nozzle")
state7 = GasState(gmodel)
v7 = flow.steady_flow_with_area_change(state6, v6, 127.0, state7)
print("#   v7=%g Mach=%g" % (v7, v7/state7.a))
print("#   state7: %s" % state7)

#------------------------------------------------------------------

print("# Part B. Nonequilibrium nozzle flow.")
from eilmer.gas import ThermochemicalReactor

# Approximation of the T4 Mach 6 nozzle profile as straight segments.
# The following are the transition points, derived from Wilson's nenzf file.
xi = [0.0, 0.150, 0.280, 0.468, 0.671, 0.984] # metres
ri = [1.0, 4.0,   6.232, 8.44,  9.92,  11.268] # ratio r/r_throat

def duct_area(x):
    "Returns A/A_throat for x im metres."
    if x < xi[0]:
        r = ri[0]
    elif x < xi[1]:
        xx = (x-xi[0])/(xi[1]-xi[0]); r = ri[0]*(1.0-xx) + ri[1]*xx
    elif x < xi[2]:
        xx = (x-xi[1])/(xi[2]-xi[1]); r = ri[1]*(1.0-xx) + ri[2]*xx
    elif x < xi[3]:
        xx = (x-xi[2])/(xi[3]-xi[2]); r = ri[2]*(1.0-xx) + ri[3]*xx
    elif x < xi[4]:
        xx = (x-xi[3])/(xi[4]-xi[3]); r = ri[3]*(1.0-xx) + ri[4]*xx
    else:
        xx = (x-xi[4])/(xi[5]-xi[4]); r = ri[4]*(1.0-xx) + ri[5]*xx
    return r*r

# Things that we'll make use of shortly.
sample_header = "# x(m) A(m**2) rho(kg/m**3) p(Pa) T(degK) e(J/kg) v(m/s) " + \
                "massf_N2 massf_O2 massf_N massf_O massf_NO dt_suggest(s) mdot(kg/s)"

def sample_data(x, area, v, gas, dt_suggest):
    mf = gas.massf_as_dict
    return "%g %g %g %g %g %g %g %g %g %g %g %g %g %g" % \
        (x, area, gas.rho, gas.p, gas.T, gas.u, v,
         mf["N2"], mf["O2"], mf["N"], mf["O"], mf["NO"],
         dt_suggest, gas.rho*v*area)

def eos_derivatives(gas0, gmodel, tol=0.0001):
    # Finite difference evaluation, assuming that gas0 is valid state already.
    gas1 = GasState(gmodel)
    gas1.copy_values(gas0)
    p0 = gas0.p; rho0 = gas0.rho; u0 = gas0.u
    #
    drho = rho0 * tol; gas1.rho = rho0 + drho
    gas1.update_thermo_from_rhou()
    dpdrho = (gas1.p - p0)/drho
    #
    gas1.rho = rho0; du = u0 * tol; gas1.u = u0 + du
    gas1.update_thermo_from_rhou()
    dpdu = (gas1.p - p0)/du
    #
    return dpdrho, dpdu

#------------------------------------------------------------------
print("# Start stepping down the nozzle expansion.")
print(sample_header)

# Gas model setup for nonequilibrium chemistry.
gmodel2 = GasModel("air-5sp-1T.lua")
nsp = gmodel2.n_species
nmodes = gmodel2.n_modes
if debug:
    print("# nsp=", nsp, " nmodes=", nmodes, " gmodel2=", gmodel2)

print("# Gas properties at the start of the expansion.")
gas0 = GasState(gmodel2)
gas0.p = state6.p
x = 0.0 # m  (location of start of expansion)
area = duct_area(x)
gas0.T = state6.T
mf = {'N2':state6.ceaSavedData['massf']['N2'],
      'O2':state6.ceaSavedData['massf']['O2'],
      'N':state6.ceaSavedData['massf']['N'],
      'O':state6.ceaSavedData['massf']['O'],
      'NO':state6.ceaSavedData['massf']['NO']}
# Normalise the mass fractions because Rowan's thermochemistry
# has tighter tolerance than the default tolerance in CEA.
mf_sum = mf['N2']+mf['O2']+mf['N']+mf['O']+mf['NO']
for k in mf.keys(): mf[k] = mf[k]/mf_sum
if debug: print("# after scaling, mf=", mf)
gas0.massf = mf # Set them all at once.
gas0.update_thermo_from_pT()
gas0.update_sound_speed()
if debug:
    print("# species=", gmodel2.species_names)
    print("# massf=", gas0.massf)
# For the same pressure and temperature, the nonequilibrium gas
# will have a different sound-speed to that for the equilibrium gas.
# The old Nenzf had a process of perturbing the flow condition such that,
# with a sufficiently large perturbation, the expanding flow will
# eventually remain supersonic.
if 0:
    # We shall lower the temperature of the gas until sound sound-speed
    # is slightly lower than the flow speed through the nozzle throat.
    v = v6
    while 1.01*gas0.a > v:
        gas0.T -= 1.0
        gas0.update_thermo_from_pT()
        gas0.update_sound_speed()
    print("# v=", v, "a=", gas0.a, "T-state6.T=", gas0.T-state6.T)
else:
    # Keep the temperature but change the velocity to give a slightly
    # supersonic condition.
    v = 1.001*gas0.a
    print("# v=", v, "a=", gas0.a, "v-v6=", v-v6)

dt_suggest = 1.0e-12  # suggested starting time-step for chemistry
print(sample_data(x, area, v, gas0, dt_suggest))

print("# Start reactions...")
reactor = ThermochemicalReactor(gmodel2, "air-5sp-1T-reactions.lua")
t = 0 # time is in seconds
t_final = 0.8e-3 # long enough to convect past exit
t_inc = 1.0e-9 # start small
t_inc_max = 1.0e-6
x_end = 0.984 # metres
while x < x_end:
    # At the start of the step...
    rho = gas0.rho; T = gas0.T; p = gas0.p; u = gas0.u
    #
    # Do the chemical increment.
    gas1 = GasState(gmodel2) # make the new one as a clone
    gas1.copy_values(gas0)
    dt_suggest = reactor.update_state(gas1, t_inc, dt_suggest)
    gas1.update_thermo_from_rhou()
    #
    du_chem = gas1.u - u
    dp_chem = gas1.p - p
    if debug: print("# du_chem=", du_chem, "dp_chem=", dp_chem)
    #
    # Update the independent variables for the end point of this step.
    # Simple, Euler update for the spatial stepper.
    t1 = t + t_inc
    x1 = x + v*t_inc
    area1 = duct_area(x1)
    #
    # Do the gas-dynamic accommodation after the chemical change.
    darea = area1 - area
    etot = u + 0.5*v*v
    dfdr, dfdu = eos_derivatives(gas1, gmodel2)
    if debug:
        print("# dfdr=", dfdr, "dfdu=", dfdu)
        print("# x=", x, "v=", v, "A=", area, "dA=", darea)
    # Linear solve to get the accommodation increments.
    #   [v*A,      rho*A,          0.0, 0.0    ]   [drho  ]   [-rho*v*dA       ]
    #   [0.0,      rho*v,          1.0, 0.0    ] * [dv    ] = [-dp_chem        ]
    #   [v*etot*A, (rho*etot+p)*A, 0.0, rho*v*A]   [dp_gda]   [-rho*v*A*du_chem]
    #   [dfdr,     0.0,           -1.0, dfdu   ]   [du_gda]   [0.0             ]
    #
    # Compute the accommodation increments using expressions from Maxima.
    denom = area*(rho*rho*v*v - dfdr*rho*rho - dfdu*p)
    drho = (area*(dp_chem - du_chem*dfdu)*rho*rho - darea*rho**3*v*v) / denom
    dv = (area*(du_chem*dfdu - dp_chem)*rho +
          darea*(dfdr*rho*rho + dfdu*p))*v / denom
    dp_gda = -((darea*dfdr*rho**3 + area*dfdu*du_chem*rho*rho + darea*dfdu*p*rho)*v*v
               - area*dfdr*dp_chem*rho*rho - area*dfdu*dp_chem*p) / denom
    du_gda = -(area*(du_chem*rho*rho*v*v - du_chem*dfdr*rho*rho - dp_chem*p)
               + darea*p*rho*v*v) / denom
    if debug:
        print("# drho=", drho, "dv=", dv, "dp_gda=", dp_gda, "du_gda=", du_gda)
        print("# residuals=", v*area*drho + rho*area*dv + rho*v*darea,
              rho*v*dv + (dp_gda + dp_chem),
              v*etot*drho*area + (rho*etot+p)*area*dv +
              rho*v*area*(du_gda + du_chem) + v*(rho*etot+p)*darea,
              dfdr*drho - dp_gda + dfdu*du_gda)
    # Add the accommodation increments.
    gas1.rho = gas0.rho + drho
    v1 = v + dv
    p1_check = gas1.p + dp_gda
    gas1.u = gas1.u + du_gda
    gas1.update_thermo_from_rhou()
    gas1.update_sound_speed()
    if debug:
        print("# At new point x1=", x1, "v1=", v1,
              ": gas1.p=", gas1.p, "p1_check=", p1_check,
              "rel_error=", abs(gas1.p-p1_check)/p1_check)
    # Have now finished the chemical and gas-dynamic update.
    print(sample_data(x1, area1, v1, gas1, dt_suggest))
    # House-keeping for the next step.
    v = v1; t = t1; x = x1; area = area1
    gas0.copy_values(gas1) # gas0 will be used in the next iteration
    t_inc = min(t_inc*1.01, t_inc_max)
    # import sys; sys.exit()

print("# Exit condition:")
print("#   temperature (K)", gas1.T)
print("#   pressure (kPa)", gas1.p/1000)
print("#   density (kg/m^3)", gas1.rho)
print("#   velocity (m/s)", v1)
print("#   Mach", v1/gas1.a)
print("#   rho*v^2 (kPa)", gas1.rho*v1*v1/1000)
print("#   area ratio", area1)
print("#   species=", gmodel2.species_names)
print("#   massf=", gas1.massf)
