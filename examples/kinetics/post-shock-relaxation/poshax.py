# poshax.py
#
# Post-shock-relaxation calculation using the general gas model.
#
# PJ, 2019-12-17, Ported from the equivalent Lua script.
#
# Started with Rowan's description in Section 4.4 of his PhD thesis
# and then built upon the gas-model API, the kinetics API,
# the quasi-1D gas flow functions, and linearized 1D flow constraints
# as used in the duct-with-hydrogen-combustion example.
# This is a proof of concept code, that you will probably adapt.
#
# Typical use:
# $ python3 poshax.py
#
#------------------------------------------------------------------
from gdtk.gas import GasModel, GasState, GasFlow, ThermochemicalReactor
debug = False

print("Initialise a gas model.")
print("air 7 species Gupta-et-al reactions")
gmodel = GasModel("air-7sp-gas-model.lua")
nsp = gmodel.n_species
nmodes = gmodel.n_modes
if debug: print("nsp=", nsp, " nmodes=", nmodes, " gmodel=", gmodel)
reactor = ThermochemicalReactor(gmodel, "air-7sp-chemistry.lua")
# The example here matches the case discussed on page 63 of the thesis.
state1 = GasState(gmodel)
state1.p = 133.3 # Pa
state1.T = 300.0 # degree K
state1.massf = {"N2":0.78, "O2":0.22}
print("Free stream conditions, before the shock.")
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("    state1: %s" % state1)
mach1 = 12.28
v1 = mach1 * state1.a
print("mach1:", mach1, "v1:", v1)

print("Stationary normal shock with chemically-frozen gas.")
flow = GasFlow(gmodel)
state2 = GasState(gmodel)
v2, vg = flow.normal_shock(state1, v1, state2)
print("    v2=", v2, "vg=", vg)
print("    state2: %s" % state2)

#------------------------------------------------------------------
print("Decide what to write to the output file.")
output_filename = "air-Mach_12.3.data"
write_massfractions = True
write_molefractions = True

# Output data format is a little like the old poshax code.
header = "# 1:x(m) 2:T(degK) 3:p(Pa) 4:rho(kg/m**3) 5:e(J/kg) 6:v(m/s)"
columnNum = 7
for name in gmodel.species_names:
    if write_massfractions:
        header += (" %d:massf_%s" % (columnNum, name)); columnNum += 1
    if write_molefractions:
        header += (" %d:molef_%s" % (columnNum, name)); columnNum += 1

header += " %d:dt_suggest(s)\n" % (columnNum)

def string_data(x, v, gas, dt_suggest):
    mystr = "%g %g %g %g %g %g" % (x, gas.T, gas.p, gas.rho, gas.u, v)
    for m,c in zip(gas.massf, gas.conc):
        if write_massfractions: mystr += (" %g" % m)
        if write_molefractions: mystr += (" %g" % c)
    mystr += " %g\n" % dt_suggest
    return mystr

#------------------------------------------------------------------
print("Relaxing flow starts here.")

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
print("Step in time, allowing the gas to react and drift along in x.")
gas0 = GasState(gmodel)
gas0.copy_values(state2) # start with post-shock conditions
x = 0 # position of shock, m
x_final = 1.0 # the marching will stop at this point, metres
v = v2 # velocity of gas post-shock, m/s
dt_suggest = 1.0e-8  # suggested starting time-step for chemistry updater
#
f = open(output_filename, 'w')
f.write(header)
f.write(string_data(x, v, gas0, dt_suggest))
#
t = 0 # time of drift is in seconds
t_final = 3.0e-3 # User-selected drift time, s
t_inc = 0.1e-6 # User-selected time steps, s
nsteps = int(t_final / t_inc)
# Notes:
#  Select t_final to be long enough to reach x_final.
#  Select t_inc small enough to capture fast changes near x=0.
#
for j in range(1, nsteps+1):
    # At the start of the step...
    rho = gas0.rho; T = gas0.T; p = gas0.p; u = gas0.u
    #
    # Do the chemical increment.
    gas1 = GasState(gmodel) # make the new one as a clone
    gas1.copy_values(gas0)
    dt_suggest = reactor.update_state(gas1, t_inc, dt_suggest)
    gas1.update_thermo_from_rhou()
    #
    du_chem = gas1.u - u
    dp_chem = gas1.p - p
    if debug: print("# du_chem=", du_chem, "dp_chem=", dp_chem)
    #
    # Do the gas-dynamic accommodation after the chemical change.
    etot = u + 0.5*v*v
    dfdr, dfdu = eos_derivatives(gas1, gmodel)
    if debug: print("# dfdr=", dfdr, "dfdu=", dfdu)
    # Linear solve to get the accommodation increments.
    #   [v,      rho,        0.0, 0.0  ]   [drho  ]   [0.0           ]
    #   [0.0,    rho*v,      1.0, 0.0  ] * [dv    ] = [-dp_chem      ]
    #   [v*etot, rho*etot+p, 0.0, rho*v]   [dp_gda]   [-rho*v*du_chem]
    #   [dfdr,   0.0,       -1.0, dfdu ]   [du_gda]   [0.0           ]
    #
    # Compute the accommodation increments using expressions from Maxima.
    denom = rho*rho*v*v - dfdr*rho*rho - dfdu*p
    drho = (dp_chem - du_chem*dfdu)*rho*rho / denom
    dv = -(dp_chem - du_chem*dfdu)*rho*v / denom
    dp_gda = -(du_chem*dfdu*rho*rho*v*v - dfdr*dp_chem*rho*rho -
               dfdu*dp_chem*p) / denom
    du_gda = -(du_chem*rho*rho*v*v - du_chem*dfdr*rho*rho - dp_chem*p) / denom
    if debug:
        print("# drho=", drho, "dv=", dv, "dp_gda=", dp_gda, "du_gda=", du_gda)
        print("# residuals=", v*drho + rho*dv, rho*v*dv + dp_gda + dp_chem,
              v*etot*drho + (rho*etot+p)*dv + rho*v*du_gda + rho*v*du_chem,
              dfdr*drho - dp_gda + dfdu*du_gda)
    # Add the accommodation increments.
    gas1.rho = gas0.rho + drho
    v1 = v + dv
    p1_check = gas1.p + dp_gda
    gas1.u = gas1.u + du_gda
    gas1.update_thermo_from_rhou()
    if debug:
        print("# At new point for step ", j, ": gas1.p=", gas1.p,
              "p1_check=", p1_check,
              "rel_error=", abs(gas1.p-p1_check)/p1_check)
    # Have now finished the chemical and gas-dynamic update.
    t = t + t_inc
    x = x + 0.5*(v + v1) * t_inc
    f.write(string_data(x, v1, gas1, dt_suggest))
    if x > x_final: break
    # House-keeping for the next step.
    v = v1
    gas0.copy_values(gas1) # gas0 will be used in the next iteration
print("Done stepping.")
