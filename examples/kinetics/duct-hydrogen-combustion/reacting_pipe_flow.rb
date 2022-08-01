# reacting_pipe_flow.rb
#
# Combustion in a supersonic stream of constant cross-section.
# This is set up to approximate the Bittker-Scullin case 3,
# as used by Fabian Zander in the hydrogen-combustion test case.
#
# PJ 2011-06-21 first version written for Eilmer3 test case
#    2016-03-19 adapted to Lua from reacting-pipe-flow.py for Eilmer3
#    2018-04-21 updated to use current gas model and reaction calls
#    2019-12-03 Ruby flavour using the CFFI gas module.
#
# Run with the commands:
# $ prep-gas combusting-species.inp h2-o2-n2-9sp.lua
# $ prep-chem h2-o2-n2-9sp.lua Bittker-Scullin.lua h2-o2-n2-9sp-18r.lua
# $ ruby reacting_pipe_flow.rb
#
#------------------------------------------------------------------
$LOAD_PATH << '~/dgdinst/lib'
require 'gdtk/gas'

# Things that we'll make use of shortly.
sample_header = "# x(m) rho(kg/m**3) p(Pa) T(degK) e(J/kg) v(m/s) "
sample_header << "massf_OH massf_H2O dt_suggest(s)"

def sample_data(x, v, gas, dt_suggest)
  return "%g %g %g %g %g %g %g %g %g " % \
         [x, gas.rho, gas.p, gas.T, gas.u, v,
          gas.massf[gas.gmodel.species_names.index("OH")],
          gas.massf[gas.gmodel.species_names.index("H2O")],
          dt_suggest]
end

def eos_derivatives(gas0, gmodel, tol=0.0001)
  # Finite difference evaluation, assuming that gas0 is valid state already.
  gas1 = GasState.new(gmodel)
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
end

#------------------------------------------------------------------
# Start the main script...
debug = false
puts "# Reacting pipe flow -- Bittker-Scullin test case 3."
puts sample_header

# Gas model setup
gmodel = GasModel.new("h2-o2-n2-9sp.lua")
nsp = gmodel.n_species
nmodes = gmodel.n_modes
if debug then
  puts "nsp=#{nsp} nmodes=#{nmodes} gmodel=#{gmodel}"
end

puts "# Gas properties at the start of the pipe."
gas0 = GasState.new(gmodel)
gas0.p = 96.87e3 # Pa
x = 0.0 # m  (inlet of pipe)
v = 4551.73 # m/s
gas0.T = 1559.0 # degree K
gas0.molef = {"O2"=>0.1480, "N2"=>0.5562, "H2"=>0.2958}
gas0.update_thermo_from_pT()
dt_suggest = 1.0e-8  # suggested starting time-step for chemistry updater
puts sample_data(x, v, gas0, dt_suggest)

puts "# Start reactions..."
reactor = ThermochemicalReactor.new(gmodel, "h2-o2-n2-9sp-18r.lua") 
t = 0 # time is in seconds
t_final = 22.0e-6
t_inc = 0.05e-6
nsteps = (t_final / t_inc).to_i
(1..nsteps).each do |j|
  # At the start of the step...
  rho = gas0.rho; myT = gas0.T; p = gas0.p; u = gas0.u
  #
  # Do the chemical increment.
  gas1 = GasState.new(gmodel) # make the new one as a clone
  gas1.copy_values(gas0)
  dt_suggest = reactor.update_state(gas1, t_inc, dt_suggest)
  gas1.update_thermo_from_rhou()
  #
  du_chem = gas1.u - u
  dp_chem = gas1.p - p
  if debug then
    puts "# du_chem=#{du_chem} dp_chem=#{dp_chem}"
  end
  #
  # Do the gas-dynamic accommodation after the chemical change.
  etot = u + 0.5*v*v
  dfdr, dfdu = eos_derivatives(gas1, gmodel)
  if debug then
    puts "# dfdr=#{dfdr} dfdu=#{dfdu}"
  end
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
  if debug then
    puts "# drho=#{drho} dv=#{dv} dp_gda=#{dp_gda} du_gda=#{du_gda}"
    puts "# residuals=%g %g %g %g" %
         [v*drho + rho*dv, rho*v*dv + dp_gda + dp_chem,
          v*etot*drho + (rho*etot+p)*dv + rho*v*du_gda + rho*v*du_chem,
          dfdr*drho - dp_gda + dfdu*du_gda]
  end
  # Add the accommodation increments.
  gas1.rho = gas0.rho + drho
  v1 = v + dv
  p1_check = gas1.p + dp_gda
  gas1.u = gas1.u + du_gda
  gas1.update_thermo_from_rhou()
  if debug then
    puts "# At new point for step #{j}: gas1.p=#{gas1.p} " \
         "p1_check=#{p1_check} rel_error=#{(gas1.p-p1_check).abs/p1_check}"
  end
  # Have now finished the chemical and gas-dynamic update.
  t = t + t_inc
  x = x + 0.5*(v + v1) * t_inc
  puts sample_data(x, v1, gas1, dt_suggest)
  # House-keeping for the next step.
  v = v1
  gas0.copy_values(gas1) # gas0 will be used in the next iteration
end
puts "# Done."
