--[[ cj-det-wav.lua

Author: Rowan J. Gollan
Date: 2020-05-10

This script defines an Arrhenius 1-D Chapman-Jouguet (C-J)
detonation wave. This test case has been used by (amongst others):

*with same parameters as this script*
Wang, Shu, Yee, Sjoegreen (2012)
J. Comput. Phys. 231:pp. 190-214

Yee, Kotov, Wang and Shu (2013)
J. Comput Phys. 241:pp. 266--291

Kotov, Yee, Panesi, Prabhu and Wray (2014)
J. Comput. Phys. 269:pp. 215--233

*with different parameters*

Helzel, LeVeque and Warnecke (2000)
SIAM J. Sci. Comput. 22(4):pp 1489--1510

Tosatto and Vigevano (2008)
J. Comput. Phys. 227:pp. 2317--2343

--]]

config.title = "Arrhenius 1-D C-J detonation wave"
config.dimensions = 2

nsp, nmodes, gm = setGasModel('ideal-gas-AB-model.lua')
config.reacting = true
config.T_frozen = 0.0 -- Set this very low since initial temperature in unburnt gas is 1.0
config.reactions_file = "ideal-gas-AB-model.lua"

-- Sketch of domain with initial condition
--[[

   +-------------------------------------+
   |           ||                        |
   | burnt gas ||      unburnt gas       |
   |           ||                        |
   +-------------------------------------+
  x=0         x=10                      x=30
               | initial position
               | of detonation front

--]]

-- Gas model parameters
dofile("ideal-gas-AB-model.lua")
q0 = IdealGasAB.q
gamma = IdealGasAB.gamma

-- Flow states
-- Unburnt gas
u_u = 0.0
Qu = GasState:new{gm}
Qu.massf["A"] = 1.0; Qu.massf["B"] = 0.0
Qu.rho = 1.0
Qu.p = 1.0
gm:updateThermoFromRHOP(Qu)
unburntState = FlowState:new{p=Qu.p, T=Qu.T, massf=Qu.massf, velx=u_u}
-- Burnt gas
p_u = Qu.p; rho_u = Qu.rho
b = -p_u - rho_u*q0*(gamma - 1)
c = p_u^2 + 2*(gamma - 1)*p_u*rho_u*q0/(gamma + 1)
p_b = -b + math.sqrt(b^2 - c)
numer = rho_u*(p_b*(gamma + 1) - p_u)
denom = gamma*p_b
rho_b = numer/denom
S_cj = (rho_u*u_u + math.sqrt(gamma*p_b*rho_b))/rho_u
u_b = S_cj - math.sqrt(gamma*p_b/rho_b)
Qb = GasState:new{gm}
Qb.massf['A'] = 0.0; Qb.massf['B'] = 1.0
Qb.p = p_b
Qb.rho = rho_b
gm:updateThermoFromRHOP(Qb)
burntState = FlowState:new{p=Qb.p, T=Qb.T, massf=Qb.massf, velx=u_b}
print("")
print("Conditions of burnt gas:")
print("u_b= ", u_b)
print("rho_b= ", rho_b)
print("p_b= ", p_b)
print("T_b= ", Qb.T)
print("")

detFrontInit = 10 -- initial position of detonation front
function fillFn(x, y, z)
   if (x < detFrontInit) then
      return burntState
   else
      return unburntState
   end
end

-- Geometry
L = 30
h = 1
quad = CoonsPatch:new{p00=Vector3:new{x=0.0, y=0.0},
                      p10=Vector3:new{x=L,   y=0.0},
                      p11=Vector3:new{x=L,   y=h  },
                      p01=Vector3:new{x=0.0, y=h  }}

-- Grid
nxcells = 200
nycells = 2
grid = StructuredGrid:new{psurface=quad, niv=nxcells+1, njv=nycells+1}

-- Block
-- Note: At the left-end of the domain we feed in burnt gas at the appropriate
-- conditions. The right-end of domain is a wall, so we intend to terminate
-- the simulation before anything interesting interacts at the wall.
FluidBlock:new{grid=grid, initialState=fillFn,
               bcList={west=InFlowBC_Supersonic:new{flowState=burntState}}}

config.flux_calculator = "ausmdv"
config.max_time = 1.8
config.max_step = 100000
config.dt_init = 0.5e-3
config.fixed_time_step = true
config.dt_max = 1.0
config.dt_plot = config.max_time/10
