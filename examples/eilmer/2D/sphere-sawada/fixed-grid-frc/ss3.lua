-- ss3.lua
--
-- Sphere in equilibrium air modelling Case 3 from
--    K. Sawada & E. Dendou (2001)
--    Validation of hypersonic chemical equilibrium flow calculations
--    using ballistic-range data.
--    Shock Waves (2001) Vol. 11, pp 43--51
--
-- Experimental shock stand-off distance is 2.59mm
-- Sawada & Dendou CFD value:               2.56mm
-- ceq value:                               2.54mm
--
-- Authors: PAJ and RJG
-- Versions:
--    Tcl version: 22-Jan-2004, derived from rbody.
--    Python version: ss3.py, 04-Apr-2005, 10-Aug-2006, 27-Nov-2006
--                    12-Nov-2008 by RJG for use in Elmer3
--    Lua version: 22-Aug-2016 by PJ for use in Eilmer4
--                 2019-05-26 shape grid with Billig correlation.
--                 2021-06-25 use finite rate chemistry

config.title = "Sphere in hypersonic flow of air in chemical equilibrium."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
config.cfl_value = 0.5

nsp, nmodes, gm = setGasModel('air-5sp.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
-- Free-stream flow definition
p_inf = 20.0e3   -- Pa
T_inf = 296.0    -- degrees K
vx_inf = 4.68e3  -- flow speed, m/s
massf_inf = {O2=0.233, N2=0.767}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=vx_inf, massf=massf_inf}
initial = FlowState:new{p=0.3*p_inf, T=T_inf, massf=massf_inf}

-- Reacting flow needs a reaction file
config.reactions_file = 'air-5sp-6r.lua'
config.reacting = true

print("Building grid.")
local billig_patch = require "billig_patch"
R = 31.8e-3  -- radius of sphere, in metres
M_inf = vx_inf/inflow.a
print("M_inf=", M_inf)
bp = billig_patch.make_patch{Minf=M_inf, R=R, scale=0.8}
cf_circum = RobertsFunction:new{end0=true, end1=false, beta=1.02}
grid = StructuredGrid:new{psurface=bp.patch, niv=61, njv=61,
			  cfList={west=cf_circum, east=cf_circum}}
blk0 = FBArray:new{grid=grid, initialState=initial, label="blk",
		       bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
			       north=OutFlowBC_Simple:new{}}, 
		       nib=2, njb=2}
-- We have left east and south as (default) slip-walls

-- Set a few more config options
config.max_time = 15.0*R/vx_inf -- allow time to settle at nose
config.max_step = 10000000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/8
