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
--
-- Authors: PAJ and RJG
-- Versions:
--    Tcl version: 22-Jan-2004, derived from rbody.
--    Python version: ss3.py, 04-Apr-2005, 10-Aug-2006, 27-Nov-2006
--                    12-Nov-2008 by RJG for use in Elmer3
--    Lua version: 22-Aug-2016 by PJ for use in Eilmer4
--                 2019-05-26 shape grid with Billig correlation.

config.title = "Sphere in hypersonic flow of air in chemical equilibrium."
print(config.title)
config.dimensions = 2
config.axisymmetric = true

nsp, nmodes, gm = setGasModel('cea-lut-air.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
-- Free-stream flow definition
p_inf = 20.0e3   -- Pa
T_inf = 296.0    -- degrees K
vx_inf = 4.68e3  -- flow speed, m/s
inflow = FlowState:new{p=p_inf, T=T_inf, velx=vx_inf}
initial = FlowState:new{p=0.3*p_inf, T=T_inf}

print("Building grid.")
require "billig_patch"
R = 31.8e-3  -- radius of sphere, in metres
M_inf = vx_inf/inflow.a
print("M_inf=", M_inf)
bp = billig_patch.make_patch{Minf=M_inf, R=R, scale=0.8}
cf_circum = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grid = StructuredGrid:new{psurface=bp.patch, niv=61, njv=61,
			  cfList={west=cf_circum, east=cf_circum}}
blk0 = FluidBlockArray{grid=grid, initialState=initial, label="blk",
		       bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
			       north=OutFlowBC_Simple:new{}}, 
		       nib=1, njb=4}
-- We have left east and south as (default) slip-walls

-- Set a few more config options
config.max_time = 10.0*R/vx_inf -- allow time to settle at nose
config.max_step = 40000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/4

dofile("sketch-domain.lua")
