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

title = "Sphere in hypersonic flow of air in chemical equilibrium."
config.dimensions = 2
config.axisymmetric = true
config.solver_mode = 'transient'
config.flux_calculator = "hanel"
config.grid_motion = "shock_fitting"
config.gasdynamic_update_scheme = "moving_grid_2_stage"

nsp, nmodes, gm = setGasModel('air-5sp-eq.lua')

-- Free-stream flow definition
p_inf = 20.0e3   -- Pa
T_inf = 296.0    -- degrees K
vx_inf = 4.68e3  -- flow speed, m/s
inflow = FlowState:new{p=p_inf, T=T_inf, velx=vx_inf}
initial = FlowState:new{p=0.3*p_inf, T=T_inf}

flowDict = {
   initial=inflow,
   inflow=inflow
}

local billig_patch = require "billig_patch"
R = 31.8e-3  -- radius of sphere, in metres
M_inf = vx_inf/inflow.a
print("M_inf=", M_inf)
bp = billig_patch.make_patch{Minf=M_inf, R=R, scale=0.8}
cf_east = RobertsFunction:new{end0=true, end1=false, beta=1.1}
cf_west = RobertsFunction:new{end0=true, end1=false, beta=1.05}

grid = StructuredGrid:new{psurface=bp.patch, niv=32, njv=128,
			  cfList={west=cf_west, east=cf_east}}

registerFluidGridArray{
   grid=grid,
   nib=1, njb=4,
   fsTag='initial',
   shock_fitting=true,
   bcTags={west='inflow_sf', north='outflow', east='wall'}
}

bcDict = {
   inflow_sf=InFlowBC_ShockFitting:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{},
   wall=WallBC_WithSlip0:new{}
}

makeFluidBlocks(bcDict, flowDict)

body_flow_time = R/vx_inf
config.cfl_value = 0.5
config.interpolation_delay = 3*body_flow_time
config.shock_fitting_delay = body_flow_time
config.max_time = 10.0*body_flow_time -- allow time to settle at nose
config.max_step = 10000000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/16
