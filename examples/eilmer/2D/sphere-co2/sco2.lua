-- sco2.lua
-- To model the experiment described in
-- D. Liao et al. (2020)
-- Measurement and numerical simulation of shock standoff distances
-- over hypersonic spheres in CO_2 in a ballistic range.
-- Shock Waves 30:131–138
--
-- PJ 2021-08-12
--
config.title = "Sphere fired into carbon dioxide, shock-fitted boundary."
print(config.title)
-- Parameters for test S-2
V_inf = 2845.0  -- m/s
p_inf = 7.43e3  -- Pa
T_inf = 293.2   -- K
radius = 0.005  -- m
-- Experimental value of normalized shock stand-off distance
-- delta/radius = 0.0705
--
nsp, nmodes, gm = setGasModel('co2-5sp-2T.gas')
config.reactions_file = 'co2-5sp-2T-chemistry.chem'
config.energy_exchange_file = 'co2-mixture-energy-exchange.kin'
initial = FlowState:new{p=p_inf/3.0, T=T_inf, T_modes={T_inf},
                        velx=0.0, massf={CO2=1.0}}
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf},
                       velx=V_inf, massf={CO2=1.0}}
--
config.dimensions = 2
config.axisymmetric = true
config.viscous = true
config.reacting = true
--
print "Building grid."
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-radius, y=0.0}
c = Vector3:new{x=0.0, y=radius}
--
d = Vector3:new{x=-1.5*radius, y=0}
e = Vector3:new{x=-1.5*radius, y=radius}
f = Vector3:new{x=-radius, y=2.0*radius}
g = Vector3:new{x=0.0, y=3.0*radius}
--
psurf = makePatch{north=Line:new{p0=g, p1=c},
		  east=Arc:new{p0=b, p1=c, centre=a},
		  south=Line:new{p0=d, p1=b},
		  west=Bezier:new{points={d, e, f, g}}}
-- Cluster toward the shock and the wall to capture both
-- the post-shock relaxation zone and the viscous boundary layer.
rcfx = RobertsFunction:new{end0=true, end1=true, beta=1.10}
grid = StructuredGrid:new{psurface=psurf, niv=41, njv=81,
                          cfList= {south=rcfx,north=rcfx}}
--
blk = FBArray:new{
   grid=grid, initialState=initial,
   bcList={west=InFlowBC_ShockFitting:new{flowState=inflow},
           east=WallBC_NoSlip_FixedT:new{Twall=T_inf, group="loads"},
           north=OutFlowBC_Simple:new{}},
   nib=2, njb=2
}
identifyBlockConnections()
--
-- Set a few more config options
bft = (radius*2)/V_inf
print("bft=", bft, "(characteristic time for flow over body)")
config.flux_calculator = "ausmdv"
-- config.gasdynamic_update_scheme = "backward_euler"
config.gasdynamic_update_scheme = "moving-grid-1-stage"
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 1*bft
config.cfl_schedule = {{0.0,0.25}, {4*bft,0.25}, {5*bft,0.5}}
config.dt_init = 1.0e-9
config.reaction_time_delay = 1*bft
config.reaction_fraction_schedule = {{1*bft, 0.01}, {2*bft, 0.1}, {3*bft, 1.0}}
-- Just to get this exercise up and running, set the following
-- viscous parameters.  Would like to eventially use the default
-- "least_squares" and "cells".
config.viscous_delay = 3*bft
config.viscous_signal_factor = 0.1
config.spatial_deriv_calc = "divergence"
config.spatial_deriv_locn = "vertices"
config.max_time = 20*bft
config.max_step = 400000
config.dt_plot = bft
config.write_loads = true
config.dt_loads = bft
config.max_invalid_cells = 10
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
