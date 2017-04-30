-- sphere.lua
-- This script is used to setup a simlulation of Rose and Stark's
-- hemispherical heat-transfer probe set in a shock tube.
--
-- PAJ and RJG, 2016-08-28
-- Parts of this file go way back to before 2010 with bits taken from
-- the mbcns2/lehr_sphere example.
--
config.title = "Rose and Stark experiment"
R = 6.6e-3 -- nose radis, metres

-- free stream conditions
-- taken from Table 1, Entry 5
nsp, nmodes, gmodel = setGasModel('cea-lut-air.lua')
p_init = 6.7 -- Pa
p_inf = 535.6 -- Pa
vx_inf = 2436.5 -- m/s
T_inf = 2573.5 -- K
inflow = FlowState:new{p=p_inf, T=T_inf, velx=vx_inf}
initial = FlowState:new{p=p_init, T=T_inf}

-- Now set some configuration options
body_flow_time = R/vx_inf
t_final = 10 * body_flow_time -- allow time to establish
ni = 32; nj = 64
config.axisymmetric = true
config.viscous = true
config.spatial_deriv_calc = "divergence"
config.spatial_deriv_locn = "vertices"
config.viscous_signal_factor = 0.1
config.viscous_delay = 2 * body_flow_time
config.flux_calc = AUSMDV
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 3 * body_flow_time  
config.max_time = t_final
config.max_step = 800000
config.dt_init = 1.0e-12
config.cfl_value = 0.4
config.dt_plot = config.max_time/10.0

-- Set up the geometry for defining the grid
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-R, y=0.0}
c = Vector3:new{x=0.0, y=R}
d = { Vector3:new{x=-1.5*R, y=0.0}, Vector3:new{x=-1.5*R, y=R},
      Vector3:new{x=-R, y=2*R}, Vector3:new{x=0.0, y=3*R} }
-- Set up surface and grid
sphere_edge = Arc:new{p0=b, p1=c, centre=a}
psurf = makePatch{north=Line:new{p0=d[#d], p1=c}, south=Line:new{p0=d[1], p1=b},
		  east=sphere_edge, west=Bezier:new{points=d}}
cf_radial = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grid = StructuredGrid:new{psurface=psurf, niv=ni+1, njv=nj+1,
			  cfList={north=cf_radial, south=cf_radial}}

blk = FluidBlockArray{grid=grid, fillCondition=initial, label='blk',
		      bcList={west=InFlowBC_ShockFitting:new{flowCondition=inflow},
			      east=WallBC_NoSlip_FixedT:new{Twall=296.0},
			      north=OutFlowBC_Simple:new{}},
		      nib=1, njb=4}

dofile("sketch-domain.lua")
