-- lehr.lua
-- Author: Rowan J. Gollan
-- Date: 2016-01-27
--
-- This script is used to setup a simlulation of Lehr's
-- hemispherical projectile fired into a detonable gas.
--
-- Reference:
--   Lehr, H. (1972)
--   Acta Astronautica, 17, pp.589--597
--
config.title = "Lehr experiment M=3.55"
R = 7.5e-3 -- nose radis, metres

-- free stream conditions
-- taken from Table 1, Entry 5
nsp, nmodes, gmodel = setGasModel('Rogers-Schexnayder-gas-model.lua')
p_inf = 186.0/760.0*101325.0 -- Pa
u_inf = 1892 -- m/s
T_inf = 292 -- K
molef_inf = {H2=2/3, O2=1/3}
massf_inf = gmodel:molef2massf(molef_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, massf=massf_inf}
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0, massf=massf_inf}

-- Now set some configuration options
body_flow_length = R/u_inf
t_final = 10 * body_flow_length -- allow time to establish
ni = 128; nj = 256
config.axisymmetric = true
config.reacting = true
config.reactions_file = 'Rogers-Schexnayder-reac-file.lua'
config.reaction_time_delay = 0.0
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 3 * body_flow_length  
config.max_time = t_final
config.max_step = 800000
config.dt_init = 1.0e-10
config.cfl_value = 0.4
config.dt_plot = config.max_time/100.0

-- Set up the geometry for defining the grid
a = {x=0.0, y=0.0}
b = {x=-R, y=0.0}
c = {x=0.0, y=R}
d = {{x=-1.5*R,y=0.0}, {x=-1.5*R,y=R}, {x=-R,y=2*R}, {x=0.0,y=3*R}}
-- Set up surface and grid
psurf = CoonsPatch:new{
   north=Line:new{p0=d[#d], p1=c}, east=Arc:new{p0=b, p1=c, centre=a},
   south=Line:new{p0=d[1], p1=b}, west=Bezier:new{points=d}
}
-- Set up some mild clustering towards body so that the heat
-- release zone can be resolved
-- cluster = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grid = StructuredGrid:new{psurface=psurf, niv=ni+1, njv=nj+1,}

blk = FBArray:new{grid=grid, fillCondition=initial, label='blk',
		  bcList={west=InFlowBC_ShockFitting:new{flowCondition=inflow},
			  north=OutFlowBC_Simple:new{}},
                  nib=1, njb=4}

----------------------------------------------------------------------------
-- History:
--   2016-01-27 : RJG updated for eilmerD Lua input
--   2015-03-15 : RJG re-worked eilmer3 example for M = 3.55
--   2010-02-27 : PJ adapted for eilmer3
--                bits taken from sphere-heat-transfer and
--                mbcns2/lehr_sphere
