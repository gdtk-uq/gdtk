-- lehr.lua
-- Author: Rowan J. Gollan & Peter J.
-- Date: 2020-01-27
--
-- This script is used to setup a simlulation of Lehr's
-- hemispherical projectile fired into a detonable gas.
--
-- Reference:
--   Lehr, H. (1972)
--   Acta Astronautica, 17, pp.589--597
--
config.title = "Lehr experiment M=4.79"
R = 7.5e-3 -- nose radis, metres

-- free stream conditions
-- sound-speed 403m/s taken from Table 1 in Lehr (1972)
-- Mach number 4.79 to match example used by Greg Wilson.
nsp, nmodes, gmodel = setGasModel('Rogers-Schexnayder-gas-model.lua')
p_inf = 320.0/760.0*101325.0 -- Pa
u_inf = 1931 -- m/s 
T_inf = 292 -- K
molef_inf = {H2=2/6.76, O2=1/6.76, N2=3.76/6.76} -- stoichiometric mix
massf_inf = gmodel:molef2massf(molef_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, massf=massf_inf}
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0, massf=massf_inf}

-- Now set some configuration options
config.axisymmetric = true
config.reacting = true
config.reactions_file = 'Rogers-Schexnayder-reac-file.lua'
config.reaction_time_delay = 0.0
config.flux_calculator = "ausmdv"
config.max_time = 100.0e-6 -- allow time to establish
config.max_step = 800000
config.dt_init = 1.0e-10
config.cfl_value = 0.4
config.dt_plot = config.max_time/100.0

-- Set up the geometry for defining the grid
a = {x=0.0, y=0.0}
b = {x=-R, y=0.0}
c = {x=0.0, y=R}
d = {{x=-1.3*R,y=0.0}, {x=-1.3*R,y=0.7*R}, {x=-0.87*R,y=1.4*R}, {x=0.0,y=2.1*R}}
-- Set up surface and grid
psurf = CoonsPatch:new{
   north=Line:new{p0=d[#d], p1=c}, east=Arc:new{p0=b, p1=c, centre=a},
   south=Line:new{p0=d[1], p1=b}, west=Bezier:new{points=d}
}
-- Set up some mild clustering towards body so that the heat
-- release zone can be resolved
-- cluster = RobertsFunction:new{end0=false, end1=true, beta=1.2}
ni = 200; nj = 300
grid = StructuredGrid:new{psurface=psurf, niv=ni+1, njv=nj+1,}

blk = FBArray:new{grid=grid, fillCondition=initial, label='blk',
		  bcList={west=InFlowBC_Supersonic:new{flowCondition=inflow},
			  north=OutFlowBC_Simple:new{}},
                  nib=2, njb=4}

----------------------------------------------------------------------------
-- History:
-- 27-Feb-2010
-- Adapted bits from sphere-heat-transfer and Rowan's mbcsn2/lehr_sphere.
-- This is the continuation of the simulation on a better fitted grid.
-- 27-01-2020
-- Eilmer4 flavour.
