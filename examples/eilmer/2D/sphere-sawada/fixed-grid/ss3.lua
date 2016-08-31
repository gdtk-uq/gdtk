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
-- The grid is a bit wasteful because the shock lies close to
-- the body for equilibrium air, however, this grid layout 
-- (as used in rbody) allows us to play with perfect-gas models
-- without hitting the inflow boundary with the shock.
--
-- Authors: PAJ and RJG
-- Versions:
--    Tcl version: 22-Jan-2004, derived from rbody.
--    Python version: ss3.py, 04-Apr-2005, 10-Aug-2006, 27-Nov-2006
--                    12-Nov-2008 by RJG for use in Elmer3
--    Lua version: 22-Aug-2016 by PJ for use in Eilmer4

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

print "Building grid."
R = 31.8e-3  -- radius of sphere, in metres
deg2rad = math.pi / 180.0
alpha1 = 135.0 * deg2rad
alpha2 = 50.8 * deg2rad
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-1.0*R, y=0.0}
c = Vector3:new{x=math.cos(alpha1)*R, y=math.sin(alpha1)*R}
d = Vector3:new{x=0.0, y=R}
e = Vector3:new{x=math.cos(alpha2)*R, y=math.sin(alpha2)*R}
f = Vector3:new{x=1.4*R, y=1.5*R}
g = Vector3:new{x=1.5*R, y=2.5*R}
h = Vector3:new{x=1.5*R, y=3.5*R}
i = Vector3:new{x=-1.5*R, y=0.0}
j = Vector3:new{x=-1.5*R, y=1.5*R}
k = Vector3:new{x=-1.0*R, y=2.8*R}

bc = Arc:new{p0=b, p1=c, centre=a}
cd = Arc:new{p0=c, p1=d, centre=a}
de = Arc:new{p0=d, p1=e, centre=a}

psurf = makePatch{north=Bezier:new{points={h, g, f, e}},
		  east=Polyline:new{segments={bc, cd, de}},
		  south=Line:new{p0=i, p1=b},
		  west=Bezier:new{points={i, j, k, h}}}
cf_circum = RobertsFunction:new{end0=true, end1=false, beta=1.1}
cf_radial = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grid = StructuredGrid:new{psurface=psurf, niv=61, njv=61,
			  cfList={west=cf_circum, east=cf_circum,
				  south=cf_radial, north=cf_radial}}
blk0 = SBlockArray{grid=grid, fillCondition=initial, label="blk",
		   bcList={west=InFlowBC_Supersonic:new{flowCondition=inflow},
			   north=OutFlowBC_Simple:new{}}, 
		   nib=1, njb=4}
-- We have left east and south as (default) slip-walls

-- Set a few more config options
config.flux_calculator = "adaptive"
config.max_time = 10.0*R/vx_inf -- allow time to settle at nose
config.max_step = 40000
config.dt_init = 1.0e-9
config.cfl_value = 0.5 
config.dt_plot = config.max_time/4

dofile("sketch-domain.lua")
