-- Author: Rowan J. Gollan
-- Date: 2020-09-11
--
-- This script simulates the Busemann inlet depicted in Fig. 1 in Moelder and Szpiro.
--
-- Reference:
-- Moelder and Szpiro (1966)
-- Busemann inlet for hypersonic speeds.
-- Journal of Spacecraft, 3(8):1303--1304
--

config.title = "Full Busemann inlet with M_inf=7.12, M_2=3.00, M_3=2.27, and theta_s=17.2deg"
config.dimensions = 2
config.axisymmetric = true

M_inf = 7.12
p_inf = 1000.0
T_inf = 350.0
setGasModel('ideal-air-gas-model.lua')
dummy = FlowState:new{p=p_inf, T=T_inf}
vx_inf = M_inf*dummy.a
inflow = FlowState:new{p=p_inf, T=T_inf, velx=vx_inf}


contour = Spline2:new{filename='bwall-spline.inp'}
cnt0 = contour(0.0)
cnt1 = contour(1.0)
h = cnt1.y
isoL = 4.0*h
thrt = Line:new{p0=cnt1, p1=Vector3:new{x=cnt1.x + isoL, y=cnt1.y}}
wall_constr = Polyline:new{segments={contour, thrt}}
wall = ArcLengthParameterizedPath:new{underlying_path=wall_constr}
a = Vector3:new{x=cnt0.x, y=0.0}
b = Vector3:new{x=wall(1.0).x, y=0.0}
inBndry = Line:new{p0=a, p1=cnt0}
outBndry = Line:new{p0=b, p1=wall(1.0)}
axis = Line:new{p0=a, p1=b}


nic = 400
njc = 50

quad0 = makePatch{north=wall, east=outBndry, south=axis, west=inBndry}
grid0 = StructuredGrid:new{psurface=quad0, niv=nic+1, njv=njc+1}
--[=[
r_grid = {{0.0, 1/3, 2/3, 1.0},
          {0.0, 1/3-0.3, 2/3+0.3, 1.0},
          {0.0, 1/3-0.3, 2/3+0.3, 1.0},
          {0.0, 1/3, 2/3, 1.0}}
s_grid = {{0.0, 0.0, 0.0, 0.0},
          {1.0/3, 1.0/3-0.3, 1.0/3-0.3, 1.0/3},
          {2.0/3, 2.0/3+0.3, 2.0/3+0.3, 2.0/3},
          {1.0, 1.0, 1.0, 1.0}}
grid_constr = StructuredGrid:new{psurface=quad0, niv=nic+1, njv=njc+1, r_grid=r_grid, s_grid=s_grid}
--r_grid, s_grid = grid_constr:determine_rs_grids(quad0)
--grid0 = StructuredGrid:new{psurface=quad0, niv=nic+1, njv=njc+1, r_grid=r_grid, s_grid=s_grid}
--]=]

blk0 = FluidBlock:new{grid=grid0, initialState=inflow,
     bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
	     east=OutFlowBC_Simple:new{}}
}

config.max_time = 5.0*(wall(1.0).x/vx_inf)
config.dt_init = 1.0e-6
config.dt_plot = config.max_time / 10.0
config.max_step = 5000


