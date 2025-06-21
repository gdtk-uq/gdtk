-- An M-flow test case.
--
-- This example should match the weak-shock example in
-- Moelder (1967). The M-flow is an internal axisymmetric
-- flow. It begins with a conical shock and has a wall shape
-- to maintain the conicity of the shock.
--
-- Author: RJG
-- Date: 2025-06-15 -- initial build
--       2025-06-21 -- change to more orthogonal domain boundaries

local rad = math.rad
local pi = math.pi
local tan = math.tan
local atan2 = math.atan2

theta_s = rad(140.0)

config.dimensions = 2
config.axisymmetric = true
config.solver_mode = 'steady'

-- Read contour points as spline, then determine other geometry
fname = 'm-flow-contour.dat'
contour = Spline2:new{filename=fname}

-- construction points
d = contour(0.0)
c = contour(1.0)
dp = contour(0.001) -- perturbed slightly, to compute slope
lip_angle = atan2(dp.y - d.y, dp.x - d.x)
cp = contour(0.999) -- perturbed slightly
tail_angle = atan2(c.y - cp.y, c.x - cp.x)
Lx = (c.x - d.x)
y_lip = 0.02*Lx
a_y = (1 - y_lip)*d.y
dx_a = (d.y - a_y)*tan(lip_angle)

alpha = (theta_s - pi)
beta = tail_angle + pi/2

x_shock = (d.x*tan(alpha) - d.y - c.x*tan(beta) + c.y)/(tan(alpha) - tan(beta))
y_shock = tan(alpha)*(x_shock - d.x) + d.y
dy_shock = c.y - y_shock
b_y = c.y - 1.1*dy_shock
b_x = ((b_y - c.y)/tan(beta)) + c.x


a = {x=d.x+dx_a, y=a_y}
b = {x=b_x, y=b_y}

-- boundaries
ab = Line:new{p0=a, p1=b}
ad = Line:new{p0=a, p1=d}
bc = Line:new{p0=b, p1=c}

-- patch
quad = ControlPointPatch:new{north=contour, south=ab, east=bc, west=ad,
                             ncpi=10, ncpj=4, guide_patch="channel"}

-- grid
cf_n = GeometricFunction:new{a=0.04, r=1.1, N=20, reverse=true}
cf_t = GeometricFunction:new{a=0.001, r=1.02, N=170}
nx = 200
ny = 20
grid = registerFluidGridArray{
   grid = StructuredGrid:new{psurface=quad, niv=nx+1, njv=ny+1,
                             cfList={north=cf_t, south=cf_t}},
   nib=8,
   njb=1,
   fsTag = "inflow",
   bcTags = {west="inflow", south="inflow", east="outflow"}
}

-- gas and flow conditions
setGasModel('ideal-air.gas')
p_inf = 1.0e5
T_inf = 300.0
M_inf = 5
dummy = FlowState:new{p=p_inf, T=T_inf}
u_inf = M_inf*dummy.a
print("a= ", dummy.a)
print("u_inf= ", u_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
flowDict = {inflow=inflow}
bcDict = {
  inflow = InFlowBC_Supersonic:new{flowState=inflow},
  outflow = OutFlowBC_Simple:new{}
}

-- now build blocks
makeFluidBlocks(bcDict, flowDict)

-- settings for the solver
config.flux_calculator= "ausmdv"
config.interpolation_order = 2
config.extrema_clipping = false

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-9,
   number_of_phases = 2,
   max_steps_in_initial_phases = { 50 },
   use_physicality_check = true,
   max_linear_solver_iterations = 50,
   total_snapshots = 3,
   steps_between_status = 1,
   steps_between_snapshots = 5,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.95,
   start_cfl = 2.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 0.9
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 20.0
}
