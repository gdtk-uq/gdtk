-- An M-flow test case.
--
-- This example should match the weak-shock example in
-- Moelder (1967). The M-flow is an internal axisymmetric
-- flow. It begins with a conical shock and has a wall shape
-- to maintain the conicity of the shock.
--
-- Author: RJG
-- Date: 2025-06-15
--

theta_s = math.rad(140.0)

config.dimensions = 2
config.axisymmetric = true
config.solver_mode = 'steady'

-- Read contour points as spline, then determine other geometry
fname = 'm-flow-contour.dat'
contour = Spline2:new{filename=fname}

-- construction points
d = contour(0.0)
c = contour(1.0)
Lx = (c.x - d.x)
y_shock = Lx*math.tan(math.pi - theta_s)
y_frac = 0.05
a = {x=d.x, y=(1 - y_frac)*d.y}
b = {x=c.x, y=d.y - 1.1*y_shock}

-- boundaries
ab = Line:new{p0=a, p1=b}
ad = Line:new{p0=a, p1=d}
bc = Line:new{p0=b, p1=c}

-- patch
quad = CoonsPatch:new{north=contour, south=ab, east=bc, west=ad}

-- grid
cf_n = RobertsFunction:new{end0=false, end1=true, beta=1.05}
cf_t = RobertsFunction:new{end0=true, end1=false, beta=1.1}
nx = 50
ny = 50
grid = registerFluidGrid{
   grid = StructuredGrid:new{psurface=quad, niv=nx+1, njv=ny+1,
                             cfList={west=cf_n, east=cf_n, south=cf_t, north=cf_t}},
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
   stop_on_relative_residual = 1.0e-6,
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
