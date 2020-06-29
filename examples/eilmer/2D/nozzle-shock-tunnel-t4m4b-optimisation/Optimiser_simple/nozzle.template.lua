-- nozzle-template.lua
-- Optimize the shape of an axisymmetric Mach 4 nozzle for T4
-- 18-Jun-2020
-- This script was built from a combination of work done by:
-- Daniel Smith, Peter Jacobs, Tamara Sopek, and Wilson Chan.

config.title = "T4 Mach4 nozzle with air in chemical equilibrium."
print(config.title)

config.dimensions = 2
config.axisymmetric = true
config.viscous = false
config.cfl_value = 0.4
config.flux_calculator = "adaptive"
config.max_time = 0.001 -- allow enough to reach steady state
config.max_step = 1000000
config.dt_init = 1.0e-11

-- To accelerate the calculation, use block-marching.
config.block_marching = true
config.nib = 31
config.njb = 2
config.propagate_inflow_data = true

-- For the simulation, we are going to use the tabulated gas behaviour
-- because it will be much faster than the CEA-backed gas model.
nsp, nmodes, gm = setGasModel('cea-lut-air.lua')
state0 = GasState:new{gm}

-- The inflow for the simulation is at the nozzle throat.
state0.p = 4.56e6
state0.T = 1861.0
throat_V = 832
gm:updateThermoFromPT(state0); gm:updateTransCoeffs(state0)
throat_mu = state0.mu
throat_rho = state0.rho
print("M_throat=", throat_V/state0.a)
inflow = FlowState:new{p=state0.p, T=state0.T, velx=throat_V}

-- Nozzle profile is built from a Bezier curve...
-- Read original shape as coordinates from the Bezier-control-pts-t4-m7.data.
local file = io.open("Bezier-control-pts-t4m4b-initial.data", "r")
local bezx = {}; bezy = {}
for line in file:lines() do
    local t = {}
    if line:sub(1,1)=="#" then -- line starts with "#"
       -- do nothing
    else
       for s in line:gmatch("%S+") do
           t[#t+1] = tonumber(s)
       end
       bezx[#bezx+1] = t[1]; bezy[#bezy+1] = t[2]
    end
end
file:close()
-- Adjust those points based on data from optimiser,
-- which are substituted for the symbolic names below.
dy = {0.0, 0.0, $dy1, $dy2, $dy3, $dy4, $dy5, $dy6, $dy7}
-- then, populate table containing Bezier points.
bpoints = {}
for i=1, #bezx do
   bpoints[#bpoints+1] = {x=bezx[i], y=bezy[i]+dy[i]}
end
-- and define the throat geometry.
L_thr = 0.0125 -- axial length of throat
x_start = bpoints[1].x  -- start of supersonic expansion
y_throat = bpoints[1].y -- throat radius

-- The throat-region is the constant-area section up to the
-- start of the conical expansion.
throat_region = CoonsPatch:new{
   p00=Vector3:new{x=-L_thr, y=0.0},
   p10=Vector3:new{x=x_start, y=0.0},
   p11=Vector3:new{x=x_start, y=y_throat},
   p01=Vector3:new{x=-L_thr, y=y_throat}}
-- Supersonic expansion is fully defined by its north edge.
nozzle_profile = ArcLengthParameterizedPath:new{underlying_path=Bezier:new{points=bpoints}}
exp_region = NozzleExpansionPatch:new{north=nozzle_profile}
-- Define structured grids for both regions.
print("Building grid.")
x_cf_throat = RobertsFunction:new{end0=true, end1=true, beta=1.48}
x_cf = RobertsFunction:new{end0=true, end1=false, beta=1.1}
y_cf = nil -- for inviscid simulations
if config.viscous then
   y_cf = RobertsFunction:new{end0=false, end1=true, beta=1.02}
end
throat_grid = StructuredGrid:new{
   psurface=throat_region, niv=12, njv=21,
   cfList={west=y_cf, east=y_cf,south=x_cf_throat, north=x_cf_throat}
}
exp_grid = StructuredGrid:new{
   psurface=exp_region, niv=151, njv=21,
   cfList={west=y_cf, east=y_cf,south=x_cf, north=x_cf}
}
-- Divide the full domain into many blocks, mainly in the axial direction.
-- We are going to run with 2 MPI tasks, so set njb=2.
throat_fba = FBArray:new{grid=throat_grid, initialState=inflow, label="throat",
                         bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                                 north=WallBC_NoSlip_FixedT:new{Twall=300.0}},
                         nib=1, njb=2}
exp_fba = FBArray:new{grid=exp_grid, initialState=inflow, label="expansion",
                      bcList={north=WallBC_NoSlip_FixedT:new{Twall=300.0},
                              east=OutFlowBC_Simple:new{}},
                      nib=30, njb=2}
identifyBlockConnections()
