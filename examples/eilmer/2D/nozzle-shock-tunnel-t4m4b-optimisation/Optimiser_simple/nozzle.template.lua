-- nozzle-template.lua
-- Optimize an axisymmetric RST nozzle for T4
-- Daniel Smith 18-Jun-2020
-- This script was built from a combination of work done by Peter Jacobs, 
-- Tamara Sopek, and Wilson Chan. 

config.title = "T4 Mach7 nozzle with air in chemical equilibrium."
print(config.title)

-- Set a few more config options
config.dimensions = 2
config.axisymmetric = true
config.viscous = false -- true
config.cfl_value = 0.4 -- default value is 0.5
config.flux_calculator = "adaptive"
config.max_time = 0.001 -- allow enough to reach steady state
config.max_step = 1000000
config.dt_init = 1.0e-11

-- To accelerate the calculation...
config.block_marching = true -- to accelerate calculations, utilise block-marching capability
config.nib = 31
config.njb = 2
config.propagate_inflow_data = true

--inflow conditions
gm = GasModel:new{'cea-air5species-gas-model.lua'}
state0 = GasState:new{gm}
state0.p = 4.56e6
state0.T = 1861.0
throat_V = 832
gm:updateThermoFromPT(state0); gm:updateTransCoeffs(state0)

throat_mu = state0.mu
throat_rho = state0.rho

print("Flow velocity at nozzle throat = ", throat_V)
print("Dynamic viscosity at nozzle throat = ", throat_mu)
print("M_throat=", throat_V/state0.a)

-- setup simulation
print("Second, set up the simulation grid and the initial flow state.")
-- For the simulation, we are going to use the tabulated gas behaviour
-- because it will be much faster than the CEA-backed gas model.
nsp, nmodes, gm = setGasModel('cea-lut-air.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)

-- The inflow for the simulation is at the nozzle throat.
inflow = FlowState:new{p=state0.p, T=state0.T, velx=throat_V}
print("Species mass fractions at throat, according to CEA gas model")
for k, v in pairs(state0.ceaSavedData.massf) do
   print(string.format("massf[%s] = %g", k, v))
end

-- Read coordinates from the Bezier-control-pts-t4-m7.data.
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

-- now we want to adjust those bezier points based on inputs from optimiser
-- reading those in here

dy1 = $dy1
dy2 = $dy2
dy3 = $dy3
dy4 = $dy4
dy5 = $dy5
dy6 = $dy6
dy7 = $dy7

dy = {0.0,0.0,dy1,dy2,dy3,dy4,dy6,dy6,dy7} -- first two bezier points are always unchanged

-- then, populate table containing Bezier points.
bez_points = {}
for npoints=1,#bezx do
   bez_points[#bez_points+1] = Vector3:new{x=bezx[npoints], y=bezy[npoints] + dy[npoints]}
end
-- and define some critical geometrical points of this nozzle.
L_thr = 0.0125             -- axial length of throat
x_start = bez_points[1].x  -- start of supersonic expansion
y_throat = bez_points[1].y -- throat radius

-- The throat-region is the constant-area section up to the
-- start of the conical expansion.
throat_region = CoonsPatch:new{
   p00=Vector3:new{x=-L_thr, y=0.0},
   p10=Vector3:new{x=x_start, y=0.0},
   p11=Vector3:new{x=x_start, y=y_throat},
   p01=Vector3:new{x=-L_thr, y=y_throat}}
-- Supersonic expansion is fully defined by its north edge.
exp_region = NozzleExpansionPatch:new{north=Bezier:new{points=bez_points}}
-- Define structured grids for both regions.
print "Building grid."
x_cf_throat = RobertsFunction:new{end0=true, end1=true, beta=1.45792369398}
-- Note that the Bezier wall has a natural clustering of points toward the throat.
-- The following weak cluster function x_cf mitigates that clustering.
x_cf = RobertsFunction:new{end0=false, end1=true, beta=1.5}
y_cf = RobertsFunction:new{end0=false, end1=true, beta=1.02}
throat_grid = StructuredGrid:new{psurface=throat_region, niv=12, njv=21}--,cfList={west=y_cf, east=y_cf,south=x_cf_throat, north=x_cf_throat}}
exp_grid = StructuredGrid:new{psurface=exp_region, niv=151, njv=21}--,cfList={west=y_cf, east=y_cf,south=x_cf, north=x_cf}}
-- Divide the full domain into many blocks, mainly in the axial direction.
throat_fba = FBArray:new{grid=throat_grid, initialState=inflow, label="throat",
                         bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                                 north=WallBC_NoSlip_FixedT:new{Twall=300.0}},
                         nib=1, njb=2}
exp_fba = FBArray:new{grid=exp_grid, initialState=inflow, label="expansion",
                      bcList={north=WallBC_NoSlip_FixedT:new{Twall=300.0},
                              east=OutFlowBC_Simple:new{}},
                      nib=30, njb=2}
identifyBlockConnections()