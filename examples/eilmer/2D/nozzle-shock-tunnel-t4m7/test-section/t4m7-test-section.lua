-- t4m7-test-section.lua
-- Tamara Sopek, 2020-21-08

config.title = "T4 test section for Mach7 nozzle flow with air in chemical equilibrium."
print(config.title)

-- Set a few more config options
config.dimensions = 2
config.axisymmetric = true
config.viscous = true
config.turbulence_model = "k_omega"
config.cfl_value = 0.4 -- default value is 0.5
config.flux_calculator = "adaptive"
config.max_time = 0.002 -- allow enough to reach steady state
config.max_step = 1000000
config.dt_init = 1.0e-11

-- To accelerate the calculation...
config.block_marching = true -- to accelerate calculations, utilise block-marching capability
config.nib = 31
config.njb = 4
config.propagate_inflow_data = true

-- Get flow conditions at the nozzle throat by performing ESTCj-like calculations,
-- using the CEA-backed gas model.
print("First, compute the flow conditions at the nozzle throat.")
print("shock-tube fill conditions")
gm = GasModel:new{'cea-air13species-gas-model.lua'}
state1 = GasState:new{gm}
state1.p = 131.0e3; state1.T = 300.0
gm:updateThermoFromPT(state1); gm:updateTransCoeffs(state1)
print("state1:"); printValues(state1)
--
print("normal shock, given shock speed")
Vs = 1644.0
state2, V2, Vg = gasflow.normal_shock(state1, Vs)
gm:updateThermoFromPT(state2); gm:updateTransCoeffs(state2)
print("V2=", V2, "Vg=", Vg)
print("state2:"); printValues(state2)
--
print("reflected shock")
state5, Vr = gasflow.reflected_shock(state2, Vg)
gm:updateThermoFromPT(state5); gm:updateTransCoeffs(state5)
print("Vr=", Vr)
print("state5:"); printValues(state5)
--
print("Expand from stagnation (with ratio of pressure to match observation)")
state5s, V5s = gasflow.expand_from_stagnation(state5, 8.32/23.61)
gm:updateThermoFromPT(state5s); gm:updateTransCoeffs(state5s)
print("V5s=", V5s, " Mach=", V5s/state5s.a)
print("state5s:"); printValues(state5s)
print("(h5s-h1)=", gm:enthalpy(state5s) - gm:enthalpy(state1)) 
--
print("Expand to throat condition (Mach 1.0001)")
state6, V6 = gasflow.expand_to_mach(state5s, 1.0001)
gm:updateThermoFromPT(state6); gm:updateTransCoeffs(state6)
print("V6=", V6, " Mach=", V6/state6.a)
print("state6:"); printValues(state6)

--   
throat_V = V6
print("Flow velocity at nozzle throat = ", throat_V)
throat_mu = state6.mu
print("Dynamic viscosity at nozzle throat = ", throat_mu)
throat_rho = state6.rho
--

--
print("Quasi-one-dimensional expansion to estimated test-flow condition.")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 27)
gm:updateThermoFromPT(state7); gm:updateTransCoeffs(state7)
print("V7=", V7, " Mach=", V7/state7.a)
print("state7:"); printValues(state7)

print("Second, set up the simulation grid and the initial flow state.")
-- Estimate turbulence quantities for the inflow. Note that for 
-- fully laminar nozzles these values are adjusted later once 
-- the grid has been defined.
throat_tke = 1.5 * (throat_V * 0.05)^2
throat_mu_t = 100.0 * throat_mu
throat_omega = throat_rho * throat_tke / throat_mu_t
print("Inflow turbulence: tke=", throat_tke, "omega=", throat_omega)
print("These will be adjusted later if nozzle flow is fully laminar.")

-- For the simulation, we are going to use the tabulated gas behaviour
-- because it will be much faster than the CEA-backed gas model.
nsp, nmodes, gm = setGasModel('cea-lut-air.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
-- The inflow for the simulation is at the nozzle throat.
inflow = FlowState:new{p=state6.p, T=state6.T, velx=V6, tke=throat_tke, omega=throat_omega}
print("Species mass fractions at throat, according to CEA gas model")
for k, v in pairs(state6.ceaSavedData.massf) do
   print(string.format("massf[%s] = %g", k, v))
end

-- Read coordinates from the Bezier-control-pts-t4-m7.data.
local file = io.open("Bezier-control-pts-t4-m7.data", "r")
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

-- then, populate table containing Bezier points.
bez_points = {}
for npoints=1,#bezx do
   bez_points[#bez_points+1] = Vector3:new{x=bezx[npoints], y=bezy[npoints]}
end
-- and define some critical geometrical points of this nozzle.
L_thr = 0.02                -- axial length of throat
x_start = bez_points[1].x   -- start of supersonic expansion
y_throat =  bez_points[1].y -- throat radius
x_exit = bez_points[#bez_points].x
y_exit = bez_points[#bez_points].y

-- Supersonic expansion
nozzle_profile = ArcLengthParameterizedPath:new{underlying_path=Bezier:new{points=bez_points}}
exp_region = NozzleExpansionPatch:new{north=nozzle_profile}

---------------------------------------------------------------
-- Additional lines to allow simulation of test section flow --
---------------------------------------------------------------
-- User-defined geometrical dimensions of the test section
ts_sim_length = 1.2  -- test section simulated length
ts_radius = 0.225  -- test section radius
config.max_time = 2.0e-3 -- allow enough to reach steady state
config.dt_plot = 0.1e-3
config.max_step = 10000000000
nozzle_soln_jobName = "t4m7"
nozzle_soln_tindx = 1
nozzle_soln_nBlocks = 124
nozzle_soln_dir = "../wall-functions/"
initial = FlowState:new{p=100.0, T=300.0, velx=0.0, tke=0.001, omega=1.0}

-- Switch off block-marching
config.block_marching = false 

-- Locate most south-eastern vertex to truncate the expansion region,
-- using the x-coordinate as the target for our secant root-finder
nozzle_soln = FlowSolution:new{jobName=nozzle_soln_jobName, tindx=nozzle_soln_tindx, 
                                nBlocks=nozzle_soln_nBlocks, dir=nozzle_soln_dir}
target = nozzle_soln:vtx{ib=116, i=0, j=0, k=0}.x

-- Use root-finder to locate r0 of expansion region that corresponds 
-- to the target x-location to start the truncated expansion region
function error_in_r(r)
    return (target - exp_region:eval(r, 0.0).x)
end
local secant = require "secant"
trunc_r0, err = secant.solve(error_in_r, 0.8, 0.99, 1.0e-12)
if err then
   print("Secant iteration failed.")
   print("error message is:", err)
end

-- Reset the extent of the expansion_region to the new r0, such that 
-- we are only re-computing a small portion of the end of the nozzle.
trunc_exp_region = SubRangedSurface:new{underlying_surface=exp_region, 
                                        r0=trunc_r0, r1=1.0, s0=0.0, s1=1.0}

-- Define structured grid for now-truncated nozzle region
x_cf_nozzle = RobertsFunction:new{end0=true, end1=true, beta=1.00}
y_cf_nozzle = RobertsFunction:new{end0=false, end1=true, beta=1.002}
exp_grid = StructuredGrid:new{psurface=trunc_exp_region, niv=51, njv=81,
			      cfList={west=y_cf_nozzle, east=y_cf_nozzle,
				      south=x_cf_nozzle, north=x_cf_nozzle}}

-- Set up fluid block for now-truncated nozzle region
exp_fba = FBArray:new{grid=exp_grid, initialState=nozzle_soln, label="truncated_expansion",
		      bcList={west=InFlowBC_StaticProfile:new{fileName="extract-slice.data", 
		                                              match="xyA-to-xyA"},
                              north=WallBC_NoSlip_FixedT:new{Twall=300.0, wall_function=true},
                              east=OutFlowBC_Simple:new{}}, 
                      nib=1, njb=1}

-- Now add on test section geometry
ts_core_region = CoonsPatch:new{
    p00=Vector3:new{x=x_exit, y=0.0},
    p10=Vector3:new{x=x_exit+ts_sim_length, y=0.0},
    p11=Vector3:new{x=x_exit+ts_sim_length, y=y_exit},
    p01=Vector3:new{x=x_exit, y=y_exit}}
ts_wall_region = CoonsPatch:new{
    p00=Vector3:new{x=x_exit, y=y_exit},
    p10=Vector3:new{x=x_exit+ts_sim_length, y=y_exit},
    p11=Vector3:new{x=x_exit+ts_sim_length, y=ts_radius},
    p01=Vector3:new{x=x_exit, y=ts_radius}}

-- Define structured grid for test section blocks
ts_niv = 646
x_cf_ts = RobertsFunction:new{end0=true, end1=true, beta=1.0}
ts_core_grid = StructuredGrid:new{psurface=ts_core_region, niv=ts_niv, njv=81,
                                  cfList={west=y_cf_nozzle,
                                          east=RobertsFunction:new{end0=false, end1=true, beta=1.15},
                                          south=x_cf_ts, north=x_cf_ts}}
ts_wall_grid = StructuredGrid:new{psurface=ts_wall_region, niv=ts_niv, njv=51,
                                  cfList={west=RobertsFunction:new{end0=true, end1=false, beta=1.002},
                                          east=RobertsFunction:new{end0=true, end1=true, beta=1.09},
                                          south=x_cf_ts, north=x_cf_ts}}

-- Define full domain into blocks
ts_core_fba = FBArray:new{grid=ts_core_grid, initialState=initial, label="ts-core",
                          bcList={east=OutFlowBC_Simple:new{}},
                          nib=7, njb=1}
ts_wall_fba = FBArray:new{grid=ts_wall_grid, initialState=initial, label="ts-wall",
                          bcList={west=WallBC_WithSlip:new{label="ts-west-wall"},
                                  north=WallBC_WithSlip:new{label="ts-north-wall"},
                                  east=OutFlowBC_Simple:new{}},
                          nib=7, njb=2}

identifyBlockConnections()
