-- t4m4b.lua : T4 Mach4B nozzle
-- Tamara Sopek, 2020-11-04

config.title = "T4 Mach4B nozzle with air in chemical equilibrium."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
config.viscous = true
config.flux_calculator = "adaptive"
--
config.max_time = 0.001 -- allow enough to reach steady state
config.max_step = 1000000
config.dt_init = 1.0e-11
config.cfl_value = 0.4 -- default value is 0.5
--
config.block_marching = true -- to accelerate calculations, utilise block-marching capability.
config.nib = 31
config.njb = 2
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
print("Quasi-one-dimensional expansion to estimated test-flow condition.")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 27.0)
gm:updateThermoFromPT(state7); gm:updateTransCoeffs(state7)
print("V7=", V7, " Mach=", V7/state7.a)
print("state7:"); printValues(state7)
--
print("Second, set up the simulation grid and the initial flow state.")

-- Set up flow state.
-- For the simulation, we are going to use the tabulated gas behaviour
-- because it will be much faster than the CEA-backed gas model.
nsp, nmodes, gm = setGasModel('cea-lut-air.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
-- The inflow for the simulation is at the nozzle throat.
inflow = FlowState:new{p=state6.p, T=state6.T, velx=V6}
print("Species mass fractions at throat, according to CEA gas model")
for k, v in pairs(state6.ceaSavedData.massf) do
   print(string.format("massf[%s] = %g", k, v))
end

-- Read coordinates from the Bezier-control-pts-t4-m4.data.
local file = io.open("Bezier-control-pts-t4-m4b.data", "r")
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
nozzle_profile = ArcLengthParameterizedPath:new{underlying_path=Bezier:new{points=bez_points}}
exp_region = NozzleExpansionPatch:new{north=nozzle_profile}
-- Define structured grids for both regions.
print "Building grid."
x_cf_throat = RobertsFunction:new{end0=true, end1=true, beta=1.46}
x_cf = RobertsFunction:new{end0=true, end1=false, beta=1.1}
y_cf = RobertsFunction:new{end0=false, end1=true, beta=1.02}
throat_grid = StructuredGrid:new{psurface=throat_region, niv=49, njv=81,
				 cfList={west=y_cf, east=y_cf,
					 south=x_cf_throat, north=x_cf_throat}}
exp_grid = StructuredGrid:new{psurface=exp_region, niv=751, njv=81,
			      cfList={west=y_cf, east=y_cf,
				      south=x_cf, north=x_cf}}
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
