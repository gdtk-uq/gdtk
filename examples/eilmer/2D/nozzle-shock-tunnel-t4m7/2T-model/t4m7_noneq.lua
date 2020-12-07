-- t4m7_noneq.lua
-- Tamara Sopek, 4-Dec-2020

config.title = "T4 Mach7 nozzle with air in chemical & thermal non-equilibrium."
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
--
config.reacting = true -- turn on reactions
config.reactions_file = 'air-5sp-6r-2T.lua'
config.energy_exchange_file = 'air-energy-exchange.lua'

-- Get flow conditions at the nozzle throat by performing ESTCj-like calculations,
-- using the CEA-backed gas model.
print("First, compute the flow conditions at the nozzle throat.")
print("shock-tube fill conditions")
gm = GasModel:new{'cea-air5species-gas-model.lua'}
state1 = GasState:new{gm}
state1.p = 200.0e3; state1.T = 300.0
gm:updateThermoFromPT(state1); gm:updateTransCoeffs(state1)
print("state1:"); printValues(state1)
--
print("normal shock, given shock speed")
Vs = 1678.62
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
state5s, V5s = gasflow.expand_from_stagnation(state5, 19.3253e6/state5.p)  
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
print("Quasi-one-dimensional expansion to estimated test-flow condition.")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 129.72)
gm:updateThermoFromPT(state7); gm:updateTransCoeffs(state7)
print("V7=", V7, " Mach=", V7/state7.a)
print("state7:"); printValues(state7)

print("Second, set up the simulation grid and the initial flow state.")
-- Estimate turbulence quantities for the inflow. 
turb_intensity = 0.05
turb_to_laminar_visc_ratio = 100.0
throat_tke = 1.5 * (throat_V * turb_intensity)^2
throat_mu_t = turb_to_laminar_visc_ratio * throat_mu
throat_omega = throat_rho * throat_tke / throat_mu_t
print("Inflow turbulence: tke=", throat_tke, "omega=", throat_omega)

-- Set up flow state.
nsp, nmodes, gm = setGasModel('air-5sp-2T.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)

-- Work out species mass fractions at throat
print("Species mass fractions at throat, according to gas model")
massf_throat = {}; massf_sum = 0.0
for k, v in pairs(state6.ceaSavedData.massf) do
   print(string.format("massf[%s] = %g", k, v))
   massf_throat[k] = v
   massf_sum = massf_sum + v
end
-- Sanity check in case the CEA calculated terms don't add to 1.0
-- Usually out by approx 1.0e-6. Just add this to the N2 term as its
-- influence will certainly be neglible, according to Will Landsberg.
if massf_sum ~= 1.0 then
  print("Mass fraction adjustment to sum mass to 1.0")
  print("Original mass sum = ", massf_sum)
  massf_throat['N2'] = massf_throat['N2'] + (1.0 - massf_sum)
end

-- The inflow for the simulation is at the nozzle throat.
inflow = FlowState:new{p=state6.p, T=state6.T, T_modes={state6.T}, velx=V6, tke=throat_tke, 
                       omega=throat_omega, massf=massf_throat}
print(inflow)

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

-- The throat-region is the constant-area section up to the
-- start of the contoured expansion.
throat_region = CoonsPatch:new{
   p00=Vector3:new{x=-L_thr, y=0.0},
   p10=Vector3:new{x=x_start, y=0.0},
   p11=Vector3:new{x=x_start, y=y_throat},
   p01=Vector3:new{x=-L_thr, y=y_throat}}
-- Supersonic expansion is fully defined by its north edge
nozzle_profile = ArcLengthParameterizedPath:new{underlying_path=Bezier:new{points=bez_points}}
exp_region = NozzleExpansionPatch:new{north=nozzle_profile}
-- Define structured grids for both regions.
print "Building grid."
x_cf_throat = RobertsFunction:new{end0=true, end1=true, beta=1.53} 
x_cf = RobertsFunction:new{end0=true, end1=false, beta=1.1} --1.1
y_cf = RobertsFunction:new{end0=false, end1=true, beta=1.01}
throat_grid = StructuredGrid:new{psurface=throat_region, niv=31, njv=81,
				 cfList={west=y_cf, east=y_cf,
					 south=x_cf_throat, north=x_cf_throat}}
exp_grid = StructuredGrid:new{psurface=exp_region, niv=601, njv=81,
			      cfList={west=y_cf, east=y_cf,
				      south=x_cf, north=x_cf}}
-- Divide the full domain into many blocks, mainly in the axial direction.
-- Wall functions used for transition in throat or at x_tr=0 to allow faster simulations from using 
-- weaker cell clustering. Since facility nozzles should have non-separated BL which wall functions
-- model well, it is appropriate to used them here. 
throat_fba = FBArray:new{grid=throat_grid, initialState=inflow, label="throat",
                         bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                                 north=WallBC_NoSlip_FixedT:new{Twall=300.0, wall_function=true}}, 
                         nib=1, njb=4}
exp_fba = FBArray:new{grid=exp_grid, initialState=inflow, label="expansion",
                      bcList={north=WallBC_NoSlip_FixedT:new{Twall=300.0, wall_function=true},
                              east=OutFlowBC_Simple:new{}}, 
                      nib=30, njb=4}
identifyBlockConnections()

-- Define approximately where boundary layer transitions to turbulence.
-- We set the location of the start of the TurbulentZone at the start of
-- the nozzle.
x_tr = 0.05; x1 = bez_points[#bez_points].x
y0 = 0.0; y1 = bez_points[#bez_points].y
turbZone_botLeftCorner = Vector3:new{x=x_tr, y=y0}
turbZone_topRightCorner = Vector3:new{x=x1, y=y1}
TurbulentZone:new{p0=turbZone_botLeftCorner, p1=turbZone_topRightCorner}
--
