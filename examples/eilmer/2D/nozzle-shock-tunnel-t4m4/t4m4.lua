-- t4m4.lua
-- A demonstration of doing a NENZFR-like calculation in Eilmer4.
-- with model of the original T4 Mach4 nozzle and operating conditions
-- from Wilson's 2015 example t4-m4-nozzle.py
-- PAJ, 2017-12-04

config.title = "T4 Mach4 nozzle with air in chemical equilibrium."
print(config.title)

print("First, compute the flow conditions at the nozzle throat.")
-- Using the CEA-backed gas model.
print("shock-tube fill conditions")
gm = GasModel:new{'cea-air13species-gas-model.lua'}
state1 = GasState:new{gm}
state1.p = 131.0e3; state1.T = 300.0
gm:updateThermoFromPT(state1)
print("state1:"); printValues(state1)
--
print("normal shock, given shock speed")
Vs = 1644.0
state2, V2, Vg = gasflow.normal_shock(state1, Vs)
print("V2=", V2, "Vg=", Vg)
print("state2:"); printValues(state2)
--
print("reflected shock")
state5, Vr = gasflow.reflected_shock(state2, Vg)
print("Vr=", Vr)
print("state5:"); printValues(state5)
--
print("Expand from stagnation (with ratio of pressure to match observation)")
state5s, V5s = gasflow.expand_from_stagnation(state5, 8.32/23.61)
print("V5s=", V5s, " Mach=", V5s/state5s.a)
print("state5s:"); printValues(state5s)
print("(h5s-h1)=", gm:enthalpy(state5s) - gm:enthalpy(state1)) 
--
print("Expand to throat condition (Mach 1.0001)")
state6, V6 = gasflow.expand_to_mach(state5s, 1.0001)
print("V6=", V6, " Mach=", V6/state6.a)
print("state6:"); printValues(state6)
--
print("Quasi-one-dimensional expansion to estimated test-flow condition.")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 27.0)
print("V7=", V7, " Mach=", V7/state7.a)
print("state7:"); printValues(state7)

print("Second, set up the simulation grid and the initial flow state.")
config.dimensions = 2
config.axisymmetric = true

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

print "Building grid."
-- Coordinates are from the nenzfr file Bezier-control-pts-t4-m4.data.
-- There is an implied straight (conical) section from x_start to the
-- first point in the Bezier contour that defines the flow-straightening
-- section of the nozzle wall.
-- Dimensions in metres.
L_thr = 0.0125
x_start = 0.0 -- start of supersonic expansion
y_throat = 1.250e-2
bez_points = {
   Vector3:new{x=9.843000e-02, y=3.343000e-02},
   Vector3:new{x=1.811060e-01, y=5.080216e-02},
   Vector3:new{x=2.637820e-01, y=6.172235e-02},
   Vector3:new{x=3.464580e-01, y=6.373297e-02},
   Vector3:new{x=4.291340e-01, y=6.660485e-02},
   Vector3:new{x=5.118100e-01, y=6.611000e-02}}

-- The throat-region is the constant-area section up to the
-- start of the xonical expansion.
throat_region = CoonsPatch:new{
   p00=Vector3:new{x=-L_thr, y=0.0},
   p10=Vector3:new{x=x_start, y=0.0},
   p11=Vector3:new{x=x_start, y=y_throat},
   p01=Vector3:new{x=-L_thr, y=y_throat}}
-- Supersonic expansion
a0 = Vector3:new{x=x_start, y=0.0}
a1 = Vector3:new{x=x_start, y=y_throat}
b0 = Vector3:new{x=bez_points[1].x, y=0.0}
c0 = Vector3:new{x=bez_points[#bez_points].x, y=0.0}
exp_region = ChannelPatch:new{
   north=Polyline:new{segments={Line:new{p0=a1, p1=bez_points[1]},
				Bezier:new{points=bez_points}}},
   south=Line:new{p0=a0, p1=c0},
   ruled=true, pure2D=true}
-- Define structured grids for both regions.
x_cf_throat = RobertsFunction:new{end0=true, end1=true, beta=2.0}
x_cf = RobertsFunction:new{end0=true, end1=false, beta=1.1}
y_cf = RobertsFunction:new{end0=false, end1=true, beta=1.02}
throat_grid = StructuredGrid:new{psurface=throat_region, niv=39, njv=41,
				 cfList={west=y_cf, east=y_cf,
					 south=x_cf_throat, north=x_cf_throat}}
exp_grid = StructuredGrid:new{psurface=exp_region, niv=601, njv=41,
			      cfList={west=y_cf, east=y_cf,
				      south=x_cf, north=x_cf}}
-- Divide the full domain into many blocks, mainly in the axial direction.
throat_blk = FluidBlockArray{grid=throat_grid, initialState=inflow, label="throat",
			     bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
				     north=WallBC_NoSlip_FixedT:new{Twall=300.0}}, 
			     nib=1, njb=2}
exp_blk = FluidBlockArray{grid=exp_grid, initialState=inflow, label="expansion",
			  bcList={north=WallBC_NoSlip_FixedT:new{Twall=300.0},
				  east=OutFlowBC_SimpleExtrapolate:new{}}, 
			  nib=30, njb=2}
identifyBlockConnections()

-- Set a few more config options
config.viscous = true
config.flux_calculator = "adaptive"
config.max_time = 0.001 -- allow enough to reach steady state
config.max_step = 800000
config.dt_init = 1.0e-9
-- To accelerate the calculation...
config.block_marching = true
config.nib = 31
config.njb = 2
config.propagate_inflow_data = true
