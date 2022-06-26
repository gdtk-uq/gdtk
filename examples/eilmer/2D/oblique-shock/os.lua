-- os.lua
-- Oblique shock that is defined by flow conditions around the box.
-- PJ 2022-06-26

config.title = "Oblique shock in a box"
print(config.title)
config.dimensions = 2
--
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
-- Compute the flowstates either side of the oblique shock.
M1 = 2.0
gs1 = GasState:new{gm}
gs1.p = 100.0e5
gs1.T = 300.0
gm:updateThermoFromPT(gs1)
gm:updateSoundSpeed(gs1)
V1 = M1 * gs1.a
fs1 = FlowState:new{p=gs1.p, T=gs1.T, velx=V1}
print("fs1=", fs1)
--
beta = math.rad(45.0)
gs2, theta, V2 = gasflow.theta_oblique(gs1, V1, beta)
print("gs2=", fs2)
print("theta=", theta, "radians")
print("V2=", V2)
fs2 = FlowState:new{p=gs2.p, T=gs2.T,
                    velx=V2*math.cos(theta), vely=V2*math.sin(theta)}
print("fs2=", fs2)
-- Set up a box in the (x,y)-plane and grid it.
quad0 = CoonsPatch:new{p00={x=0,y=0}, p10={x=1,y=0},
                       p11={x=1,y=1}, p01={x=0,y=1}}
grid0 = StructuredGrid:new{psurface=quad0, niv=41, njv=41}
-- Set up the flow domain as a single block with boundary conditions.
bc = InOutFlowBC_DualState:new{flowState1=fs1, flowState2=fs2,
                               p={x=0,y=0}, n={x=-1.0,y=1.0}}
blk0 = FluidBlock:new{grid=grid0, initialState=fs1,
                      bcList={east=bc, north=bc, south=bc, west=bc}}
--
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5
