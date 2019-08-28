-- exptube.lua
-- An expansion tube demonstration, MPI flavour.
-- Although the set-up is more complicated than the shared-memory variant,
-- this MPI example should scale to very large simulations.
--
-- PJ 2019-06-14
--
config.title = "MPI Expansion tube demo with UDF-bc diaphragm."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
--
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
-- We will set up three FluidBlocks to represent 
-- [0] the driver tube
-- [1] shock tube
-- [2] acceleration tube
fill_gases = {[0]=FlowState:new{p=1.0e6, T=3000.0},
	      FlowState:new{p=5.0e3, T=300.0},
	      FlowState:new{p=100.0, T=300.0}}
Ls = {[0]=0.0, 0.5, 1.0, 2.5} -- break-points for tubes 
R = 0.010
--
-- Set up three quads in the (x,y)-plane by first defining
-- the corner nodes, then define the patches with those corners.
quads = {}
for i = 0, 2 do
   quads[i] = CoonsPatch:new{p00=Vector3:new{x=Ls[i], y=0.0},
			     p10=Vector3:new{x=Ls[i+1], y=0.0},
			     p11=Vector3:new{x=Ls[i+1], y=R},
			     p01=Vector3:new{x=Ls[i], y=R}}
end
--
-- Mesh the patches with a uniform discretisation based on a cell size.
grids = {}; blks = {}
local dx = 0.010
for i = 0, 2 do
   local nx = math.floor((Ls[i+1]-Ls[i])/dx); local ny = 2
   grids[i] = StructuredGrid:new{psurface=quads[i], niv=nx+1, njv=ny+1}
   blks[i] = FluidBlock:new{grid=grids[i], initialState=fill_gases[i]}
end
identifyBlockConnections()
-- The diaphragm boundary condition is between FluidBlocks [1] and [2]
blks[1].bcList[east] = ExchangeBC_FullFacePlusUDF:new{
   otherBlock=2, otherFace=west,
   fileName="diaphragm.lua"
}
blks[2].bcList[west] = ExchangeBC_FullFacePlusUDF:new{
   otherBlock=1, otherFace=east,
   fileName="diaphragm.lua"
}
-- A place to record the state of the diaphragm.
config.user_pad_length = 1
user_pad_data = {0}
-- We want the shock-tube block [1] that sets the rupture state
-- to be on the MPI master task 0.  Its user_pad_data is broadcast.
mpiDistributeBlocks{ntasks=3, dist="load-balance", preassign={[1]=0}}
-- The function that sets the diaphragm state is also a user-defined function.
config.udf_supervisor_file='supervisor.lua'
--
-- Do a little more setting of global data.
config.max_time = 1.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 10.0e-6
config.dt_history = 1.0e-6
