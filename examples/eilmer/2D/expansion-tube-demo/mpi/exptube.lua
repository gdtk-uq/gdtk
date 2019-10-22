-- exptube.lua
-- An expansion tube demonstration, MPI flavour.
-- Although the set-up is more complicated than the shared-memory variant,
-- this MPI example should scale to very large simulations.
--
-- PJ 2019-06-14, 2019-10-17
--
config.title = "MPI Expansion tube demo with UDF-bc diaphragm."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
config.max_time = 1.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 10.0e-6
config.dt_history = 1.0e-6
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
   quads[i] = CoonsPatch:new{p00={x=Ls[i], y=0.0}, p10={x=Ls[i+1], y=0.0},
			     p11={x=Ls[i+1], y=R}, p01={x=Ls[i], y=R}}
end
--
-- Mesh the patches with a uniform discretisation based on a cell size.
grids = {}
local dx = 0.010
for i = 0, 2 do
   local nx = math.floor((Ls[i+1]-Ls[i])/dx); local ny = 2
   grids[i] = StructuredGrid:new{psurface=quads[i], niv=nx+1, njv=ny+1}
end
--
-- Construct several FluidBlocks per tube.
nibList = {[0]=2, 3, 3}
fba = {}
for i = 0, 2 do
   fba[i] = FBArray:new{grid=grids[i], nib=nibList[i], initialState=fill_gases[i]}
end
identifyBlockConnections()
--
-- The diaphragm boundary condition is between final block in shock tube
-- and first block in acceleration tube.
upstream_blk = fba[1].blockArray[3][1]
downstream_blk = fba[2].blockArray[1][1]
print("For diaphragm, block ids are:")
print("   upstream:", upstream_blk.id)
print("   downstream:", downstream_blk.id)
upstream_blk.bcList[east] = ExchangeBC_FullFacePlusUDF:new{
   otherBlock=downstream_blk.id, otherFace=west, fileName="diaphragm.lua"
}
downstream_blk.bcList[west] = ExchangeBC_FullFacePlusUDF:new{
   otherBlock=upstream_blk.id, otherFace=east, fileName="diaphragm.lua"
}
-- Final block in the acceleration tube is open at the downstream end.
fba[2].blockArray[3][1].bcList[east] = OutFlowBC_SimpleExtrapolate:new{}
--
-- A place to record the state of the diaphragm.
-- 0 == diaphragm closed
-- 1 == diaphragm open (ruptured)
config.user_pad_length = 1
user_pad_data = {0}
-- We want the shock-tube right-most block that sets the rupture state
-- to be on the MPI master task 0.  Its user_pad_data is broadcast.
mpiDistributeBlocks{ntasks=8, dist="load-balance",
                    preassign={[upstream_blk.id]=0}}
-- The function that sets the diaphragm state is also a user-defined function.
config.udf_supervisor_file='supervisor.lua'
