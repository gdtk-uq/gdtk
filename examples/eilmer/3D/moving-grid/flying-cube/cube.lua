-- Structured, spherical mesh around a cube.
-- 2021-09-26: PJ, adapted from the sphere example.
-- 2021-10-01: FZ, converted to flying cube with moving mesh

-- Setting up the gas parameters and values
config.dimensions = 3
config.axisymmetric = false
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=500.0, T=300.0}
inflow = FlowState:new{p=760.0, T=71.0, velx=1005.0}
config.viscous = true

-- Set up the moving grid components
config.gasdynamic_update_scheme = "moving_grid_1_stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "move_grid.lua"

-- Domain is constructed by 6 blocks attached to the outside of the cube.
-- The blocks are labeled in accordance with the face they are attached to.
a_cube = 0.0127
R_mesh = 0.075
vols = {}
faces = {"east", "west", "north", "south", "top", "bottom"}
for _,f in ipairs(faces) do
   cube_face = CubePatch:new{a=a_cube, centre={0,0,0}, face_name=f}
   outer_face = SpherePatch:new{radius=R_mesh, centre={0,0,0}, face_name=f}
   if f == "east" or f == "south" or f == "top" then
      vols[f] = TwoSurfaceVolume:new{face0=cube_face, face1=outer_face, ruled_direction="k"}
   else
      -- west, north, bottom
      vols[f] = TwoSurfaceVolume:new{face0=outer_face, face1=cube_face, ruled_direction="k"}
   end
end

-- Define the 6 grids.
nFactor = 2
N_edge = 20*nFactor  -- cells along edge of cube faces that are mapped onto the sphere
N_normal = 20*nFactor  -- cells in sphere normal direction
grids = {}
for _,f in ipairs(faces) do
   if f == "east" or f == "south" or f == "top" then
      cf = RobertsFunction:new{end0=true, end1=false, beta=1.05}
   else
      cf = RobertsFunction:new{end0=false, end1=true, beta=1.05}
   end
   cfList = {edge04=cf, edge15=cf, edge26=cf, edge37=cf}
   grids[f] = StructuredGrid:new{pvolume=vols[f], cfList=cfList,
                                 niv=N_edge, njv=N_edge, nkv=N_normal}
end

-- and, finally, the FluidBlocks
blks = {}
nib = 3; njb = 3; nkb = 3
wall_bc = WallBC_NoSlip_FixedT:new{Twall=300.0, group='walls'}
inflow_bc = InOutFlowBC_Ambient:new{flowState=inflow}
for _,f in ipairs(faces) do
   if f == "east" or f == "south" or f == "top" then
      bcList = {top=inflow_bc, bottom=wall_bc}
   else
      bcList = {bottom=inflow_bc, top=wall_bc}
   end
   blks[f] = FBArray:new{grid=grids[f], initialState=initial,
                         bcList=bcList, nib=nib, njb=njb, nkb=nkb}
end
identifyBlockConnections()

-- Set up the run-time loads for computing the cube movement
run_time_loads={
   {group="walls", moment_centre=Vector3:new{x= 0.0, y= 0.0, z=0.0}}
       }

-- Setting up the 'standard' simulation parameters
config.max_time = 50.0e-3  -- seconds
config.max_step = 30000000
config.cfl_value = 0.5
config.stringent_cfl = true
config.dt_plot = 0.5e-3

-- These calculations are pretty tough, so allow some bad cells
config.adjust_invalid_cell_data = true
config.max_invalid_cells = 15

-- We need to do some extra work in a controlling file
config.udf_supervisor_file='udf-process.lua'

-- Run time loads are used to calculate the aerodynamic loading
config.compute_run_time_loads = true
config.run_time_loads_count = 1

-- We use a 'user_pad' to store all our information 
config.user_pad_length = 8
-- 1=angle, 2=angular_velocity, x, xdot, y, ydot, z, zdot
ang_vel = 7111 * 2 * math.pi / 60 -- 7500rpm
user_pad_data = {0.0, ang_vel, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

-- Distribute the blocks for running MPI jobs
mpiTasks = mpiDistributeBlocks{ntasks=100, dist="load-balance"}