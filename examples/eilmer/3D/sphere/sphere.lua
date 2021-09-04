-- Structured mesh around a sphere.
--
-- Author: Ingo Jahn & Peter J.
-- Created: 22/3/2021
-- Modified:
--   2021-04-03 simplified to use the SpherePatch and TwoSurfaceVolume classes.
--   2021-09-04 get the orientation correct for half of the blocks.
--
config.dimensions = 3
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
--
-- Domain is constructed by 6 blocks attached to the outside of the sphere.
-- The blocks are labeled in accordance with the face they are attached to.
R_sphere = 1.0
R_mesh = 4.0
vols = {}
faces = {"east", "west", "north", "south", "top", "bottom"}
for _,f in ipairs(faces) do
   sphere_face = SpherePatch:new{radius=R_sphere, centre={0,0,0}, face_name=f}
   outer_face = SpherePatch:new{radius=R_mesh, centre={0,0,0}, face_name=f}
   if f == "east" or f == "south" or f == "top" then
      vols[f] = TwoSurfaceVolume:new{face0=sphere_face, face1=outer_face, ruled_direction="k"}
   else
      -- west, north, bottom
      vols[f] = TwoSurfaceVolume:new{face0=outer_face, face1=sphere_face, ruled_direction="k"}
   end
end
-- Define the 6 grids.
N_edge = 20  -- cells along edge of cube faces that are mapped onto the sphere
N_normal = 20  -- cells in sphere normal direction
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
nib = 1; njb = 1; nkb = 1
wall_bc = WallBC_WithSlip:new{}
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
--
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
-- config.dt_init = 1.0e-6
-- config.cfl_value = 0.5
