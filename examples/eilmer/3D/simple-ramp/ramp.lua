-- ramp.lua
-- A simple 3D simulation of flow over a ramp with 10-degree deflection.
-- PJ & RG
-- 2015-02-24 -- adapted from the Lua version of cone20 and 
--               the Python version simple_ramp
config.title = "Mach 1.5 flow over a 10-degree ramp."
print(config.title)
config.dimensions = 3
config.viscous = false
config.max_time = 5.0e-3
config.max_step = 1000
config.gasdynamic_update = "euler"
config.interpolation_order = 2
config.dt_plot = 1.0e-3
config.dt_history = 1.0e-5

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0, vely=0.0}

function simpleBoxCorners(values)
   -- Creates a corner coordinate list for a simple box.
   local xPos = values.xPos or 0.0
   local yPos = values.yPos or 0.0
   local zPos = values.zPos or 0.0
   local xSize = values.xSize or 1.0
   local ySize = values.ySize or 1.0
   local zSize = values.zSize or 1.0
   p0 = Vector3:new{xPos,       yPos,       zPos}
   p1 = Vector3:new{xPos+xSize, yPos,       zPos}
   p2 = Vector3:new{xPos+xSize, yPos+ySize, zPos}
   p3 = Vector3:new{xPos,       yPos+ySize, zPos}
   p4 = Vector3:new{xPos,       yPos,       zPos+zSize}
   p5 = Vector3:new{xPos+xSize, yPos,       zPos+zSize}
   p6 = Vector3:new{xPos+xSize, yPos+ySize, zPos+zSize}
   p7 = Vector3:new{xPos,       yPos+ySize, zPos+zSize}
   return {p0, p1, p2, p3, p4, p5, p6, p7}
end

-- First block is the region in front of the ramp. 10x40(x4)
-- cluster down, toward the wedge surface
cluster_k = RobertsFunction:new{end0=true, end1=false, beta=1.2}
nocf = LinearFunction:new{}
cflist = {nocf, nocf, nocf, nocf, nocf, nocf, nocf, nocf,
	  cluster_k, cluster_k, cluster_k, cluster_k}
vol0 = TFIVolume:new{vertices=simpleBoxCorners{xSize=0.2,ySize=0.1}}
grid0 = StructuredGrid:new{pvolume=vol0, niv=11, njv=5, nkv=41, cfList=cflist}
-- For the grid over the ramp, start with a regular box... 30x40(x4)
blk1Corners = simpleBoxCorners{xPos=0.2,xSize=0.8,ySize=0.1}
-- Now, raise the end of the ramp.
-- Remember that Lua indexing starts at 1.
blk1Corners[2]:z(0.8 * math.tan(math.pi * 10.0/180.0))
blk1Corners[3]:z(blk1Corners[2]:z())
grid1 = StructuredGrid:new{pvolume=TFIVolume:new{vertices=blk1Corners}, 
			   niv=31, njv=5, nkv=41, cfList=cflist}

blk0 = SBlock:new{label="first-block", grid=grid0, fillCondition=initial}
blk0.bcList[west] = SupInBC:new{flowCondition=inflow}
blk1 = SBlock:new{label="second-block", grid=grid1, fillCondition=initial,
		  hcellList={{1,1,2},{20,1,1}}}
blk1.bcList[east] = ExtrapolateOutBC:new{}
identifyBlockConnections()
