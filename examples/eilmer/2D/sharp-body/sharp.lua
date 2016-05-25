-- sharp.lua
config.title = "Mach 3 flow over a sharp 2D body"
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=2000.0, vely=0.0}

-- Geometry of flow domain.
function y(x)
   -- (x,y)-space path for x>=0
   if x <= 3.291 then
      return -0.008333 + 0.609425*x - 0.092593*x*x
   else
      return 1.0
   end
end

function xypath(t)
   -- Parametric path with 0<=t<=1.
   local x = 10.0 * t
   local yval = y(x)
   if yval < 0.0 then
      yval = 0.0
   end
   return {x=x, y=yval}
end

a = Vector3:new{x=-1.0, y=0.0}; b = Vector3:new{ x=0.0, y=0.0}
c = Vector3:new{x=10.0, y=1.0}; d = Vector3:new{x=10.0, y=7.0}
e = Vector3:new{ x=0.0, y=7.0}; f = Vector3:new{x=-1.0, y=7.0}
-- lower boundary including body surface
ab = Line:new{p0=a, p1=b}; bc = LuaFnPath:new{luaFnName="xypath"} 
-- upper boundary
fe = Line:new{p0=f, p1=e}; ed = Line:new{p0=e, p1=d}
-- vertical lines
af = Line:new{p0=a, p1=f}; be = Line:new{p0=b, p1=e}
cd = Line:new{p0=c, p1=d} 
-- Mesh the patches, with particular discretisation.
ny = 60
clustery = RobertsFunction:new{end0=true, end1=false, beta=1.3}
clusterx = RobertsFunction:new{end0=true, end1=false, beta=1.2}
grid0 = StructuredGrid:new{psurface=makePatch{north=fe, east=be, south=ab, west=af},
			   cfList={east=clustery, west=clustery},
			   niv=17, njv=ny+1}
grid1 = StructuredGrid:new{psurface=makePatch{north=ed, east=cd, south=bc, west=be},
			   cfList={north=clusterx,south=clusterx,west=clustery},
			   niv=81, njv=ny+1}
-- Define the flow-solution blocks.
blk0 = SBlock:new{grid=grid0, fillCondition=inflow}
blk1 = SBlock:new{grid=grid1, fillCondition=initial}
-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_Supersonic:new{flowCondition=inflow}
blk1.bcList[east] = OutFlowBC_Simple:new{}

config.max_time = 15.0e-3  -- seconds
config.max_step = 2500
config.dt_init = 1.0e-6
