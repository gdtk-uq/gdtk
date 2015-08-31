-- rmi.lua
config.title = "Richtmyer-Meshkov Instability"
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
state1 = FlowState:new{p=100.0e3, T=300.0, velx=0.0, vely=0.0} -- pre-shock
state2 = FlowState:new{p=497.8e3, T=531.1, velx=469.9, vely=0.0} -- post-shock
state3 = FlowState:new{p=100.0e3, T=75.0, velx=0.0, vely=0.0} -- dense gas


-- Geometry of flow domain.
L = 15.0 -- length of downstream part of duct
H = 1.0  -- half-height of duct
delta = 0.1 -- magnitude of interface perturbation
-- We're going to simulate both sides of the symmetry plane.
a0 = Vector3:new{-L/5, -H}; a1 = Vector3:new{-L/5,H}
b0 = Vector3:new{delta, -H}; b1 = Vector3:new{delta, H}
c0 = Vector3:new{L, -H}; c1 = Vector3:new{L, H}
-- Shape of the disturbed interface.
function xypath(t)
   -- Parametric path with 0<=t<=1.
   local y = -1.0 * (1.0-t) + 1.0*t
   local x = delta*math.cos(2*math.pi*t)
   return {x, y, 0.0}
end

-- On the left of the interface
region1 = AOPatch:new{south=Line:new{a0,b0}, north=Line:new{a1,b1},
		      west=Line:new{a0,a1}, east=LuaFnPath:new{"xypath"}}
-- On the right of the interface
region2 = AOPatch:new{south=Line:new{b0,c0}, north=Line:new{b1,c1},
		      west=LuaFnPath:new{"xypath"}, east=Line:new{c0,c1}}
		      
nx = 80
grid1 = StructuredGrid:new{psurface=region1, niv=2*nx+1, njv=2*nx+1}
grid2 = StructuredGrid:new{psurface=region2, niv=10*nx+1, njv=2*nx+1}
-- Define the boundary conditions that we care about 
-- and the flow-solution blocks.
bcList1 = {north=nil, east=nil, south=nil,
	   west=SupInBC:new{flowCondition=state2, label="inflow-boundary"}}
blk1 = SBlockArray{grid=grid1, nib=2, njb=2,
		   fillCondition=state1, bcList=bcList1, label="BLOCK-1"}
bcList2 = {north=nil, east=ExtrapolateOutBC:new{label="outflow-boundary"},
	   south=nil, west=nil}
blk2 = SBlockArray{grid=grid2, nib=10, njb=2,
		   fillCondition=state3, bcList=bcList2, label="BLOCK-2"}
identifyBlockConnections()

config.max_time = 40.0e-3  -- seconds
config.max_step = 10000
config.dt_init = 5.0e-6
