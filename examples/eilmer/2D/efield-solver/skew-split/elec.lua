-- elec.lua
-- Test file for electric field solver
-- NNG
-- 2021-01-20

config.title = "Testing Electric Field Solver"
print(config.title)
config.dimensions = 2
config.solve_electric_field = true
config.electric_field_count = 1
config.field_conductivity_model = "constant"
config.new_flow_format = true -- We use the new flow format to write out the flow data
config.flow_format = "eilmer4binary"

-- The gas model is defined via a gas-model file.
nsp, nmodes, gmodel = setGasModel('ideal-air.lua')

gs = GasState:new{gmodel}
gs.p = 101.325e3 -- Pa
T = 300.0
gs.T = T

initial = FlowState:new{p=gs.p, T=gs.T, velx=0.0}
print("Initial: ")
print(initial)

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
D = 0.2

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=1.0, y=0.0}

c = Vector3:new{x=0.0+D, y=1.0}
d = Vector3:new{x=1.0+D, y=1.0}

ab = Line:new{p0=a, p1=b}
bd = Line:new{p0=b, p1=d}
ac = Line:new{p0=a, p1=c}
cd = Line:new{p0=c, p1=d}

quad = makePatch{north=cd, east=bd, south=ab, west=ac}

-- Mesh the patches, with particular discretisation.
nx = 32; ny = 32;
grid = StructuredGrid:new{psurface=quad, niv=nx+1, njv=ny+1}
blks = FBArray:new{nib=1, njb=2, grid=grid,
                     initialState=initial,
                     bcList={
                             north = WallBC_WithSlip:new{field_bc = FixedField:new{value=2.0}},
                             east  = WallBC_WithSlip:new{field_bc = ZeroNormalGradient:new{}},
                             west  = WallBC_WithSlip:new{field_bc = ZeroNormalGradient:new{}},
                             south = WallBC_WithSlip:new{field_bc = FixedField:new{value=1.0}},
                     }
}

-- global data.
config.max_step = 1
