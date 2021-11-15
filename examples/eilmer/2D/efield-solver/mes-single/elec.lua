-- elec.lua
-- Test file for electric field solver
-- NNG
-- 2021-01-20

config.title = "Testing Electric Field Solver"
print(config.title)
config.dimensions = 2
config.solve_electric_field = true
config.field_conductivity_model = "test"
config.new_flow_format = true -- We use the new flow format to write out the flow data
config.flow_format = "eilmer4text"

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=101.35e3, T=300.0, velx=0.0}

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=1.0, y=0.0}
c = Vector3:new{x=0.0, y=1.0}
d = Vector3:new{x=1.0, y=1.0}
ab = Line:new{p0=a, p1=b}
bd = Line:new{p0=b, p1=d}
ac = Line:new{p0=a, p1=c}
cd = Line:new{p0=c, p1=d}
quad = makePatch{north=cd, east=bd, south=ab, west=ac}

-- Mesh the patches, with particular discretisation.
nx = 16; ny = 16;
grid = StructuredGrid:new{psurface=quad, niv=nx+1, njv=ny+1}
blk = FluidBlock:new{grid=grid,
                     initialState=initial,
                     bcList={north = WallBC_WithSlip:new{field_bc=FixedGradient_Test:new{}},
                             east  = WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
                             south = WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
                             west  = WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
                     }
}

-- No boundary conditions needed, just yet
identifyBlockConnections()

-- global data.
config.max_step = 1
