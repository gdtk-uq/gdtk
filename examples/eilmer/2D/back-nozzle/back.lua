-- back.py
-- Conical nozzle from Back, Massier and Gier (1965)
-- Peter J. 2015-10-21 adpated from the Python version.
config.title = "Flow through a conical nozzle."
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
-- The stagnation gas represents a reservoir condition.
stagnation_gas = FlowState:new{p=500.0e3, T=300.0}
low_pressure_gas = FlowState:new{p=30.0, T=300.0}

-- Define geometry.
-- The original paper specifies sizes in inches, Eilmer3 works in metres.
inch = 0.0254 -- metres
L_subsonic = 3.0 * inch
L_nozzle = 3.0 * inch
R_tube = 1.5955 * inch
R_throat = 0.775 * inch
R_curve = 1.55 * inch -- radius of curvature of throat profile
theta = 15.0 * math.pi / 180.0 -- radians

-- Compute the centres of curvature for the contraction profile.
height = R_throat + R_curve
hypot = R_tube + R_curve
base = math.sqrt(hypot*hypot - height*height)
centre_A = Vector3:new{0.0, height}
centre_B = Vector3:new{-base, 0.0}
fraction = R_tube/hypot
intersect_point = centre_B + Vector3:new{fraction*base, fraction*height}

-- Assemble nodes from coordinates.
z0 = Vector3:new{-L_subsonic, 0.0}
p0 = Vector3:new{-L_subsonic, R_tube}
z1 = Vector3:new{centre_B:x(), centre_B:y()} -- initialize from a previously defined Node
p1 = centre_B + Vector3:new{0.0, R_tube}
p2 = Vector3:new{intersect_point:x(), intersect_point:y()}
z2 = Vector3:new{p2:x(), 0.0}  -- on the axis, below p2
z3 = Vector3:new{0.0, 0.0}
p3 = Vector3:new{0.0, R_throat}
-- Compute the details of the conical nozzle
p4 = Vector3:new{R_curve*math.sin(theta), height - R_curve*math.cos(theta)}
z4 = Vector3:new{p4:x(), 0.0}
L_cone = L_nozzle - p4:x()
p5 = p4 + Vector3:new{L_cone, L_cone*math.tan(theta)}
z5 = Vector3:new{p5:x(), 0.0}

north0 = Polyline:new{Line:new{p0,p1},Arc:new{p1,p2,centre_B},Arc:new{p2,p3,centre_A}}
east0west1 = Line:new{z3, p3}
south0 = Line:new{z0, z3}
west0 = Line:new{z0, p0}
north1 = Polyline:new{Arc:new{p3,p4,centre_A}, Line:new{p4,p5}}
east1 = Line:new{z5, p5}
south1 = Line:new{z3, z5}

-- Define the blocks, boundary conditions and set the discretisation.
nx0 = 50; nx1 = 60; ny = 30
grid0 = StructuredGrid:new{psurface=makePatch{north0, east0west1, south0, west0},
			   niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=makePatch{north1, east1, south1, east0west1},
			   niv=nx1+1, njv=ny+1}
subsonic_region = SBlock:new{grid=grid0, fillCondition=stagnation_gas,
			     label="subsonic-region"}
supersonic_region = SBlock:new{grid=grid1, fillCondition=low_pressure_gas,
			       label="supersonic-region", 
			       hcellList={{1,1},{nx1-1,1}}}
-- History locations near throat and exit
identifyBlockConnections()
subsonic_region.bcList[west] = InFlowBC_FromStagnation:new{stagCondition=stagnation_gas,
							   label="inflow-boundary"}
supersonic_region.bcList[east] = OutFlowBC_Simple:new{label="outflow-boundary"}

-- Do a little more setting of global data.
config.axisymmetric = 1
config.flux_calculator = "adaptive"
config.max_time = 4.0e-3  -- seconds
config.max_step = 50000
config.dt_init = 1.0e-7
config.dt_plot = 0.2e-3
config.dt_history = 10.0e-6

