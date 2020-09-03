-- nozzle-template.lua
-- Optimize an axisymmetric bell nozzle for ENGG7601 assignment.
-- Peter J. 2017-09-06 adpated from the Back nozzle simulation,
-- hence the dimensions in inches.
-- PJ, 2018-03-04, turn into a template file
--
config.title = "Flow through a rocket nozzle."
print(config.title)
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')

-- The stagnation gas represents a reservoir condition inside the rocket motor.
-- The low_pressure_gas is an arbitrary fill condition for two of the blocks.
-- It will be swept away.
-- The external_stream will provide an environment for the rocket's exhaust gas.
stagnation_gas = FlowState:new{p=5.0e6, T=3000.0}
low_pressure_gas = FlowState:new{p=30.0, T=300.0}
external_stream = FlowState:new{p=8.0e3, T=300.0, velx=2.0e3}

-- Define geometry of our rocket motor and nozzle.
-- The original paper by Back etal (for a lab experiment on supersonic nozzles)
-- specifies sizes in inches, Eilmer works in metres.
inch = 0.0254 -- metres
L_subsonic = 3.0 * inch
L_nozzle = 6.0 * inch
R_tube = 1.5955 * inch
R_throat = 0.775 * inch
R_curve = 1.55 * inch -- radius of curvature of throat profile
--
-- The following three angles set the shape of the supersonic expansion
-- part of the nozzle.
-- The profile is defined by a circular arc, followed by a Bezier-curve
-- with 5 defining points {b0, b1, b2, b3, b4} whose positions are set
-- by the angles theta_init, alpha, beta, theta_cone.
-- With theta_init=theta_cone defining the nominally-straight conical nozzle.
-- You may vary alpha and beta away from zero, to generate a curve
-- to replace the straight profile of the nominal cone.
-- The values alpha=0 and beta=0 will give you a Bezier curve that
-- happens to be a straight line.
-- Set theta_init > theta_cone to get a rapidly expanding thrust surface.
--
theta_cone = math.rad($theta_cone) -- nominal straight-cone angle
theta_init = math.rad($theta_init) -- starting angle for thrust nozzle
alpha = math.rad($alpha)  -- angle for setting b2 in Bezier curve
beta = math.rad($beta)  -- angle for setting b3 in Bezier curve

-- Compute the centres of curvature for the contraction profile.
height = R_throat + R_curve
hypot = R_tube + R_curve
base = math.sqrt(hypot*hypot - height*height)
centre_A = Vector3:new{x=0.0, y=height}
centre_B = Vector3:new{x=-base, y=0.0}
fraction = R_tube/hypot
intersect_point = centre_B + Vector3:new{x=fraction*base, y=fraction*height}

-- Assemble nodes from coordinates.
z0 = Vector3:new{x=-L_subsonic, y=0.0}
p0 = Vector3:new{x=-L_subsonic, y=R_tube}
z1 = Vector3:new{centre_B} -- initialize from a previously defined Node
p1 = centre_B + Vector3:new{x=0.0, y=R_tube}
p2 = Vector3:new{intersect_point}
z2 = Vector3:new{x=p2.x, y=0.0}  -- on the axis, below p2
z3 = Vector3:new{x=0.0, y=0.0}
p3 = Vector3:new{x=0.0, y=R_throat}
-- Compute the details of the conical nozzle.
-- Circular arc to p4, followed by straight line at angle theta to p5.
p4 = Vector3:new{x=R_curve*math.sin(theta_init),
		 y=height-R_curve*math.cos(theta_init)}
z4 = Vector3:new{x=p4.x, y=0.0}
L_cone = L_nozzle - p4.x
R_exit = p4.y + L_cone*math.tan(theta_cone)
p5 = Vector3:new{x=p4.x+L_cone, y=R_exit}
z5 = Vector3:new{x=p5.x, y=0.0}
-- Final nodes define the Bezier curve.
b0 = p4
b1 = p4 + 0.2*L_cone*Vector3:new{x=1.0, y=math.tan(theta_init)}
b2 = p4 + 0.4*L_cone*Vector3:new{x=1.0, y=math.tan(theta_init+alpha)}
b3 = p5 - 0.3*L_cone*Vector3:new{x=1.0, y=math.tan(theta_cone-beta)}
b4 = p5
-- Some space downstream of the nozzle exit
z6 = Vector3:new{x=z5.x+L_nozzle, y=0.0}
p6 = Vector3:new{x=z6.x, y=R_exit}
q5 = Vector3:new{x=z5.x, y=2*R_exit}
q6 = Vector3:new{x=z6.x, y=q5.y}

north0 = Polyline:new{segments={Line:new{p0=p0,p1=p1},
				Arc:new{p0=p1,p1=p2,centre=centre_B},
				Arc:new{p0=p2,p1=p3,centre=centre_A}}}
east0west1 = Line:new{p0=z3, p1=p3}
south0 = Line:new{p0=z0, p1=z3}
west0 = Line:new{p0=z0, p1=p0}
north1 = Polyline:new{segments={Arc:new{p0=p3,p1=p4,centre=centre_A},
				Bezier:new{points={b0, b1, b2, b3, b4}}}}
east1 = Line:new{p0=z5, p1=p5}
south1 = Line:new{p0=z3, p1=z5}
-- The subsonic and supersonic parts of the nozzle have complicated edges.
patch0 = CoonsPatch:new{north=north0, east=east0west1, south=south0, west=west0}
patch1 = CoonsPatch:new{north=north1, east=east1, south=south1, west=east0west1}
-- The downstream region is just two rectangular boxes.
patch2 = CoonsPatch:new{p00=z5, p10=z6, p11=p6, p01=p5}
patch3 = CoonsPatch:new{p00=p5, p10=p6, p11=q6, p01=q5}

-- Define the blocks, boundary conditions and
-- set the discretisation to join cells consistently.
nx0 = 50; nx1 = 100; nx2 = 80; ny = 30
grid0 = StructuredGrid:new{psurface=patch0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=patch1, niv=nx1+1, njv=ny+1}
grid2 = StructuredGrid:new{psurface=patch2, niv=nx2+1, njv=ny+1}
grid3 = StructuredGrid:new{psurface=patch3, niv=nx2+1, njv=ny+1}
subsonic_region = FluidBlock:new{grid=grid0, initialState=stagnation_gas}
supersonic_region = FluidBlock:new{grid=grid1, initialState=low_pressure_gas}
downstream_region = FluidBlock:new{grid=grid2, initialState=low_pressure_gas}
external_region = FluidBlock:new{grid=grid3, initialState=external_stream}

-- History locations near throat and exit
setHistoryPoint{ib=1, i=1, j=1}
setHistoryPoint{ib=1, i=nx1-1, j=1}

-- Boundary conditions for all of the blocks.
-- First stitch together adjoining blocks,
identifyBlockConnections()
-- then, directly specify the stagnation conditions for the subsonic inflow.
subsonic_region.bcList['west'] = InFlowBC_FromStagnation:new{stagnationState=stagnation_gas}
-- to get loads on thrust surface, add that boundary condition to the group
supersonic_region.bcList['north'] = WallBC_WithSlip:new{group="loads"}
downstream_region.bcList['east'] = OutFlowBC_Simple:new{}
external_region.bcList['east'] = OutFlowBC_Simple:new{}
external_region.bcList['west'] = InFlowBC_Supersonic:new{flowState=external_stream}

-- Do a little more setting of global data.
config.axisymmetric = 1
config.flux_calculator = "adaptive"
config.max_time = 1.0e-3  -- seconds
config.max_step = 50000
config.dt_init = 1.0e-7
config.dt_plot = 0.1e-3
config.dt_history = 10.0e-6
config.dt_loads = 1.0e-3
config.write_loads = true
