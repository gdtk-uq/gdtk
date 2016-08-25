-- back.py
-- Conical nozzle from Back, Massier and Gier (1965)
-- Peter J. 2015-10-21 adpated from the Python version.
--          
config.title = "Flow through a conical nozzle."
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
fixed_pressure = false
if fixed_pressure then
   -- The stagnation gas represents a reservoir condition.
   stagnation_gas = FlowState:new{p=500.0e3, T=300.0}
else
   -- We'll specify a mass_flux and let the pressure be adjusted.
   -- We still need a stagnation condition and the temperature is fixed.
   stagnation_gas = FlowState:new{p=400.0e3, T=300.0}
   -- Note that we start with a known wrong stagnation pressure.
end
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
-- Compute the details of the conical nozzle
p4 = Vector3:new{x=R_curve*math.sin(theta), y=height-R_curve*math.cos(theta)}
z4 = Vector3:new{x=p4.x, y=0.0}
L_cone = L_nozzle - p4.x
p5 = p4 + Vector3:new{x=L_cone, y=L_cone*math.tan(theta)}
z5 = Vector3:new{x=p5.x, y=0.0}

north0 = Polyline:new{segments={Line:new{p0=p0,p1=p1},
				Arc:new{p0=p1,p1=p2,centre=centre_B},
				Arc:new{p0=p2,p1=p3,centre=centre_A}}}
east0west1 = Line:new{p0=z3, p1=p3}
south0 = Line:new{p0=z0, p1=z3}
west0 = Line:new{p0=z0, p1=p0}
north1 = Polyline:new{segments={Arc:new{p0=p3,p1=p4,centre=centre_A},
				Line:new{p0=p4,p1=p5}}}
east1 = Line:new{p0=z5, p1=p5}
south1 = Line:new{p0=z3, p1=z5}
patch0 = makePatch{north=north0, east=east0west1, south=south0, west=west0}
patch1 = makePatch{north=north1, east=east1, south=south1, west=east0west1}

-- Define the blocks, boundary conditions and set the discretisation.
nx0 = 50; nx1 = 60; ny = 30
grid0 = StructuredGrid:new{psurface=patch0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=patch1, niv=nx1+1, njv=ny+1}
subsonic_region = SBlock:new{grid=grid0, fillCondition=stagnation_gas, label="subsonic-region"}
supersonic_region = SBlock:new{grid=grid1, fillCondition=low_pressure_gas, label="supersonic-region"} 
-- History locations near throat and exit
setHistoryPoint{ib=1, i=1, j=1}
setHistoryPoint{ib=1, i=nx1-1, j=1}
identifyBlockConnections()
if fixed_pressure then
   -- Directly specify the stagnation conditions for the subsonic inflow.
   subsonic_region.bcList[west] = InFlowBC_FromStagnation:new{stagCondition=stagnation_gas,
							      label="inflow-boundary"}
else
   -- Specify the inflow mass_flux (kg/s/m^^2) across inlet and guess the stagnation condition.
   subsonic_region.bcList[west] = InFlowBC_FromStagnation:new{stagCondition=stagnation_gas,
							      mass_flux=275.16, relax_factor=0.2,
							      label="inflow-boundary"}
end
supersonic_region.bcList[east] = OutFlowBC_Simple:new{label="outflow-boundary"}

-- Do a little more setting of global data.
config.axisymmetric = 1
config.flux_calculator = "adaptive"
config.max_time = 4.0e-3  -- seconds
config.max_step = 50000
config.dt_init = 1.0e-7
config.dt_plot = 0.2e-3
config.dt_history = 10.0e-6

