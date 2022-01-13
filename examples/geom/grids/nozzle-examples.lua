-- Example structured grids on converging-diverging nozzle.
-- Author: Rowan J. Gollan
-- Date: 2022-01-13
--
-- To generate this domain, I have used Peter's nozzle-optimize-bell example.
-- The chosen parameters are the ones that are an optimum as given
-- in technical note describing this example.
--
-- To run these examples and produce VTK files:
--
--   $ e4shared --custom-script --script-file=nozzle-examples.lua
--

niv = 100
njv = 20

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
theta_cone = math.rad(30.0) -- nominal straight-cone angle
theta_init = math.rad(30.0) -- starting angle for thrust nozzle
alpha = math.rad(0.0)  -- angle for setting b2 in Bezier curve
beta = math.rad(0.0)  -- angle for setting b3 in Bezier curve

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
theta_cone = math.rad(18.493) -- nominal straight-cone angle
theta_init = math.rad(29.679) -- starting angle for thrust nozzle
alpha = math.rad(-6.233)  -- angle for setting b2 in Bezier curve
beta = math.rad(5.956)  -- angle for setting b3 in Bezier curve

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

north = Polyline:new{segments={Line:new{p0=p0,p1=p1},
                               Arc:new{p0=p1,p1=p2,centre=centre_B},
                               Arc:new{p0=p2,p1=p3,centre=centre_A},
                               Arc:new{p0=p3,p1=p4,centre=centre_A},
                               Bezier:new{points={b0, b1, b2, b3, b4}}}
                    }
south = Line:new{p0=z0, p1=z5}
west = Line:new{p0=z0, p1=p0}
east = Line:new{p0=z5, p1=p5}


-- Example 1: Coons patch grid
coonsPatch = CoonsPatch:new{north=north, east=east, south=south, west=west}
coonsGrid = StructuredGrid:new{psurface=coonsPatch, niv=niv, njv=njv}
coonsGrid:write_to_vtk_file('nozzle-coons-patch-grid.vtk')

-- Example 2: AO-patch grid
aoPatch = AOPatch:new{north=north, east=east, south=south, west=west}
aoGrid = StructuredGrid:new{psurface=aoPatch, niv=niv, njv=njv}
aoGrid:write_to_vtk_file('nozzle-AO-patch-grid.vtk')

-- Example 2a: AO-patch grid with non-default background grid
aoPatch2 = AOPatch:new{north=north, east=east, south=south, west=west, nx=40, ny=10}
aoGrid2 = StructuredGrid:new{psurface=aoPatch2, niv=niv, njv=njv}
aoGrid2:write_to_vtk_file('nozzle-AO-patch-grid-2.vtk')

-- Example 3: Channel patch
channelPatch = ChannelPatch:new{south=south, north=north}
channelGrid = StructuredGrid:new{psurface=channelPatch, niv=niv, njv=njv}
channelGrid:write_to_vtk_file('nozzle-channel-patch-grid.vtk')

-- Example 4: Ruled surface
ruledSurface = RuledSurface:new{edge0=south, edge1=north}
ruledGrid = StructuredGrid:new{psurface=ruledSurface, niv=niv, njv=njv}
ruledGrid:write_to_vtk_file('nozzle-ruled-surface-grid.vtk')

-- Example 5: Control-point patch
-- Use channel patch as a means to generate control points
-- that are orthogonal to surface
ctrlPts = {}
ncpi = 11
ncpj = 4
dr = 1/(ncpi - 2)
ds = 1/(ncpj - 2)
r_pos = {0}
for i=1,ncpi-2 do
   r_pos[#r_pos+1] = dr/2 + (i-1)*dr
end
r_pos[#r_pos+1] = 1
s_pos = {0}
for j=1,ncpj-2 do
   if j == ncpj-2 then ds = 0.9*ds end
   s_pos[#s_pos+1] = ds/2 + (j-1)*ds
end
s_pos[#s_pos+1] = 1

for i=1,ncpi-1 do -- npci-1: don't use channel path on final i
   ctrlPts[i] = {}
   for j=1,ncpj do
      ctrlPts[i][j] = channelPatch(r_pos[i], s_pos[j])
      -- Ensure that second row of control points are orthogonal
      if j == 2 then ctrlPts[i][j].x = ctrlPts[i][1].x end
   end
end
-- at imax
i = ncpi
ctrlPts[i] = {}
for j=1,ncpj do
   ctrlPts[i][j] = coonsPatch(r_pos[i], s_pos[j])
end

dofile('write_ctrl_points.lua')
writeCtrlPts2VTK(ctrlPts, 'nozzle-ctrl-pts.vtk')

ctrlPtPatch = ControlPointPatch:new{north=north, south=south, east=east, west=west,
                                    control_points=ctrlPts}
ctrlPtGrid = StructuredGrid:new{psurface=ctrlPtPatch, niv=niv, njv=njv}
ctrlPtGrid:write_to_vtk_file('nozzle-ctrl-pt-patch-grid.vtk')









