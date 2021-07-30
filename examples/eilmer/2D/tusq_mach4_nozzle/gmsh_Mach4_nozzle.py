import gmsh
import sys
import numpy as np

nFactor = 1

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("mach4nozzle")

# Let's read in our nozzle wall profile
xNoz, yNoz = np.loadtxt('03_mach4_RDE_nozzle_bl_corr.dat', unpack=True)

# Add a 'cluster' for the walls
dCell = 1e-3/nFactor

# Go through and create points for each of these
nozProf = []
for x,y in zip(xNoz, yNoz):
    tmp = gmsh.model.geo.addPoint(x, y, 0.0, dCell)
    nozProf.append(tmp)

nozzle = gmsh.model.geo.addSpline(nozProf)

# Let's try using my previous geometry code from lua to calculate the points
bg = [xNoz[0], yNoz[0]]
lg = [0.0, 0.0]
R = 12 # From Sivells nozzle design
throat_r = yNoz[0] # From Sivell output
radius1 = throat_r * R # From Sivells definitions
radius2 = 20e-3
cg = [0.0, bg[1] + radius1] # Main arc origin
dx_offset = np.sqrt(radius1**2 - (radius1 - 0.005)**2)
dg = [bg[0] - dx_offset, bg[1] + 0.005]
theta2 = np.arccos(dx_offset / radius1)
rad2x = dg[0] + radius2 * np.cos(theta2)
rad2y = dg[1] + radius2 * np.sin(theta2)
eg = [rad2x, rad2y] # Small arc origin
fg = [eg[0] - radius2, eg[1]]

'''
h____________g      c
|            |                 __/|
|             \f  e         __/   |
|              \d        __/      |
|               \___b____/        |
|                                 |
|                                 |
|i________________________________|j

c = centre of large radius (b --> d)
e = centre of small radius (d --> f)
'''

# Now convert these to points
b = nozProf[0]
c = gmsh.model.geo.addPoint(cg[0], cg[1], 0.0)
d = gmsh.model.geo.addPoint(dg[0], dg[1], 0.0)
e = gmsh.model.geo.addPoint(eg[0], eg[1], 0.0)
f = gmsh.model.geo.addPoint(fg[0], fg[1], 0.0)
# And add the extras
g = gmsh.model.geo.addPoint(fg[0], 0.065, 0.0, dCell)
h = gmsh.model.geo.addPoint(fg[0] - 0.13, 0.065, 0.0, dCell)
i = gmsh.model.geo.addPoint(fg[0] - 0.13, 0.0, 0.0, dCell*5)
j = gmsh.model.geo.addPoint(xNoz[-1], 0.0, 0.0, dCell*5)
# Set the mesh size for b, d and f as much smaller, but the same
gmsh.model.geo.mesh.setSize([(0, b), (0, d), (0, f)], dCell)
# Adding in some new points for my unstructured b.l. attempt
blOffset = 0.005
# Offset the large radius
radius3 = radius1 + blOffset
theta3 = np.pi/2 - theta2
BLd_x = cg[0] - radius3 * np.sin(theta3)
BLd_y = cg[1] - radius3 * np.cos(theta3)
radius4 = radius2 + blOffset
# Once we have the geometric components we can add the new points
BLb = gmsh.model.geo.addPoint(xNoz[0], yNoz[0]-blOffset, 0.0, dCell)
BLd = gmsh.model.geo.addPoint(BLd_x, BLd_y, 0.0, dCell)
BLf = gmsh.model.geo.addPoint(eg[0] - radius4, eg[1], 0.0, dCell)
BLg = gmsh.model.geo.addPoint(fg[0] - blOffset, 0.065, 0.0, dCell)
BLk = gmsh.model.geo.addPoint(xNoz[-1], yNoz[-1]-blOffset, 0.0, dCell)
# Add in some new test section points
l0 = gmsh.model.geo.addPoint(xNoz[-1]+0.01, yNoz[-1]+0.01, 0.0, dCell)
l = gmsh.model.geo.addPoint(xNoz[-1]+0.01, yNoz[-1]+0.1, 0.0, dCell*10)
m = gmsh.model.geo.addPoint(xNoz[-1]+0.2, 0.0, 0.0, dCell*10)
n = gmsh.model.geo.addPoint(xNoz[-1]+0.2, yNoz[-1]+0.1, 0.0, dCell*10)
# Now we can compute all the lines we need
blk0n = gmsh.model.geo.addLine(BLg, g)
blk0e = gmsh.model.geo.addLine(g, f)
blk0s = gmsh.model.geo.addLine(BLf, f)
blk0w = gmsh.model.geo.addLine(BLg, BLf)
blk1e = gmsh.model.geo.addCircleArc(f, e, d)
blk1s = gmsh.model.geo.addLine(BLd, d)
blk1w = gmsh.model.geo.addCircleArc(BLf, e, BLd)
blk2e = gmsh.model.geo.addCircleArc(d, c, b)
blk2s = gmsh.model.geo.addLine(BLb, b)
blk2w = gmsh.model.geo.addCircleArc(BLd, c, BLb)
# Now offset the nozzle profile
BL0 = gmsh.model.geo.copy([(1, nozzle)])
gmsh.model.geo.translate(BL0, 0.0, -blOffset, 0)
blk3e = gmsh.model.geo.addLine(BLk, nozProf[-1])
# Now the unstructured block 4 boundaries (note weird nomenclature)
blk4n = gmsh.model.geo.addLine(h, BLg)
blk4e = gmsh.model.geo.addLine(j, BLk)
blk4s = gmsh.model.geo.addLine(i, j)
blk4w = gmsh.model.geo.addLine(i, h)
# Add in the test section
blk5wn0 = gmsh.model.geo.addLine(nozProf[-1], l0)
blk5wn = gmsh.model.geo.addLine(l0, l)
blk5n = gmsh.model.geo.addLine(l,n)
blk5e = gmsh.model.geo.addLine(m,n)
blk5s = gmsh.model.geo.addLine(j,m)
# Set up my boundary conditions
wall = gmsh.model.addPhysicalGroup(1, [blk4n, blk0n, blk0e, blk1e, blk2e, nozzle, blk5wn0, blk5wn])
gmsh.model.setPhysicalName(1, wall, "WallFixedT")

inflow = gmsh.model.addPhysicalGroup(1, [blk4w])
gmsh.model.setPhysicalName(1, inflow, "Inflow")

axis = gmsh.model.addPhysicalGroup(1, [blk4s, blk5s])
gmsh.model.setPhysicalName(1, axis, "Axis")

outflow = gmsh.model.addPhysicalGroup(1, [blk5n, blk5e])
gmsh.model.setPhysicalName(1, outflow, "Outflow")

# Make sure gmsh saves all elements!
gmsh.option.setNumber("Mesh.SaveAll", 1)

# Now make the loop (taken from t1 example)
curve0 = gmsh.model.geo.addCurveLoop([blk0n, blk0e, blk1e, blk2e, nozzle, -blk3e, -int(BL0[0][1]), -blk2w, -blk1w, -blk0w])
curve1 = gmsh.model.geo.addCurveLoop([blk4n, blk0w, blk1w, blk2w, int(BL0[0][1]), -blk4e, -blk4s, blk4w])
curve2 = gmsh.model.geo.addCurveLoop([blk4e, blk3e, blk5wn0, blk5wn, blk5n, -blk5e, -blk5s])

# And the surface (t1 example)
surf0 = gmsh.model.geo.addPlaneSurface([curve0])
surf1 = gmsh.model.geo.addPlaneSurface([curve1])
surf2 = gmsh.model.geo.addPlaneSurface([curve2])

gmsh.model.geo.synchronize()

# Convert to unstructured
gmsh.model.mesh.setTransfiniteCurve(blk0n, 10*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk3e, 10*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk0e, 50*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk0w, 50*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk1e, 50*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk1w, 50*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk2e, 50*nFactor)
gmsh.model.mesh.setTransfiniteCurve(blk2w, 50*nFactor)
gmsh.model.mesh.setTransfiniteCurve(nozzle, 250*nFactor)
gmsh.model.mesh.setTransfiniteCurve(BL0[0][1], 250*nFactor)

gmsh.model.mesh.setTransfiniteSurface(surf0, "Left", [BLg, g, nozProf[-1], BLk])

gmsh.model.geo.mesh.setRecombine(2, surf0)

gmsh.model.mesh.generate(2)
gmsh.write('mach4nozzleStage1.su2')

# Launch the GUI to see the results:
# Uncomment next two lines for GUI
#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.finalize()