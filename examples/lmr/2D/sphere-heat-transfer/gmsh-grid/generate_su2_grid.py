import gmsh

gmsh.initialize()

# general settings
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.Algorithm", 5) # Delaunay algorithm
gmsh.option.setNumber("Mesh.SaveAll", 1)

# geometric properties
size = 1.0e-4
R = 6.6e-03 # sphere radius in meters

# define points
centre = gmsh.model.geo.add_point(0.0, 0.0, 0.0, size)
a = gmsh.model.geo.add_point(-R, 0.0, 0.0, size)
b = gmsh.model.geo.add_point(0.0, R, 0.0, size)
c = []
c.append(gmsh.model.geo.add_point(-1.5*R, 0.0, 0.0, size))
c.append(gmsh.model.geo.add_point(-1.5*R, R,   0.0, size))
c.append(gmsh.model.geo.add_point(-R,     2*R, 0.0, size))
c.append(gmsh.model.geo.add_point(0.0,    3*R, 0.0, size))

# define paths
gmsh.model.geo.add_line(c[-1], b)
gmsh.model.geo.add_circle_arc(b, centre, a)
gmsh.model.geo.add_line(a, c[0])
gmsh.model.geo.add_bezier(c)

# define curve loop
gmsh.model.geo.add_curve_loop([1,2,3,4])

# define surface
gmsh.model.geo.add_plane_surface([1])

# define boundary names
gmsh.model.geo.add_physical_group(1, [1], name="outflow")
gmsh.model.geo.add_physical_group(1, [2], name="wall")
gmsh.model.geo.add_physical_group(1, [3], name="slip_wall")
gmsh.model.geo.add_physical_group(1, [4], name="inflow")

gmsh.model.geo.synchronize()

# define prism layer
boundary_layer = gmsh.model.mesh.field.add("BoundaryLayer")
gmsh.model.mesh.field.set_numbers(boundary_layer, "CurvesList", [2])
gmsh.model.mesh.field.set_numbers(boundary_layer, "PointsList", [a,b])
gmsh.model.mesh.field.set_number(boundary_layer, "SizeFar", size)
gmsh.model.mesh.field.set_number(boundary_layer, "Size", 1e-6)
gmsh.model.mesh.field.set_number(boundary_layer, "Quads", 1)
gmsh.model.mesh.field.set_number(boundary_layer, "Ratio", 1.1)
gmsh.model.mesh.field.set_number(boundary_layer, "Thickness", 3.0e-4)
gmsh.model.mesh.field.setAsBoundaryLayer(boundary_layer)

# generate 2D grid
gmsh.model.mesh.generate(2)

# write grid to SU2 format
gmsh.write("sphere.su2")

# uncomment the line below to visualize the mesh
# gmsh.fltk.run()

gmsh.finalize()
