-- Example for reference manual
-- RJG, 2025-02-12
line = Line:new{p0={x=0.0}, p1={x=1.0}}
npts = 15

geol = GeometricFunction:new{a=0.002, r=1.1, N=60}
g0 = StructuredGrid:new{path=line, niv=npts, cf=geol}
g0:write_to_vtk_file("geometric-cluster-left.vtk")

geol2 = GeometricFunction:new{a=0.004, r=1.2, N=60}
g1 = StructuredGrid:new{path=line, niv=npts, cf=geol2}
g1:write_to_vtk_file("geometric-cluster-left2.vtk")

geor = GeometricFunction:new{a=0.004, r=1.2, N=60, reverse=true}
g2 = StructuredGrid:new{path=line, niv=npts, cf=geor}
g2:write_to_vtk_file("geometric-cluster-right.vtk")
