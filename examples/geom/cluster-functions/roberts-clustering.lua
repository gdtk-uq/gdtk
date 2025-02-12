-- Example for reference manual
-- RJG, 2025-02-12
line = Line:new{p0={x=0.0}, p1={x=1.0}}
npts = 15

rbl = RobertsFunction:new{end0=true, end1=false, beta=1.05}
g0 = StructuredGrid:new{path=line, niv=npts, cf=rbl}
g0:write_to_vtk_file("roberts-cluster-left.vtk")

rbr = RobertsFunction:new{end0=false, end1=true, beta=1.05}
g0 = StructuredGrid:new{path=line, niv=npts, cf=rbr}
g0:write_to_vtk_file("roberts-cluster-right.vtk")

rbb = RobertsFunction:new{end0=true, end1=true, beta=1.05}
g0 = StructuredGrid:new{path=line, niv=npts, cf=rbb}
g0:write_to_vtk_file("roberts-cluster-both.vtk")

