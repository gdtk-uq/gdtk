-- Example for reference manual
-- RJG, 2025-02-12
line = Line:new{p0={x=0.0}, p1={x=1.0}}
npts = 15

gaussl = GaussianFunction:new{m=0.1, s=0.1, ratio=0.1}
g0 = StructuredGrid:new{path=line, niv=npts, cf=gaussl}
g0:write_to_vtk_file("gaussian-cluster-left.vtk")

gaussm = GaussianFunction:new{m=0.5, s=0.2, ratio=0.1}
g1 = StructuredGrid:new{path=line, niv=npts, cf=gaussm}
g1:write_to_vtk_file("gaussian-cluster-mid.vtk")

gaussr = GaussianFunction:new{m=0.75, s=0.1, ratio=0.2}
g2 = StructuredGrid:new{path=line, niv=npts, cf=gaussr}
g2:write_to_vtk_file("gaussian-cluster-right.vtk")
