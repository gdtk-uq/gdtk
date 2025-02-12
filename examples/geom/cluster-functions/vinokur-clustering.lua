-- Example for reference manual
-- RJG, 2025-02-12
line = Line:new{p0={x=0.0}, p1={x=1.0}}
npts = 15

vin0 = VinokurFunction:new{n=15, s1=0.01, sn=0.02}
g0 = StructuredGrid:new{path=line, niv=npts, cf=vin0}
g0:write_to_vtk_file("vinokur-cluster-ex0.vtk")

vin1 = VinokurFunction:new{n=15, s1=0.005, sn=0.02}
g1 = StructuredGrid:new{path=line, niv=npts, cf=vin1}
g1:write_to_vtk_file("vinokur-cluster-ex1.vtk")

vin2 = VinokurFunction:new{n=15, s1=0.005, sn=0.03}
g2 = StructuredGrid:new{path=line, niv=npts, cf=vin2}
g2:write_to_vtk_file("vinokur-cluster-ex2.vtk")

