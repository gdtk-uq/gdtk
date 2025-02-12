-- Example for reference manual
-- RJG, 2025-02-12
line = Line:new{p0={x=0.0}, p1={x=1.0}}
npts = 15

vin0 = VinokurGeomHybridFunction:new{n=15, s1=0.01, n1=3, r1=1.1, sn=0.02, nn=5, rn=1.1}
g0 = StructuredGrid:new{path=line, niv=npts, cf=vin0}
g0:write_to_vtk_file("vin-geom-cluster-ex0.vtk")

vin1 = VinokurGeomHybridFunction:new{n=15, s1=0.01, n1=9, r1=1.3, sn=0.02, nn=5, rn=1.05}
g1 = StructuredGrid:new{path=line, niv=npts, cf=vin1}
g1:write_to_vtk_file("vin-geom-cluster-ex1.vtk")

vin2 = VinokurGeomHybridFunction:new{n=15, s1=0.03, n1=5, r1=1.3, sn=0.03, nn=5, rn=1.05}
g2 = StructuredGrid:new{path=line, niv=npts, cf=vin2}
g2:write_to_vtk_file("vin-geom-cluster-ex2.vtk")


