-- 3D Flat Plate pair for testing corner cell issues.
-- 2023-02-27: NNG

-- dimensions in (m)
Lx = 0.15
Ly = 0.05
Lz = 0.05
il = 0.01

nicell = 70; nicell0 = 10; njcell = 48; nkcell = 48

clusterx0 = RobertsFunction:new{end0=false, end1=true, beta=1.1}
clusterx1 = RobertsFunction:new{end0=true, end1=false, beta=1.05}
clustery = GeometricFunction:new{a=0.001, r=1.2, N=njcell}
clusterz = GeometricFunction:new{a=0.001, r=1.2, N=nkcell}

cfList0 = {edge01=clusterx0, edge12=clustery, edge32=clusterx0, edge03=clustery,
	       edge45=clusterx0, edge56=clustery, edge76=clusterx0, edge47=clustery,
	       edge04=clusterz, edge15=clusterz, edge26=clusterz, edge37=clusterz}

p000 = Vector3:new{x=-il, y=0.0, z=0.0}
p100 = Vector3:new{x=0.0,  y=0.0, z=0.0}
p010 = Vector3:new{x=-il, y=Ly,  z=0.0}
p110 = Vector3:new{x=0.0,  y=Ly,  z=0.0}
p001 = Vector3:new{x=-il, y=0.0, z=Lz}
p101 = Vector3:new{x=0.0,  y=0.0, z=Lz}
p011 = Vector3:new{x=-il, y=Ly,  z=Lz}
p111 = Vector3:new{x=0.0,  y=Ly,  z=Lz}
grid0 = StructuredGrid:new{pvolume=TFIVolume:new{vertices={p000, p100, p110, p010, p001, p101, p111, p011}}, niv=nicell0+1, njv=njcell+1, nkv=nkcell+1, cfList=cfList0}
ugrid0 = UnstructuredGrid:new{sgrid=grid0}
ugrid0:set_boundaryset_tag(3, "inflow")

cfList1 = {edge01=clusterx1, edge12=clustery, edge32=clusterx1, edge03=clustery,
	       edge45=clusterx1, edge56=clustery, edge76=clusterx1, edge47=clustery,
	       edge04=clusterz, edge15=clusterz, edge26=clusterz, edge37=clusterz}
p000 = Vector3:new{x=0.0, y=0.0, z=0.0}
p100 = Vector3:new{x=Lx,  y=0.0, z=0.0}
p010 = Vector3:new{x=0.0, y=Ly,  z=0.0}
p110 = Vector3:new{x=Lx,  y=Ly,  z=0.0}
p001 = Vector3:new{x=0.0, y=0.0, z=Lz}
p101 = Vector3:new{x=Lx,  y=0.0, z=Lz}
p011 = Vector3:new{x=0.0, y=Ly,  z=Lz}
p111 = Vector3:new{x=Lx,  y=Ly,  z=Lz}
grid1 = StructuredGrid:new{pvolume=TFIVolume:new{vertices={p000, p100, p110, p010, p001, p101, p111, p011}}, niv=nicell+1, njv=njcell+1, nkv=nkcell+1, cfList=cfList1}
ugrid1 = UnstructuredGrid:new{sgrid=grid1}
ugrid1:set_boundaryset_tag(1, "outflow")
ugrid1:set_boundaryset_tag(2, "wall")
ugrid1:set_boundaryset_tag(5, "wall")

ugrid1:joinGrid(ugrid0)
ugrid1:write_to_su2_file('grid.su2')
