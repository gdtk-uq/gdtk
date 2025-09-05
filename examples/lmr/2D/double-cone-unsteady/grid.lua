-- ==========================================================
-- SU2 grid generation script
-- ==========================================================
-- adapted from a script provided by H. G. Hornung 28-Jan-2020 

-- define some grid refinement paramters
ni0 = 80; nj0 = 80 -- We'll scale discretization off these values
factor = 2.0
ni0 = ni0*factor; nj0 = nj0*factor

-- define some geometric parameters
L1 = 0.1 -- m
L2 = 0.05 -- m
th1 = 10.0 -- deg
th2 = 70.0 -- deg
th1r = math.rad(th1)
th2r = math.rad(th2)

-- define nodes
a00 = Vector3:new{x=-L1/2.5, y=0.0}
a0 = Vector3:new{x=0.0, y=0.0}
b0 = Vector3:new{x=L1*math.cos(th1r), y=L1*math.sin(th1r)} -- junction between cones
b1 = Vector3:new{x=b0.x-L1/1.0*math.sin((th1r+th2r)/2), y=b0.y+L1/1.0*math.cos((th1r+th2r)/2)} 
a11 = Vector3:new{x=a00.x, y=b1.y}
a1 = Vector3:new{x=0.0-L1/3, y=a11.y}
c0 = Vector3:new{x=b0.x+L2*math.cos(th2r), y=b0.y+L2*math.sin(th2r)} -- downstream-edge of second cone
c1 = Vector3:new{x=c0.x-L2/1.2,y=c0.y+L1/1.5} -- out in the free stream
d0 = Vector3:new{x=c0.x+L1/4, y=c0.y} -- down-stream edge of domain
d1 = Vector3:new{x=d0.x,y=c1.y}

-- define clustering
rcfx = RobertsFunction:new{end0=true,end1=true,beta=1.2}
rcfy = RobertsFunction:new{end0=true,end1=false,beta=1.005}

-- define structured grids and convert to unstructured grids
-- Note:  0,1,2,3 == west,east,south,north
grid0 = StructuredGrid:new{psurface=CoonsPatch:new{p00=a00,p10=a0,p11=a1,p01=a11},
                            cfList={east=rcfy,west=rcfy}, niv=81, njv= nj0+1}
ugrid0 = UnstructuredGrid:new{sgrid=grid0}
ugrid0:set_boundaryset_tag(3, "inflow")
ugrid0:set_boundaryset_tag(2, "slip")
ugrid0:set_boundaryset_tag(0, "inflow")

grid1 = StructuredGrid:new{psurface=CoonsPatch:new{p00=a0,p10=b0,p11=b1,p01=a1},
			    cfList= {north=rcfx,east=rcfy,south=rcfx,west=rcfy},
			    niv=ni0+1, njv=nj0+1}
ugrid1 = UnstructuredGrid:new{sgrid=grid1}
ugrid1:set_boundaryset_tag(3, "inflow")
ugrid1:set_boundaryset_tag(2, "wall")
ugrid1:set_boundaryset_tag(0, "inflow")

grid2 = StructuredGrid:new{psurface=CoonsPatch:new{p00=b0,p10=c0,p11=c1,p01=b1},
			    cfList= {north=None,east=rcfy,south=None,west=rcfy},
			    niv=ni0+1, njv=nj0+1}   
ugrid2 = UnstructuredGrid:new{sgrid=grid2}
ugrid2:set_boundaryset_tag(3, "inflow")
ugrid2:set_boundaryset_tag(2, "wall")

grid3 = StructuredGrid:new{psurface=CoonsPatch:new{p00=c0,p10=d0,p11=d1,p01=c1},
			    cfList= {north=None,east=rcfy,south=None,west=rcfy},
			    niv=tonumber(ni0/2.0)+1, njv=nj0+1}   
ugrid3 = UnstructuredGrid:new{sgrid=grid3}
ugrid3:set_boundaryset_tag(3, "inflow")
ugrid3:set_boundaryset_tag(1, "outflow")
ugrid3:set_boundaryset_tag(2, "wall")

-- join all structured grids into a single grid
ugrid = ugrid0
ugrid:joinGrid(ugrid1)
ugrid:joinGrid(ugrid2)
ugrid:joinGrid(ugrid3)

-- write out unstructured grid in .su2 format
ugrid:write_to_su2_file("grid.su2")
