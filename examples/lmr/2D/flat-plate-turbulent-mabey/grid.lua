-- ==========================================================
-- SU2 grid generation script
-- ==========================================================
-- adapted from the eilmer 4.0 script
--
--         wall
--        c---------b
-- flow=> |         |
--        d         |
--          -\-     |
--    flow=>    -\- |
--        0         a ----> x
--

-- define some grid refinement paramters
ni0 = 200; nj0 = 100 -- We'll scale discretization off these values
factor = 2.0
ni0 = ni0*factor; nj0 = nj0*factor

-- define some geometric parameters
L = 1.40 -- metres
H = 0.40 * L

-- define nodes
a = Vector3:new{x=L, y=0.0}
b = Vector3:new{x=L, y=H}
c = Vector3:new{x=0.0, y=H}
d = Vector3:new{x=0.0, y=3.0*H/4.0}

-- define clustering
rcfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
rcfyW = RobertsFunction:new{end0=false,end1=true,beta=1.00074}
rcfyE = RobertsFunction:new{end0=false,end1=true,beta=1.00014}

-- define structured grids and convert to unstructured grids
-- Note:  0,1,2,3 == west,east,south,north
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=d, p10=a, p11=b, p01=c},
                          cfList={north=rcfx, east=rcfyE, south=rcfx, west=rcfyW},
                          niv=ni0, njv= nj0+1}
ugrid = UnstructuredGrid:new{sgrid=grid}
ugrid:set_boundaryset_tag(0, "inflow")
ugrid:set_boundaryset_tag(1, "outflow")
ugrid:set_boundaryset_tag(2, "inflow")
ugrid:set_boundaryset_tag(3, "wall")

-- write out unstructured grid in .su2 format
ugrid:write_to_su2_file("grid.su2")
