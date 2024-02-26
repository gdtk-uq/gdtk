
-- Geometry of the flow domain
L = 600.0e-3 -- metres
H = 0.20 * L
--
--         wall
--        c---------b
-- flow=> |         |
--        d         |
--          -\-     |
--    flow=>    -\- |
--        0         a ----> x
--
a = Vector3:new{x=L, y=0.0}; b = Vector3:new{x=L, y=H};
c = Vector3:new{x=0.0, y=H}; d = Vector3:new{x=0.0, y=3.0*H/4.0}
patch = CoonsPatch:new{p00=d, p10=a, p11=b, p01=c}

-- For testing
--niv = 640+1; njv = 384+1;
--cflist = {north=cfx, east=GeometricFunction:new{a=0.0005, r=1.2, N=njv, reverse=true},
--          south=cfx, west=GeometricFunction:new{a=0.0010, r=1.2, N=njv, reverse=true}}

-- For actual converged simulations
--niv = 86+1; njv = 46+1;
--cflist = {north=cfx, east=GeometricFunction:new{a=0.001, r=1.2, N=njv, reverse=true},
--          south=cfx, west=GeometricFunction:new{a=0.002, r=1.2, N=njv, reverse=true}}
--cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
-- For actual converged simulations
niv = 164+1; njv = 86+1;
cflist = {north=cfx, east=GeometricFunction:new{a=0.0002, r=1.2, N=njv, reverse=true},
          south=cfx, west=GeometricFunction:new{a=0.0004, r=1.2, N=njv, reverse=true}}
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}


grd = StructuredGrid:new{psurface=patch, niv=niv, njv=njv, cfList=cflist}
ugrid = UnstructuredGrid:new{sgrid=grd}

ugrid:set_boundaryset_tag(0, "inflow")
ugrid:set_boundaryset_tag(1, "outflow")
ugrid:set_boundaryset_tag(2, "inflow")
ugrid:set_boundaryset_tag(3, "wall")

ugrid:write_to_su2_file('grid.su2')
