-- mabey.lua : Turbulent flow over a flat plate
-- Dimir Y.X. Pot, Wilson Y.K. Chan, 2018-03-14
-- Ported from Eilmer3
--  Mabey test case (AGARDograph 223 - Test series 7402)
--  (Referenced from Fernholz & Finley (1977),
--  AGARDograph No. 223, "A critical compilation of
--  compressible turbulent boundary layer data.")
--

-- Geometry of the flow domain
L = 1.40 -- metres
H = 0.40 * L
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
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
cflist = {north=cfx, east=RobertsFunction:new{end0=false,end1=true,beta=1.0014},
	  south=cfx, west=RobertsFunction:new{end0=false,end1=true,beta=1.0074}}
sgrid = StructuredGrid:new{psurface=patch, niv=129/3, njv=97/3, cfList=cflist}
ugrid = UnstructuredGrid:new{sgrid=sgrid}
ugrid:write_to_su2_file('grid.su2')
