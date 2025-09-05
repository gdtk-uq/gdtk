
-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=1.0, y=0.0}
c = Vector3:new{x=0.0, y=1.0}
d = Vector3:new{x=1.0, y=1.0}
ab = Line:new{p0=a, p1=b}
bd = Line:new{p0=b, p1=d}
ac = Line:new{p0=a, p1=c}
cd = Line:new{p0=c, p1=d}
quad = makePatch{north=cd, east=bd, south=ab, west=ac}

-- Mesh the patches, with particular discretisation.
nx = 16; ny = 16;
ni1 = nx + 1; nj1 = ny + 1

-- 2. Grids
grid1 = registerFluidGrid{
   grid=StructuredGrid:new{psurface=quad, niv=ni1, njv=nj1},
   fsTag="initial",
   bcTags={
      north="north",
      south="south",
      east="east",
      west="west"
   }
}

identifyGridConnections()

