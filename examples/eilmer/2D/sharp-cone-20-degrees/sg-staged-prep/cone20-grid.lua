-- cone20-grid.lua
-- Simple grid-specification file for Eilmer4.
-- PJ 2021-03-16 adapted from the cone20.lua example.
-- Process with the command:
-- $ e4shared --prep-grid --job=cone20
--
config.dimensions = 2
--
-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.2, y=0.0}
c = Vector3:new{x=1.0, y=0.29118}
d = Vector3:new{x=1.0, y=1.0}
e = Vector3:new{x=0.2, y=1.0}
f = Vector3:new{x=0.0, y=1.0}
ab = Line:new{p0=a, p1=b} -- lower boundary, axis
bc = Line:new{p0=b, p1=c} -- lower boundary, cone surface
fe = Line:new{p0=f, p1=e}; ed = Line:new{p0=e, p1=d} -- upper boundary
af = Line:new{p0=a, p1=f} -- vertical line, inflow
be = Line:new{p0=b, p1=e} -- vertical line, between quads
cd = Line:new{p0=c, p1=d} -- vertical line, outflow
quad0 = makePatch{north=fe, east=be, south=ab, west=af}
quad1 = makePatch{north=ed, east=cd, south=bc, west=be, gridType="ao"}
--
-- Mesh the patches, with particular discretisation and
-- register the grids for use in the simulation setup phase.
nx0 = 10; nx1 = 30; ny = 40
grid0 = registerGrid{
      grid=StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1},
      tag="inflow-gas",
      bcTags={west="inflow"}
}
grid1 = registerGrid{
   grid=StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1},
   tag="initial-gas",
   bcTags={east="outflow"}
}
identifyGridConnections()
--
dofile("sketch-domain.lua")
