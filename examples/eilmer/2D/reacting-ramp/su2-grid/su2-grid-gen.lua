-- SU2 grid generation script
-- e4shared --custom-script --script-file=su2-grid-gen.lua

theta = 15 * math.pi/180.0 -- ramp angle, radians
L = 1.0 -- ramp length, meters

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=-0.5, y=0.0}
b = Vector3:new{x=0.0, y=0.0}
c = Vector3:new{x=1.0, y=L*math.tan(theta)}
d = Vector3:new{x=1.0, y=1.0}
e = Vector3:new{x=0.0, y=1.0}
f = Vector3:new{x=-0.5, y=1.0}
ab = Line:new{p0=a, p1=b}; be = Line:new{p0=b, p1=e} -- lower boundary including cone surface
fe = Line:new{p0=f, p1=e}; ed = Line:new{p0=e, p1=d} -- upper boundary
af = Line:new{p0=a, p1=f}; be = Line:new{p0=b, p1=e} -- vertical lines
bc = Line:new{p0=b, p1=c}; cd = Line:new{p0=c, p1=d} -- vertical lines
quad0 = makePatch{north=fe, east=be, south=ab, west=af}
quad1 = makePatch{north=ed, east=cd, south=bc, west=be} --, gridType="ao"}
-- Mesh the patches, with particular discretisation.
nx0 = 20; nx1 = 40; ny = 40
grid0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1}

-- join all structured grids into a single grid
singleBlockGrid = grid0
singleBlockGrid:joinGrid(grid1, "east")

-- convert structured grid to unstructured grid
usgrid=UnstructuredGrid:new{sgrid=singleBlockGrid}

-- write out unstructured grid in .su2 format
usgrid:write_to_su2_file("ramp15.su2")
