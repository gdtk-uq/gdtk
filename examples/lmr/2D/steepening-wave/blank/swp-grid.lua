
-- Define the domain
N=128
nProcesses = 4
xmin = 0.0; xmax = 2.0

a = Vector3:new{x = xmin, y = -0.01}
b = Vector3:new{x = xmax, y = -0.01}
c = Vector3:new{x = a.x, y = 0.01}
d = Vector3:new{x = b.x, y = 0.01}

south0 = Line:new{p0 = a, p1 = b}
west0 = Line:new{p0 = a, p1 = c}
north0 = Line:new{p0 = c, p1 = d}
east0 = Line:new{p0 = b, p1 = d}

psurf0 = makePatch{north = north0, east = east0, south = south0, west = west0}

grid0 = registerFluidGridArray{
   nib=nProcesses, njb=1,
   grid=StructuredGrid:new{psurface=psurf0, niv=N+1, njv=4},
   fsTag="initial",
}
identifyGridConnections()
