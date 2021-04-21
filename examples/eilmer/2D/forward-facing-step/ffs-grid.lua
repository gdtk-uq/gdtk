-- ffs-grid.lua
-- Prepare grids for forward-facing-step.
-- $ e4shared --job=ffs --prep-grids
-- PJ 2021-04-21
--
config.dimensions = 2 -- Not really needed, because 2 is the default.
--
-- Geometry of the flow domain.
--
--    a2-----b2-----------------c2  O
-- I  |      |                  |   U
-- N  |   1  |        2         |   T
-- F  |      |                  |   F
-- L  a1-----b1-----------------c1  L
-- O  |   0  |  SLIP-WALLS          O
-- W  a0-----b0                     W
--
a0 = {x=0.0, y=0.0}; a1 = {x=0.0, y=0.2}; a2 = {x=0.0, y=1.0}
b0 = {x=0.6, y=0.0}; b1 = {x=0.6, y=0.2}; b2 = {x=0.6, y=1.0}
c1 = {x=3.0, y=0.2}; c2 = {x=3.0, y=1.0}
surf0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
surf1 = CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2}
surf2 = CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2}
--
-- Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = math.floor(0.6/dx); nbc = math.floor(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = math.floor(0.2/dx); n12 = math.floor(0.8/dx)
print("n01=", n01, "n12=", n12)
grid0 = StructuredGrid:new{psurface=surf0, niv=nab+1, njv=n01+1}
grid1 = StructuredGrid:new{psurface=surf1, niv=nab+1, njv=n12+1}
grid2 = StructuredGrid:new{psurface=surf2, niv=nbc+1, njv=n12+1}
--
-- Set up boundary condition tags for the ones that we care about.
bcTags01 = {west="INFLOW"}
bcTags2 = {east="OUTFLOW"}
--
-- Define the individual grids and stitch them together.
ga0 = registerGridArray{grid=grid0, nib=1, njb=1, fsTag="INFLOW", bcTags=bcTags01}
ga1 = registerGridArray{grid=grid1, nib=1, njb=4, fsTag="INFLOW", bcTags=bcTags01}
ga2 = registerGridArray{grid=grid2, nib=4, njb=4, fsTag="INFLOW", bcTags=bcTags2}
identifyGridConnections()
