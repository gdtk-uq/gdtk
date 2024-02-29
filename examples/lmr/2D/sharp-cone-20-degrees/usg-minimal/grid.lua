-- grid.lua
print("Set up unstructured grid for a Mach 1.5 flow over a 20 degree cone.")
--
-- 1. Geometry
a0 = {x=0.0, y=0.0};     a1 = {x=0.0, y=1.0}
b0 = {x=0.2, y=0.0};     b1 = {x=0.2, y=1.0}
c0 = {x=1.0, y=0.29118}; c1 = {x=1.0, y=1.0}
--
quad0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
quad1 = AOPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
--
-- 2. Grids
sgrid0 = StructuredGrid:new{psurface=quad0, niv=11, njv=41}
sgrid1 = StructuredGrid:new{psurface=quad1, niv=31, njv=41}
-- We make a single structured grid and then
-- construct an unstructured grid from it.
sgrid0:joinGrid(sgrid1, "east")
grid0 = registerFluidGrid{
   grid=UnstructuredGrid:new{sgrid=sgrid0},
   fsTag="inflow",
   bcTags={[Face.west]="inflow", [Face.east]="outflow",
           [Face.south]="wall", [Face.north]="outflow"}
}
