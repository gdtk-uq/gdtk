print("Periodic shear layer -- construct the grids and make connections.")
-- PJ 2024-03-26: adapt from Eilmer4 example

config.dimensions = 2

H = 0.010 -- y-direction layer thickness in metres
L = 0.100 -- x-direction wavelength in metres
ymin = -20.0*H; ymax = 20.0*H
xmin = -L; xmax = L
domain = CoonsPatch:new{p00={x=xmin,y=ymin}, p10={x=xmax,y=ymin},
                        p01={x=xmin,y=ymax}, p11={x=xmax,y=ymax}}
factor = 2
niv = math.floor(60*factor)+1
njv = math.floor(120*factor)+1
grid0 = StructuredGrid:new{psurface=domain, niv=niv, njv=njv}

nib=2; njb = 3
iga = registerFluidGridArray{grid=grid0, nib=nib, njb=njb, fsTag="initial"}

-- At this point, the blocks within the GridArray will be connected but
-- we want the overall domain to be periodic in the x-direction, as well.
-- We manually connect the east faces of the blocks on the east end of the GridArray
-- to the west faces of the corresponding blocks on the west-end of the same array.
-- To get convenient access to the individual grids within the GridArray,
-- we first get the actual GridArray object from the global list and then
-- access the individual grids within that object.
ga = gridArraysList[iga+1]
for j=1, njb do
   connectGrids(ga.myGrids[1][j].id, "west", ga.myGrids[nib][j].id, "east")
end
