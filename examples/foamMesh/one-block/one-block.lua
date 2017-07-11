-- This example builds a unit square and converts it to a wedge.
-- Corners of blocks
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=1.0, y=0.0}
c = Vector3:new{x=0.0, y=1.0}
d = Vector3:new{x=1.0, y=1.0}

-- Lines connecting blocks.
ab = Line:new{p0=a, p1=b}
ac = Line:new{p0=a, p1=c}
bd = Line:new{p0=b, p1=d}
cd = Line:new{p0=c, p1=d}

-- Define patch (which are parametric surfaces, no discretisation at this point.)
quad0 = CoonsPatch:new{north=cd, east=bd, south=ab, west=ac}

-- Define 2D grid on patch, clustering can be added if desired
nxcells = 10; nycells = 5
grid0 = StructuredGrid:new{psurface=quad0, niv=nxcells+1, njv=nycells+1}

-- Create a 3D wedge grid from the 2D structured grid
wedge0 = grid0:makeWedgeGrid{dtheta=0.2, symmetric=true, label="wedge0"}

-- Define OpenFoam block (a "grid" with labels)
blk0 = FoamBlock:new{grid=wedge0,
		     bndry_labels={west="i-00", south="w-00", east="o-00"}}








