-- An example of a clipped cavity, as found in
-- Sections 2.1.9 and 2.1.10 of the OpenFOAM manual.
--
--             BC=w-01
--          f-----g-----h
--          |  b3 |  b4 |
--  BC=w-00 |     |     | BC=w-00  
--          c-----d-----e
--          |  b0 |
--          |     |
--          a-----b
--          BC=w-00
--
-- Authors: IJ and RJG
-- Date: 2017-06-29

-- Global settings go first
axisymmetric = false
turbulence_model = "S-A" -- other option is: "k-epsilon"

-- Corners of blocks
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.6, y=0.0}
c = Vector3:new{x=0.0, y=0.4}
d = Vector3:new{x=0.6, y=0.4}
e = Vector3:new{x=1.0, y=0.4}
f = Vector3:new{x=0.0, y=1.0}
g = Vector3:new{x=0.6, y=1.0}
h = Vector3:new{x=1.0, y=1.0}

-- Lines connecting blocks.
ab = Line:new{p0=a, p1=b} -- horizontal line (lowest level)
cd = Line:new{p0=c, p1=d}; de = Line:new{p0=d, p1=e} -- horizontal lines (mid level)
fg = Line:new{p0=f, p1=g}; gh = Line:new{p0=g, p1=h} -- horizontal lines (top level)
ac = Line:new{p0=a, p1=c}; cf = Line:new{p0=c, p1=f} -- vertical lines (left)
bd = Line:new{p0=b, p1=d}; dg = Line:new{p0=d, p1=g} -- vertical lines (mid)
eh = Line:new{p0=e, p1=h} -- vertical line (right)

-- Define patches (which are parametric surfaces, no discretisation at this point.)
quad0 = CoonsPatch:new{north=cd, east=bd, south=ab, west=ac}
quad1 = CoonsPatch:new{north=fg, east=dg, south=cd, west=cf}
quad2 = CoonsPatch:new{north=gh, east=eh, south=de, west=dg}

-- Define grids. Here's where discretisation is added to a Patch
nx0cells = 12; nx1cells = 8;
ny0cells = 8; ny1cells = 12
grid0 = StructuredGrid:new{psurface=quad0, niv=nx0cells+1, njv=ny0cells+1}
grid1 = StructuredGrid:new{psurface=quad1, niv=nx0cells+1, njv=ny1cells+1}
grid2 = StructuredGrid:new{psurface=quad2, niv=nx1cells+1, njv=ny1cells+1}

-- Lastly, define the blocks.
blk0 = FoamBlock:new{grid=grid0, 
		     bndry_labels={west="w-01", south="w-01", east="w-01"}}
blk1 = FoamBlock:new{grid=grid1,
		     bndry_labels={west="w-01", north="w-00"}}
blk2 = FoamBlock:new{grid=grid2,
		     bndry_labels={south="w-01", east="w-01", north="w-00"}}




