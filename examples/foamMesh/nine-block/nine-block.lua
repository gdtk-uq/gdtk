-- An example showing the usage of foamMesh.
-- Here a grid is generated consisting of 9 blocks arranged in a rectangle.
--
--                BC=w-00
--          d0----d1----d2----d3
--          |  b6 |  b7 |  b8 |
--          |     |     |     |
--          c0----c1----c2----c3
--          |  b3 |  b4 |  b5 |
--  BC=i-00 |     |     |     |   BC=o-00  
--          b0----b1----b2----b3
--          |  b0 |  b1 |  b3 |
--          |     |     |     |
--          a0----a1----a2----a3
--                BC=w-00
-- 
-- Test case to show that linking of all faces works.
--
-- Authors: IJ and RJG
-- Date: 2017-06-29

-- Global settings go first
axisymmetric = false
turbulence_model = "S-A" -- other option is: "k-epsilon"

-- Corners of blocks
a0 = Vector3:new{x=0.00, y=0.0}
b0 = Vector3:new{x=0.33, y=0.0}
c0 = Vector3:new{x=0.66, y=0.0}
d0 = Vector3:new{x=1.00, y=0.0}
a1 = Vector3:new{x=0.00, y=0.33}
b1 = Vector3:new{x=0.33, y=0.33}
c1 = Vector3:new{x=0.66, y=0.33}
d1 = Vector3:new{x=1.00, y=0.33}
a2 = Vector3:new{x=0.00, y=0.66}
b2 = Vector3:new{x=0.33, y=0.66}
c2 = Vector3:new{x=0.66, y=0.66}
d2 = Vector3:new{x=1.00, y=0.66}
a3 = Vector3:new{x=0.00, y=1.0}
b3 = Vector3:new{x=0.33, y=1.0}
c3 = Vector3:new{x=0.66, y=1.0}
d3 = Vector3:new{x=1.00, y=1.0}

-- Define patches (which are parametric surfaces, no discretisation at this point.)
surf = {}
surf[0] = CoonsPatch:new{p00=a0, p10=b0, p01=a1, p11=b1}
surf[1] = CoonsPatch:new{p00=b0, p10=c0, p01=b1, p11=c1}
surf[2] = CoonsPatch:new{p00=c0, p10=d0, p01=c1, p11=d1}
surf[3] = CoonsPatch:new{p00=a1, p10=b1, p01=a2, p11=b2}
surf[4] = CoonsPatch:new{p00=b1, p10=c1, p01=b2, p11=c2}
surf[5] = CoonsPatch:new{p00=c1, p10=d1, p01=c2, p11=d2}
surf[6] = CoonsPatch:new{p00=a2, p10=b2, p01=a3, p11=b3}
surf[7] = CoonsPatch:new{p00=b2, p10=c2, p01=b3, p11=c3}
surf[8] = CoonsPatch:new{p00=c2, p10=d2, p01=c3, p11=d3}

-- Define grids. Here's where discretisation is added to a Patch
nxcells = 5; nycells = 8
grid={}
grid[0] = StructuredGrid:new{psurface=surf[0], niv=nxcells+1, njv=nycells+1}
grid[1] = StructuredGrid:new{psurface=surf[1], niv=nxcells+1, njv=nycells+1}
grid[2] = StructuredGrid:new{psurface=surf[2], niv=nxcells+1, njv=nycells+1}
grid[3] = StructuredGrid:new{psurface=surf[3], niv=nxcells+1, njv=nycells+1}
grid[4] = StructuredGrid:new{psurface=surf[4], niv=nxcells+1, njv=nycells+1}
grid[5] = StructuredGrid:new{psurface=surf[5], niv=nxcells+1, njv=nycells+1}
grid[6] = StructuredGrid:new{psurface=surf[6], niv=nxcells+1, njv=nycells+1}
grid[7] = StructuredGrid:new{psurface=surf[7], niv=nxcells+1, njv=nycells+1}
grid[8] = StructuredGrid:new{psurface=surf[8], niv=nxcells+1, njv=nycells+1}


-- Lastly, define the blocks.
blk={}
blk[0] = FoamBlock:new{grid=grid[0], 
		     bndry_labels={west="i-00", south="w-00"}}
blk[1] = FoamBlock:new{grid=grid[1], 
		     bndry_labels={south="w-00"}}
blk[2] = FoamBlock:new{grid=grid[2], 
		     bndry_labels={south="w-00", east="o-00", }}
blk[3] = FoamBlock:new{grid=grid[3], 
		     bndry_labels={west="i-00"}}
blk[4] = FoamBlock:new{grid=grid[4]}
--                     bndry_labels={} }
blk[5] = FoamBlock:new{grid=grid[5], 
		     bndry_labels={east="o-00"}}
blk[6] = FoamBlock:new{grid=grid[6], 
		     bndry_labels={west="i-00", north="w-00"}}
blk[7] = FoamBlock:new{grid=grid[7], 
		     bndry_labels={north="w-00"}}
blk[8] = FoamBlock:new{grid=grid[8], 
		     bndry_labels={north="w-00", east="o-00", }}


