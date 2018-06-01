-- An example showing the usage of foamMesh.
-- Here a grid consisting of two 3-D blocks is generated.
-- The grid consists of a cube and a second block that flares outwards in all 
-- directions. The diagram below shows the frontal view. The rear wall has the 
-- same shape, but point labels are of form a1, b1, ...
--
--             BC=w-01   -f0
--                    --   |
--          d0----e0-      |
--          |  b0 |    b1  | 
--  BC=i-00 |     |        |   BC=o-00  
--          a0----b0-      |
--                    --   |
--            BC=w-00    -c0
-- 
-- Test case to show that linking of all faces works correctly in 3-D
--
-- Authors: IJ
-- Date: 2018-04-01

-- Corners of blocks
a0 = Vector3:new{x=0.0, y=0.0, z=0.0}
b0 = Vector3:new{x=1.0, y=0.0, z=0.0}
c0 = Vector3:new{x=2.0, y=-0.25, z=-0.25}
d0 = Vector3:new{x=0.0, y=1.0, z=0.0}
e0 = Vector3:new{x=1.0, y=1.0, z=0.0}
f0 = Vector3:new{x=2.0, y=1.25, z=-0.25}
a1 = Vector3:new{x=0.0, y=0.0, z=0.5}
b1 = Vector3:new{x=1.0, y=0.0, z=0.5}
c1 = Vector3:new{x=2.0, y=-0.25, z=0.75}
d1 = Vector3:new{x=0.0, y=1.0, z=0.5}
e1 = Vector3:new{x=1.0, y=1.0, z=0.5}
f1 = Vector3:new{x=2.0, y=1.25, z=0.75}

-- Lines connecting blocks.
ab = Line:new{p0=a0, p1=b0}
ac = Line:new{p0=a0, p1=c0}
bd = Line:new{p0=b0, p1=d0}
cd = Line:new{p0=c0, p1=d0}
ac = Line:new{p0=a0, p1=c0}
bd = Line:new{p0=b0, p1=d0}
cd = Line:new{p0=c0, p1=d0}

-- Define patch (which are parametric surfaces, no discretisation at this point.)
--quad0 = CoonsPatch:new{north=cd, east=bd, south=ab, west=ac}

vol0 = TFIVolume:new{vertices={a0,b0,e0,d0,a1,b1,e1,d1}}
vol1 = TFIVolume:new{vertices={b0,c0,f0,e0,b1,c1,f1,e1}}

-- Define 2D grid on patch, clustering can be added if desired
nxcells = 10; nycells = 10; nzcells = 3

grid0 = StructuredGrid:new{pvolume=vol0, niv=nxcells+1, njv=nycells+1, nkv=nzcells+1}
grid1 = StructuredGrid:new{pvolume=vol1, niv=nxcells+1, njv=nycells+1, nkv=nzcells+1}

-- Define OpenFoam block (a "grid" with labels)
blk0 = FoamBlock:new{grid=grid0,
		     bndry_labels={west="i-00", south="w-00",              north="w-01", top="w-02", bottom="w-02"}}
blk1 = FoamBlock:new{grid=grid1,
		     bndry_labels={             south="w-00", east="o-00", north="w-01", top="w-02", bottom="w-02"}}







