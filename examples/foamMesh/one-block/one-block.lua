-- An example showing the usage of foamMesh.
-- Here an axi-symmetric mesh consisting of a single axis is generated 
--
--
--
--             BC=w-01
--           c--------d
--           |        |
--  BC=i-00  |   b0   |  BC=o-00
--           |        |
--      -----a--------b----- Axis of Rotation
--             BC=w-00
--
-- Authors: IJ and RJG
-- Date: 2017-06-29

axisymmetric = true
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

-- Define OpenFoam block (a "grid" with labels)
blk0 = FoamBlock:new{grid=grid0,
		     bndry_labels={west="i-00", south="w-00", east="o-00", north="w-01"}}








