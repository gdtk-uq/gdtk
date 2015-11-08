-- Author: Rowan J. Gollan
-- Date: 2015-11-08
--[[

This script is a template script that can
be used as part of a testing routine to
exercise every possible connection in 3D
of right-handed blocks in a structured grid
code. The script is specialised by reading
in a particular connection to test in a
"case.txt" file.

This connection test script requires some explanation,
mostly so I'll remember how it works in the future.
In this preamble to the file, I'll document main idea
of how this test works. The comments that appear
later in the file will attempt to described the
implementatation.

Test case overview
------------------
The test involves two blocks placed side-by-side
with a face in common providing a connection.
The initial gas state in both blocks varies
spatially such that there is a non-uniform
distribution of pressure and density across
the connecting face. This non-uniformity has
no symmetry across the connceting plane.
This is important for the testing
(to be explained shortly).
The arrangement of the blocks is shown here
(except that face [BCGF] of the left block
physically coincides with face [ADHE] of
the right block):

      H------------G          H------------G
     /|           /|         /|           /|
    / |          / |        / |          / |
   D--+---------C  |       D--+---------C  |
   |  | LFT BLK |  |       |  |RGHT BLK |  |
   |  |         |  |       |  |         |  |
   |  E---------+--F       |  E---------+--F
   | /          | /        | /          | /
   |/           |/         |/           |/
   A------------B          A------------B

The test is run by letting the flow solver
iterate the state of the system forward
in time by one timestep. The non-uniform
nature of the flow properties gives a
unique flow field after one time step.
We build a reference case of how that
flow state should look after one step.
We will then swap the arrangement of the block
connections and test that the flow state
after one step matches our reference case.
If if doesn't, it will be a strong indication
that a copying error has occurred during
the exchange routines for that particular
block connection type.

Implementation of the test
--------------------------
The "trick" to this test is that the left
and right blocks are always located in the
same physical position, no matter which connection
is being tested. By keeping the blocks in the
same physical position, it means that the expected
state of the system after one step is always the
same in physical space. This simplifies the
extraction of the data. We can use the "--probe"
option of the post-processor to probe data at
the same physical location, and then check that
the output matches our reference case.

Most of the complexity in the script below is about
rotating logical orientations of the blocks so that
we can test all of the connections. When reading the
script, keep in mind that the alphabetic labels
always refer to the physical locations of block
corners as shown in the figure above. The numeric labels
refer to the logical labels used internally in the code
to designate the corners of 3D block.

--]]

-- Read the case 
file = io.open("case.txt", "r")
faceA = tostring(file:read("*line"))
faceB = tostring(file:read("*line"))
orientation = tonumber(file:read("*number"))
file:close()

-- Now find the connecting vertex pairs by looking
-- up the tabled 'vtxPairs3D' proviede in 'blk_conn.lua'
vtxPairs = vtxPairs3D[faceA..faceB..orientation]

-- Proceed with some general simulation setup
config.title = "3D Connection tester"
config.dimensions = 3

setGasModel('ideal-air.lua')
config.max_time = 0.6e-3
config.max_step = 1
config.dt_init = 1.0e-6
config.flux_calculator = "ausmdv"

nn = 4
p_base = 1.0e5
T_base = 348.3

-- Set up the physical aspects of the test
-- The next two functions provide the initial state of
-- the gas in the LFT and RGHT cubes. This initial state
-- varies in 3D space.
function lftCubeFill(x, y, z)
   p_inc_factor = x * (5*y^2) * (10*z^3)
   p = p_base*(1.0 + p_inc_factor)
   return FlowState:new{p=p, T=T_base}
end

function rghtCubeFill(x, y, z)
   p_inc_factor = x * (5*y^2) * (10*z^3)
   p = p_base*(0.8 + p_inc_factor)
   return FlowState:new{p=p, T=T_base}
end

-- Set up a function to give vertex positions of a cube
--[[
          H------------G
         /|           /|
        / |          / |
       D--+---------C  |
       |  |         |  |
       |  |         |  |
       |  E---------+--F
       | /          | /
       |/           |/
       A------------B

--]]
function cube(origin, length)
   x0 = origin:x(); y0 = origin:y(); z0 = origin:z()
   A = Vector3:new{x0, y0, z0}
   B = Vector3:new{x0+length, y0, z0}
   C = Vector3:new{x0+length, y0+length, z0}
   D = Vector3:new{x0, y0+length, z0}
   E = Vector3:new{x0, y0, z0+length}
   F = Vector3:new{x0+length, y0, z0+length}
   G = Vector3:new{x0+length, y0+length, z0+length}
   H = Vector3:new{x0, y0+length, z0+length}
   return {A=A, B=B, C=C, D=D, E=E, F=F, G=G, H=H}
end

-- Build cubes (left and right)
-- We store the cube corners based on their
-- physical location (alphabetic labels) 
-- in 'cL' (corners of left cube)
-- and 'cR' (corners of right cube)
cL = cube(Vector3:new{0.0, 0.0, 0.0}, 1.0)
cR = cube(Vector3:new{1.0, 0.0, 0.0}, 1.0)

-- Now we consider the logical arrangement of the corners (vertices) 
-- depending on which face of the left block is being used in the
-- connection faceA.
-- lftVtxs: stores the corners in the correct order (0, 1, 2, ...)
--          such that the correct face coincides with face [BCGF] of the
--          left block
-- lftFace: store the index locations associated with the corners [BCGF]
local lftVtxs
local lftFace
if faceA == north then
   lftVtxs = {cL.D, cL.A, cL.B, cL.C, cL.H, cL.E, cL.F, cL.G}
   lftFace = {B=2, C=3, F=6, G=7}
elseif faceA == south then
   lftVtxs = {cL.B, cL.C, cL.D, cL.A, cL.F, cL.G, cL.H, cL.E}
   lftFace = {B=0, C=1, F=4, G=5}
elseif faceA == east then
   lftVtxs = {cL.A, cL.B, cL.C, cL.D, cL.E, cL.F, cL.G, cL.H}
   lftFace = {B=1, C=2, F=5, G=6}
elseif faceA == west then
   lftVtxs = {cL.F, cL.E, cL.H, cL.G, cL.B, cL.A, cL.D, cL.C}
   lftFace = {B=4, C=7, F=0, G=3}
elseif faceA == top then
   lftVtxs = {cL.E, cL.A, cL.D, cL.H, cL.F, cL.B, cL.C, cL.G}
   lftFace = {B=5, C=6, F=4, G=7}
elseif faceA == bottom then
   lftVtxs = {cL.B, cL.F, cL.G, cL.C, cL.A, cL.E, cL.H, cL.D}
   lftFace = {B=0, C=3, F=1, G=2}
end

-- On a given face, there is always a set of corners on the
-- opposite face across the block. These sets of opposites
-- occur in the pairs: east/west, north/south, and top/bottom.
-- We store these connecting sets for later use. 

ewPair = {[0]=1, [1]=0, [3]=2, [2]=3, [7]=6, [6]=7, [4]=5, [5]=4}
nsPair = {[3]=0, [0]=3, [2]=1, [1]=2, [6]=5, [5]=6, [7]=4, [4]=7}
tbPair = {[4]=0, [0]=4, [5]=1, [1]=5, [6]=2, [2]=6, [7]=3, [3]=7}

-- Our left face is made of corners: B, C, F and G.
-- From the figure above, we know the connections to the right
-- face:
--       B (of left) +--+ A (of right)
--       C (of left) +--+ D (of right)
--       F (of left) +--+ E (of right)
--       G (of left) +--+ H (of right)
-- Now we already have a list of all pairs of vertices involved
-- in the connection. So we'll look through the list of vertices
-- and find the numeric value that matches with corner B of
-- the left face. When we've found that, we know that the matching
-- vertex is corner B of the right face. Continuing in that manner
-- for corners C, F and G of the left face, we can find the numeric
-- labels for the matching corners D, E and H on the right face.
-- (It would be possible to condense the code that does this search
--  but it gets very obsucre, very quickly.)
local rghtPts = {}
pB = lftFace.B
for _,vp in ipairs(vtxPairs) do
   if pB == vp[1] then rghtPts.A = vp[2] end
end
pC = lftFace.C
for _,vp in ipairs(vtxPairs) do
   if pC == vp[1] then rghtPts.D = vp[2] end
end
pF = lftFace.F
for _,vp in ipairs(vtxPairs) do
   if pF == vp[1] then rghtPts.E = vp[2] end
end
pG = lftFace.G
for _,vp in ipairs(vtxPairs) do
   if pG == vp[1] then rghtPts.H = vp[2] end
end

-- Now we get to use the connecting sets we used above.
-- Since we know the numeric points at locations 
-- A, E, H and D on the right face, AND we know what
-- type of connecting face it is, then we can instantly
-- fill out the opposing corners of the right cube
-- (B, F, G, C) using the connection mapping. 
if faceB == north or faceB == south then
   rghtPts.B = nsPair[rghtPts.A]
   rghtPts.F = nsPair[rghtPts.E]
   rghtPts.G = nsPair[rghtPts.H]
   rghtPts.C = nsPair[rghtPts.D]
elseif faceB == east or faceB == west then
   rghtPts.B = ewPair[rghtPts.A]
   rghtPts.F = ewPair[rghtPts.E]
   rghtPts.G = ewPair[rghtPts.H]
   rghtPts.C = ewPair[rghtPts.D]
elseif faceB == top or faceB == bottom then
   rghtPts.B = tbPair[rghtPts.A]
   rghtPts.F = tbPair[rghtPts.E]
   rghtPts.G = tbPair[rghtPts.H]
   rghtPts.C = tbPair[rghtPts.D]
end

-- We need to list the right points in order (0, 1, 2, etc. but
-- counting from 1 [hence i+1]).
-- So we basically build a reverse map of what we've collected
-- in rghtPts.
rghtVtxs = {}
for cnr,i in pairs(rghtPts) do
   rghtVtxs[i+1] = cR[cnr]
end

-- In the last step, we build the left and right blocks from
-- our collections of vertices.
lftGrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices=lftVtxs},
			     niv=nn+1, njv=nn+1, nkv=nn+1}
lftBlk = SBlock:new{grid=lftGrid, fillCondition=lftCubeFill}

rghtGrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices=rghtVtxs},
			      niv=nn+1, njv=nn+1, nkv=nn+1}
rghtBlk = SBlock:new{grid=rghtGrid, fillCondition=rghtCubeFill}

identifyBlockConnections()


   
   




