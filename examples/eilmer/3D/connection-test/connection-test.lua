-- read the case and find the vertex pairs
file = io.open("case.txt", "r")
faceA = tostring(file:read("*line"))
faceB = tostring(file:read("*line"))
orientation = tonumber(file:read("*number"))
file:close()

vtxPairs = vtxPairs3D[faceA..faceB..orientation]

print("faceA= ", faceA, " faceB= ", faceB, " ornt= ", orientation)

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

-- build cubes (left and right)
cL = cube(Vector3:new{0.0, 0.0, 0.0}, 1.0)
cR = cube(Vector3:new{1.0, 0.0, 0.0}, 1.0)
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

ewPair = {[0]=1, [1]=0, [3]=2, [2]=3, [7]=6, [6]=7, [4]=5, [5]=4}
nsPair = {[3]=0, [0]=3, [2]=1, [1]=2, [6]=5, [5]=6, [7]=4, [4]=7}
tbPair = {[4]=0, [0]=4, [5]=1, [1]=5, [6]=2, [2]=6, [7]=3, [3]=7}

-- find the vertices on the connecting face of the right cube
-- We could condense this next bit of code, but it gets hard
-- to understand later
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

for k,v in pairs(rghtPts) do print(k,v) end

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

rghtVtxs = {}
for cnr,i in pairs(rghtPts) do
   rghtVtxs[i+1] = cR[cnr]
end

lftGrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices=lftVtxs},
			     niv=nn+1, njv=nn+1, nkv=nn+1}
lftBlk = SBlock:new{grid=lftGrid, fillCondition=lftCubeFill}

rghtGrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices=rghtVtxs},
			      niv=nn+1, njv=nn+1, nkv=nn+1}
rghtBlk = SBlock:new{grid=rghtGrid, fillCondition=rghtCubeFill}

identifyBlockConnections()


   
   




