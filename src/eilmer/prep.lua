-- prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- 
print("Loading prep.lua...")

-- As an experiment to see if we can avoid SEGFAULTs, we tried turning off
-- the garbage collector in the Lua interpreter.  It made no difference.
-- Disabling the garbage collector on the D side did make a difference.
-- collectgarbage("stop")

-- Storage for later definitions of Block objects
blocks = {}

-- Storgage for later definitions of SolidBlock objects
solidBlocks = {}

-- Symbolic names for identifying boundaries
north = "north"; NORTH = "north"
east = "east"; EAST = "east"
south = "south"; SOUTH = "south"
west = "west"; WEST = "west"
top = "top"; TOP = "top"
bottom = "bottom"; BOTTOM = "bottom"

function faceList(dimensions)
   local myList = {north, east, south, west}
   if dimensions == 3 then 
      table.insert(myList, top)
      table.insert(myList, bottom)
   end
   return myList
end

-- -----------------------------------------------------------------------

-- Connections between 2D blocks, described as vertex-pairs.
-- When specifying a connection between 2D blocks,
-- we actually specify the names of the connecting faces.
-- The connection orientation is always 0.
tabulatedData = {
   {{{3,0},{2,1}}, {north, south}},
   {{{3,3},{2,0}}, {north, west}},
   {{{3,2},{2,3}}, {north, north}},
   {{{3,1},{2,2}}, {north, east}},
   {{{0,0},{3,1}}, {west,  south}},
   {{{0,3},{3,0}}, {west,  west}},
   {{{0,2},{3,3}}, {west,  north}},
   {{{0,1},{3,2}}, {west,  east}},
   {{{1,0},{0,1}}, {south, south}},
   {{{1,3},{0,0}}, {south, west}},
   {{{1,2},{0,3}}, {south, north}},
   {{{1,1},{0,2}}, {south, east}},
   {{{2,0},{1,1}}, {east,  south}},
   {{{2,3},{1,0}}, {east,  west}},
   {{{2,2},{1,3}}, {east,  north}},
   {{{2,1},{1,2}}, {east,  east}},
} -- end tabulatedData

-- Since there are a number of permutations for each valid vertex-vertex
-- connection set, define a function that can help sort the pairs into
-- a standard order.
function cmpVtxPair(a, b)
   return (a[1] < b[1]) or (a[1] == b[1] and a[2] < b[2])
end

-- Table unpacking function from Section 5.1 in Programming in Lua
function unpack(t, i, n)
   i = i or 1
   n = n or #t
   if i <= n then
      return t[i], unpack(t, i+1, n)
   end
end

-- A couple of convenience functions
function tostringVtxPair(t)
   return string.format("{%d,%d}", t[1], t[2])
end
function tostringVtxPairList(t)
   str = "{" .. tostringVtxPair(t[1])
   for i=2, #t do
      str = str .. "," .. tostringVtxPair(t[i])
   end
   str = str .. "}"
   return str
end
function tostringConnection(t)
   local faceA, faceB, orientation = unpack(t)
   orientation = orientation or 0
   return string.format("{%s, %s, %d}", faceA, faceB, orientation)
end

connections2D = {}
vtxPairs2D = {}
for _,v in ipairs(tabulatedData) do
   local vtxPairs, connection = unpack(v)
   table.sort(vtxPairs, cmpVtxPair)
   connections2D[vtxPairs] = connection
   local this_face, other_face = unpack(connection)
   vtxPairs2D[{this_face, other_face}] = vtxPairs
end

-- Connections between 3D blocks, described as sets of vertex pairs.
-- When specifying a connection, we specify the set of paired vertices.
-- The following table provides a map from (A,B) inter-block connections
-- that are defined by (A-vtx,B-vtx) pairs to inter-block connections that
-- are defined by (A-face, B-face, orientation, axis-map) tuples.
-- The orientation is used in the C-code to decide where to transfer data
-- at the block boundaries and the axis-map is used to interpret GridPro
-- connectivity data.  The i,j,k axes of *this* block are aligned with the
-- specified axes of the *other* block. 
tabulatedData = {
   {{{3,2},{7,6},{6,7},{2,3}}, {north, north, 0, '-i-j+k'}},
   {{{3,3},{7,2},{6,6},{2,7}}, {north, north, 1, '+k-j+i'}},
   {{{3,7},{7,3},{6,2},{2,6}}, {north, north, 2, '+i-j-k'}},
   {{{3,6},{7,7},{6,3},{2,2}}, {north, north, 3, '-k-j-i'}},

   {{{3,0},{7,4},{6,5},{2,1}}, {north, south, 0, '+i+j+k'}},
   {{{3,1},{7,0},{6,4},{2,5}}, {north, south, 1, '+k+j-i'}},
   {{{3,5},{7,1},{6,0},{2,4}}, {north, south, 2, '-i+j-k'}},
   {{{3,4},{7,5},{6,1},{2,0}}, {north, south, 3, '-k+j+i'}},

   {{{3,1},{7,5},{6,6},{2,2}}, {north, east, 0, '+j-i+k'}},
   {{{3,2},{7,1},{6,5},{2,6}}, {north, east, 1, '+k-i-j'}},
   {{{3,6},{7,2},{6,1},{2,5}}, {north, east, 2, '-j-i-k'}},
   {{{3,5},{7,6},{6,2},{2,1}}, {north, east, 3, '-k-i+j'}},

   {{{3,3},{7,7},{6,4},{2,0}}, {north, west, 0, '-j+i+k'}},
   {{{3,0},{7,3},{6,7},{2,4}}, {north, west, 1, '+k+i+j'}},
   {{{3,4},{7,0},{6,3},{2,7}}, {north, west, 2, '+j+i-k'}},
   {{{3,7},{7,4},{6,0},{2,3}}, {north, west, 3, '-k+i-j'}},

   {{{3,4},{7,7},{6,6},{2,5}}, {north, top, 0, '+i-k+j'}},
   {{{3,5},{7,4},{6,7},{2,6}}, {north, top, 1, '+j-k-i'}},
   {{{3,6},{7,5},{6,4},{2,7}}, {north, top, 2, '-i-k-j'}},
   {{{3,7},{7,6},{6,5},{2,4}}, {north, top, 3, '-j-k+i'}},

   {{{3,1},{7,2},{6,3},{2,0}}, {north, bottom, 0, '-i+k+j'}},
   {{{3,0},{7,1},{6,2},{2,3}}, {north, bottom, 1, '+j+k+i'}},
   {{{3,3},{7,0},{6,1},{2,2}}, {north, bottom, 2, '+i+k-j'}},
   {{{3,2},{7,3},{6,0},{2,1}}, {north, bottom, 3, '-j+k-i'}},

   {{{1,2},{5,6},{4,7},{0,3}}, {south, north, 0, '+i+j+k'}},
   {{{1,3},{5,2},{4,6},{0,7}}, {south, north, 1, '-k+j+i'}},
   {{{1,7},{5,3},{4,2},{0,6}}, {south, north, 2, '-i+j-k'}},
   {{{1,6},{5,7},{4,3},{0,2}}, {south, north, 3, '+k+j-i'}},

   {{{1,0},{5,4},{4,5},{0,1}}, {south, south, 0, '-i-j+k'}},
   {{{1,1},{5,0},{4,4},{0,5}}, {south, south, 1, '-k-j-i'}},
   {{{1,5},{5,1},{4,0},{0,4}}, {south, south, 2, '+i-j-k'}},
   {{{1,4},{5,5},{4,1},{0,0}}, {south, south, 3, '+k-j+i'}},

   {{{1,1},{5,5},{4,6},{0,2}}, {south, east, 0, '-j+i+k'}},
   {{{1,2},{5,1},{4,5},{0,6}}, {south, east, 1, '-k+i-j'}},
   {{{1,6},{5,2},{4,1},{0,5}}, {south, east, 2, '+j+i-k'}},
   {{{1,5},{5,6},{4,2},{0,1}}, {south, east, 3, '+k+i+j'}},

   {{{1,3},{5,7},{4,4},{0,0}}, {south, west, 0, '+j-i+k'}},
   {{{1,0},{5,3},{4,7},{0,4}}, {south, west, 1, '-k-i+j'}},
   {{{1,4},{5,0},{4,3},{0,7}}, {south, west, 2, '-j-i-k'}},
   {{{1,7},{5,4},{4,0},{0,3}}, {south, west, 3, '+k-i-j'}},

   {{{1,4},{5,7},{4,6},{0,5}}, {south, top, 0, '-i+k+j'}},
   {{{1,5},{5,4},{4,7},{0,6}}, {south, top, 1, '-j+k-i'}},
   {{{1,6},{5,5},{4,4},{0,7}}, {south, top, 2, '+i+k-j'}},
   {{{1,7},{5,6},{4,5},{0,4}}, {south, top, 3, '+j+k+i'}},

   {{{1,1},{5,2},{4,3},{0,0}}, {south, bottom, 0, '+i-k+j'}},
   {{{1,0},{5,1},{4,2},{0,3}}, {south, bottom, 1, '-j-k+i'}},
   {{{1,3},{5,0},{4,1},{0,2}}, {south, bottom, 2, '-i-k-j'}},
   {{{1,2},{5,3},{4,0},{0,1}}, {south, bottom , 3, '+j-k-i'}},

   {{{2,2},{6,6},{5,7},{1,3}}, {east, north, 0, '-j+i+k'}},
   {{{2,3},{6,2},{5,6},{1,7}}, {east, north, 1, '-j-k+i'}},
   {{{2,7},{6,3},{5,2},{1,6}}, {east, north, 2, '-j-i-k'}},
   {{{2,6},{6,7},{5,3},{1,2}}, {east, north, 3, '-j+k-i'}},

   {{{2,0},{6,4},{5,5},{1,1}}, {east, south, 0, '+j-i+k'}},
   {{{2,1},{6,0},{5,4},{1,5}}, {east, south, 1, '+j-k-i'}},
   {{{2,5},{6,1},{5,0},{1,4}}, {east, south, 2, '+j+i-k'}},
   {{{2,4},{6,5},{5,1},{1,0}}, {east, south, 3, '+j+k+i'}},

   {{{2,1},{6,5},{5,6},{1,2}}, {east, east, 0, '-i-j+k'}},
   {{{2,2},{6,1},{5,5},{1,6}}, {east, east, 1, '-i-k-j'}},
   {{{2,6},{6,2},{5,1},{1,5}}, {east, east, 2, '-i+j-k'}},
   {{{2,5},{6,6},{5,2},{1,1}}, {east, east, 3, '-i+k+j'}},

   {{{2,3},{6,7},{5,4},{1,0}}, {east, west, 0, '+i+j+k'}},
   {{{2,0},{6,3},{5,7},{1,4}}, {east, west, 1, '+i-k+j'}},
   {{{2,4},{6,0},{5,3},{1,7}}, {east, west, 2, '+i-j-k'}},
   {{{2,7},{6,4},{5,0},{1,3}}, {east, west, 3, '+i+k-j'}},

   {{{2,4},{6,7},{5,6},{1,5}}, {east, top, 0, '-k-i+j'}},
   {{{2,5},{6,4},{5,7},{1,6}}, {east, top, 1, '-k-j-i'}},
   {{{2,6},{6,5},{5,4},{1,7}}, {east, top, 2, '-k+i-j'}},
   {{{2,7},{6,6},{5,5},{1,4}}, {east, top, 3, '-k+j+i'}},

   {{{2,1},{6,2},{5,3},{1,0}}, {east, bottom, 0, '+k+i+j'}},
   {{{2,0},{6,1},{5,2},{1,3}}, {east, bottom, 1, '+k-j+i'}},
   {{{2,3},{6,0},{5,1},{1,2}}, {east, bottom, 2, '+k-i-j'}},
   {{{2,2},{6,3},{5,0},{1,1}}, {east, bottom, 3, '+k+j-i'}},

   {{{0,2},{4,6},{7,7},{3,3}}, {west, north, 0, '+j-i+k'}},
   {{{0,3},{4,2},{7,6},{3,7}}, {west, north, 1, '+j+k+i'}},
   {{{0,7},{4,3},{7,2},{3,6}}, {west, north, 2, '+j+i-k'}},
   {{{0,6},{4,7},{7,3},{3,2}}, {west, north, 3, '+j-k-i'}},

   {{{0,0},{4,4},{7,5},{3,1}}, {west, south, 0, '-j+i+k'}},
   {{{0,1},{4,0},{7,4},{3,5}}, {west, south, 1, '-j+k-i'}},
   {{{0,5},{4,1},{7,0},{3,4}}, {west, south, 2, '-j-i-k'}},
   {{{0,4},{4,5},{7,1},{3,0}}, {west, south, 3, '-j-k+i'}},

   {{{0,1},{4,5},{7,6},{3,2}}, {west, east, 0, '+i+j+k'}},
   {{{0,2},{4,1},{7,5},{3,6}}, {west, east, 1, '+i+k-j'}},
   {{{0,6},{4,2},{7,1},{3,5}}, {west, east, 2, '+i-j-k'}},
   {{{0,5},{4,6},{7,2},{3,1}}, {west, east, 3, '+i-k+j'}},

   {{{0,3},{4,7},{7,4},{3,0}}, {west, west, 0, '-i-j+k'}},
   {{{0,0},{4,3},{7,7},{3,4}}, {west, west, 1, '-i+k+j'}},
   {{{0,4},{4,0},{7,3},{3,7}}, {west, west, 2, '-i+j-k'}},
   {{{0,7},{4,4},{7,0},{3,3}}, {west, west, 3, '-i-k-j'}},

   {{{0,4},{4,7},{7,6},{3,5}}, {west, top, 0, '+k+i+j'}},
   {{{0,5},{4,4},{7,7},{3,6}}, {west, top, 1, '+k+j-i'}},
   {{{0,6},{4,5},{7,4},{3,7}}, {west, top, 2, '+k-i-j'}},
   {{{0,7},{4,6},{7,5},{3,4}}, {west, top, 3, '+k-j+i'}},

   {{{0,1},{4,2},{7,3},{3,0}}, {west, bottom, 0, '-k-i+j'}},
   {{{0,0},{4,1},{7,2},{3,3}}, {west, bottom, 1, '-k+j+i'}},
   {{{0,3},{4,0},{7,1},{3,2}}, {west, bottom, 2, '-k+i-j'}},
   {{{0,2},{4,3},{7,0},{3,1}}, {west, bottom, 3, '-k-j-i'}},

   {{{5,2},{6,6},{7,7},{4,3}}, {top, north, 0, '+i+k-j'}},
   {{{5,3},{6,2},{7,6},{4,7}}, {top, north, 1, '-k+i-j'}},
   {{{5,7},{6,3},{7,2},{4,6}}, {top, north, 2, '-i-k-j'}},
   {{{5,6},{6,7},{7,3},{4,2}}, {top, north, 3, '+k-i-j'}},

   {{{5,0},{6,4},{7,5},{4,1}}, {top, south, 0, '-i+k+j'}},
   {{{5,1},{6,0},{7,4},{4,5}}, {top, south, 1, '-k-i+j'}},
   {{{5,5},{6,1},{7,0},{4,4}}, {top, south, 2, '+i-k+j'}},
   {{{5,4},{6,5},{7,1},{4,0}}, {top, south, 3, '+k+i+j'}},

   {{{5,1},{6,5},{7,6},{4,2}}, {top, east, 0, '-j+k-i'}},
   {{{5,2},{6,1},{7,5},{4,6}}, {top, east, 1, '-k-j-i'}},
   {{{5,6},{6,2},{7,1},{4,5}}, {top, east, 2, '+j-k-i'}},
   {{{5,5},{6,6},{7,2},{4,1}}, {top, east, 3, '+k+j-i'}},

   {{{5,3},{6,7},{7,4},{4,0}}, {top, west, 0, '+j+k+i'}},
   {{{5,0},{6,3},{7,7},{4,4}}, {top, west, 1, '-k+j+i'}},
   {{{5,4},{6,0},{7,3},{4,7}}, {top, west, 2, '-j-k+i'}},
   {{{5,7},{6,4},{7,0},{4,3}}, {top, west, 3, '+k-j+i'}},

   {{{5,4},{6,7},{7,6},{4,5}}, {top, top, 0, '-i+j-k'}},
   {{{5,5},{6,4},{7,7},{4,6}}, {top, top, 1, '-j-i-k'}},
   {{{5,6},{6,5},{7,4},{4,7}}, {top, top, 2, '+i-j-k'}},
   {{{5,7},{6,6},{7,5},{4,4}}, {top, top, 3, '+j+i-k'}},

   {{{5,1},{6,2},{7,3},{4,0}}, {top, bottom, 0, '+i+j+k'}},
   {{{5,0},{6,1},{7,2},{4,3}}, {top, bottom, 1, '-j+i+k'}},
   {{{5,3},{6,0},{7,1},{4,2}}, {top, bottom, 2, '-i-j+k'}},
   {{{5,2},{6,3},{7,0},{4,1}}, {top, bottom, 3, '+j-i+k'}},

   {{{0,2},{3,6},{2,7},{1,3}}, {bottom, north, 0, '-i+k+j'}},
   {{{0,3},{3,2},{2,6},{1,7}}, {bottom, north, 1, '+k+i+j'}},
   {{{0,7},{3,3},{2,2},{1,6}}, {bottom, north, 2, '+i-k+j'}},
   {{{0,6},{3,7},{2,3},{1,2}}, {bottom, north, 3, '-k-i+j'}},

   {{{0,0},{3,4},{2,5},{1,1}}, {bottom, south, 0, '+i+k-j'}},
   {{{0,1},{3,0},{2,4},{1,5}}, {bottom, south, 1, '+k-i-j'}},
   {{{0,5},{3,1},{2,0},{1,4}}, {bottom, south, 2, '-i-k-j'}},
   {{{0,4},{3,5},{2,1},{1,0}}, {bottom, south, 3, '-k+i-j'}},

   {{{0,1},{3,5},{2,6},{1,2}}, {bottom, east, 0, '+j+k+i'}},
   {{{0,2},{3,1},{2,5},{1,6}}, {bottom, east, 1, '+k-j+i'}},
   {{{0,6},{3,2},{2,1},{1,5}}, {bottom, east, 2, '-j-k+i'}},
   {{{0,5},{3,6},{2,2},{1,1}}, {bottom, east, 3, '-k+j+i'}},

   {{{0,3},{3,7},{2,4},{1,0}}, {bottom, west, 0, '-j+k-i'}},
   {{{0,0},{3,3},{2,7},{1,4}}, {bottom, west, 1, '+k+j-i'}},
   {{{0,4},{3,0},{2,3},{1,7}}, {bottom, west, 2, '+j-k-i'}},
   {{{0,7},{3,4},{2,0},{1,3}}, {bottom, west, 3, '-k-j-i'}},

   {{{0,4},{3,7},{2,6},{1,5}}, {bottom, top, 0, '+i+j+k'}},
   {{{0,5},{3,4},{2,7},{1,6}}, {bottom, top, 1, '+j-i+k'}},
   {{{0,6},{3,5},{2,4},{1,7}}, {bottom, top, 2, '-i-j+k'}},
   {{{0,7},{3,6},{2,5},{1,4}}, {bottom, top, 3, '-j+i+k'}},

   {{{0,1},{3,2},{2,3},{1,0}}, {bottom, bottom, 0, '-i+j-k'}},
   {{{0,0},{3,1},{2,2},{1,3}}, {bottom, bottom, 1, '+j+i-k'}},
   {{{0,3},{3,0},{2,1},{1,2}}, {bottom, bottom, 2, '+i-j-k'}},
   {{{0,2},{3,3},{2,0},{1,1}}, {bottom, bottom, 3, '-j-i-k'}}
} -- end of tabulatedData

connections3D = {}
vtxPairs3D = {}
-- When reading GridPro block connectivity file, 
-- we need to look up Eilmer notation for connection orientations.
eilmer_orientation = {}
for _,v in ipairs(tabulatedData) do
   local vtxPairs, connection = unpack(v)
   table.sort(vtxPairs, cmpVtxPair)
   connections3D[vtxPairs] = connection
   local this_face, other_face, orientation, axis_map = unpack(connection)
   vtxPairs3D[{this_face, other_face, orientation}] = vtxPairs
   eilmer_orientation[{this_face, other_face, axis_map}] = orientation
end

function to_eilmer_axis_map(gridpro_ijk)
   -- Convert from GridPro axis_map string to Eilmer3 axis_map string.
   -- From GridPro manual, Section 7.3.2 Connectivity Information.
   -- Example, 123 --> '+i+j+k'
   local axis_map = {[0]='xx', [1]='+i', [2]='+j', [3]='+k',
		     [4]='-i', [5]='-j', [6]='-k'}
   if type(gridpro_ijk) == "number" then
      gridpro_ijk = string.format("%03d", gridpro_ijk)
   end
   if type(gridpro_ijk) ~= "string" then
      print("Expected a string or integer of three digits but got:", gridpro_ijk)
      os.exit(-1)
   end
   eilmer_ijk = axis_map[tonumber(string.sub(gridpro_ijk, 1, 1))] ..
      axis_axis_map[tonumber(string.sub(gridpro_ijk, 2, 2))] ..
      map[tonumber(string.sub(gridpro_ijk, 3, 3))]
   return other_ijk
end

-- -----------------------------------------------------------------------
-- Classes for constructing boundary conditions.
-- Each boundary condition is composed of lists of actions to do
-- at specific points in the superloop of the main simulation code.

-- For the classes below, we just follow the prototype pattern
-- given in Ierusalimchy's book "Programming in Lua"

-- Base class and subclasses for GhostCellEffect
GhostCellEffect = {
   type = ""
}
function GhostCellEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

InternalCopyThenReflect = GhostCellEffect:new()
InternalCopyThenReflect.type = "internal_copy_then_reflect"
function InternalCopyThenReflect:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

FlowStateCopy = GhostCellEffect:new{flowCondition=nil}
FlowStateCopy.type = "flowstate_copy"
function FlowStateCopy:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "flowstate": %s', self.flowCondition:toJSONString())
   str = str .. '}'
   return str
end

ExtrapolateCopy = GhostCellEffect:new{xOrder=0}
ExtrapolateCopy.type = "extrapolate_copy"
function ExtrapolateCopy:tojson()
   local str = string.format('          {"type": "%s", "x_order": %d}', self.type, self.xOrder)
   return str
end

FixedPT = GhostCellEffect:new{p_out=1.0e5, T_out=300.0}
FixedPT.type = "fixed_pressure_temperature"
function FixedPT:tojson()
   local str = string.format('          {"type": "%s", "p_out": %f, "T_out": %f}',
			     self.type, self.p_out, self.T_out)
   return str
end

FullFaceExchangeCopy = GhostCellEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
FullFaceExchangeCopy.type = "full_face_exchange_copy"
function FullFaceExchangeCopy:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

UserDefinedGhostCell = GhostCellEffect:new{fileName='user-defined-bc.lua'}
UserDefinedGhostCell.type = "user_defined"
function UserDefinedGhostCell:tojson()
   local str = string.format('         {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s"', self.fileName)
   str = str .. '}'
   return str
end

-- Base class and subclasses for BoundaryInterfaceEffect
BoundaryInterfaceEffect = {
   type = ""
}
function BoundaryInterfaceEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

CopyCellData = BoundaryInterfaceEffect:new()
CopyCellData.type = "copy_cell_data"
function CopyCellData:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

ZeroVelocity = BoundaryInterfaceEffect:new()
ZeroVelocity.type = "zero_velocity"
function ZeroVelocity:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

FixedT = BoundaryInterfaceEffect:new{Twall=nil}
FixedT.type = "fixed_temperature"
function FixedT:tojson()
   local str = string.format('          {"type": "%s",', self.type)
   str = str .. string.format(' "Twall": %f', self.Twall)
   str = str .. '}'
   return str
end

UpdateThermoTransCoeffs = BoundaryInterfaceEffect:new()
UpdateThermoTransCoeffs.type = "update_thermo_trans_coeffs"
function UpdateThermoTransCoeffs:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

WallKOmega = BoundaryInterfaceEffect:new()
WallKOmega.type = "wall_k_omega"
function WallKOmega:tojson()
   local str = string.format('          {"type" : "%s"}', self.type)
   return str
end

TemperatureFromGasSolidInterface = BoundaryInterfaceEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
TemperatureFromGasSolidInterface.type = "temperature_from_gas_solid_interface"
function TemperatureFromGasSolidInterface:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

UserDefinedInterface = BoundaryInterfaceEffect:new{fileName='user-defined-bc.lua'}
UserDefinedInterface.type = "user_defined"
function UserDefinedInterface:tojson()
   local str = string.format('         {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s"', self.fileName)
   str = str .. '}'
   return str
end

-- Base class and subclasses for BoundaryFluxEffect
BoundaryFluxEffect = {
   type = ""
}
function BoundaryFluxEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

EnergyFluxFromAdjacentSolid = BoundaryFluxEffect:new{otherBlock=nil, otherFace=nil, orientation=-1}
EnergyFluxFromAdjacentSolid.type = "energy_flux_from_adjacent_solid"
function EnergyFluxFromAdjacentSolid:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

-- Class for BoundaryCondition

BoundaryCondition = {
   label = "",
   myType = "",
   preReconAction = {},
   preSpatialDerivAction = {},
   postDiffFluxAction = {}
}
function BoundaryCondition:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end
function BoundaryCondition:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", \n', self.label)
   str = str .. '        "pre_recon_action": [\n'
   for i,effect in ipairs(self.preReconAction) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preReconAction then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "pre_spatial_deriv_action": [\n'
   for i,effect in ipairs(self.preSpatialDerivAction) do
      str = str .. effect:tojson()
      if i ~= #self.preSpatialDerivAction then str = str .. "," end
   end
   str = str .. '\n        ],\n'
   str = str .. '        "post_diff_flux_action": [\n'
   for i,effect in ipairs(self.postDiffFluxAction) do
      str = str .. effect:tojson()
      if i ~= #self.postDiffFluxAction then str = str .. "," end
   end
   str = str .. '\n        ]\n'
   str = str .. '    }'
   return str
end

SlipWallBC = BoundaryCondition:new()
SlipWallBC.myType = "SlipWall"
function SlipWallBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

FixedTWallBC = BoundaryCondition:new()
FixedTWallBC.myType = "FixedTWall"
function FixedTWallBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(), ZeroVelocity:new(),
			       FixedT:new{Twall=o.Twall},
			       UpdateThermoTransCoeffs:new(),
			       WallKOmega:new() }
   return o
end

AdiabaticWallBC = BoundaryCondition:new()
AdiabaticWallBC.myType = "FixedTWall"
function AdiabaticWallBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(), ZeroVelocity:new(),
			       WallKOmega:new() }
   return o
end

SupInBC = BoundaryCondition:new()
SupInBC.myType = "SupIn"
function SupInBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { FlowStateCopy:new{flowCondition=o.flowCondition} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

ExtrapolateOutBC = BoundaryCondition:new()
ExtrapolateOutBC.myType = "ExtrapolateOut"
function ExtrapolateOutBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

FixedPTOutBC = BoundaryCondition:new()
FixedPTOutBC.myType = "FixedPTOut"
function FixedPTOutBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { ExtrapolateCopy:new{xOrder = o.xOrder},
			FixedPT:new{p_out=o.p_out, T_out=o.T_out} }
   o.preSpatialDerivAction = { CopyCellData:new() }
   return o
end

FullFaceExchangeBC = BoundaryCondition:new()
FullFaceExchangeBC.myType = "FullFaceExchange"
function FullFaceExchangeBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { FullFaceExchangeCopy:new{otherBlock=o.otherBlock,
						 otherFace=o.otherFace,
						 orientation=o.orientation} }
   o.preSpatialDerivAction = { UpdateThermoTransCoeffs:new() }
   return o
end

UserDefinedBC = BoundaryCondition:new()
UserDefinedBC.myType = "UserDefined"
function UserDefinedBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { UserDefinedGhostCell:new{fileName=o.fileName} }
   o.preSpatialDerivAction = { UserDefinedInterface:new{fileName=o.fileName} } 
   return o
end

AdjacentToSolidBC = BoundaryCondition:new()
AdjacentToSolidBC.myType = "AdjacentToSolid"
function AdjacentToSolidBC:new(o)
   o = BoundaryCondition.new(self, o)
   o.preReconAction = { InternalCopyThenReflect:new() }
   o.preSpatialDerivAction = { CopyCellData:new(), ZeroVelocity:new(),
			       TemperatureFromGasSolidInterface:new{otherBlock=o.otherBlock,
							    otherFace=o.otherFace,
							    orientation=o.orientation},
			       WallKOmega:new() }
   o.postDiffFluxAction = { EnergyFluxFromAdjacentSolid:new{otherBlock=o.otherBlock,
							    otherFace=o.otherFace,
							    orientation=o.orientation }
   }
   return o
end

-- -----------------------------------------------------------------------

-- Class for Block construction (based on a StructuredGrid).
SBlock = {
   myType = "SBlock"
} -- end Block

function SBlock:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block, for later construction of the config file.
   -- Note that we want block id to start at zero for the D code.
   o.id = #(blocks)
   blocks[#(blocks)+1] = o
   -- Must have a grid and fillCondition
   assert(o.grid, "need to supply a grid")
   assert(o.fillCondition, "need to supply a fillCondition")
   -- Fill in default values, if already not set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("BLOCK-%d", o.id)
   o.omegaz = o.omegaz or 0.0
   o.bcList = o.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or SlipWallBC:new()
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells.
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if config.dimensions == 3 then
      o.nkc = o.grid:get_nkv() - 1
   else
      o.nkc = 1
   end
   -- The following table p for the corner locations,
   -- is to be used later for testing for block connections.
   o.p = {}
   if config.dimensions == 3 then
      o.p[0] = o.grid:get_vtx(0, 0, 0)
      o.p[1] = o.grid:get_vtx(o.nic, 0, 0)
      o.p[2] = o.grid:get_vtx(o.nic, o.njc, 0)
      o.p[3] = o.grid:get_vtx(0, o.njc, 0)
      o.p[4] = o.grid:get_vtx(0, 0, o.nkc)
      o.p[5] = o.grid:get_vtx(o.nic, 0, o.nkc)
      o.p[6] = o.grid:get_vtx(o.nic, o.njc, o.nkc)
      o.p[7] = o.grid:get_vtx(0, o.njc, o.nkc)
   else
      o.p[0] = o.grid:get_vtx(0, 0)
      o.p[1] = o.grid:get_vtx(o.nic, 0)
      o.p[2] = o.grid:get_vtx(o.nic, o.njc)
      o.p[3] = o.grid:get_vtx(0, o.njc)
   end
   -- print("Block id=", o.id, "p0=", tostring(o.p[0]), "p1=", tostring(o.p[1]),
   --       "p2=", tostring(o.p[2]), "p3=", tostring(o.p[3]))
   return o
end

function SBlock:tojson()
   local str = string.format('"block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   str = str .. string.format('    "omegaz": %f,\n', self.omegaz)
   str = str .. string.format('    "nhcell": %d,\n', #(self.hcellList))
   for i = 1, #(self.hcellList) do
      local hcell = self.hcellList[i]
      if config.dimensions == 3 then
	 str = str .. string.format('    "history-cell-%d": [%d, %d, %d],\n', 
				    i-1, hcell[1], hcell[2], hcell[3])
      else
	 str = str .. string.format('    "history-cell-%d": [%d, %d],\n',
				    i-1, hcell[1], hcell[2])
      end
   end
   -- Boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

-- ---------------------------------------------------------------------------

function closeEnough(vA, vB, tolerance)
   -- Decide if two Vector quantities are close enough to being equal.
   -- This will be used to test that the block corners coincide.
   tolerance = tolerance or 1.0e-4
   return (vabs(vA - vB)/(vabs(vA + vB)+1.0)) <= tolerance
end

function connectBlocks(blkA, faceA, blkB, faceB, orientation)
   print("connectBlocks: blkA=", blkA.id, "faceA=", faceA, 
	 "blkB=", blkB.id, "faceB=", faceB, "orientation=", orientation)
   blkA.bcList[faceA] = FullFaceExchangeBC:new{otherBlock=blkB.id, otherFace=faceB,
					       orientation=orientation}
   blkB.bcList[faceB] = FullFaceExchangeBC:new{otherBlock=blkA.id, otherFace=faceA,
					       orientation=orientation}
   -- [TODO] need to test for matching corner locations and 
   -- consistent numbers of cells
end

function isPairInList(targetPair, pairList)
   local count = 0
   for _,v in ipairs(pairList) do
      if (v[1] == targetPair[1] and v[2] == targetPair[2]) or
	 (v[2] == targetPair[1] and v[1] == targetPair[2]) 
      then
	 count = count + 1
      end
   end
   return count > 0
end

function verticesAreCoincident(A, B, vtxPairs, tolerance)
   tolerance = tolerance or 1.0e-6
   local allVerticesAreClose = true
   for _,v in ipairs(vtxPairs) do
      -- print("A.id=", A.id, "B.id=", B.id, "vtxPair=", tostringVtxPair(v))
      -- print("  A.p=", tostring(A.p[v[1]]), "B.p=", tostring(B.p[v[2]]))
      if vabs(A.p[v[1]] - B.p[v[2]]) > tolerance then
	 allVerticesAreClose = false
      end
   end
   return allVerticesAreClose
end

function identifyBlockConnections(blockList, excludeList, tolerance)
   -- Identify block connections by trying to match corner points.
   -- Parameters:
   -- blockList: the list of SBlock objects to be included in the search.
   --    If nil, the whole collection is searched.
   -- excludeList: list of pairs of SBlock objects that should not be
   --    included in the search for connections.
   -- tolerance: spatial tolerance for the colocation of vertices
   blockList = blockList or blocks
   excludeList = excludeList or {}
   tolerance = tolerance or 1.0e-6
   for _,A in ipairs(blockList) do
      for _,B in ipairs(blockList) do
	 if (A ~= B) and (not isPairInList({A, B}, excludeList)) then
	    -- print("Proceed with test for coincident vertices.")
	    local connectionCount = 0
	    if config.dimensions == 2 then
	       -- print("2D test")
	       for vtxPairs,connection in pairs(connections2D) do
		  -- print("vtxPairs=", tostringVtxPairList(vtxPairs),
		  --       "connection=", tostringConnection(connection))
		  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
		     local faceA, faceB = unpack(connection)
		     connectBlocks(A, faceA, B, faceB, 0)
		     connectionCount = connectionCount + 1
		  end
	       end
	    else
	       -- print("   3D test")
	       for vtxPairs,connection in pairs(connections3D) do
		  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
		     local faceA, faceB, orientation = unpack(connection)
		     connectBlocks(A, faceA, B, faceB, orientation)
		     connectionCount = connectionCount + 1
		  end
	       end
	    end
	    if connectionCount > 0 then
	       -- So we don't double-up on connections.
	       excludeList[#excludeList+1] = {A,B}
	    end
	 end -- if (A ~= B...
      end -- for _,B
   end -- for _,A
end

-- ---------------------------------------------------------------------------

function SBlockArray(t)
   -- Expect one table as argument, with named fields.
   -- Returns an array of blocks defined over a single region.
   assert(t.grid, "need to supply a grid")
   assert(t.fillCondition, "need to supply a fillCondition")
   t.omegaz = t.omegaz or 0.0
   t.bcList = t.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      t.bcList[face] = t.bcList[face] or SlipWallBC:new()
   end
   t.xforceList = t.xforceList or {}
   -- Numbers of subblocks in each coordinate direction
   t.nib = t.nib or 1
   t.njb = t.njb or 1
   t.nkb = t.nkb or 1
   if config.dimensions == 2 then
      t.nkb = 1
   end
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells in the D-domain.
   local nic_total = t.grid:get_niv() - 1
   local dnic = math.floor(nic_total/t.nib)
   local njc_total = t.grid:get_njv() - 1
   local dnjc = math.floor(njc_total/t.njb)
   local nkc_total = t.grid:get_nkv() - 1
   local dnkc = math.floor(nkc_total/t.nkb)
   if config.dimensions == 2 then
      nkc_total = 1
      dnkc = 1
   end
   local blockArray = {} -- will be a multi-dimensional array indexed as [i][j][k]
   local blockCollection = {} -- will be a single-dimensional array
   for ib = 1, t.nib do
      blockArray[ib] = {}
      local i0 = (ib-1) * dnic
      if (ib == t.nib) then
	 -- Last block has to pick up remaining cells.
	 dnic = nic_total - i0
      end
      for jb = 1, t.njb do
	 local j0 = (jb-1) * dnjc
	 if (jb == t.njb) then
	    dnjc = njc_total - j0
	 end
	 if config.dimensions == 2 then
	    -- 2D flow
	    local subgrid = t.grid:subgrid(i0,dnic+1,j0,dnjc+1)
	    local bcList = {north=SlipWallBC:new(), east=SlipWallBC:new(),
			    south=SlipWallBC:new(), west=SlipWallBC:new()}
	    if ib == 1 then
	       bcList[west] = t.bcList[west]
	    end
	    if ib == t.nib then
	       bcList[east] = t.bcList[east]
	    end
	    if jb == 1 then
	       bcList[south] = t.bcList[south]
	    end
	    if jb == t.njb then
	       bcList[north] = t.bcList[north]
	    end
	    new_block = SBlock:new{grid=subgrid, omegaz=t.omegaz,
				   fillCondition=t.fillCondition, bcList=bcList}
	    blockArray[ib][jb] = new_block
	    blockCollection[#blockCollection+1] = new_block
	 else
	    -- 3D flow, need one more level in the array
	    blockArray[ib][jb] = {}
	    for kb = 1, t.nkb do
	       local k0 = (kb-1) * dnkc
	       if (kb == t.nkb) then
		  dnkc = nkc_total - k0
	       end
	       local subgrid = t.grid:subgrid(i0,dnic+1,j0,dnjc+1,k0,dnkc+1)
	       local bcList = {north=SlipWallBC:new(), east=SlipWallBC:new(),
			       south=SlipWallBC:new(), west=SlipWallBC:new(),
			       top=SlipWallBC:new(), bottom=SlipWallBC:new()}
	       if ib == 1 then
		  bcList[west] = t.bcList[west]
	       end
	       if ib == t.nib then
		  bcList[east] = t.bcList[east]
	       end
	       if jb == 1 then
		  bcList[south] = t.bcList[south]
	       end
	       if jb == t.njb then
		  bcList[north] = t.bcList[north]
	       end
	       if kb == 1 then
		  bcList[bottom] = t.bcList[bottom]
	       end
	       if kb == t.nkb then
		  bcList[top] = t.bcList[top]
	       end
	       new_block = SBlock:new{grid=subgrid, omegaz=t.omegaz,
				      fillCondition=t.fillCondition,
				      bcList=bcList}
	       blockArray[ib][jb][kb] = new_block
	       blockCollection[#blockCollection+1] = new_block
	    end -- kb loop
	 end -- dimensions
      end -- jb loop
   end -- ib loop
   -- Make the inter-subblock connections
   if #blockCollection > 1 then
      identifyBlockConnections(blockCollection)
   end
   return blockArray
end

-- ---------------------------------------------------------------------------
-- Classes related to Solid blocks and boundary conditions
SolidBoundaryInterfaceEffect = {
   type = ""
}
function SolidBoundaryInterfaceEffect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

SolidBIE_FixedT = SolidBoundaryInterfaceEffect:new{Twall=300.0}
SolidBIE_FixedT.type = "fixed_temperature"
function SolidBIE_FixedT:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"Twall": %12.6e }', self.Twall)
   return str
end

SolidBIE_UserDefined = SolidBoundaryInterfaceEffect:new{fileName='user-defined-solid-bc.lua'}
SolidBIE_UserDefined.type = "user_defined"
function SolidBIE_UserDefined:tojson()
   local str = string.format('          {"type": "%s", ', self.type)
   str = str .. string.format('"filename": "%s" }', self.fileName)
   return str
end

-- Class for SolidBoundaryCondition
-- This class is a convenience class: it translates a high-level
-- user name for the boundary condition into a sequence of
-- lower-level operators.
SolidBoundaryCondition = {
   label = "",
   myType = "",
   setsFluxDirectly = false,
   preSpatialDerivAction = {}
}
function SolidBoundaryCondition:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end
function SolidBoundaryCondition:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", \n', self.label)
   str = str .. string.format('        "sets_flux_directly": %s,\n', tostring(self.setsFluxDirectly))
   str = str .. '        "pre_spatial_deriv_action": [\n'
   for i,effect in ipairs(self.preSpatialDerivAction) do
      str = str .. effect:tojson()
      -- Extra code to deal with annoying JSON trailing comma deficiency
      if i ~= #self.preSpatialDerivAction then str = str .. "," end
   end
   str = str .. '\n        ]\n    }'
   return str
end

SolidFixedTBC = SolidBoundaryCondition:new()
SolidFixedTBC.myType = "SolidFixedT"
function SolidFixedTBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivAction = { SolidBIE_FixedT:new{Twall=o.Twall} }
   return o
end

SolidUserDefinedBC = SolidBoundaryCondition:new()
SolidUserDefinedBC.myType = "SolidUserDefined"
function SolidUserDefinedBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.preSpatialDerivAction = { SolidBIE_UserDefined:new{fileName=o.fileName} }
   return o
end

SolidAdjacentToGasBC = SolidBoundaryCondition:new()
SolidAdjacentToGasBC.myType = "SolidAdjacentToGas"
function SolidAdjacentToGasBC:new(o)
   o = SolidBoundaryCondition.new(self, o)
   o.setsFluxDirectly = true
   return o
end

-- Class for SolidBlock construction (based on a StructuredGrid)
SSolidBlock = {
   myType = "SSolidBlock"
} -- end SSolidBlock

function SSolidBlock:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block for later construction of the config file.
   -- Note that we want the block id to start at zero for the D code.
   o.id = #solidBlocks
   solidBlocks[#solidBlocks+1] = o
   -- Must have a grid and initial temperature
   assert(o.grid, "need to supply a grid")
   assert(o.initTemperature, "need to supply an initTemperature")
   assert(o.properties, "need to supply physical properties for the block")
   -- Fill in some defaults, if not already set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("SOLIDBLOCK-%d", o.id)
   o.bcList = o.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or SolidFixedTBC:new{Twall=300}
   end
   -- Extract some information from the StructuredGrid
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if config.dimensions == 3 then
      o.nkc = o.grid.get_nkv() - 1
   else
      o.nkc = 1
   end
   -- The following table p for the corner locations,
   -- is to be used later for testing for block connections.
   o.p = {}
   if config.dimensions == 3 then
      o.p[0] = o.grid:get_vtx(0, 0, 0)
      o.p[1] = o.grid:get_vtx(o.nic, 0, 0)
      o.p[2] = o.grid:get_vtx(o.nic, o.njc, 0)
      o.p[3] = o.grid:get_vtx(0, o.njc, 0)
      o.p[4] = o.grid:get_vtx(0, 0, o.nkc)
      o.p[5] = o.grid:get_vtx(o.nic, 0, o.nkc)
      o.p[6] = o.grid:get_vtx(o.nic, o.njc, o.nkc)
      o.p[7] = o.grid:get_vtx(0, o.njc, o.nkc)
   else
      o.p[0] = o.grid:get_vtx(0, 0)
      o.p[1] = o.grid:get_vtx(o.nic, 0)
      o.p[2] = o.grid:get_vtx(o.nic, o.njc)
      o.p[3] = o.grid:get_vtx(0, o.njc)
   end
   return o
end

function SSolidBlock:tojson()
   local str = string.format('"solid_block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   str = str .. '    "properties": {\n'
   str = str .. string.format('       "rho": %.6e,\n', self.properties.rho)
   str = str .. string.format('       "k": %.6e,\n', self.properties.k)
   str = str .. string.format('       "Cp": %.6e\n', self.properties.Cp)
   str = str .. '    },\n'
   -- Boundary conditions
      for _,face in ipairs(faceList(config.dimensions)) do
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

-- --------------------------------------------------------------------

function write_control_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"dt_init": %e,\n', config.dt_init))
   f:write(string.format('"dt_max": %e,\n', config.dt_max))
   f:write(string.format('"cfl_value": %e,\n', config.cfl_value))
   f:write(string.format('"stringent_cfl": %s,\n', tostring(config.stringent_cfl)))
   f:write(string.format('"fixed_time_step": %s,\n', tostring(config.fixed_time_step)))
   f:write(string.format('"dt_reduction_factor": %e,\n', config.dt_reduction_factor))
   f:write(string.format('"print_count": %d,\n', config.print_count))
   f:write(string.format('"cfl_count": %d,\n', config.cfl_count))
   f:write(string.format('"max_time": %e,\n', config.max_time))
   f:write(string.format('"max_step": %d,\n', config.max_step))
   f:write(string.format('"dt_plot": %e,\n', config.dt_plot))
   f:write(string.format('"dt_history": %e,\n', config.dt_history))
   f:write(string.format('"write_at_step": %d,\n', config.write_at_step))
   f:write(string.format('"halt_now": %d\n', config.halt_now))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write("}\n")
   f:close()
end

function write_config_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"title": "%s",\n', config.title))
   f:write(string.format('"gas_model_file": "%s",\n', config.gas_model_file))
   f:write(string.format('"dimensions": %d,\n', config.dimensions))
   f:write(string.format('"axisymmetric": %s,\n',
			 tostring(config.axisymmetric)))
   f:write(string.format('"interpolation_order": %d,\n', config.interpolation_order))
   f:write(string.format('"gasdynamic_update_scheme": "%s",\n',
			 config.gasdynamic_update_scheme))
   f:write(string.format('"MHD": %s,\n', tostring(config.MHD)))
   f:write(string.format('"separate_update_for_viscous_terms": %s,\n',
			 tostring(config.separate_update_for_viscous_terms)))
   f:write(string.format('"separate_update_for_k_omega_source": %s,\n', 
			 tostring(config.separate_update_for_k_omega_source)))
   f:write(string.format('"apply_bcs_in_parallel": %s,\n',
			 tostring(config.apply_bcs_in_parallel)))
   f:write(string.format('"max_invalid_cells": %d,\n', config.max_invalid_cells))
   f:write(string.format('"thermo_interpolator": "%s",\n', 
			 string.lower(config.thermo_interpolator)))
   f:write(string.format('"interpolate_in_local_frame": %s,\n', 
			 tostring(config.interpolate_in_local_frame)))
   f:write(string.format('"apply_limiter": %s,\n', tostring(config.apply_limiter)))
   f:write(string.format('"extrema_clipping": %s,\n', tostring(config.extrema_clipping)))

   f:write(string.format('"flux_calculator": "%s",\n', config.flux_calculator))
   f:write(string.format('"compression_tolerance": %e,\n', config.compression_tolerance))
   f:write(string.format('"shear_tolerance": %e,\n', config.shear_tolerance))
   f:write(string.format('"M_inf": %e,\n', config.M_inf))

   f:write(string.format('"viscous": %s,\n', tostring(config.viscous)))
   f:write(string.format('"viscous_signal_factor": %e,\n', config.viscous_signal_factor))

   f:write(string.format('"turbulence_model": "%s",\n',
			 string.lower(config.turbulence_model)))
   f:write(string.format('"turbulence_prandtl_number": %g,\n',
			 config.turbulence_prandtl_number))
   f:write(string.format('"turbulence_schmidt_number": %g,\n',
			 config.turbulence_schmidt_number))
   f:write(string.format('"max_mu_t_factor": %e,\n', config.max_mu_t_factor))
   f:write(string.format('"transient_mu_t_factor": %e,\n', config.transient_mu_t_factor))

   f:write(string.format('"udf_source_terms_file": "%s",\n', config.udf_source_terms_file))
   f:write(string.format('"udf_source_terms": %s,\n', tostring(config.udf_source_terms)))

   f:write(string.format('"reacting": %s,\n', tostring(config.reacting)))
   f:write(string.format('"reactions_file": "%s",\n', config.reactions_file))

   f:write(string.format('"control_count": %d,\n', config.control_count))
   f:write(string.format('"nblock": %d,\n', #(blocks)))

   f:write(string.format('"block_marching": %s,\n',
			 tostring(config.block_marching)))
   f:write(string.format('"nib": %d,\n', config.nib))
   f:write(string.format('"njb": %d,\n', config.njb))
   f:write(string.format('"nkb": %d,\n', config.nkb))
   f:write(string.format('"propagate_inflow_data": %s,\n',
			 tostring(config.propagate_inflow_data)))
 
   f:write(string.format('"udf_solid_source_terms_file": "%s",\n', config.udf_solid_source_terms_file))
   f:write(string.format('"udf_solid_source_terms": %s,\n', tostring(config.udf_solid_source_terms)))
   f:write(string.format('"nsolidblock": %d,\n', #solidBlocks))

   for i = 1, #blocks do
      f:write(blocks[i]:tojson())
   end
   for i = 1, #solidBlocks do
      f:write(solidBlocks[i]:tojson())
   end
  
   f:write('"dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
   f:write("}\n")

   f:close()
end

function write_times_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# tindx sim_time dt_global\n");
   f:write(string.format("%04d %12.6e %12.6e\n", 0, 0.0, config.dt_init))
   f:close()
end

function write_block_list_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# indx label\n")
   for i = 1, #(blocks) do
      f:write(string.format("%4d %s\n", blocks[i].id, blocks[i].label))
   end
   f:close()
end

function build_job_files(job)
   print("Build job files for ", job)
   write_config_file(job .. ".config")
   write_control_file(job .. ".control")
   write_times_file(job .. ".times")
   write_block_list_file(job .. ".list")
   os.execute("mkdir -p grid/t0000")
   os.execute("mkdir -p flow/t0000")
   if #solidBlocks >= 1 then
      os.execute("mkdir -p solid-grid/t0000")
      os.execute("mkdir -p solid/t0000")
   end
   for i = 1, #blocks do
      local id = blocks[i].id
      print("Block id=", id)
      local fileName = "grid/t0000/" .. job .. string.format(".grid.b%04d.t0000", id)
      blocks[i].grid:write_to_text_file(fileName)
      os.execute("gzip -f " .. fileName)
      local fileName = "flow/t0000/" .. job .. string.format(".flow.b%04d.t0000", id)
      write_initial_flow_file(fileName, blocks[i].grid, blocks[i].fillCondition, 0.0)
      os.execute("gzip -f " .. fileName)
   end
   for i = 1, #solidBlocks do
      local id = solidBlocks[i].id
      print("SolidBlock id=", id)
      local fileName = "solid-grid/t0000/" .. job .. string.format(".solid-grid.b%04d.t0000", id)
      solidBlocks[i].grid:write_to_text_file(fileName)
      os.execute("gzip -f " .. fileName)
      local fileName = "solid/t0000/" .. job .. string.format(".solid.b%04d.t0000", id)
      writeInitialSolidFile(fileName, solidBlocks[i].grid,
			    solidBlocks[i].initTemperature, solidBlocks[i].properties, 0.0)
      os.execute("gzip -f " .. fileName)
   end

   print("Done building files.")
end

print("Done loading e4prep.lua")
