-- A module for sticking the pieces related to block connections
-- in 2D and 3D
--
-- Authors: PJ and RJG
-- Date: 2015-10-01
--          Extracted from prep.lua

-- Symbolic names for identifying boundaries.
-- Although we no longer use these names in our Lua code,
-- there may be user-written scripts that do use these definitions.
north = "north"; NORTH = "north"
east = "east"; EAST = "east"
south = "south"; SOUTH = "south"
west = "west"; WEST = "west"
top = "top"; TOP = "top"
bottom = "bottom"; BOTTOM = "bottom"
west_face_indx = 0
east_face_indx = 1
south_face_indx = 2
north_face_indx = 3
bottom_face_indx = 4
top_face_indx = 5

function faceList(dimensions)
   local myList = {'west', 'east', 'south', 'north'}
   if dimensions == 3 then
      table.insert(myList, 'bottom')
      table.insert(myList, 'top')
   end
   return myList
end

-- -----------------------------------------------------------------------

-- Connections between 2D blocks, described as vertex-pairs.
-- When specifying a connection between 2D blocks,
-- we actually specify the names of the connecting faces.
-- The connection orientation is always 0.
local tabulatedData = {
   {{{3,0},{2,1}}, {'north', 'south'}},
   {{{3,3},{2,0}}, {'north', 'west'}},
   {{{3,2},{2,3}}, {'north', 'north'}},
   {{{3,1},{2,2}}, {'north', 'east'}},
   {{{0,0},{3,1}}, {'west',  'south'}},
   {{{0,3},{3,0}}, {'west',  'west'}},
   {{{0,2},{3,3}}, {'west',  'north'}},
   {{{0,1},{3,2}}, {'west',  'east'}},
   {{{1,0},{0,1}}, {'south', 'south'}},
   {{{1,3},{0,0}}, {'south', 'west'}},
   {{{1,2},{0,3}}, {'south', 'north'}},
   {{{1,1},{0,2}}, {'south', 'east'}},
   {{{2,0},{1,1}}, {'east',  'south'}},
   {{{2,3},{1,0}}, {'east',  'west'}},
   {{{2,2},{1,3}}, {'east',  'north'}},
   {{{2,1},{1,2}}, {'east',  'east'}},
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
local tabulatedData = {
   {{{3,2},{7,6},{6,7},{2,3}}, {'north', 'north', 0, '-i-j+k'}},
   {{{3,3},{7,2},{6,6},{2,7}}, {'north', 'north', 1, '+k-j+i'}},
   {{{3,7},{7,3},{6,2},{2,6}}, {'north', 'north', 2, '+i-j-k'}},
   {{{3,6},{7,7},{6,3},{2,2}}, {'north', 'north', 3, '-k-j-i'}},

   {{{3,0},{7,4},{6,5},{2,1}}, {'north', 'south', 0, '+i+j+k'}},
   {{{3,1},{7,0},{6,4},{2,5}}, {'north', 'south', 1, '+k+j-i'}},
   {{{3,5},{7,1},{6,0},{2,4}}, {'north', 'south', 2, '-i+j-k'}},
   {{{3,4},{7,5},{6,1},{2,0}}, {'north', 'south', 3, '-k+j+i'}},

   {{{3,1},{7,5},{6,6},{2,2}}, {'north', 'east', 0, '+j-i+k'}},
   {{{3,2},{7,1},{6,5},{2,6}}, {'north', 'east', 1, '+k-i-j'}},
   {{{3,6},{7,2},{6,1},{2,5}}, {'north', 'east', 2, '-j-i-k'}},
   {{{3,5},{7,6},{6,2},{2,1}}, {'north', 'east', 3, '-k-i+j'}},

   {{{3,3},{7,7},{6,4},{2,0}}, {'north', 'west', 0, '-j+i+k'}},
   {{{3,0},{7,3},{6,7},{2,4}}, {'north', 'west', 1, '+k+i+j'}},
   {{{3,4},{7,0},{6,3},{2,7}}, {'north', 'west', 2, '+j+i-k'}},
   {{{3,7},{7,4},{6,0},{2,3}}, {'north', 'west', 3, '-k+i-j'}},

   {{{3,4},{7,7},{6,6},{2,5}}, {'north', 'top', 0, '+i-k+j'}},
   {{{3,5},{7,4},{6,7},{2,6}}, {'north', 'top', 1, '+j-k-i'}},
   {{{3,6},{7,5},{6,4},{2,7}}, {'north', 'top', 2, '-i-k-j'}},
   {{{3,7},{7,6},{6,5},{2,4}}, {'north', 'top', 3, '-j-k+i'}},

   {{{3,1},{7,2},{6,3},{2,0}}, {'north', 'bottom', 0, '-i+k+j'}},
   {{{3,0},{7,1},{6,2},{2,3}}, {'north', 'bottom', 1, '+j+k+i'}},
   {{{3,3},{7,0},{6,1},{2,2}}, {'north', 'bottom', 2, '+i+k-j'}},
   {{{3,2},{7,3},{6,0},{2,1}}, {'north', 'bottom', 3, '-j+k-i'}},

   {{{1,2},{5,6},{4,7},{0,3}}, {'south', 'north', 0, '+i+j+k'}},
   {{{1,3},{5,2},{4,6},{0,7}}, {'south', 'north', 1, '-k+j+i'}},
   {{{1,7},{5,3},{4,2},{0,6}}, {'south', 'north', 2, '-i+j-k'}},
   {{{1,6},{5,7},{4,3},{0,2}}, {'south', 'north', 3, '+k+j-i'}},

   {{{1,0},{5,4},{4,5},{0,1}}, {'south', 'south', 0, '-i-j+k'}},
   {{{1,1},{5,0},{4,4},{0,5}}, {'south', 'south', 1, '-k-j-i'}},
   {{{1,5},{5,1},{4,0},{0,4}}, {'south', 'south', 2, '+i-j-k'}},
   {{{1,4},{5,5},{4,1},{0,0}}, {'south', 'south', 3, '+k-j+i'}},

   {{{1,1},{5,5},{4,6},{0,2}}, {'south', 'east', 0, '-j+i+k'}},
   {{{1,2},{5,1},{4,5},{0,6}}, {'south', 'east', 1, '-k+i-j'}},
   {{{1,6},{5,2},{4,1},{0,5}}, {'south', 'east', 2, '+j+i-k'}},
   {{{1,5},{5,6},{4,2},{0,1}}, {'south', 'east', 3, '+k+i+j'}},

   {{{1,3},{5,7},{4,4},{0,0}}, {'south', 'west', 0, '+j-i+k'}},
   {{{1,0},{5,3},{4,7},{0,4}}, {'south', 'west', 1, '-k-i+j'}},
   {{{1,4},{5,0},{4,3},{0,7}}, {'south', 'west', 2, '-j-i-k'}},
   {{{1,7},{5,4},{4,0},{0,3}}, {'south', 'west', 3, '+k-i-j'}},

   {{{1,4},{5,7},{4,6},{0,5}}, {'south', 'top', 0, '-i+k+j'}},
   {{{1,5},{5,4},{4,7},{0,6}}, {'south', 'top', 1, '-j+k-i'}},
   {{{1,6},{5,5},{4,4},{0,7}}, {'south', 'top', 2, '+i+k-j'}},
   {{{1,7},{5,6},{4,5},{0,4}}, {'south', 'top', 3, '+j+k+i'}},

   {{{1,1},{5,2},{4,3},{0,0}}, {'south', 'bottom', 0, '+i-k+j'}},
   {{{1,0},{5,1},{4,2},{0,3}}, {'south', 'bottom', 1, '-j-k+i'}},
   {{{1,3},{5,0},{4,1},{0,2}}, {'south', 'bottom', 2, '-i-k-j'}},
   {{{1,2},{5,3},{4,0},{0,1}}, {'south', 'bottom' , 3, '+j-k-i'}},

   {{{2,2},{6,6},{5,7},{1,3}}, {'east', 'north', 0, '-j+i+k'}},
   {{{2,3},{6,2},{5,6},{1,7}}, {'east', 'north', 1, '-j-k+i'}},
   {{{2,7},{6,3},{5,2},{1,6}}, {'east', 'north', 2, '-j-i-k'}},
   {{{2,6},{6,7},{5,3},{1,2}}, {'east', 'north', 3, '-j+k-i'}},

   {{{2,0},{6,4},{5,5},{1,1}}, {'east', 'south', 0, '+j-i+k'}},
   {{{2,1},{6,0},{5,4},{1,5}}, {'east', 'south', 1, '+j-k-i'}},
   {{{2,5},{6,1},{5,0},{1,4}}, {'east', 'south', 2, '+j+i-k'}},
   {{{2,4},{6,5},{5,1},{1,0}}, {'east', 'south', 3, '+j+k+i'}},

   {{{2,1},{6,5},{5,6},{1,2}}, {'east', 'east', 0, '-i-j+k'}},
   {{{2,2},{6,1},{5,5},{1,6}}, {'east', 'east', 1, '-i-k-j'}},
   {{{2,6},{6,2},{5,1},{1,5}}, {'east', 'east', 2, '-i+j-k'}},
   {{{2,5},{6,6},{5,2},{1,1}}, {'east', 'east', 3, '-i+k+j'}},

   {{{2,3},{6,7},{5,4},{1,0}}, {'east', 'west', 0, '+i+j+k'}},
   {{{2,0},{6,3},{5,7},{1,4}}, {'east', 'west', 1, '+i-k+j'}},
   {{{2,4},{6,0},{5,3},{1,7}}, {'east', 'west', 2, '+i-j-k'}},
   {{{2,7},{6,4},{5,0},{1,3}}, {'east', 'west', 3, '+i+k-j'}},

   {{{2,4},{6,7},{5,6},{1,5}}, {'east', 'top', 0, '-k-i+j'}},
   {{{2,5},{6,4},{5,7},{1,6}}, {'east', 'top', 1, '-k-j-i'}},
   {{{2,6},{6,5},{5,4},{1,7}}, {'east', 'top', 2, '-k+i-j'}},
   {{{2,7},{6,6},{5,5},{1,4}}, {'east', 'top', 3, '-k+j+i'}},

   {{{2,1},{6,2},{5,3},{1,0}}, {'east', 'bottom', 0, '+k+i+j'}},
   {{{2,0},{6,1},{5,2},{1,3}}, {'east', 'bottom', 1, '+k-j+i'}},
   {{{2,3},{6,0},{5,1},{1,2}}, {'east', 'bottom', 2, '+k-i-j'}},
   {{{2,2},{6,3},{5,0},{1,1}}, {'east', 'bottom', 3, '+k+j-i'}},

   {{{0,2},{4,6},{7,7},{3,3}}, {'west', 'north', 0, '+j-i+k'}},
   {{{0,3},{4,2},{7,6},{3,7}}, {'west', 'north', 1, '+j+k+i'}},
   {{{0,7},{4,3},{7,2},{3,6}}, {'west', 'north', 2, '+j+i-k'}},
   {{{0,6},{4,7},{7,3},{3,2}}, {'west', 'north', 3, '+j-k-i'}},

   {{{0,0},{4,4},{7,5},{3,1}}, {'west', 'south', 0, '-j+i+k'}},
   {{{0,1},{4,0},{7,4},{3,5}}, {'west', 'south', 1, '-j+k-i'}},
   {{{0,5},{4,1},{7,0},{3,4}}, {'west', 'south', 2, '-j-i-k'}},
   {{{0,4},{4,5},{7,1},{3,0}}, {'west', 'south', 3, '-j-k+i'}},

   {{{0,1},{4,5},{7,6},{3,2}}, {'west', 'east', 0, '+i+j+k'}},
   {{{0,2},{4,1},{7,5},{3,6}}, {'west', 'east', 1, '+i+k-j'}},
   {{{0,6},{4,2},{7,1},{3,5}}, {'west', 'east', 2, '+i-j-k'}},
   {{{0,5},{4,6},{7,2},{3,1}}, {'west', 'east', 3, '+i-k+j'}},

   {{{0,3},{4,7},{7,4},{3,0}}, {'west', 'west', 0, '-i-j+k'}},
   {{{0,0},{4,3},{7,7},{3,4}}, {'west', 'west', 1, '-i+k+j'}},
   {{{0,4},{4,0},{7,3},{3,7}}, {'west', 'west', 2, '-i+j-k'}},
   {{{0,7},{4,4},{7,0},{3,3}}, {'west', 'west', 3, '-i-k-j'}},

   {{{0,4},{4,7},{7,6},{3,5}}, {'west', 'top', 0, '+k+i+j'}},
   {{{0,5},{4,4},{7,7},{3,6}}, {'west', 'top', 1, '+k+j-i'}},
   {{{0,6},{4,5},{7,4},{3,7}}, {'west', 'top', 2, '+k-i-j'}},
   {{{0,7},{4,6},{7,5},{3,4}}, {'west', 'top', 3, '+k-j+i'}},

   {{{0,1},{4,2},{7,3},{3,0}}, {'west', 'bottom', 0, '-k-i+j'}},
   {{{0,0},{4,1},{7,2},{3,3}}, {'west', 'bottom', 1, '-k+j+i'}},
   {{{0,3},{4,0},{7,1},{3,2}}, {'west', 'bottom', 2, '-k+i-j'}},
   {{{0,2},{4,3},{7,0},{3,1}}, {'west', 'bottom', 3, '-k-j-i'}},

   {{{5,2},{6,6},{7,7},{4,3}}, {'top', 'north', 0, '+i+k-j'}},
   {{{5,3},{6,2},{7,6},{4,7}}, {'top', 'north', 1, '-k+i-j'}},
   {{{5,7},{6,3},{7,2},{4,6}}, {'top', 'north', 2, '-i-k-j'}},
   {{{5,6},{6,7},{7,3},{4,2}}, {'top', 'north', 3, '+k-i-j'}},

   {{{5,0},{6,4},{7,5},{4,1}}, {'top', 'south', 0, '-i+k+j'}},
   {{{5,1},{6,0},{7,4},{4,5}}, {'top', 'south', 1, '-k-i+j'}},
   {{{5,5},{6,1},{7,0},{4,4}}, {'top', 'south', 2, '+i-k+j'}},
   {{{5,4},{6,5},{7,1},{4,0}}, {'top', 'south', 3, '+k+i+j'}},

   {{{5,1},{6,5},{7,6},{4,2}}, {'top', 'east', 0, '-j+k-i'}},
   {{{5,2},{6,1},{7,5},{4,6}}, {'top', 'east', 1, '-k-j-i'}},
   {{{5,6},{6,2},{7,1},{4,5}}, {'top', 'east', 2, '+j-k-i'}},
   {{{5,5},{6,6},{7,2},{4,1}}, {'top', 'east', 3, '+k+j-i'}},

   {{{5,3},{6,7},{7,4},{4,0}}, {'top', 'west', 0, '+j+k+i'}},
   {{{5,0},{6,3},{7,7},{4,4}}, {'top', 'west', 1, '-k+j+i'}},
   {{{5,4},{6,0},{7,3},{4,7}}, {'top', 'west', 2, '-j-k+i'}},
   {{{5,7},{6,4},{7,0},{4,3}}, {'top', 'west', 3, '+k-j+i'}},

   {{{5,4},{6,7},{7,6},{4,5}}, {'top', 'top', 0, '-i+j-k'}},
   {{{5,5},{6,4},{7,7},{4,6}}, {'top', 'top', 1, '-j-i-k'}},
   {{{5,6},{6,5},{7,4},{4,7}}, {'top', 'top', 2, '+i-j-k'}},
   {{{5,7},{6,6},{7,5},{4,4}}, {'top', 'top', 3, '+j+i-k'}},

   {{{5,1},{6,2},{7,3},{4,0}}, {'top', 'bottom', 0, '+i+j+k'}},
   {{{5,0},{6,1},{7,2},{4,3}}, {'top', 'bottom', 1, '-j+i+k'}},
   {{{5,3},{6,0},{7,1},{4,2}}, {'top', 'bottom', 2, '-i-j+k'}},
   {{{5,2},{6,3},{7,0},{4,1}}, {'top', 'bottom', 3, '+j-i+k'}},

   {{{0,2},{3,6},{2,7},{1,3}}, {'bottom', 'north', 0, '-i+k+j'}},
   {{{0,3},{3,2},{2,6},{1,7}}, {'bottom', 'north', 1, '+k+i+j'}},
   {{{0,7},{3,3},{2,2},{1,6}}, {'bottom', 'north', 2, '+i-k+j'}},
   {{{0,6},{3,7},{2,3},{1,2}}, {'bottom', 'north', 3, '-k-i+j'}},

   {{{0,0},{3,4},{2,5},{1,1}}, {'bottom', 'south', 0, '+i+k-j'}},
   {{{0,1},{3,0},{2,4},{1,5}}, {'bottom', 'south', 1, '+k-i-j'}},
   {{{0,5},{3,1},{2,0},{1,4}}, {'bottom', 'south', 2, '-i-k-j'}},
   {{{0,4},{3,5},{2,1},{1,0}}, {'bottom', 'south', 3, '-k+i-j'}},

   {{{0,1},{3,5},{2,6},{1,2}}, {'bottom', 'east', 0, '+j+k+i'}},
   {{{0,2},{3,1},{2,5},{1,6}}, {'bottom', 'east', 1, '+k-j+i'}},
   {{{0,6},{3,2},{2,1},{1,5}}, {'bottom', 'east', 2, '-j-k+i'}},
   {{{0,5},{3,6},{2,2},{1,1}}, {'bottom', 'east', 3, '-k+j+i'}},

   {{{0,3},{3,7},{2,4},{1,0}}, {'bottom', 'west', 0, '-j+k-i'}},
   {{{0,0},{3,3},{2,7},{1,4}}, {'bottom', 'west', 1, '+k+j-i'}},
   {{{0,4},{3,0},{2,3},{1,7}}, {'bottom', 'west', 2, '+j-k-i'}},
   {{{0,7},{3,4},{2,0},{1,3}}, {'bottom', 'west', 3, '-k-j-i'}},

   {{{0,4},{3,7},{2,6},{1,5}}, {'bottom', 'top', 0, '+i+j+k'}},
   {{{0,5},{3,4},{2,7},{1,6}}, {'bottom', 'top', 1, '+j-i+k'}},
   {{{0,6},{3,5},{2,4},{1,7}}, {'bottom', 'top', 2, '-i-j+k'}},
   {{{0,7},{3,6},{2,5},{1,4}}, {'bottom', 'top', 3, '-j+i+k'}},

   {{{0,1},{3,2},{2,3},{1,0}}, {'bottom', 'bottom', 0, '-i+j-k'}},
   {{{0,0},{3,1},{2,2},{1,3}}, {'bottom', 'bottom', 1, '+j+i-k'}},
   {{{0,3},{3,0},{2,1},{1,2}}, {'bottom', 'bottom', 2, '+i-j-k'}},
   {{{0,2},{3,3},{2,0},{1,1}}, {'bottom', 'bottom', 3, '-j-i-k'}}
} -- end of tabulatedData

connections3D = {}
vtxPairs3D = {}
vtxPairs3DwTableKey = {}
-- When reading GridPro block connectivity file,
-- we need to look up Eilmer notation for connection orientations.
eilmer_orientation = {}
for _,v in ipairs(tabulatedData) do
   local vtxPairs, connection = unpack(v)
   table.sort(vtxPairs, cmpVtxPair)
   connections3D[vtxPairs] = connection
   local this_face, other_face, orientation, axis_map = unpack(connection)
   vtxPairs3D[this_face..other_face..orientation] = vtxPairs
   -- Also add key-as-table form
   vtxPairs3DwTableKey[{this_face, other_face, orientation}] = vtxPairs
   eilmer_orientation[this_face..other_face..axis_map] = orientation
end

----------------------------------------------------------------

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

function closeEnough(vA, vB, tolerance)
   -- Decide if two Vector quantities are close enough to being equal.
   -- This will be used to test that the block or grid corners coincide.
   tolerance = tolerance or 1.0e-4
   return (vabs(vA - vB)/(vabs(vA + vB)+1.0)) <= tolerance
end

function verticesAreCoincident(A, B, vtxPairs, tolerance)
   if false then -- debug
      print("allVerticesAreClose:")
   end
   tolerance = tolerance or 1.0e-6
   local allVerticesAreClose = true
   for _,v in ipairs(vtxPairs) do
      if false then -- debug
         print("  test vertex pair:")
         print("  A.id=", A.id, "B.id=", B.id, "vtxPair=", tostringVtxPair(v))
         print("    A.p=", tostring(A.p[v[1]]), "B.p=", tostring(B.p[v[2]]))
      end
      if not closeEnough(A.p[v[1]], B.p[v[2]], tolerance) then
	 allVerticesAreClose = false
      end
   end
   if false then -- debug
      print("  result:", allVerticesAreClose)
   end
   return allVerticesAreClose
end
