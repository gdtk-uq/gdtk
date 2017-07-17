-- foam-mesh.lua
-- Lua helper definitions for the D program: foam-mesh
--
-- Authors: Rowan G., Ingo J., and Peter J.
-- Date: 2017-07-03

faceMap = {
   north=0,
   east=1,
   south=2,
   west=3,
   top=4,
   bottom=5
}

function checkAllowedNames(myTable, allowedNames)
   local setOfNames = {}
   local namesOk = true
   for i,name in ipairs(allowedNames) do
      setOfNames[name] = true
   end
   for k,v in pairs(myTable) do
      if not setOfNames[k] then
	 print("Warning: Invalid name: ", k)
	 namesOk = false
      end
   end
   return namesOk
end

-- Allowed boundary label prefixes
bndryLabelPrefixes = {"w-", -- for walls
                      "i-", -- for inlets
		      "o-", -- for outlets
		      "s-", -- for symmetry
		      "p-", -- for patches
}

function checkBndryLabels(bndryList)
   for k,v in pairs(bndryList) do
      local labelOK = false
      for _,allowedPrefix in ipairs(bndryLabelPrefixes) do
	 pattern = string.gsub(allowedPrefix, "%-", "%%%-")
	 i, j = string.find(v, pattern)
	 if (i == 1) then
	    labelOK = true
	 end
      end
      -- If labelOK is still false at end, then this particular
      -- label was badly formed.
      if not labelOK then
	 print(string.format("The boundary label '%s' is not allowed.", v))
	 print("Allowed label names start with the following prefixes:")
	 for _,allowedPrefix in ipairs(bndryLabelPrefixes) do
	    print(allowedPrefix)
	 end
	 os.exit(1)
      end
   end
end
	      
-- Storage for global collection of boundary labels
globalBndryLabels = {}

-- Storage for FoamBlock objects
blks = {}

-- Class definition
FoamBlock = {}
function FoamBlock:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"grid", "bndry_labels"})
   assert(flag, "Invalid name for item supplied to FoamBlock:new().")
   setmetatable(o, self)
   self.__index = self
   -- Make a record of this block for later use when writing out.
   o.id = #(blks)
   blks[#(blks)+1] = o
   if (o.grid == nil) then
      error("A 'grid' object must be supplied to FoamBlock:new().")
   end
   if (o.grid:get_dimensions() ~= 3) then
      errMsg = "The 'grid' object supplied to FoamBlock:new() must be a 3D grid.\n"
      errMsg = errMsg .. "You can convert a 2D grid using 'makeSlabGrid' or 'makeWedgeGrid' functions."
      error(errMsg)
   end
   if (o.grid:get_type() == "structured_grid") then
      -- We'll need to convert to unstructured before we can proceed
      o.grid = UnstructuredGrid:new{sgrid=o.grid}
   end
   checkBndryLabels(o.bndry_labels)
   -- Populate the unset bndry_labels with the internal defaults
   for face, idx in pairs(faceMap) do
      o.bndry_labels[face] = o.bndry_labels[face] or o.grid:get_boundaryset_tag(idx)
   end
   -- Add the unique boundary labels to the global collection
   for _,bl in pairs(o.bndry_labels) do
      globalBndryLabels[bl] = true
   end
   return o
end

function writeMeshes()
   if #blks > 1 then
      print("foamMesh only works on the first block presently.")
   end
   if (verbosityLevel >= 1) then
      print("Writing out grid into 'polyMesh/'")
   end
   blks[1].grid:writeOpenFoamPolyMesh("constant")
end

function main(verbosityLevel)
   vrbLvl = verbosityLevel
   writeMeshes()
end   

