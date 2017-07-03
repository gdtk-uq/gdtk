-- foam-mesh.lua
-- Lua helper definitions for the D program: foam-mesh
--
-- Authors: Rowan G., Ingo J., and Peter J.
-- Date: 2017-07-03

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
	 i,j = string.find(allowedPrefix, v)
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
	      

-- Storage for FoamMesh objects
meshes = {}

FoamMesh = {}

function FoamMesh:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"psurface", "pvolume", "project3D",
				      "niv", "njv", "nkv",
				      "bndryList"})
   assert(flag, "Invalid name for item supplied to FoamMesh:new().")
   setmetatable(o, self)
   self.__index = self
   -- Make a record of this mesh for later use when writing out.
   o.id = #(meshes)
   meshes[#(meshes)+1] = o
   -- Need either a psurface or a pvolume
   if (o.psurface == nil and o.pvolume == nil) then
      print("Either a 'psurface' or 'pvolume' must be supplied to FoamMesh:new().")
      os.exit(1)
   end
   -- Let's do some more checks before going on with grid generation
   checkBndryLabels(o.bndryList)
   if o.psurface then
      assert(o.niv, "need to provide 'niv' in FoamMesh:new().")
      assert(o.njv, "need to provide 'njv' in FoamMesh:new().")
   end
   if o.pvolume then
      -- Check we have been given discretisation settings.
      assert(o.niv, "need to provide 'niv' in FoamMesh:new().")
      assert(o.njv, "need to provide 'njv' in FoamMesh:new().")
      assert(o.nkv, "need to provide 'nkv' in FoamMesh:new().")
   end
   -- Let's work on psurface, if present
   if o.psurface then
      -- We'll set the nkv to 2 for a one-cell slice
      o.nkv = 2
      -- We need a project3D type
      if o.project3D == nil then
	 print("A 'project3D' type of 'slab' or 'wedge' must be supplied when a 'psurface' is given to FoamMesh:new().")
	 os.exit(1)
      end
      if o.project3D == 'slab' then
	 o.pvolume = SlabVolume:new{face0123=o.psurface, dz=Vector3:new{x=0, y=0, z=1}}
      elseif o.project3D == 'wedge' then
	 o.pvolume = WedgeVolume:new{face0123=o.psurface, dtheta=1}
      else
	 print("The 'project3D' type ", o.pvolume, " is unknown.")
	 print("This error occurred in FoamMesh:new().")
	 os.exit(1)
      end
   end
   -- So at this point, we either have our handed volume,
   -- or we created one from a surface.
   -- Now it's time to build a grid.
   -- We build a structured grid, then directly convert it
   -- to an unstructured grid.
   -- TODO: Handle clustering.
   local sgrid = StructuredGrid:new{pvolume=o.pvolume,
				    niv=o.niv, njv=o.njv, nkv=o.nkv}
   o.usgrid = UnstructuredGrid:new{sgrid=sgrid, new_label="mesh-"~tostring(o.id)}
   return o
end
   

