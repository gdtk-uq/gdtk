-- openfoam.lua
-- Lua helper definitions for OpenFOAM.
--
-- Authors: Rowan G., Ingo J., and Peter J.
-- Date: 2017-07-03
--
-- Reece O. 2025-02-07
-- These definitions were moved from foam-mesh.lua to facilitate grid 
-- generation for PATO.

-- Global settings
axisymmetric = false -- Options are: "true" and "false".
collapseEdges = false -- Options are: "true" and "false".
openFoamDimensions = 2
dtheta = 0.01 -- reduced to approx 1 degree
dz = 0.2

faceMap = {
   west=0,
   east=1,
   south=2,
   north=3,
   bottom=4,
   top=5
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
bndryLabelPrefixes = {
   "w-", -- for walls
   "i-", -- for inlets
   "o-", -- for outlets
   "s-", -- for symmetry
   "p-", -- for patches
}

patoBcPrefixes = {
   "w-", -- for walls
   "a-", -- for adiabatic boundaries
   "s-", -- for symmetry
}

function checkBndryLabels(bndryList, bcPrefixes)
   for k,v in pairs(bndryList) do
      local labelOK = false
      for _,allowedPrefix in ipairs(bcPrefixes) do
         pattern = string.gsub(allowedPrefix, "%%-", "%%%%-")
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
         for _,allowedPrefix in ipairs(bcPrefixes) do
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
   --if (o.grid:get_dimensions() ~= 2) then
   --   errMsg = "The 'grid' object supplied to FoamBlock:new() must be a 2D grid.\n"
   --   error(errMsg)
   --end
   if (o.bndry_labels == nil) then
      o.bndry_labels = {}  -- create empty list for cases where block doesn't have outwards facing faces.
   end
   if (o.grid:get_type() ~= "structured_grid") then
      error("The 'grid' object supplied to FoamBlock:new() must be a structured grid.\n")
   end
   if (o.grid:get_dimensions() == 2) then
      -- Construct a slab or wedge, as appropriate.
      -- First, test that the 2D grid is "face-up" in the z-direction.
      local cell_0_properties = quadProperties{p0=o.grid:get_vtx(0, 0), p1=o.grid:get_vtx(1, 0),
                                               p2=o.grid:get_vtx(1, 1), p3=o.grid:get_vtx(0, 1)}
      if cell_0_properties.n.z < 0.0 then
         error("The 'grid' object supplied to FoamBlock:new() is facing in the negative z-direction.\n")
      end
      if (axisymmetric) then
         newGrid = o.grid:makeWedgeGrid{dtheta=dtheta, symmetric=true}
      else
         newGrid = o.grid:makeSlabGrid{dz=dz}
      end
   elseif (o.grid:get_dimensions() == 3) then
         newGrid = o.grid
         openFoamDimensions = 3
   else
      error("The 'grid' object supplied to FoamBlock:new() must be a 2D or 3D grid.\n")
   end
   -- and then convert to unstructured
   o.ugrid = UnstructuredGrid:new{sgrid=newGrid}

   -- Now look over the boundary labels.
   if o.pato then
      bcPrefixes = patoBcPrefixes
   else
      bcPrefixes = bndryLabelPrefixes
   end
   checkBndryLabels(o.bndry_labels, bcPrefixes)
   -- Add "top", "bottom" labels for 2-D
   if (o.grid:get_dimensions() == 2) then
      if (axisymmetric) then
         o.bndry_labels.top = "wedge-front"
         o.bndry_labels.bottom = "wedge-rear"
      else
         o.bndry_labels.top = "FrontAndBack"
         o.bndry_labels.bottom = "FrontAndBack"
      end
   end
   -- Populate the unset bndry_labels with the defaults
   for _,face in ipairs({"north", "east", "south", "west","top","bottom"}) do
      o.bndry_labels[face] = o.bndry_labels[face] or "unassigned"
   end
   -- Add the unique boundary labels to the global collection
   for _,bl in pairs(o.bndry_labels) do
      globalBndryLabels[bl] = true
   end
   return o
end

function markInternalBoundaries(grid, blks)
   for ib, blk in ipairs(blks) do
      for bndry,_ in pairs(blk.bndry_labels) do
         -- Convert to boundaryset in master grid
         iBndry = (ib-1)*6 + faceMap[bndry]
         if grid:is_boundaryset_empty(iBndry) then
            blk.bndry_labels[bndry] = "internal"
         end
      end
   end
end

function amendTags(grid, blks)
   if (vrbLvl >= 2) then
      print("Amending tags in master mesh.")
   end
   nBoundaries = grid:get_nboundaries()
   for iBndry=0,nBoundaries-1 do
      origTag = grid:get_boundaryset_tag(iBndry)
      newTag = string.format("%s-%04d", origTag, math.floor(iBndry/6))
      grid:set_boundaryset_tag(iBndry, newTag)
   end
   if (vrbLvl >= 2) then
      print("   DONE: Amending tags in master mesh.")
   end
end

function runCollapseEdges()
   if (vrbLvl >= 1) then
      print("Running OpenFOAM command: collapseEdges.")
   end
   -- A 2-step process:
   -- 1. Place the collapseDict file in place.
   dgdDir = os.getenv("DGD")
   collapseDictFile = string.format("%s/share/foamMesh-templates/collapseDict", dgdDir)
   retVal = os.execute("test -d system")
   if retVal ~= 0 then
      os.execute("mkdir system")
   end
   cmd = string.format("cp %s system/", collapseDictFile)
   os.execute(cmd)
   -- 2. Run the collapeEdges command
   cmd = "collapseEdges -overwrite -noZero"
   os.execute(cmd)
   if (vrbLvl >= 1) then
      print("   DONE: Running OpenFOAM command: collapseEdges.")
   end
end

 function runCreatePatchDict()
   if (vrbLvl >= 1) then
      print("Running OpenFOAM command: createPatchDict.")
   end
   local flag = os.execute("createPatch -overwrite")
   assert(flag, "Cannot find command createPatch, check that OpenFOAM environment has been loaded.")
   if (vrbLvl >= 1) then
      print("   DONE: Running OpenFOAM command: createPatch.")
   end
end

function runCreatePatchDict_empty()
   if (vrbLvl >= 1) then
      print("Creating template.")
   end
   -- Now copy required template files in place.
   foamTmpltDir = os.getenv("DGD").."/share/foamMesh-templates"
   filesToCopy = {"createPatchDict"}
   for _,f in ipairs(filesToCopy) do
      cmd = string.format("cp %s/%s %s/", foamTmpltDir, f, "system")
      os.execute(cmd)
   end
   if (vrbLvl >= 1) then
      print("   DONE: Creating templates.")
   end
   if (vrbLvl >= 1) then
      print("Running OpenFOAM command: createPatchDict.")
   end
   local flag = os.execute("createPatch -overwrite")
   assert(flag, "Cannot find command createPatch, check that OpenFOAM environment has been loaded.")
   if (vrbLvl >= 1) then
      print("   DONE: Running OpenFOAM command: createPatch.")
   end
end

function writeFoamHeader(f)
   f:write(string.format("// Auto-generated by foamMesh on %s\n", os.date("%d-%b-%Y at %X")))
   f:write("/*--------------------------------*- C++ -*----------------------------------*| \n")
   f:write("| =========                 |                                                 | \n")
   f:write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n")
   f:write("|  \\    /   O peration     | Version:  2.4.0                                 | \n")
   f:write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n")
   f:write("|    \\/     M anipulation  |                                                 | \n")
   f:write("|*---------------------------------------------------------------------------*/ \n")
   f:write("FoamFile \n")
   f:write("{\n")
   f:write("    version     2.0;\n")
   f:write("    format      ascii;\n")
end

function writeStandardBC(f)
      if label == "FrontAndBack" then
         f:write("    frontAndBack \n")
         f:write("    { \n")
         f:write("    empty,\n")
         f:write("    }\n")
      end
      if label == "wedge-front" then
         f:write("    wedge-front \n")
         f:write("    { \n")
         f:write("    symmetry;\n")
         f:write("    }\n")
      end
      if label == "wedge-rear" then
         f:write("    wedge-rear \n")
         f:write("    { \n")
         f:write("    symmetry;\n")
         f:write("    }\n")
      end
      if label == "unassigned" then
         -- do nothing
      end
end

function writeMesh()
   if (vrbLvl >= 1) then
      print("Writing out grid into 'constant/polyMesh/'")
   end
   myMesh:writeOpenFoamPolyMesh("constant")
   if (vrbLvl >= 1) then
      print("   DONE: Writing out grid into 'constant/polyMesh/'")
   end
end

function runRenumberMesh()
   if (vrbLvl >= 1) then
      print("Running OpenFOAM command: renumberMesh.")
   end
   -- Check if 0 exists.
   retVal = os.execute("test -d 0")
   if retVal == 0 then
      -- 0/ already exists. We don't want renumberMesh
      if (vrbLvl >= 1) then
         print("   SKIPPED: Running OpenFOAM command: renumberMesh.")
         print("   Run manually once /0 is set-up.")
      end
   else
      local flag = os.execute("renumberMesh -overwrite -noZero")
      assert(flag, "Cannot find command renumberMesh, check that OpenFOAM environment has been loaded.")
      if (vrbLvl >= 1) then
         print("   DONE: Running OpenFOAM command: renumberMesh.")
      end
   end
end

function runCheckMesh()
   if (vrbLvl >= 1) then
      print("Running OpenFOAM command: checkMesh.")
   end
   local flag = os.execute("checkMesh")
   assert(flag, "Cannot find command checkMesh, check that OpenFOAM environment has been loaded.")
   if (vrbLvl >= 1) then
      print("   DONE: Running OpenFOAM command: checkMesh.")
   end
end

function clearPolyMesh()
   if (vrbLvl >= 1) then
      print("Clearing PolyMesh folder.")
   end
   -- Check if constant/polyMesh exists.
   retVal = os.execute("test -d constant/polyMesh")
   if retVal == 0 then
       os.execute("rm -r constant/polyMesh")
   end
   if (vrbLvl >= 1) then
      print("   DONE: existing polyMesh folder deleted.")
   end
end
