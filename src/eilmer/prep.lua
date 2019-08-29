-- prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- 
print("Loading prep.lua...")

-- This module is not required by prep itself, but we load it here
-- to make it available in the user's script.
require 'billig'

require 'lua_helper'
local deepclone = lua_helper.deepclone

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

NGHOST = 2

require 'blk_conn'
-- Let's pull the symbols out of the blk_conn module
-- and make them global in this namespace
for k,v in pairs(blk_conn) do
   _G[k] = v
end

require 'bc'
for k,v in pairs(bc) do
   _G[k] = v
end

require 'gridpro'
-- Make these functions global so that users may call them
-- in the configuration script
applyGridproConnectivity = gridpro.applyGridproConnectivity
applyGridproBoundaryConditions = gridpro.applyGridproBoundaryConditions

-- Storage for steady-state solver settings
sssOptionsHidden = { -- hidden from user
   -- set defaults here
   use_preconditioner = true,
   frozen_preconditioner_count = 1,
   start_preconditioning = 1,
   ilu_fill = 0,
   precondition_matrix_type = "block_diagonal",
   use_scaling = true,
   use_complex_matvec_eval = false,
   number_pre_steps = 10,
   number_total_steps = 100,
   max_number_attempts = 3,
   stop_on_relative_global_residual = 1.0e-99,
   stop_on_absolute_global_residual = 1.0e-99,
   -- Restarted preconditioned FGMRES settings
   max_outer_iterations = 10,
   max_restarts = 10,
   number_inner_iterations = 5,
   -- Options for start-up phase
   number_start_up_steps = 5,
   cfl0 = 1.0,
   eta0 = 0.5,
   tau0 = 0.1,
   sigma0 = 1.0e-8,
   p0 = 0.75,
   -- Options for inexact Newton phase
   cfl1 = 10.0,
   tau1 = 0.1,
   sigma1 = 1.0e-8,
   p1 = 1.0,
   eta_strategy = "constant",
   eta1 = 0.5,
   eta1_max = 0.9,
   eta1_min = 0.01,
   eta_ratio_per_step = 0.9,
   gamma = 0.9,
   alpha = 2.0,
   -- Options related to writing out snapshots and diagnostics
   snapshots_count = 10,
   number_total_snapshots = 5,
   write_diagnostics_count = 20,
   write_loads_count = 20,

   __index = function (t, k) 
      return sssOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if sssOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'SteadyStateSolver' table.", k))
      else
	 sssOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 sssOptionsHidden.__newindex(t, k, v)
      end
   end
}

SteadyStateSolver = {}
setmetatable(SteadyStateSolver, sssOptionsHidden)

-- Storage for shape sensitivity calculator settings
sscOptionsHidden = { -- hidden from user
   -- set defaults here
   pseudotime = false,
   pseudotime_lhs_jacobian_order = 1,
   adjoint_precondition_matrix_order = 0,
   read_frozen_limiter_values_from_file = false,
   -- sensitivity parameters
   epsilon = 1.0e-30,
   -- GMRES parameters
   maxOuterIterations = 10,
   maxRestarts = 10,
   cfl0=1.0,
   eta = 0.1, 
   stop_on_relative_global_residual = 1.0e-99,
   -- Bezier curve fit parameters
   tol_bezier_curve_fit = 1.0e-06,
   max_steps_bezier_curve_fit = 10000,
   -- user-defined file
   user_defined_objective_file = "dummy-obj-file.lua",
   
   __index = function (t, k) 
      return sscOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if sscOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'ShapeSensitivityCalculator' table.", k))
      else
	 sscOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 sscOptionsHidden.__newindex(t, k, v)
      end
   end
}

ShapeSensitivityCalculator = {}
setmetatable(ShapeSensitivityCalculator, sscOptionsHidden)

-- Storage for solid domain loose update settings
sdluOptionsHidden = { -- hidden from user
   -- set defaults here
   max_newton_iterations = 10,
   tolerance_newton_update = 1.0e-2,
   max_gmres_iterations = 10,
   tolerance_gmres_solve = 1.0e-3,
   perturbation_size = 1.0e-2,

   __index = function (t, k) 
      return sdluOptionsHidden[k]
   end,
   __newindex = function (t, k, v)
      if sdluOptionsHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'SolidDomainLooseUpdate' table.", k))
      else
	 sdluOptionsHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 sdluOptionsHidden.__newindex(t, k, v)
      end
   end
}

SolidDomainLooseUpdate = {}
setmetatable(SolidDomainLooseUpdate, sdluOptionsHidden)

-- ---------------------------------------------------------------------------------------
-- A class to build FlowState objects as Lua tables.

FlowState_defaults = {
   p = 100.0e3, -- pressure, Pa
   T = 300.0, -- static (trans-rotational) temperature, K
   p_e = 0.0, -- electron pressure
   quality = 1.0,
   --
   -- Once we have an initialized GasModel, we will patch
   -- the following couple of attributes
   massf = {["air"]=1.0, }, -- mass fractions
   T_modes = {}, -- temperatures for other internal energy modes, K
   --
   mu = 0.0,
   k = 0.0,
   --
   velx = 0.0, -- velocity component in x-direction, m/s
   vely = 0.0, -- velocity component in y-direction, m/s
   velz = 0.0, -- velocity component in z-direction, m/s
   Bx = 0.0, -- magnetic field component in x-direction
   By = 0.0,
   Bz = 0.0,
   psi = 0.0, -- magnetic-field parameter
   divB = 0.0, -- magnetic field divergence parameter
   tke = 0.0, -- turbulent kinetic energy, (m/s)^^2
   omega = 1.0, -- turbulence frequency, 1/s
   mu_t = 0.0, -- turbulence viscosity, Pa.s
   k_t = 0.0, -- turbulence conductivity
   S = 0, -- shock indicator
}

-- Prototype for the FlowState class.
FlowState = {
   myType = "FlowState",
}

function FlowState:new(o)
   local flag = type(self)=='table' and self.myType=='FlowState'
   if not flag then
      error("Make sure that you are using FlowState:new{} and not FlowState.new{}", 2)
   end
   o = o or {}
   -- We limit the names of the entries so that the user does not supply items
   -- (such as density or energy) thinking that we will be setting the FlowState
   -- object using those values.
   flag =  checkAllowedNames(o, {"p", "T", "T_modes", "p_e",
                                 "quality", "massf",
                                 "velx", "vely", "velz",
                                 "Bx", "By", "Bz", "psi", "divB",
                                 "tke", "omega", "mu_t", "k_t", "S"})
   if not flag then
      error("Invalid name for item supplied to FlowState constructor.", 2)
   end
   -- We want to make consistent values/calculations in the context
   -- of an initialized gas model, but we only want one such gas model and gas state.
   if FlowState.gm == nil then
      local gm = getGasModel()
      local nsp = gm:nSpecies()
      local names = {}
      for i = 0, nsp-1 do
         local spName = gm:speciesName(i)
         spName = string.gsub(spName, '_plus', '+')
         spName = string.gsub(spName, '_minus', '-')
         names[#names+1] = spName
      end
      local nmodes = gm:nModes()
      FlowState.gm = gm
      FlowState.nSpecies = nsp
      FlowState.speciesNames = names
      FlowState.nModes = nmodes
      -- Patch a couple of the default values, to be consistent with the GasModel
      local massf = {[names[1]]=1.0, }
      for i = 2, nsp do massf[names[i]] = 0.0 end
      FlowState_defaults.massf = massf
      local T_modes = {}
      for i = 1, nmodes do T_modes[#T_modes+1] = FlowState_defaults.T end
      FlowState_defaults.T_modes = T_modes
   end
   -- Now, fill in default values for the FlowState object being constructed.
   -- If an item is not already present, copy the default value into the object,
   -- being careful to make new tables for massf and T_modes.
   for k, v in pairs(FlowState_defaults) do
      if o[k] == nil then
	 if k == "massf" then
	    o.massf = {}
	    for species, fraction in pairs(v) do o.massf[species] = fraction end
	 elseif k == "T_modes" then
	    o.T_modes = {}
	    for i,temperature in ipairs(v) do o.T_modes[i] = temperature end
	 else
	    o[k] = v
	 end
      end
   end
   -- On first use, we will not yet have a background gas state object.
   if FlowState.Q == nil then FlowState.Q = GasState:new{FlowState.gm} end
   -- Use the background GasState to compute some other variables.
   local Q = FlowState.Q
   Q.p = o.p
   Q.T = o.T
   Q.quality = o.quality;
   local massf = {}
   if FlowState.nSpecies > 1 then
      for k, v in pairs(o.massf) do massf[k] = v end
   else
      massf[FlowState.speciesNames[1]] = 1.0
   end
   Q.massf = massf
   if FlowState.nModes > 0 then
      if o.T_modes == nil then
	 -- We did not receive any T_modes, assume in equilibrium with T.
	 local T_modes = {}
	 for i = 1, FlowState.nModes do T_modes[#T_modes] = o.T end
	 o.T_modes = T_modes
      else
	 -- We have been given something, which is expected
	 -- to be either an array or a single number.
	 if type(o.T_modes) == "number" then
	    local my_array = {}
	    for i = 1, FlowState.nModes do my_array[#my_array+1] = o.T_modes end
	    o.T_modes = my_array
	 end
      end
      if not(#o.T_modes == FlowState.nModes) then
	 error(string.format("Wrong number of T_modes values: %d.", #o.T_modes), 2)
      end
      Q.T_modes = o.T_modes
   end
   local gm = FlowState.gm
   gm:updateThermoFromPT(Q)
   gm:updateSoundSpeed(Q)
   gm:updateTransCoeffs(Q)
   -- If we add entries here, be sure to keep in sync with list
   -- of validFlowStateFields in luaflowstate.d.
   o.rho = Q.rho
   o.a = Q.a
   o.mu = Q.mu
   o.k = Q.k
   setmetatable(o, self)
   self.__index = self
   return o
end

function FlowState:toJSONString()
   -- Since writing the JSON data is all about getting the values
   -- into the Dlang domain, we'll just delegate this work to the
   -- wrapped Dlang FlowState class.
   local fs = _FlowState:new{p=self.p, T=self.T, T_modes=self.T_modes, p_e=self.p_e,
			     quality=self.quality, massf=self.massf,
			     velx=self.velx, vely=self.vely, velz=self.velz,
			     Bx=self.Bx, By=self.By, Bz=self.Bz, psi=self.psi, divB=self.divB,
			     tke=self.tke, omega=self.omega, mu_t=self.mu_t, k_t=self.k_t,
			     S=self.S}
   return fs:toJSONString()
end

function FlowState:__tostring()
   local str = "FlowState{p=" .. tostring(self.p)
   str = str .. ", T=" .. tostring(self.T)
   str = str .. ", quality=" .. tostring(self.quality)
   str = str .. ", massf={"
   for k, v in pairs(self.massf) do
      str = str .. '["' ..tostring(k) .. '"]=' .. tostring(v) .. ', '
   end
   str = str .. "}"
   str = str .. ", T_modes={"
   for i, v in ipairs(self.T_modes) do
      str = str .. tostring(v) .. ", "
   end
   str = str .. "}"
   str = str .. ", p_e=" .. tostring(self.p_e)
   str = str .. ", velx=" .. tostring(self.velx)
   str = str .. ", vely=" .. tostring(self.vely)
   str = str .. ", velz=" .. tostring(self.velz)
   str = str .. ", Bx=" .. tostring(self.Bx)
   str = str .. ", By=" .. tostring(self.By)
   str = str .. ", Bz=" .. tostring(self.Bz)
   str = str .. ", psi=" .. tostring(self.psi)
   str = str .. ", divB=" .. tostring(self.divB)
   str = str .. ", tke=" .. tostring(self.tke)
   str = str .. ", omega=" .. tostring(self.omega)
   str = str .. ", mu_t=" .. tostring(self.mu_t)
   str = str .. ", k_t=" .. tostring(self.k_t)
   str = str .. ", S=" .. tostring(self.S)
   str = str .. "}"
   return str
end

-- ---------------------------------------------------------------------------------------

-- Storage for later definitions of FluidBlock objects.
-- Note that the index for this array starts at 1, in the Lua way.
-- The block ids start at 0 to be like the indexing inside the D code.
-- Yes, this is confusing...
fluidBlocks = {}
-- Storage for later definitions of FluidBlockArray objects.
fluidBlockArrays = {}
-- We may want to look up the blocks via labels rather than numerical id
-- in user-defined procedures.
-- The following dictionaries store the connections.
fluidBlocksDict = {}
fluidBlockArraysDict = {}

-- The user may assign the MPI task id for eack block manually
-- but, if they don't, a default distribution will be made.
mpiTasks = nil

-- Storgage for later definitions of SolidBlock objects
solidBlocks = {}

-- Storage for history cells
historyCells = {}
solidHistoryCells = {}

-- Storage for special zones
ignitionZones = {}
reactionZones = {}
turbulentZones = {}
suppressReconstructionZones = {}

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
      error("Expected a string or integer of three digits but got:"..tostring(gridpro_ijk))
   end
   local eilmer_ijk = axis_map[tonumber(string.sub(gridpro_ijk, 1, 1))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 2, 2))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 3, 3))]
   return eilmer_ijk
end

-- -----------------------------------------------------------------------

-- Class for gas dynamics FluidBlock construction.
FluidBlock = {
   myType = "FluidBlock",
} -- end FluidBlock

function FluidBlock:new(o)
   local flag = type(self)=='table' and self.myType=='FluidBlock'
   if not flag then
      error("Make sure that you are using FluidBlock:new{} and not FluidBlock.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"grid", "initialState", "fillCondition", "active",
                                "label", "omegaz", "bcList", "bcDict",
                                "hcellList", "xforceList", "fluidBlockArrayId"})
   if not flag then
      error("Invalid name for item supplied to FluidBlock constructor.", 2)
   end
   if o.initialState == nil then o.initialState = o.fillCondition end -- try old name
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block, for later construction of the config file.
   -- Note that we want block id to start at zero for the D code.
   o.id = #(fluidBlocks)
   fluidBlocks[#(fluidBlocks)+1] = o
   o.label = o.label or string.format("FluidBlock-%d", o.id)
   if fluidBlocksDict[o.label] then
      error('Have previously defined a FluidBlock with label "' .. o.label .. '"', 2)
   end
   fluidBlocksDict[o.label] = o.id
   -- Set to -1 if NOT part of a fluid-block-array, otherwise use supplied value
   o.fluidBlockArrayId = o.fluidBlockArrayId or -1
   -- Must have a grid and initialState
   assert(o.grid, "need to supply a grid")
   assert(o.initialState, "need to supply a initialState")
   if getmetatable(o.initialState) == FlowSolution then
      -- Let's build a initialState function here from a FlowSolution.
      o.initialState = makeFlowStateFn(o.initialState)
   end
   -- Fill in default values, if already not set
   if o.active == nil then
      o.active = true
   end
   o.omegaz = o.omegaz or 0.0
   if o.bcList then
      o.bcList = deepclone(o.bcList, false)
   else
      o.bcList = {}
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
   -- Check the grid information.
   if config.dimensions ~= o.grid:get_dimensions() then
      local msg = string.format("Mismatch in dimensions, config %d grid %d.",
				config.dimensions, o.grid.get_dimensions())
      error(msg)
   end
   if o.grid:get_type() == "structured_grid" then
      -- Extract some information from the StructuredGrid
      -- Note 0-based indexing for vertices and cells.
      o.nic = o.grid:get_niv() - 1
      o.njc = o.grid:get_njv() - 1
      if config.dimensions == 3 then
	 o.nkc = o.grid:get_nkv() - 1
      else
	 o.nkc = 1
      end
      o.ncells = o.nic * o.njc * o.nkc
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
      -- print("FluidBlock id=", o.id, "p0=", tostring(o.p[0]), "p1=", tostring(o.p[1]),
      --       "p2=", tostring(o.p[2]), "p3=", tostring(o.p[3]))
      -- Attach default boundary conditions for those not specified.
      for _,face in ipairs(faceList(config.dimensions)) do
	 o.bcList[face] = o.bcList[face] or WallBC_WithSlip:new()
      end
   end
   if o.grid:get_type() == "unstructured_grid" then
      -- Extract some information from the UnstructuredGrid
      o.ncells = o.grid:get_ncells()
      o.nvertices = o.grid:get_nvertices()
      o.nfaces = o.grid:get_nfaces()
      o.nboundaries = o.grid:get_nboundaries()
      -- Attach boundary conditions from list or from the dictionary of conditions.
      for i = 0, o.nboundaries-1 do
	 local mybc = o.bcList[i]
	 if (mybc == nil) and o.bcDict then
	    local tag = o.grid:get_boundaryset_tag(i)
	    mybc = o.bcDict[tag]
	 end
	 mybc = mybc or WallBC_WithSlip:new() -- default boundary condition
	 o.bcList[i] = mybc
      end
   end
   return o
end

function FluidBlock:tojson()
   -- [TODO] Should refactor the boundary condition checking and error messages.
   -- [TODO] Only the loops are different and we should error() rather than exit(), I think. PJ 2017-04-29
   local str = string.format('"block_%d": {\n', self.id)
   str = str .. string.format('    "type": "%s",\n', self.myType)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "fluidBlockArrayId": %d,\n', self.fluidBlockArrayId)
   str = str .. string.format('    "omegaz": %.18e,\n', self.omegaz)
   str = str .. string.format('    "grid_type": "%s",\n', self.grid:get_type())
   if self.grid:get_type() == "structured_grid" then
      str = str .. string.format('    "nic": %d,\n', self.nic)
      str = str .. string.format('    "njc": %d,\n', self.njc)
      str = str .. string.format('    "nkc": %d,\n', self.nkc)
      -- Boundary conditions for structured grid.
      for _,face in ipairs(faceList(config.dimensions)) do
	 if not self.bcList[face].is_gas_domain_bc then
	    local msg = string.format("Boundary condition problem for block:%d, face:%s\n", self.id, face)
	    msg = msg.."   This boundary condition should be a gas domain b.c.\n"
	    msg = msg.."   The preparation stage cannot complete successfully.\n"
	    error(msg)
	 end
	 if not self.bcList[face].is_configured then
	    local msg = string.format("Boundary condition problem for block:%d, face:%s\n", self.id, face)
	    msg = msg.."   This boundary condition was not configured correctly.\n"
	    msg = msg.."   If you used one of the standard boundary conditions,\n"
	    msg = msg.."   did you remember to call the b.c constructor as bcName:new{}?\n"
	    msg = msg.."   If you have custom configured the boundary condition,\n"
	    msg = msg.."   did you remember to set the 'is_configured' flag to true?\n"
	    error(msg)
	 end
	 str = str .. string.format('    "boundary_%s": ', face) ..
	    self.bcList[face]:tojson() .. ',\n'
      end
   end
   if self.grid:get_type() == "unstructured_grid" then
      str = str .. string.format('    "ncells": %d,\n', self.ncells)
      str = str .. string.format('    "nvertices": %d,\n', self.nvertices)
      str = str .. string.format('    "nfaces": %d,\n', self.nfaces)
      str = str .. string.format('    "nboundaries": %d,\n', self.nboundaries)
      -- Boundary conditions for the unstructured grid
      for i = 0, self.nboundaries-1 do
	 if not self.bcList[i].is_gas_domain_bc then
	    local msg = string.format("Boundary condition problem for block:%d, boundary:%d\n", self.id, i)
	    msg = msg.."   This boundary condition should be a gas domain b.c.\n"
	    msg = msg.."   The preparation stage cannot complete successfully.\n"
	    error(msg)
	 end
	 if not self.bcList[i].is_configured then
	    local msg = string.format("Boundary condition problem for block:%d, boundary:%d\n", self.id, i)
	    msg = msg.."   This boundary condition was not configured correctly.\n"
	    msg = msg.."   If you used one of the standard boundary conditions,\n"
	    msg = msg.."   did you remember to call the b.c constructor as bcName:new{}?\n"
	    msg = msg.."   If you have custom configured the boundary condition,\n"
	    msg = msg.."   did you remember to set the 'is_configured' flag to true?\n"
	    error(msg)
	 end
	 str = str .. string.format('    "boundary_%d": ', i) ..
	    self.bcList[i]:tojson() .. ',\n'
      end
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

-- ---------------------------------------------------------------------------
function SBlock2UBlock(blk)
   local origId = blk.id
   local origLabel = blk.label
   -- Let's swap out any exchange_over_full_face BCs and replace
   -- with exchange BCs.
   local bcList = {}
   for i,face in ipairs(faceList(config.dimensions)) do
      if blk.bcList[face].type == "exchange_over_full_face" then
	 -- We'll convert any exchange_over_full_face BCs
	 bcList[i-1] = ExchangeBC_MappedCell:new{}
      else
	 -- For all other BCs, we directly copy.
	 bcList[i-1] = blk.bcList[face]
      end
   end
   ublk = FluidBlock:new{grid=UnstructuredGrid:new{sgrid=blk.grid},
			 initialState=blk.initialState,
			 active=blk.active,
			 label=nil,
			 omegaz=blk.omegaz,
			 bcList=bcList}
   local newId = ublk.id
   local newLabel = ublk.label 
   -- Swap blocks in global list
   fluidBlocks[origId+1], fluidBlocks[newId+1] = fluidBlocks[newId+1], fluidBlocks[origId+1]
   -- Fix id and label of ublk
   fluidBlocks[origId+1].id = origId
   fluidBlocks[origId+1].label = origLabel
   -- Now remove original SFluidBlock, which is presently in pos ublk.id+1 and
   -- remove the latest dictionary entry, since the ublk no longer has that label.
   table.remove(fluidBlocks, newId+1)
   fluidBlocksDict[newLabel] = nil
end

function closeEnough(vA, vB, tolerance)
   -- Decide if two Vector quantities are close enough to being equal.
   -- This will be used to test that the block corners coincide.
   tolerance = tolerance or 1.0e-4
   return (vabs(vA - vB)/(vabs(vA + vB)+1.0)) <= tolerance
end

function connectBlocks(blkA, faceA, blkB, faceB, orientation)
   print("connectBlocks: blkA=", blkA.id, "faceA=", faceA, 
	 "blkB=", blkB.id, "faceB=", faceB, "orientation=", orientation)
   if blkA.grid:get_type() ~= "structured_grid" or blkB.grid:get_type() ~= "structured_grid" then
      error("connectBlocks() Works only for structured-grid blocks.", 2)
   end
   if blkA.myType == "FluidBlock" and blkB.myType == "FluidBlock" then
      blkA.bcList[faceA] = ExchangeBC_FullFace:new{otherBlock=blkB.id, otherFace=faceB,
						   orientation=orientation}
      blkB.bcList[faceB] = ExchangeBC_FullFace:new{otherBlock=blkA.id, otherFace=faceA,
						   orientation=orientation}
      -- [TODO] need to test for matching corner locations and consistent numbers of cells
   elseif blkA.myType == "FluidBlock" and blkB.myType == "SolidBlock" then
      -- Presently, only handle faceA == NORTH, faceB == SOUTH
      if faceA == north and faceB == south then
	 blkA.bcList[faceA] = WallBC_AdjacentToSolid:new{otherBlock=blkB.id,
							 otherFace=faceB,
							 orientation=orientation}
	 blkB.bcList[faceB] = SolidAdjacentToGasBC:new{otherBlock=blkA.id,
						       otherFace=faceA,
						       orientation=orientation}
      else
	 -- [TODO] Implement and handle other connection types.
	 local msg = "The requested FluidBlock to SolidBlock connection is not available.\n"
	 msg = msg .."FluidBlock-"..tostring(faceA).." :: SolidBlock-"..tostring(faceB)
	 error(msg, 2)
      end
   elseif blkA.myType == "SolidBlock" and blkB.myType == "FluidBlock" then
       -- Presently, only handle faceA == SOUTH, faceB == NORTH
      if faceA == south and faceB == north then
	 blkA.bcList[faceA] = SolidAdjacentToGasBC:new{otherBlock=blkB.id,
						       otherFace=faceB,
						       orientation=orientation}
	 blkB.bcList[faceB] = WallBC_AdjacentToSolid:new{otherBlock=blkA.id,
							 otherFace=faceA,
							 orientation=orientation}
      else
	 -- [TODO] Implement and handle other connection types.
	 local msg = "The requested SolidBlock to FluidBlock connection is not available.\n"
	 msg = msg.."SolidBlock-"..tostring(faceA).." :: FluidBlock-"..tostring(faceB)
	 error(msg, 2)
      end
   elseif blkA.myType == "SolidBlock" and blkB.myType == "SolidBlock" then
      -- Presently only handle EAST-WEST and WEST-EAST connections
      if ( (faceA == east and faceB == west) or ( faceA == west and faceB == east) ) then
	 blkA.bcList[faceA] = SolidConnectionBoundaryBC:new{otherBlock=blkB.id,
							    otherFace=faceB,
							    orientation=orientation}
	 blkB.bcList[faceB] = SolidConnectionBoundaryBC:new{otherBlock=blkA.id,
							    otherFace=faceA,
							    orientation=orientation}
      else
	 -- [TODO] Implement and handle other connection types for solid domains.
	 local msg = "The requested SolidBlock to SolidBlock connection is not available.\n"
	 msg = msg.."SolidBlock-"..tostring(faceA).." :: SolidBlock-"..tostring(faceB)
	 error(msg, 2)
      end
   end
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
      -- print("A.id=", A.id, "B.id=", B.id, "vtxPair=", tostringVtxPair(v)) -- DEBUG
      -- print("  A.p=", tostring(A.p[v[1]]), "B.p=", tostring(B.p[v[2]])) -- DEBUG
      if vabs(A.p[v[1]] - B.p[v[2]]) > tolerance then
	 allVerticesAreClose = false
      end
   end
   return allVerticesAreClose
end

function identifyBlockConnections(blockList, excludeList, tolerance)
   -- Identify block connections by trying to match corner points.
   -- Parameters:
   -- blockList: the list of SFluidBlock objects to be included in the search.
   --    If nil, the whole collection is searched.
   -- excludeList: list of pairs of SFluidBlock objects that should not be
   --    included in the search for connections.
   -- tolerance: spatial tolerance for the colocation of vertices
   local myBlockList = {}
   if ( blockList ) then
      for k,v in pairs(blockList) do myBlockList[k] = v end
   else -- Use the global gas blocks list
      for k,v in pairs(fluidBlocks) do myBlockList[k] = v end
   end
   -- Add solid domain block to search
   for _,blk in ipairs(solidBlocks) do myBlockList[#myBlockList+1] = blk end
   excludeList = excludeList or {}
   -- Put UFluidBlock objects into the exclude list because they don't
   -- have a simple topology that can always be matched to an SFluidBlock.
   for _,A in ipairs(myBlockList) do
      if A.grid:get_type() == "unstructured_grid" then excludeList[#excludeList+1] = A end
   end
   tolerance = tolerance or 1.0e-6
   for _,A in ipairs(myBlockList) do
      for _,B in ipairs(myBlockList) do
	 if (A ~= B) and (not isPairInList({A, B}, excludeList)) then
	    -- print("Proceed with test for coincident vertices.") -- DEBUG
	    local connectionCount = 0
	    if config.dimensions == 2 then
	       -- print("2D test A.id=", A.id, " B.id=", B.id) -- DEBUG
	       for vtxPairs,connection in pairs(connections2D) do
		  -- print("vtxPairs=", tostringVtxPairList(vtxPairs),
		  --       "connection=", tostringConnection(connection)) -- DEBUG
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

-- Class for FluidBlock-Array objects.
FBArray = {
   myType = "FBArray"
}

function FBArray:new(o)
   local flag = type(self)=='table' and self.myType=='FBArray'
   if not flag then
      error("Make sure that you are using FBArray:new{} and not FBArray.new{}", 2)
   end
   o = o or {}
   local flag = checkAllowedNames(o, {"grid", "initialState", "fillCondition",
				      "active", "label", "omegaz", "bcList",
				      "nib", "njb", "nkb"})
   if not flag then
      error("Invalid name for item supplied to FBArray constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- We will embed the FBArray identity in the individual blocks
   -- and we would like that identity to start from 0 for the D code.
   o.id = #(fluidBlockArrays)
   --
   o.label = o.label or string.format("FluidBlockArray-%d", o.id)
   if fluidBlockArraysDict[o.label] then
      error('Have previously defined a FBArray with label "' .. o.label .. '"', 2)
   end
   fluidBlockArraysDict[o.label] = o.id
   if not o.grid then
      error("You need to supply a grid to FBArray constructor.", 2)
   end
   if (not o.grid.get_type) or o.grid:get_type() ~= "structured_grid" then
      error("You need to supply a structured_grid to FBArray constructor.", 2)
   end
   if not o.initialState then
      -- try old name
      o.initialState = o.fillCondition
   end
   if not o.initialState then
      error("You need supply an initialState to FBArray constructor.", 2)
   end
   o.omegaz = o.omegaz or 0.0
   o.bcList = o.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or WallBC_WithSlip:new()
   end
   o.xforceList = o.xforceList or {}
   -- Numbers of subblocks in each coordinate direction
   o.nib = o.nib or 1
   o.njb = o.njb or 1
   o.nkb = o.nkb or 1
   if config.dimensions == 2 then
      o.nkb = 1
   end
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells in the D-domain.
   local nic_total = o.grid:get_niv() - 1
   local dnic = math.floor(nic_total/o.nib)
   local njc_total = o.grid:get_njv() - 1
   local dnjc = math.floor(njc_total/o.njb)
   local nkc_total = o.grid:get_nkv() - 1
   local dnkc = math.floor(nkc_total/o.nkb)
   if config.dimensions == 2 then
      nkc_total = 1
      dnkc = 1
   end
   o.blockArray = {} -- will be a multi-dimensional array indexed as [i][j][k]
   o.blockCollection = {} -- will be a single-dimensional array
   local nic_remaining = nic_total
   local i0 = 0
   for ib = 1, o.nib do
      o.blockArray[ib] = {}
      local nic = math.floor(nic_remaining/(o.nib-ib+1))
      if (ib == o.nib) then
         -- On last block, just use what's left
         nic = nic_remaining
      end
      nic_remaining = nic_remaining - nic
      local njc_remaining = njc_total
      local j0 = 0
      for jb = 1, o.njb do
         local njc = math.floor(njc_remaining/(o.njb-jb+1))
	 if (jb == o.njb) then
	    njc = njc_remaining
	 end
         njc_remaining = njc_remaining - njc
	 if config.dimensions == 2 then
	    -- 2D flow
            print("ib=", ib, "jb= ", jb)
            print("i0= ", i0, " nic= ", nic, " j0= ", j0, " njc= ", njc)
	    local subgrid = o.grid:subgrid(i0,nic+1,j0,njc+1)
	    local bcList = {north=WallBC_WithSlip:new(), east=WallBC_WithSlip:new(),
			    south=WallBC_WithSlip:new(), west=WallBC_WithSlip:new()}
	    if ib == 1 then
	       bcList[west] = o.bcList[west]
	    end
	    if ib == o.nib then
	       bcList[east] = o.bcList[east]
	    end
	    if jb == 1 then
	       bcList[south] = o.bcList[south]
	    end
	    if jb == o.njb then
	       bcList[north] = o.bcList[north]
	    end
	    local new_block = FluidBlock:new{grid=subgrid, omegaz=o.omegaz,
                                             initialState=o.initialState,
                                             bcList=bcList,
                                             fluidBlockArrayId=o.id}
	    o.blockArray[ib][jb] = new_block
	    o.blockCollection[#o.blockCollection+1] = new_block
	 else
	    -- 3D flow, need one more level in the array
	    o.blockArray[ib][jb] = {}
            local nkc_remaining = nkc_total
            local k0 = 0
	    for kb = 1, o.nkb do
               local nkc = math.floor(nkc_remaining/(o.nkb-kb+1))
               if (kb == o.nkb) then
                  nkc = nkc_remaining
               end
               nkc_remaining = nkc_remaining - nkc
	       local subgrid = o.grid:subgrid(i0,nic+1,j0,njc+1,k0,nkc+1)
	       local bcList = {north=WallBC_WithSlip:new(), east=WallBC_WithSlip:new(),
			       south=WallBC_WithSlip:new(), west=WallBC_WithSlip:new(),
			       top=WallBC_WithSlip:new(), bottom=WallBC_WithSlip:new()}
	       if ib == 1 then
		  bcList[west] = o.bcList[west]
	       end
	       if ib == o.nib then
		  bcList[east] = o.bcList[east]
	       end
	       if jb == 1 then
		  bcList[south] = o.bcList[south]
	       end
	       if jb == o.njb then
		  bcList[north] = o.bcList[north]
	       end
	       if kb == 1 then
		  bcList[bottom] = o.bcList[bottom]
	       end
	       if kb == o.nkb then
		  bcList[top] = o.bcList[top]
	       end
	       local new_block = FluidBlock:new{grid=subgrid, omegaz=o.omegaz,
                                                initialState=o.initialState,
                                                bcList=bcList,
                                                fluidBlockArrayId=o.id,}
	       o.blockArray[ib][jb][kb] = new_block
	       o.blockCollection[#o.blockCollection+1] = new_block
               -- Prepare k0 at end of loop, ready for next iteration
               k0 = k0 + nkc
	    end -- kb loop
	 end -- dimensions
         -- Prepare j0 at end of loop, ready for next iteration
         j0 = j0 + njc
      end -- jb loop
      -- Prepare i0 at end of loop, ready for next iteration
      i0 = i0 + nic
   end -- ib loop
   -- Make the inter-subblock connections
   if #o.blockCollection > 1 then
      identifyBlockConnections(o.blockCollection)
   end
   --
   -- Retain meta-information about the new FluidBlockArray
   -- for use later in the user-defined functions, during simulation.
   -- Note that the index of this array starts at 1 (in the Lua way).
   fluidBlockArrays[#fluidBlockArrays+1] = o
   --
   return o
end -- FBArray:new

-- Retain the original behaviour.
function FluidBlockArray(t)
   o = FBArray:new(t)
   return o.blockArray
end

function FBArray:tojson()
   local str = string.format('"fluid_block_array_%d": {\n', self.id)
   str = str .. string.format('    "nib": %d,\n', self.nib)
   str = str .. string.format('    "njb": %d,\n', self.njb)
   str = str .. string.format('    "nkb": %d,\n', self.nkb)
   str = str .. string.format('    "blockIds": [ ')
   for ib=1,#(self.blockCollection)-1 do
      str = str .. string.format('%d, ', self.blockCollection[ib].id)
   end
   str = str .. string.format('%d ]\n', self.blockCollection[#self.blockCollection].id)
   str = str .. '},\n'
   return str
end

-- Class for SolidBlock construction
SolidBlock = {
   myType = "SolidBlock",
} -- end SSolidBlock

function SolidBlock:new(o)
   local flag = type(self)=='table' and self.myType=='SolidBlock'
   if not flag then
      error("Make sure that you are using SolidBlock:new{} and not SolidBlock.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"grid", "initTemperature", "active",
                                "label", "bcList", "properties"})
   if not flag then
      error("Invalid name for item supplied to SolidBlock constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block for later construction of the config file.
   -- Note that we want the block id to start at zero for the D code.
   print("Adding new solid block.")
   o.id = #solidBlocks
   solidBlocks[#solidBlocks+1] = o
   -- Must have a grid and initial temperature
   if not o.grid then
      error("You need to supply a grid to SolidBlock constructor.", 2)
   end
   if (not o.grid.get_type) or o.grid:get_type() ~= "structured_grid" then
      error("You need to supply a structured_grid to SolidBlock constructor.", 2)
   end
   if not o.initTemperature then
      error("You need to supply an initTemperature to SolidBlock constructor.", 2)
   end
   if not o.properties then
      error("You need to supply physical properties for the block.", 2)
   end
   if (type(o.properties) == 'table') then
      local flag2 = checkAllowedNames(o.properties, {"rho", "k", "Cp",
						     "k11", "k12", "k13",
						     "k21", "k22", "k23",
						     "k31", "k32", "k33"})
      if not flag2 then
         error("Invalid name for item supplied in SolidBlock properties table.", 2)
      end
      -- Fill in the k values as 0.0 if not set.
      local kProps = {"k", "k11", "k12", "k13", "k21", "k22", "k23", "k31", "k32", "k33"}
      for _,kName in ipairs(kProps) do
	 o.properties[kName] = o.properties[kName] or 0.0
      end
   end
   -- Fill in some defaults, if not already set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("SOLIDBLOCK-%d", o.id)
   if o.bcList then
      o.bcList = deepclone(o.bcList, false)
   else
      o.bcList =  {} -- boundary conditions
   end
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or SolidAdiabaticBC:new{}
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

function SolidBlock:tojson()
   local str = string.format('"solid_block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   -- Boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      if not self.bcList[face].is_solid_domain_bc then
	 local errMsg = string.format("Boundary condition problem for solid block:%d, face:%s\n", self.id, face)
	 errMsg = errMsg .. "       This boundary condition should be a solid domain b.c.\n"
	 errMsg = errMsg .. "       The preparation stage cannot complete successfully.\n"
	 error(errMsg)
      end
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

function SolidBlockArray(t)
   -- Expect one table as argument, with named fields.
   -- Returns an array of blocks defined over a single region.
   local flag = checkAllowedNames(t, {"grid", "initTemperature", "active",
				      "label", "bcList", "properties",
				      "nib", "njb", "nkb"})
   if not flag then
      error("Invalid name for item supplied to SolidBlockArray.", 2)
   end
   if not t.grid then
      error("You need to supply a grid to SolidBlockArray.", 2)
   end
   if (not t.grid.get_type) or t.grid:get_type() ~= "structured_grid" then
      error("You need to supply a structured_grid to SolidBlockArray.", 2)
   end
   if not t.initTemperature then
      error("You need to supply an 'initTemperature' to SolidBlockArray.", 2)
   end
   if not t.properties then
      error("You need to supply 'properties' to SolidBlockArray.", 2)
   end
   t.bcList = t.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      t.bcList[face] = t.bcList[face] or SolidAdiabaticBC:new{}
   end
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
	    local bcList = {north=SolidAdiabaticBC:new{}, east=SolidAdiabaticBC:new{},
			    south=SolidAdiabaticBC:new{}, west=SolidAdiabaticBC:new{}}
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
	    local new_block = SolidBlock:new{grid=subgrid, properties=t.properties,
                                             initTemperature=t.initTemperature,
                                             bcList=bcList}
	    blockArray[ib][jb] = new_block
	    blockCollection[#blockCollection+1] = new_block
	 else
	    error("SolidBlockArray not implemented for 3D.")
	 end -- dimensions
      end -- jb loop
   end -- ib loop
   -- Make the inter-subblock connections
   if #blockCollection > 1 then
      identifyBlockConnections(blockCollection)
   end
   return blockArray
end -- SolidBlockArray


function setHistoryPoint(args)
   -- Accepts a variety of arguments:
   --  1. x, y, z coordinates
   --  setHistoryPoint{x=7.9, y=8.2, z=0.0}
   --  2. block and single-index for cell
   --  setHistoryPoint{ib=2, i=102}
   --  3. block and structured grid indices
   --  setHistoryPoint{ib=0, i=20, j=10, k=0}
   --
   local flag = checkAllowedNames(args, {"x", "y", "z", "ib", "i", "j", "k"})
   if not flag then
      error("Invalid name for item supplied to setHistoryPoint.", 2)
   end
   -- First look for x,y,z
   local found = false
   if args.x then
      local x = args.x
      local y = args.y
      local z = args.z or 0.0
      local minDist = 1.0e9 -- something very large
      local blkId = 0
      local cellId = 0
      for ib,blk in ipairs(fluidBlocks) do
	 local indx, dist = blk.grid:find_nearest_cell_centre{x=x, y=y, z=z}
	 if (dist < minDist) then
	    minDist = dist
	    blkId = ib
	    cellId = indx
	 end
      end
      -- Convert blkId to 0-offset
      blkId = blkId - 1
      historyCells[#historyCells+1] = {ib=blkId, i=cellId}
      found = true
   end
   -- Still trying; look for integer indices.
   if (not found) and args.j and args.ib then
      local ib = args.ib
      local i = args.i
      local j = args.j
      local k = args.k or 0
      -- Convert back to single_index
      local nic = fluidBlocks[ib+1].nic
      local njc = fluidBlocks[ib+1].njc
      local cellId = k * (njc * nic) + j * nic + i
      historyCells[#historyCells+1] = {ib=args.ib, i=cellId}
      found = true
   end
   if not found then
      -- Final option is to directly set the identity block and cell.
      if args.ib and args.i then
         historyCells[#historyCells+1] = {ib=args.ib, i=args.i}
      else
         error("Could not identify cell for setHistoryPoint.", 2)
      end
   end
   -- If we arrive here, we have successfully set the cell identity.
   -- Print a summary of its identity and location.
   local n = #historyCells
   local ib = historyCells[n].ib -- zero-based block index
   local i = historyCells[n].i
   local pos = fluidBlocks[ib+1].grid:cellCentroid(i) -- one-based block array
   print(string.format("Fluid History point [%d] ib=%d i=%d x=%g y=%g z=%g",
                       n, ib, i, pos.x, pos.y, pos.z))
   return
end

function setSolidHistoryPoint(args)
   -- Accepts a variety of arguments:
   --  1. x, y, z coordinates
   --  setSolidHistoryPoint{x=7.9, y=8.2, z=0.0}
   --  2. block and single-index for cell
   --  setSolidHistoryPoint{ib=2, i=102}
   --  3. block and structured grid indices
   --  setSolidHistoryPoint{ib=0, i=20, j=10, k=0}
   --
   local flag = checkAllowedNames(args, {"x", "y", "z", "ib", "i", "j", "k"})
   if not flag then
      error("Invalid name for item supplied to setSolidHistoryPoint.", 2)
   end
   -- First look for x,y,z
   local found = false
   if args.x then
      local x = args.x
      local y = args.y
      local z = args.z or 0.0
      local minDist = 1.0e9 -- something very large
      local blkId = 0
      local cellId = 0
      for ib,blk in ipairs(solidBlocks) do
	 local indx, dist = blk.grid:find_nearest_cell_centre{x=x, y=y, z=z}
	 if (dist < minDist) then
	    minDist = dist
	    blkId = ib
	    cellId = indx
	 end
      end
      -- Convert blkId to 0-offset
      blkId = blkId - 1
      solidHistoryCells[#solidHistoryCells+1] = {ib=blkId, i=cellId}
      found = true
   end
   -- Still trying; look for integer indices.
   if (not found) and args.j and args.ib then
      local ib = args.ib
      local i = args.i
      local j = args.j
      local k = args.k or 0
      -- Convert back to single_index
      local nic = solidBlocks[ib+1].nic
      local njc = solidBlocks[ib+1].njc
      local cellId = k * (njc * nic) + j * nic + i
      solidHistoryCells[#solidHistoryCells+1] = {ib=args.ib, i=cellId}
      found = true
   end
   if not found then
      -- Final option is to directly set the identity block and cell.
      if args.ib and args.i then
         solidHistoryCells[#solidHistoryCells+1] = {ib=args.ib, i=args.i}
      else
         error("Could not identify cell for setSolidHistoryPoint.", 2)
      end
   end
   -- If we arrive here, we have successfully set the cell identity.
   -- Print a summary of its identity and location.
   local n = #solidHistoryCells
   local ib = solidHistoryCells[n].ib -- zero-based block index
   local i = solidHistoryCells[n].i
   local pos = solidBlocks[ib+1].grid:cellCentroid(i) -- one-based block array
   print(string.format("Solid History point [%d] ib=%d i=%d x=%g y=%g z=%g",
                       n, ib, i, pos.x, pos.y, pos.z))
   return
end


function makeFlowStateFn(flowSol)
   function flowFn(x, y, z)
      -- We try to find an enclosing cell.
      local cell = flowSol:find_enclosing_cell{x=x, y=y, z=z}
      -- If that fails, we'll just grab a 'nearest' cell.
      -- This should never fail.
      if not cell.ib then
	 cell = flowSol:find_nearest_cell_centre{x=x, y=y, z=z}
      end
      cell.fmt = "FlowState"
      -- The table for a cell's data should be enough to set the FlowState.
      return flowSol:get_cell_data(cell) 
   end
   return flowFn
end


function checkCellVolumes(t)
   if not t then
      t = {}
   end
   -- Stop reporting cells after this limit
   local badCellLimit = t.badCellLimit or 20
   local badCells = {}
   local badCellCount = 0
   for ib,blk in ipairs(fluidBlocks) do
      local grid = blk.grid
      for idx=0,grid:get_ncells()-1 do
	 local vol = grid:cellVolume(idx)
	 if vol <= 0 then
	    badCellCount = badCellCount + 1
	    badCells[#badCells+1] = {ib, idx}
	    if badCellCount >= badCellLimit then
	       return false, badCells
	    end
	 end
      end
   end
   if #badCells > 0 then
      return false, badCells
   else
      return true, badCells
   end
end


function mpiDistributeBlocks(args)
   -- Assign blocks to MPI tasks,
   -- keeping a record of the MPI rank (or task) for every block.
   -- This record is stored in the global variable mpiTasks.
   --
   if args and not(type(args) == "table") then
      error("mpiDistributeBlocks expects its arguments in single table with named fields", 2);
   end
   args = args or {}
   local flag = checkAllowedNames(args, {"ntasks", "dist", "preassign"})
   if not flag then
      error("Invalid name for item supplied to mpiDistributeBlocks.", 2)
   end
   --
   local nBlocks = #fluidBlocks
   -- If nTasks is not given, assume that we want one per block.
   local nTasks = (args and args.ntasks) or nBlocks
   nTasks = math.min(nTasks, nBlocks)
   --
   local option = (args and args.dist) or "round-robin"
   --
   -- The following list will eventually hold the MPI-rank for each FluidBlock.
   local mpiTaskList = {}; for i=1,nBlocks do mpiTaskList[i] = -1 end
   --
   -- Preassigned blocks.
   -- Entries in the table are of the form blockId:mpiRank
   -- Remember that blockId and mpiRank values start at zero.
   if args and (type(args.preassign) == "table") then
      for blkId,mpirank in pairs(args.preassign) do mpiTaskList[blkId+1] = mpirank end
   end
   --
   -- Distribute the rest of the blocks.
   if option == "uniform" or option == "round-robin" or option == "roundrobin" then
      -- For each block, assign to an MPI task, round-robin order.
      local mpirank = 0
      for i=1,nBlocks do
         if mpiTaskList[i] < 0 then
            -- This block not previously assigned a rank.
            mpiTaskList[i] = mpirank
            mpirank = mpirank + 1
            if mpirank == nTasks then mpirank = 0 end
         end
      end
   elseif option == "loadbalance" or option == "load-balance" then
      -- Load-balance procedure first sorts the blocks by size...
      local blksNcells = {}
      local totalCells = 0
      for i=1,nBlocks do
         local blk = fluidBlocks[i]
         blksNcells[i] = {i, blk.ncells}
         totalCells = totalCells + blk.ncells
      end
      table.sort(blksNcells, function (a,b) return a[2] > b[2] end)
      -- ...then distributes the blocks to the tasks,
      -- biggest block first into the task with the smallest load.
      -- We shall tally the loads, in number of cells, for each MPI task.
      local taskLoads = {}; for i=1,nTasks do taskLoads[i] = 0 end
      -- Account for the preassigned blocks.
      for ib=1,nBlocks do
         mpirank = mpiTaskList[ib]
         if mpirank >= 0 then
            taskLoads[mpirank+1] = taskLoads[mpirank+1] + fluidBlocks[ib].ncells
         end
      end
      -- Distribute remaining blocks.
      for _,v in pairs(blksNcells) do
         local ib = v[1]; local ncells = v[2]
         if mpiTaskList[ib] < 0 then
            -- Add the so-far-unassigned block to the MPI task with smallest load.
            local indxSmallest = 1; local smallestLoad = taskLoads[1]
            for i=2,#taskLoads do
               if taskLoads[i] < smallestLoad then
                  indxSmallest = i; smallestLoad = taskLoads[i]
               end
            end
            mpiTaskList[ib] = indxSmallest-1 -- MPI task ids start from zero
            taskLoads[indxSmallest] = taskLoads[indxSmallest] + ncells
         end
      end
      --
      local maxmpiLoads = (math.max(unpack(taskLoads)))
      local minmpiLoads = (math.min(unpack(taskLoads)))
      local mpiProcessors = math.max(unpack(mpiTaskList)) + 1
      print("Load balancing - Distribute blocks to CPUs") 
      print(string.format("Number of processors   \t \t = %d", mpiProcessors))
      print(string.format("Ideal cell partitioning   \t = %d cells/proc", totalCells/mpiProcessors))
      print(string.format("Smallest partition factor \t = %.3f", minmpiLoads/(totalCells/mpiProcessors)))
      print(string.format("Largest partition factor  \t = %.3f", maxmpiLoads/(totalCells/mpiProcessors)))
      print(string.format("Largest processor load    \t = %.0f cells", maxmpiLoads))

   else
      error('Did not select one of "round-robin" or "load-balance". for mpiDistributeBlocks', 2) 
   end
   -- Assign the newly-constructed list to the global variable
   -- for later use in writing the job.mpimap file.
   mpiTasks = mpiTaskList
   -- Finally, return the list as we have always done, however,
   -- we expect that the caller will ignore this return value.
   return mpiTaskList
end
   
-- -----------------------------------------------------------------------

-- Classes for construction of zones.

ReactionZone = {
   p0 = nil,
   p1 = nil
}

function ReactionZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to ReactionZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(reactionZones)
   reactionZones[#(reactionZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

IgnitionZone = {
   p0 = nil,
   p1 = nil,
   T = nil -- degrees K
}

function IgnitionZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1", "T"})
   if not flag then
      error("Invalid name for item supplied to IgnitionZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(ignitionZones)
   ignitionZones[#(ignitionZones)+1] = o
   -- Must have corners and temperature
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   if not o.T then
      error("You need to supply ignition temperature T", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

TurbulentZone = {
   p0 = nil,
   p1 = nil,
}

function TurbulentZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to TurbulentZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(turbulentZones)
   turbulentZones[#(turbulentZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

SuppressReconstructionZone = {
   p0 = nil,
   p1 = nil,
}

function SuppressReconstructionZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to SuppressReconstructionZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(suppressReconstructionZones)
   suppressReconstructionZones[#(suppressReconstructionZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

-- --------------------------------------------------------------------

function write_control_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"dt_init": %.18e,\n', config.dt_init))
   f:write(string.format('"dt_max": %.18e,\n', config.dt_max))
   f:write(string.format('"cfl_value": %.18e,\n', config.cfl_value))
   f:write(string.format('"stringent_cfl": %s,\n', tostring(config.stringent_cfl)))
   f:write(string.format('"viscous_signal_factor": %.18e,\n', config.viscous_signal_factor))
   f:write(string.format('"turbulent_signal_factor": %.18e,\n', config.turbulent_signal_factor))
   f:write(string.format('"residual_smoothing_weight": %.18e,\n', config.residual_smoothing_weight))
   f:write(string.format('"residual_smoothing_iterations": %d,\n', config.residual_smoothing_iterations))
   f:write(string.format('"residual_smoothing_type": "%s",\n',
			 config.residual_smoothing_type))
   f:write(string.format('"fixed_time_step": %s,\n', tostring(config.fixed_time_step)))
   f:write(string.format('"print_count": %d,\n', config.print_count))
   f:write(string.format('"cfl_count": %d,\n', config.cfl_count))
   f:write(string.format('"max_time": %.18e,\n', config.max_time))
   f:write(string.format('"max_step": %d,\n', config.max_step))
   f:write(string.format('"dt_plot": %.18e,\n', config.dt_plot))
   f:write(string.format('"dt_history": %.18e,\n', config.dt_history))
   f:write(string.format('"dt_loads": %.18e,\n', config.dt_loads))
   f:write(string.format('"halt_now": %d,\n', config.halt_now))
   f:write('"steady_state_solver_options" : {\n')
   f:write(string.format('   "use_preconditioner": %s,\n', tostring(SteadyStateSolver.use_preconditioner)))
   f:write(string.format('   "frozen_preconditioner_count": %d,\n', SteadyStateSolver.frozen_preconditioner_count))
   f:write(string.format('   "start_preconditioning": %d,\n', SteadyStateSolver.start_preconditioning))
   f:write(string.format('   "ilu_fill": %d,\n', SteadyStateSolver.ilu_fill))
   f:write(string.format('   "precondition_matrix_type": "%s",\n', SteadyStateSolver.precondition_matrix_type))
   f:write(string.format('   "use_scaling": %s,\n', tostring(SteadyStateSolver.use_scaling)))
   f:write(string.format('   "use_complex_matvec_eval": %s,\n', tostring(SteadyStateSolver.use_complex_matvec_eval)))
   f:write(string.format('   "number_pre_steps": %d,\n', SteadyStateSolver.number_pre_steps))
   f:write(string.format('   "number_total_steps": %d,\n', SteadyStateSolver.number_total_steps))
   f:write(string.format('   "max_number_attempts": %d,\n', SteadyStateSolver.max_number_attempts))
   f:write(string.format('   "stop_on_relative_global_residual": %.18e,\n', SteadyStateSolver.stop_on_relative_global_residual))
   f:write(string.format('   "stop_on_absolute_global_residual": %.18e,\n', SteadyStateSolver.stop_on_absolute_global_residual))
   f:write(string.format('   "max_outer_iterations": %d,\n', SteadyStateSolver.max_outer_iterations))
   f:write(string.format('   "max_restarts": %d,\n', SteadyStateSolver.max_restarts))
   f:write(string.format('   "number_inner_iterations": %d,\n', SteadyStateSolver.number_inner_iterations))
   f:write(string.format('   "number_start_up_steps": %d,\n', SteadyStateSolver.number_start_up_steps))
   f:write(string.format('   "cfl0": %.18e,\n', SteadyStateSolver.cfl0))
   f:write(string.format('   "eta0": %.18e,\n', SteadyStateSolver.eta0))
   f:write(string.format('   "tau0": %.18e,\n', SteadyStateSolver.tau0))
   f:write(string.format('   "sigma0": %.18e,\n', SteadyStateSolver.sigma0))
   f:write(string.format('   "p0": %.18e,\n', SteadyStateSolver.p0))
   f:write(string.format('   "cfl1": %.18e,\n', SteadyStateSolver.cfl1))
   f:write(string.format('   "tau1": %.18e,\n', SteadyStateSolver.tau1))
   f:write(string.format('   "sigma1": %.18e,\n', SteadyStateSolver.sigma1))
   f:write(string.format('   "p1": %.18e,\n', SteadyStateSolver.p1))
   f:write(string.format('   "eta_strategy": "%s",\n', SteadyStateSolver.eta_strategy))
   f:write(string.format('   "eta1": %.18e,\n', SteadyStateSolver.eta1))
   f:write(string.format('   "eta1_max": %.18e,\n', SteadyStateSolver.eta1_max))
   f:write(string.format('   "eta1_min": %.18e,\n', SteadyStateSolver.eta1_min))
   f:write(string.format('   "eta_ratio_per_step": %.18e,\n', SteadyStateSolver.eta_ratio_per_step))
   f:write(string.format('   "gamma": %.18e,\n', SteadyStateSolver.gamma))
   f:write(string.format('   "alpha": %.18e,\n', SteadyStateSolver.alpha))
   f:write(string.format('   "snapshots_count": %d,\n', SteadyStateSolver.snapshots_count))
   f:write(string.format('   "number_total_snapshots": %d,\n', SteadyStateSolver.number_total_snapshots))
   f:write(string.format('   "write_diagnostics_count": %d,\n', SteadyStateSolver.write_diagnostics_count))
   f:write(string.format('   "write_loads_count": %d\n', SteadyStateSolver.write_loads_count))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write('    }\n')
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write("}\n")
   f:close()
end

function write_config_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"title": "%s",\n', config.title))
   f:write(string.format('"start_time": %.18e,\n', config.start_time))
   f:write(string.format('"grid_format": "%s",\n', config.grid_format))
   f:write(string.format('"flow_format": "%s",\n', config.flow_format))
   f:write(string.format('"gas_model_file": "%s",\n', config.gas_model_file))
   f:write(string.format('"udf_supervisor_file": "%s",\n', tostring(config.udf_supervisor_file)))
   if type(user_pad_data) == 'table' then
      if config.user_pad_length < #user_pad_data then
         config.user_pad_length = #user_pad_data
      end
   end
   f:write(string.format('"user_pad_length": %d,\n', config.user_pad_length))
   f:write('"user_pad_data": [')
   if type(user_pad_data) == 'table' then
      for i,e in ipairs(user_pad_data) do
         f:write(string.format('%.18e', e))
         if i < #user_pad_data then f:write(', ') end
      end
   end
   f:write('],\n')
   f:write(string.format('"sticky_electrons": %s,\n',
			 tostring(config.sticky_electrons)))
   f:write(string.format('"include_quality": %s,\n',
			 tostring(config.include_quality)))
   f:write(string.format('"dimensions": %d,\n', config.dimensions))
   f:write(string.format('"axisymmetric": %s,\n',
			 tostring(config.axisymmetric)))
   f:write(string.format('"strang_splitting": "%s",\n', config.strang_splitting))
   f:write(string.format('"gasdynamic_update_scheme": "%s",\n',
			 config.gasdynamic_update_scheme))
   f:write(string.format('"residual_smoothing": %s,\n', tostring(config.residual_smoothing)))
   f:write(string.format('"with_local_time_stepping": %s,\n', tostring(config.with_local_time_stepping)))
   f:write(string.format('"local_time_stepping_limit_factor": %d,\n', tostring(config.local_time_stepping_limit_factor)))
   f:write(string.format('"with_super_time_stepping": %s,\n', tostring(config.with_super_time_stepping)))
   f:write(string.format('"coupling_with_solid_domains": "%s",\n',
			 config.coupling_with_solid_domains))
   f:write(string.format('"solid_has_isotropic_properties": %s,\n', tostring(config.solid_has_isotropic_properties)))
   f:write(string.format('"solid_has_homogeneous_properties": %s,\n', tostring(config.solid_has_homogeneous_properties)))
   f:write('"solid_domain_loose_update_options" : {\n')
   f:write(string.format('   "max_newton_iterations" : %d,\n', SolidDomainLooseUpdate.max_newton_iterations))
   f:write(string.format('   "tolerance_newton_update" : %.18e,\n', SolidDomainLooseUpdate.tolerance_newton_update))
   f:write(string.format('   "max_gmres_iterations" : %d,\n', SolidDomainLooseUpdate.max_gmres_iterations))
   f:write(string.format('   "tolerance_gmres_solve" : %.18e,\n', SolidDomainLooseUpdate.tolerance_gmres_solve))
   f:write(string.format('   "perturbation_size" : %.18e \n', SolidDomainLooseUpdate.perturbation_size))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write('},\n')
   f:write(string.format('"MHD": %s,\n', tostring(config.MHD)))
   f:write(string.format('"MHD_static_field": %s,\n', tostring(config.MHD_static_field)))
   f:write(string.format('"MHD_resistive": %s,\n', tostring(config.MHD_resistive)))
   f:write(string.format('"divergence_cleaning": %s,\n', tostring(config.divergence_cleaning)))
   f:write(string.format('"divB_damping_length": %.18e,\n', config.divB_damping_length))
   f:write(string.format('"apply_bcs_in_parallel": %s,\n',
			 tostring(config.apply_bcs_in_parallel)))
   f:write(string.format('"flowstate_limits_max_velocity": %.18e,\n', config.flowstate_limits_max_velocity))
   f:write(string.format('"flowstate_limits_max_tke": %.18e,\n', config.flowstate_limits_max_tke))
   f:write(string.format('"flowstate_limits_min_tke": %.18e,\n', config.flowstate_limits_min_tke))
   f:write(string.format('"flowstate_limits_max_temp": %.18e,\n', config.flowstate_limits_max_temp))
   f:write(string.format('"flowstate_limits_min_temp": %.18e,\n', config.flowstate_limits_min_temp))
   f:write(string.format('"max_invalid_cells": %d,\n', config.max_invalid_cells))
   f:write(string.format('"adjust_invalid_cell_data": %s,\n', tostring(config.adjust_invalid_cell_data)))
   f:write(string.format('"report_invalid_cells": %s,\n', tostring(config.report_invalid_cells)))

   f:write(string.format('"high_order_flux_calculator": %s,\n', tostring(config.high_order_flux_calculator)))
   f:write(string.format('"flux_calculator": "%s",\n', config.flux_calculator))
   f:write(string.format('"interpolation_order": %d,\n', config.interpolation_order))
   f:write(string.format('"interpolation_delay": %.18e,\n', config.interpolation_delay))
   f:write(string.format('"suppress_radial_reconstruction_at_xaxis": %s,\n',
                         tostring(config.suppress_radial_reconstruction_at_xaxis)))
   f:write(string.format('"thermo_interpolator": "%s",\n', 
			 string.lower(config.thermo_interpolator)))
   f:write(string.format('"allow_reconstruction_for_energy_modes": %s,\n', 
			 tostring(config.allow_reconstruction_for_energy_modes)))
   f:write(string.format('"interpolate_in_local_frame": %s,\n', 
			 tostring(config.interpolate_in_local_frame)))
   f:write(string.format('"apply_limiter": %s,\n', tostring(config.apply_limiter)))
   f:write(string.format('"extrema_clipping": %s,\n', tostring(config.extrema_clipping)))
   f:write(string.format('"unstructured_limiter": "%s",\n', config.unstructured_limiter))
   f:write(string.format('"freeze_limiter_on_step": %d,\n', config.freeze_limiter_on_step))
   f:write(string.format('"use_extended_stencil": %s,\n', tostring(config.use_extended_stencil)))
   f:write(string.format('"venkat_K_value": %.18e,\n', config.venkat_K_value))
   f:write(string.format('"compression_tolerance": %.18e,\n', config.compression_tolerance))
   f:write(string.format('"shear_tolerance": %.18e,\n', config.shear_tolerance))
   f:write(string.format('"M_inf": %.18e,\n', config.M_inf))
   f:write(string.format('"artificial_compressibility": %s,\n', tostring(config.artificial_compressibility)))
   f:write(string.format('"ac_alpha": %.18e,\n', config.ac_alpha))

   f:write(string.format('"grid_motion": "%s",\n', tostring(config.grid_motion)))
   f:write(string.format('"write_vertex_velocities": %s,\n', tostring(config.write_vertex_velocities)))
   f:write(string.format('"udf_grid_motion_file": "%s",\n', tostring(config.udf_grid_motion_file)))
   
   f:write(string.format('"shock_fitting_delay": %.18e,\n', config.shock_fitting_delay))
   f:write(string.format('"shock_fitting_interpolation_order": %d,\n', config.shock_fitting_interpolation_order))
   f:write(string.format('"shock_fitting_scale_factor": %.18e,\n', config.shock_fitting_scale_factor))

   f:write(string.format('"viscous": %s,\n', tostring(config.viscous)))
   f:write(string.format('"use_viscosity_from_cells": %s,\n', tostring(config.use_viscosity_from_cells)))
   f:write(string.format('"spatial_deriv_from_many_points": %s,\n', tostring(config.spatial_deriv_from_many_points)))
   f:write(string.format('"spatial_deriv_calc": "%s",\n', config.spatial_deriv_calc))
   f:write(string.format('"spatial_deriv_locn": "%s",\n', config.spatial_deriv_locn))
   f:write(string.format('"include_ghost_cells_in_spatial_deriv_clouds": %s,\n',
			 tostring(config.include_ghost_cells_in_spatial_deriv_clouds)))
   f:write(string.format('"suppress_reconstruction_at_boundaries": %s,\n',
			 tostring(config.suppress_reconstruction_at_boundaries)))
   f:write(string.format('"suppress_reconstruction_at_captured_shocks": %s,\n',
			 tostring(config.suppress_reconstruction_at_captured_shocks)))
   f:write(string.format('"viscous_delay": %.18e,\n', config.viscous_delay))
   f:write(string.format('"shear_stress_relative_limit": %.18e,\n', config.shear_stress_relative_limit))
   f:write(string.format('"apply_shear_stress_relative_limit": %s,\n', tostring(config.apply_shear_stress_relative_limit)))
   f:write(string.format('"mass_diffusion_model": "%s",\n',
			 string.lower(config.mass_diffusion_model)))
   f:write(string.format('"constant_lewis_number": %s,\n', tostring(config.constant_lewis_number)))
   f:write(string.format('"species_specific_lewis_numbers": %s,\n', tostring(config.species_specific_lewis_numbers)))
   f:write(string.format('"lewis_number": %.18e,\n', config.lewis_number))

   f:write(string.format('"separate_update_for_viscous_terms": %s,\n',
			 tostring(config.separate_update_for_viscous_terms)))
   f:write(string.format('"separate_update_for_k_omega_source": %s,\n', 
			 tostring(config.separate_update_for_k_omega_source)))

   f:write(string.format('"turbulence_model": "%s",\n',
			 string.lower(config.turbulence_model)))
   f:write(string.format('"turbulence_prandtl_number": %.18e,\n',
			 config.turbulence_prandtl_number))
   f:write(string.format('"turbulence_schmidt_number": %.18e,\n',
			 config.turbulence_schmidt_number))
   f:write(string.format('"max_mu_t_factor": %.18e,\n', config.max_mu_t_factor))
   f:write(string.format('"transient_mu_t_factor": %.18e,\n', config.transient_mu_t_factor))
   f:write(string.format('"limit_tke_production": %s,\n', tostring(config.limit_tke_production)))
   f:write(string.format('"tke_production_limit_in_kelvins": %.18e,\n', config.tke_production_limit_in_kelvins))

   f:write(string.format('"udf_source_terms_file": "%s",\n', config.udf_source_terms_file))
   f:write(string.format('"udf_source_terms": %s,\n', tostring(config.udf_source_terms)))

   f:write(string.format('"reacting": %s,\n', tostring(config.reacting)))
   f:write(string.format('"reactions_file": "%s",\n', config.reactions_file))
   f:write(string.format('"reaction_time_delay": %.18e,\n', config.reaction_time_delay))
   f:write(string.format('"T_frozen": %.18e,\n', config.T_frozen))
   f:write(string.format('"T_frozen_energy": %.18e,\n', config.T_frozen_energy))
   f:write(string.format('"tci_model": "%s",\n', string.lower(config.tci_model)))
   f:write(string.format('"ignition_time_start": %.18e,\n', config.ignition_time_start))
   f:write(string.format('"ignition_time_stop": %.18e,\n', config.ignition_time_stop))
   f:write(string.format('"energy_exchange_file": "%s",\n', config.energy_exchange_file))

   f:write(string.format('"control_count": %d,\n', config.control_count))
   f:write(string.format('"nfluidblock": %d,\n', #(fluidBlocks)))
   f:write(string.format('"nfluidblockarrays": %d,\n', #(fluidBlockArrays)))

   f:write(string.format('"diffuse_wall_bcs_on_init": %s,\n', tostring(config.diffuse_wall_bcs_on_init)))
   f:write(string.format('"number_init_passes": %d,\n', config.number_init_passes))
   f:write(string.format('"wall_temperature_on_init": %.18e,\n', config.wall_temperature_on_init));

   f:write(string.format('"thermionic_emission_bc_time_delay": %.18e,\n', config.thermionic_emission_bc_time_delay))

   f:write('"shape_sensitivity_calculator_options" : {\n')
   f:write(string.format('   "pseudotime": %s,\n', tostring(ShapeSensitivityCalculator.pseudotime)))
   f:write(string.format('   "pseudotime_lhs_jacobian_order": %d,\n', ShapeSensitivityCalculator.pseudotime_lhs_jacobian_order))
   f:write(string.format('   "adjoint_precondition_matrix_order": %d,\n', ShapeSensitivityCalculator.adjoint_precondition_matrix_order))
   f:write(string.format('   "read_frozen_limiter_values_from_file": %s,\n', tostring(ShapeSensitivityCalculator.read_frozen_limiter_values_from_file)))
   f:write(string.format('   "epsilon": %.18e,\n', ShapeSensitivityCalculator.epsilon))
   f:write(string.format('   "maxOuterIterations": %d,\n', ShapeSensitivityCalculator.maxOuterIterations))
   f:write(string.format('   "maxRestarts": %d,\n', ShapeSensitivityCalculator.maxRestarts))
   f:write(string.format('   "cfl0": %.18e,\n', ShapeSensitivityCalculator.cfl0))
   f:write(string.format('   "eta": %.18e,\n', ShapeSensitivityCalculator.eta))
   f:write(string.format('   "stop_on_relative_global_residual": %.18e,\n', ShapeSensitivityCalculator.stop_on_relative_global_residual))
   f:write(string.format('   "tol_bezier_curve_fit": %.18e,\n', ShapeSensitivityCalculator.tol_bezier_curve_fit))
   f:write(string.format('   "max_steps_bezier_curve_fit": %d,\n', ShapeSensitivityCalculator.max_steps_bezier_curve_fit))
   f:write(string.format('   "user_defined_objective_file": "%s"\n', tostring(ShapeSensitivityCalculator.user_defined_objective_file)))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write('},\n')
   
   f:write(string.format('"block_marching": %s,\n',
			 tostring(config.block_marching)))
   f:write(string.format('"nib": %d,\n', config.nib))
   f:write(string.format('"njb": %d,\n', config.njb))
   f:write(string.format('"nkb": %d,\n', config.nkb))
   f:write(string.format('"propagate_inflow_data": %s,\n',
			 tostring(config.propagate_inflow_data)))
   f:write(string.format('"save_intermediate_results": %s,\n',
			 tostring(config.save_intermediate_results)))
   f:write(string.format('"boundary_group_for_loads": "%s",\n',
			 config.boundary_group_for_loads))
   f:write(string.format('"write_loads": %s,\n',
			 tostring(config.write_loads)))
   f:write(string.format('"compute_run_time_loads": %s,\n',
			 tostring(config.compute_run_time_loads)))
   f:write(string.format('"run_time_loads_count": %d,\n', config.run_time_loads_count))
   if run_time_loads and (#run_time_loads > 0) then
      f:write(string.format('"run_time_loads": { "ngroups" : %d, \n', #run_time_loads))
      for i,t in ipairs(run_time_loads) do
         f:write(string.format('   "group-%d" : { "groupLabel" : "%s",\n',  i-1, run_time_loads[i].group))
         f:write(string.format('       "momentCtr_x" : %.18e,\n',  run_time_loads[i].moment_centre.x))
         f:write(string.format('       "momentCtr_y" : %.18e,\n',  run_time_loads[i].moment_centre.y))
         f:write(string.format('       "momentCtr_z" : %.18e },\n',  run_time_loads[i].moment_centre.z))
      end
      f:write('   "dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
      f:write('},\n')
   end
   f:write(string.format('"nhcell": %d,\n', #historyCells))
   for i,hcell in ipairs(historyCells) do
      f:write(string.format('"history-cell-%d": [%d, %d],\n', i-1, hcell.ib, hcell.i))
   end
   f:write(string.format('"nsolidhcell": %d,\n', #solidHistoryCells))
   for i,hcell in ipairs(solidHistoryCells) do
      f:write(string.format('"solid-history-cell-%d": [%d, %d],\n', i-1, hcell.ib, hcell.i))
   end

   f:write(string.format('"n-reaction-zones": %d,\n', #reactionZones))
   for i,zone in ipairs(reactionZones) do
      f:write(string.format('"reaction-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end
   f:write(string.format('"n-ignition-zones": %d,\n', #ignitionZones))
   for i,zone in ipairs(ignitionZones) do
      f:write(string.format('"ignition-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z, zone.T))
   end
   f:write(string.format('"n-turbulent-zones": %d,\n', #turbulentZones))
   for i,zone in ipairs(turbulentZones) do
      f:write(string.format('"turbulent-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end
   f:write(string.format('"n-suppress-reconstruction-zones": %d,\n', #suppressReconstructionZones))
   for i,zone in ipairs(suppressReconstructionZones) do
      f:write(string.format('"suppress-reconstruction-zone-%d": [%.18e, %.18e, %.18e, %.18e, %.18e, %.18e],\n',
			    i-1, zone.p0.x, zone.p0.y, zone.p0.z,
			    zone.p1.x, zone.p1.y, zone.p1.z))
   end

   f:write(string.format('"udf_solid_source_terms_file": "%s",\n', config.udf_solid_source_terms_file))
   f:write(string.format('"udf_solid_source_terms": %s,\n', tostring(config.udf_solid_source_terms)))
   f:write(string.format('"nsolidblock": %d,\n', #solidBlocks))

   for i = 1, #fluidBlockArrays do
      f:write(fluidBlockArrays[i]:tojson())
   end
   for i = 1, #fluidBlocks do
      f:write(fluidBlocks[i]:tojson())
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
   f:write(string.format("%04d %.18e %.18e\n", 0, config.start_time, config.dt_init))
   f:close()
end

function write_block_list_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# indx type label ncells\n")
   for i = 1, #(fluidBlocks) do
      local blk = fluidBlocks[i]
      f:write(string.format("%4d %s %s %d\n", blk.id,
                            blk.grid:get_type(), blk.label, blk.ncells))
   end
   f:close()
end

function write_mpimap_file(fileName)
   if not mpiTasks then
      mpiDistributeBlocks()
   end
   local f = assert(io.open(fileName, "w"))
   f:write("# indx mpiTask\n")
   for i = 1, #(fluidBlocks) do
      f:write(string.format("%4d %4d\n", fluidBlocks[i].id, mpiTasks[i]))
   end
   f:close()
end

function write_fluidBlockArrays_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("-- A description of the fluidBlockArrays in Lua code.\n")
   f:write("-- Use dofile() to get the content into your interpreter.\n")
   f:write("fluidBlockArrays = {\n")
   for i = 1, #(fluidBlockArrays) do
      local fba = fluidBlockArrays[i]
      f:write(string.format("  [%d]={label=\"%s\",\n", fba.id, fba.label))
      f:write(string.format("    nib=%d, njb=%d, nkb=%d,\n", fba.nib, fba.njb, fba.nkb))
      local blkId = 0
      if config.dimensions == 3 then
         blkId = fba.blockArray[1][1][1].id
      else
         blkId = fba.blockArray[1][1].id
      end
      f:write(string.format("    p00=function() return infoFluidBlock(%d).p00 end,\n", blkId));
      if config.dimensions == 3 then
         blkId = fba.blockArray[fba.nib][1][1].id
      else
         blkId = fba.blockArray[fba.nib][1].id
      end
      f:write(string.format("    p10=function() return infoFluidBlock(%d).p10 end,\n", blkId));
      if config.dimensions == 3 then
         blkId = fba.blockArray[fba.nib][fba.njb][1].id
      else
         blkId = fba.blockArray[fba.nib][fba.njb].id
      end
      f:write(string.format("    p11=function() return infoFluidBlock(%d).p11 end,\n", blkId));
      if config.dimensions == 3 then
         blkId = fba.blockArray[1][fba.njb][1].id
      else
         blkId = fba.blockArray[1][fba.njb].id
      end
      f:write(string.format("    p01=function() return infoFluidBlock(%d).p01 end,\n", blkId));
      --
      f:write("    blockCollection={")
      for _,blk in ipairs(fba.blockCollection) do
         f:write(string.format("%d, ", blk.id))
      end
      f:write("}, -- end blockCollection\n")
      --
      f:write("    blockArray={\n")
      for ib,itable in pairs(fba.blockArray) do
         f:write(string.format("      [%d]={", ib))
         for jb,jitem in pairs(itable) do
            f:write(string.format("[%d]={", jb))
            if config.dimensions == 3 then
               for kb,blk in pairs(jitem) do
                  f:write(string.format("[%d]=%d, ", kb, blk.id))
               end
            else
               -- For 2D
               f:write(string.format("[1]=%d" , jitem.id))
            end
            f:write("}, ") -- end [ib][jb]
         end
         f:write("}, -- end [ib]\n")
      end
      f:write("    }, -- end blockArray\n")
      f:write("  },\n")
   end
   f:write("} -- end fluidBlockArrays\n")
   --
   f:write("-- Dictionary of FluidBlockArrays\n")
   f:write("fluidBlockArraysDict = {\n")
   for label,id in pairs(fluidBlockArraysDict) do
      f:write(string.format(' ["%s"]=%d,\n', label, id))
   end
   f:write("} -- end fluidBlockArraysDict\n")
   --
   f:write("-- Map from FluidBlock id to FluidBlockArray id\n")
   f:write("whichFluidBlockArray = {\n")
   for i,blk in ipairs(fluidBlocks) do
      if blk.fluidBlockArrayId >= 0 then
         f:write(string.format(" [%d]=%d,", blk.id, blk.fluidBlockArrayId))
      end
      if (i > 1 and i%5 == 0) or i == #fluidBlocks then f:write("\n") end
   end
   f:write("} -- end whichFluidBlockArray\n")
   --
   f:write("-- Map from FluidBlock id to FluidBlockArray label\n")
   f:write("whichFluidBlockArrayLabel = {\n")
   for i,blk in ipairs(fluidBlocks) do
      if blk.fluidBlockArrayId >= 0 then
         fba = fluidBlockArrays[blk.fluidBlockArrayId+1]
         f:write(string.format(' [%d]="%s",', blk.id, fba.label))
      end
      if (i > 1 and i%5 == 0) or i == #fluidBlocks then f:write("\n") end
   end
   f:write("} -- end whichFluidBlockArrayLabel\n")
   --
   f:write("-- Dictionary of FluidBlocks\n")
   f:write("fluidBlocksDict = {\n")
   for label,id in pairs(fluidBlocksDict) do
      f:write(string.format(' ["%s"]=%d,\n', label, id))
   end
   f:write("} -- end fluidBlocksDict\n")
   --
   f:write([[
function is_in_FluidBlockArray(blkId, label)
   -- Returns true if a block is in a particular FluidBlockArray.
   local myfba = whichFluidBlockArrayLabel[blkId]
   if not myfba then
      return false
   end
   if label == myfba then
      return true
   else
      return false
   end
end
]])
   --
   f:close()
end -- function write_fluidBlockArrays_file

function perform_spatial_gradient_consistency_check()
   -- Not all spatial gradient options are available, depending on the type of grid.
   -- First, search for any unstructured grids, since these are the most restricted.
   unstructuredGridsPresent = false
   for _,blk in ipairs(fluidBlocks) do
      if blk.grid:get_type() == "unstructured_grid" then
	 unstructuredGridsPresent = true
	 break
      end
   end
   if unstructuredGridsPresent then
      if config.spatial_deriv_calc == "divergence" then
	 print("NOTE: config.spatial_deriv_calc is being set to 'least_squares' because unstructured grids detected.")
	 config.spatial_deriv_calc = "least_squares"
      end
      if config.spatial_deriv_locn == "vertices" then
	 print("NOTE: config.spatial_deriv_locn is being set to 'cells' when using least squares.")
	 config.spatial_deriv_locn = "cells"
      end
   else
      -- Only structured grids are present.
      -- 2D structured grids have all options available.
      if config.dimensions == 3 then
         if config.spatial_deriv_calc == "divergence" then
            print("NOTE: config.spatial_deriv_calc is being set to 'least_squares' for 3D simulations.")
            config.spatial_deriv_calc = "least_squares"
         end
      end
   end
   if not config.spatial_deriv_from_many_points then
      if config.spatial_deriv_locn == "vertices" then
         print("NOTE: config.spatial_deriv_location is being set to 'faces' for 2-point derivatives.")
         config.spatial_deriv_locn = "faces"
      end
   end
end

function build_job_files(job)
   if #fluidBlocksForPrep == 0 then
      -- We'll set *all* blocks for processing.
      for i=1,#fluidBlocks do
         fluidBlocksForPrep[i] = fluidBlocks[i].id
      end
   end
   perform_spatial_gradient_consistency_check()
   if buildMasterFiles then
      print("Build job files for ", job)
      os.execute("mkdir -p config")
      write_config_file("config/" .. job .. ".config")
      write_control_file("config/" .. job .. ".control")
      write_times_file("config/" .. job .. ".times")
      write_block_list_file("config/" .. job .. ".list")
      write_mpimap_file("config/" .. job .. ".mpimap")
      write_fluidBlockArrays_file("config/" .. job .. ".fluidBlockArrays")
   end
   os.execute("mkdir -p grid/t0000")
   os.execute("mkdir -p flow/t0000")
   if #solidBlocks >= 1 then
      os.execute("mkdir -p solid-grid/t0000")
      os.execute("mkdir -p solid/t0000")
   end
   for i, id in ipairs(fluidBlocksForPrep) do
      print("FluidBlock id=", id)
      local idx = id+1
      local fileName = "grid/t0000/" .. job .. string.format(".grid.b%04d.t0000", id)
      if (config.grid_format == "gziptext") then
	 fluidBlocks[idx].grid:write_to_gzip_file(fileName .. ".gz")
      elseif (config.grid_format == "rawbinary") then
	 fluidBlocks[idx].grid:write_to_raw_binary_file(fileName .. ".bin")
      else
	 error(string.format("Oops, invalid grid_format: %s", config.grid_format))
      end
      local fileName = "flow/t0000/" .. job .. string.format(".flow.b%04d.t0000", id)
      if (config.flow_format == "gziptext") then
	 fileName = fileName .. ".gz"
      elseif (config.flow_format == "rawbinary") then
	 fileName = fileName .. ".bin"
      else
	 error(string.format("Oops, invalid flow_format: %s", config.flow_format))
      end
      local ifs = fluidBlocks[idx].initialState
      if type(ifs) == "table" and ifs.myType == "FlowState" then
	 -- We have one of the pure-Lua FlowState objects and we convert it to
	 -- a wrapped-D-language _FlowState object.
	 ifs = _FlowState:new(ifs)
      elseif type(ifs) == "function" then
	 -- leave alone
      elseif type(ifs) == "userdata" then
	 -- presume to be a wrapped-D-language _FlowState object already
      elseif type(ifs) == "string" then
         -- We are given the name of a flow file and we'll copy that in place directly
         existingFlowFile = ifs
         -- Presume file exists and let 'cp' command complain if it doesn't
         cmd = "cp " .. existingFlowFile .. " " .. fileName
         returnCode = os.execute(cmd)
         if returnCode ~= 0 then
            errMsg = "Error while trying to copy an existing flow solution as initial flow solution.\n"
            errMsg = errMsg .. "FluidBlock id= " .. id .. "\n"
            errMsg = errMsg .. "Specified existing flow file: " .. existingFlowFile .. "\n"
            errMsg = errMsg .. "Check this file exists and is readable.\n"
            errMsg = errMsg .. "Bailing out!\n"
            error(errMsg)
         end
         -- Otherwise succesful.
         str = string.format("Initialised FluidBlock id= %d with existing flow solution: \n\t%s", id, existingFlowFile)
         print(str)
      else
	 error("Unexpected type for initial flow state in block.")
      end
      if type(ifs) ~= "string" then
         local grid = fluidBlocks[idx].grid
         if grid:get_type() == "structured_grid" then
            write_initial_sg_flow_file(fileName, grid, ifs, config.start_time)
         else
            write_initial_usg_flow_file(fileName, grid, ifs, config.start_time)
         end
      end
   end
   for i = 1, #solidBlocks do
      local id = solidBlocks[i].id
      print("SolidBlock id=", id)
      local fileName = "solid-grid/t0000/" .. job .. string.format(".solid-grid.b%04d.t0000.gz", id)
      solidBlocks[i].grid:write_to_gzip_file(fileName)
      local fileName = "solid/t0000/" .. job .. string.format(".solid.b%04d.t0000", id)
      writeInitialSolidFile(fileName, solidBlocks[i].grid,
			    solidBlocks[i].initTemperature, solidBlocks[i].properties, config.start_time)
      os.execute("gzip -f " .. fileName)
   end
   --
   if #fluidBlocks == 0 then print("Warning: number of FluidBlocks is zero.") end
   print("Done building files.")
end


if false then
   -- Keep old names available, for now.
   -- Once we purge all of the old names from the examples,
   -- we should delete this code block
   SBlock = FluidBlock
   UBlock = FluidBlock
   SBlockArray = FluidBlockArray
   SSolidBlock = SolidBlock
   SSolidBlockArray = SolidBlockArray
end

print("Done loading prep.lua")
