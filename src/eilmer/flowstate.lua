-- Module for constructing FlowState objects, required by prep.lua.
--
-- Authors: PJ and RJG
--

-- A class to build FlowState objects as Lua tables.

local FlowState_defaults = {
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
   nuhat = 0.0, -- SA turbulent variable, (m^^2/s)
   mu_t = 0.0, -- turbulence viscosity, Pa.s
   k_t = 0.0, -- turbulence conductivity
   S = 0, -- shock indicator
}

-- Prototype for the FlowState class.
local FlowState = {
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
                                 "tke", "omega", "nuhat", "mu_t", "k_t", "S"})
   if not flag then
      error("Invalid name for item supplied to FlowState constructor.", 2)
   end
   -- We want to make consistent values/calculations in the context
   -- of an initialized gas model, so go fetch that model from GlobalConfig.
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
      -- Patch the massf default values, to be consistent with the GasModel
      local massf = {[names[1]]=1.0, }
      for i = 2, nsp do massf[names[i]] = 0.0 end
      FlowState_defaults.massf = massf
      -- Set T_modes to nil deliberately in defaults as we want user
      -- to set this.
      FlowState_defaults.T_modes = nil
   end
   -- Now, fill in default values for the FlowState object being constructed.
   -- If an item is not already present, copy the default value into the object,
   -- being careful to make new table for massf.
   -- RJG, 2019-11-03: We do not make defaults for T_modes.
   -- We'll require the user to set this explicitly.
   for k, v in pairs(FlowState_defaults) do
      if o[k] == nil then
	 if k == "massf" then
	    o.massf = {}
	    for species, fraction in pairs(v) do o.massf[species] = fraction end
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
         error(string.format("No values for T_modes supplied, but n_modes= %d.", FlowState.nModes))
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
   o.turb = initTurbulence(o, config.turbulence_model)
   setmetatable(o, self)
   self.__index = self
   return o
end

function FlowState:toJSONString()
   -- Since writing the JSON data is all about getting the values
   -- into the Dlang domain, we'll just delegate this work to the
   -- wrapped Dlang FlowState class.
   local fs = _FlowState.new{p=self.p, T=self.T, T_modes=self.T_modes, p_e=self.p_e,
			     quality=self.quality, massf=self.massf,
			     velx=self.velx, vely=self.vely, velz=self.velz,
			     Bx=self.Bx, By=self.By, Bz=self.Bz, psi=self.psi, divB=self.divB,
			     turb=self.turb, mu_t=self.mu_t, k_t=self.k_t,
			     S=self.S}
   return _FlowState.toJSONString(fs)
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
   if self.T_modes then
      str = str .. ", T_modes={"
      for i, v in ipairs(self.T_modes) do
         str = str .. tostring(v) .. ", "
      end
      str = str .. "}"
   end
   str = str .. ", p_e=" .. tostring(self.p_e)
   str = str .. ", velx=" .. tostring(self.velx)
   str = str .. ", vely=" .. tostring(self.vely)
   str = str .. ", velz=" .. tostring(self.velz)
   str = str .. ", Bx=" .. tostring(self.Bx)
   str = str .. ", By=" .. tostring(self.By)
   str = str .. ", Bz=" .. tostring(self.Bz)
   str = str .. ", psi=" .. tostring(self.psi)
   str = str .. ", divB=" .. tostring(self.divB)
   str = str .. ", turb={"
   for k, v in pairs(self.turb) do
      str = str .. tostring(v) .. ', '
   end
   str = str .. "}"
   str = str .. ", mu_t=" .. tostring(self.mu_t)
   str = str .. ", k_t=" .. tostring(self.k_t)
   str = str .. ", S=" .. tostring(self.S)
   str = str .. "}"
   return str
end

local function makeFlowStateFn(flowSol)
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

return {
   FlowState = FlowState,
   makeFlowStateFn = makeFlowStateFn
}
