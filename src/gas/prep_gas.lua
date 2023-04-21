#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Date: 08-Mar-2015
--
-- This program is used to prepare a detailed
-- gas model file for use by the eilmer4 gas
-- dynamics simulation code. It requires a fairly
-- minimal input from the user for commonly used
-- gas models. However, it also allows customisation
-- for advanced users.
--

CIDBFileName = {
   gupta = "gupta_etal_1990_CI_data.lua",
   wright = "wright_etal_CI_data.lua"
}


function split_string(str)
   tokens = {}
   for tk in string.gmatch(str, "%S+") do
      tokens[#tokens+1] = tk
   end
   return tokens
end

function writeIdealGas(f, sp, db, optsTable)
   -- We can only deal with one species for an ideal gas.
   if ( #sp > 1 ) then
      print("WARNING: More than one species is listed while trying to prepare")
      print("WARNING: an ideal gas model.")
      print("WARNING: Presently, the ideal gas model is limited to a single species")
      print("WARNING: We will build a file with the first species listed and ignore the rest.")
   end
   local s = sp[1]
   f:write("IdealGas = {\n")
   f:write(string.format("   speciesName = '%s',\n", s))
   f:write(string.format("   mMass = %.8f,\n", db[s].M.value))
   f:write(string.format("   gamma = %.8f,\n", db[s].gamma.value))
   if (db[s].entropyRefValues) then
      f:write("   entropyRefValues = {\n")
      f:write(string.format("      s1 = %.8e,\n", db[s].entropyRefValues.s1))
      f:write(string.format("      T1 = %.8f,\n", db[s].entropyRefValues.T1))
      f:write(string.format("      p1 = %.8e,\n", db[s].entropyRefValues.p1))
      f:write("   },\n")
   else
      -- Use some default value.
      -- This should NOT be an issue. We only need correct enthalpy
      -- values for computing equilibrium constants in a reacting
      -- gas. We should not be doing reacting flows with a mixture
      -- of ideal gases.
      f:write("   entropyRefValues = {\n")
      f:write(string.format("      s1 = %.8e,\n", db.default.entropyRefValues.s1))
      f:write(string.format("      T1 = %.8f,\n", db.default.entropyRefValues.T1))
      f:write(string.format("      p1 = %.8e,\n", db.default.entropyRefValues.p1))
      f:write("   },\n")
   end
   -- Our default choice for an ideal gas is a simple Sutherland viscosity
   -- model, but there are two reasons we might NOT use that.
   -- 1. We override that choice using an option in the options table.
   -- 2. The Sutherland viscosity model is not available for the particular
   --    species, so we fall back to CEA coefficients.
   if (optsTable and optsTable.viscosity == "CEA") or (db[s].sutherlandVisc == nil) then
      f:write("   viscosity = {\n")
      f:write("      model = 'CEA',\n")
      secName = "ceaViscosity"
      t = db[s][secName]
      f:write(string.format("      nsegments = %d,\n", t.nsegments))
      for i=0,t.nsegments-1 do
	 seg = "segment"..i
	 f:write(string.format("      segment%d = {\n", i))
	 f:write(string.format("         T_lower = %.1f,\n", t[seg].T_lower))
	 f:write(string.format("         T_upper = %.1f,\n", t[seg].T_upper))
	 f:write(string.format("         A = % -10.7e,\n", t[seg].A))
	 f:write(string.format("         B = % -10.7e,\n", t[seg].B))
	 f:write(string.format("         C = % -10.7e,\n", t[seg].C))
	 f:write(string.format("         D = % -10.7e,\n", t[seg].D))
	 f:write("      },\n")
      end
      f:write("   },\n")
   else
      -- We will use Sutherland viscosity
      f:write("   viscosity = {\n")
      f:write("      model = 'Sutherland',\n");
      f:write(string.format("      mu_ref = %.8e,\n", db[s].sutherlandVisc.mu_ref))
      f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandVisc.T_ref))
      f:write(string.format("      S = %.8f,\n", db[s].sutherlandVisc.S))
      f:write("   },\n")
   end
   -- As with viscosity, we'll prefer Sutherland over CEA, unless we can't find Sutherland
   -- values for this species.
   f:write("   thermCondModel = {\n")
   if optsTable and optsTable.Prandtl then
      f:write("      model = 'constPrandtl',\n")
      f:write(string.format("      Prandtl = %.8f\n", optsTable.Prandtl))
   elseif (optsTable and optsTable.thermal_conductivity == "CEA") or (db[s].sutherlandThermCond == nil) then
      f:write("      model = 'CEA',\n")
      secName = "ceaThermCond"
      t = db[s][secName]
      f:write(string.format("      nsegments = %d,\n", t.nsegments))
      for i=0,t.nsegments-1 do
	 seg = "segment"..i
	 f:write(string.format("      segment%d = {\n", i))
	 f:write(string.format("         T_lower = %.1f,\n", t[seg].T_lower))
	 f:write(string.format("         T_upper = %.1f,\n", t[seg].T_upper))
	 f:write(string.format("         A = % -10.7e,\n", t[seg].A))
	 f:write(string.format("         B = % -10.7e,\n", t[seg].B))
	 f:write(string.format("         C = % -10.7e,\n", t[seg].C))
	 f:write(string.format("         D = % -10.7e,\n", t[seg].D))
	 f:write("      },\n")
      end
   else
      -- We will use Sutherland expression
      f:write("      model = 'Sutherland',\n")
      f:write(string.format("      k_ref = %.8e,\n", db[s].sutherlandThermCond.k_ref))
      f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandThermCond.T_ref))
      f:write(string.format("      S = %.8f,\n", db[s].sutherlandThermCond.S))
   end
   f:write("   }\n")
   f:write("}\n")
end

function writeCeaThermoCoeffs(f, sp, db, optsTable)
   local t
   if optsTable and optsTable.database == "prefer-grimech" then
      if ( not db[sp].grimechThermoCoeffs ) then
	 print("WARNING: GRIMECH thermo coefficients have been selected as preferred,")
	 print("WARNING: but they could not be found for species: ", sp)
	 print("WARNING: We default back to using the CEA thermo coefficients (if available).")
         print("")
      else
	 t = db[sp].grimechThermoCoeffs
	 t.origin = "GRIMECH"
      end
   end

   if not t then
      if not db[sp].ceaThermoCoeffs  then
         print("ERROR: The table of CEA coefficients to compute thermodynamic properties")
         print("ERROR: could not be found for species: ", sp)
         print("ERROR: Bailing out!")
         os.exit(1)
      else
         t = db[sp].ceaThermoCoeffs
         t.origin = "CEA"
      end
   end
   f:write(string.format("db['%s'].thermoCoeffs = {\n", sp))
   f:write(string.format("  origin = '%s',\n", t.origin))
   f:write(string.format("  nsegments = %d, \n", t.nsegments))
   f:write("  T_break_points = { ")
   for _,T in ipairs(t.T_break_points) do
      f:write(string.format("%.2f, ", T))
   end
   f:write("},\n")
   f:write("  T_blend_ranges = { ")
   for _,T in ipairs(t.T_blend_ranges) do
      f:write(string.format("%.1f, ", T))
   end
   f:write("},\n")
   for i=0,t.nsegments-1 do
      seg = "segment"..i
      f:write(string.format("  segment%d = {\n", i))
      for _,c in ipairs(t[seg]) do
	 f:write(string.format("   % -12.9e,\n", c))
      end
      f:write("  },\n")
   end
   f:write("}\n")
end

function writeCeaTransCoeffs(f, sp, db, name, key)
   secName = "cea"..name
   if ( not db[sp][secName] ) then
      print("")
      print("------------------------------------------------------------------------------------------")
      print(string.format("WARNING: The table of CEA coefficients to compute diffusion property: %s", secName))
      print("WARNING: could not be found for species: ", sp)
      print("WARNING: As a substitute the values from 'defaults.lua' will be used.")
      print("------------------------------------------------------------------------------------------")
      db[sp][secName] = db.default[secName]
   end
   t = db[sp][secName]
   f:write(string.format("db['%s'].%s = {\n", sp, key))
   f:write("   model = 'CEA',\n")
   f:write(string.format("   nsegments = %d,\n", t.nsegments))
   for i=0,t.nsegments-1 do
      seg = "segment"..i
      f:write(string.format("   segment%d = {\n", i))
      f:write(string.format("      T_lower = %.1f,\n", t[seg].T_lower))
      f:write(string.format("      T_upper = %.1f,\n", t[seg].T_upper))
      f:write(string.format("      A = % -10.7e,\n", t[seg].A))
      f:write(string.format("      B = % -10.7e,\n", t[seg].B))
      f:write(string.format("      C = % -10.7e,\n", t[seg].C))
      f:write(string.format("      D = % -10.7e,\n", t[seg].D))
      f:write("   },\n")
   end
   f:write("}\n")
end

function writeChemkinTransCoeffs(f, sp, db, name, key)
   secName = "chemkin"..name
   if ( not db[sp][secName] ) then
      print("")
      print("------------------------------------------------------------------------------------------")
      print(string.format("WARNING: The table of Chemkin coefficients to compute diffusion property: %s", secName))
      print("WARNING: could not be found for species: ", sp)
      print("WARNING: As a substitute the CEA model and values from 'defaults.lua' will be used.")
      print("------------------------------------------------------------------------------------------")
      db[sp][secName] = db.default[secName]
   end
   t = db[sp][secName]
   f:write(string.format("db['%s'].%s = {\n", sp, key))

   if ( not db[sp][secName] ) then
      f:write("   model = 'CEA',\n")
   else
      f:write("   model = 'Chemkin',\n")
   end

   f:write(string.format("   nsegments = %d,\n", t.nsegments))
   for i=0,t.nsegments-1 do
      seg = "segment"..i
      f:write(string.format("   segment%d = {\n", i))
      f:write(string.format("      T_lower = %.1f,\n", t[seg].T_lower))
      f:write(string.format("      T_upper = %.1f,\n", t[seg].T_upper))
      f:write(string.format("      A = % -10.7e,\n", t[seg].A))
      f:write(string.format("      B = % -10.7e,\n", t[seg].B))
      f:write(string.format("      C = % -10.7e,\n", t[seg].C))
      f:write(string.format("      D = % -10.7e,\n", t[seg].D))
      f:write("   },\n")
   end
   f:write("}\n")
end

function writeThermPerfGas(f, species, db, optsTable)
   f:write("species = {")
   for _,sp in ipairs(species) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n\n")

   f:write("physical_model = 'thermally-perfect-gas'\n")
      

   f:write("energyModes = {'equilibrium'}\n")
   if (optsTable and optsTable.gas_giant_trans_props) then
      f:write("use_gas_giant_transport_properties = true\n")
      f:write(string.format("gas_giant_visc_db = '%s'\n", optsTable.gas_giant_visc_db))
      f:write(string.format("gas_giant_tc_db = '%s'\n", optsTable.gas_giant_tc_db))
      f:write(string.format("gas_giant_tc_model = '%s'\n", optsTable.gas_giant_tc_model))
   end
   f:write("db = {}\n")

   for _,sp in ipairs(species) do
      f:write(string.format("db['%s'] = {}\n", sp))
      f:write(string.format("db['%s'].atomicConstituents = { ", sp))
      for k,v in pairs(db[sp].atomicConstituents) do
         f:write(string.format("%s=%d, ", k, v))
      end
      f:write("}\n")
      f:write(string.format("db['%s'].charge = %d\n", sp, db[sp].charge))
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))

      diffusion_info_missing = false
      if db[sp].sigma then
         sigma = db[sp].sigma.value
      else
         diffusion_info_missing = true
         sigma = db.default.sigma.value
      end
      f:write(string.format("db['%s'].sigma = %.8f\n", sp, sigma))
      if db[sp].epsilon then
         epsilon = db[sp].epsilon.value
      else
          diffusion_info_missing = true
         epsilon = db.default.epsilon.value
      end
      f:write(string.format("db['%s'].epsilon = %.8f\n", sp, epsilon))
      -- Ionised species have a different potentials to LJ, so we don't mind them being missing (NNG)
      if ((db[sp].charge == 0) and diffusion_info_missing) then
          print("------------------------------------------------------------------------------------------")
          print("WARNING: Lennard-Jones potential data could not be found for species: ", sp)
          print("WARNING: As a substitute the values from 'defaults.lua' will be used.")
          print("WARNING: This may affect a multi-species calculation with mass diffusion.")
          print("------------------------------------------------------------------------------------------")
      end

      Lewis = db.default.Lewis.value
      if db[sp].Lewis then
         Lewis = db[sp].Lewis.value
      end
      f:write(string.format("db['%s'].Lewis = %.8f\n", sp, Lewis))
      writeCeaThermoCoeffs(f, sp, db, optsTable)
      -- Our preference is to use CEA transport properties where available,
      -- however we intercept on 'air' as a special case and use the
      -- the Sutherland model because CEA do not have a curve fit
      -- form for the transport properties of air.
      if ( sp == 'air' ) then
	 f:write("db['air'].viscosity = {\n")
	 f:write("   model = 'Sutherland',\n");
	 f:write(string.format("   mu_ref = %.8e,\n", db['air'].sutherlandVisc.mu_ref))
	 f:write(string.format("   T_ref = %.8f,\n", db['air'].sutherlandVisc.T_ref))
	 f:write(string.format("   S = %.8f,\n", db['air'].sutherlandVisc.S))
	 f:write("}\n")
      elseif (optsTable and optsTable.gas_giant_trans_props) then
         -- do nothing
      else
	   if optsTable and optsTable.transport_database == "prefer-chemkin" then
	      if ( not db[sp].chemkinViscosity ) then
		 print("WARNING: Chemkin viscosity coefficients have been selected as preferred,")
		 print("WARNING: but they could not be found for species: ", sp)
		 print("WARNING: We default back to using the CEA viscosity coefficients (if available).")
		 print("")
		 writeCeaTransCoeffs(f, sp, db, "Viscosity", "viscosity")
	      else
		 writeChemkinTransCoeffs(f, sp, db, "Viscosity", "viscosity")
	      end
	   else
		writeCeaTransCoeffs(f, sp, db, "Viscosity", "viscosity")
	   end
      end
      if ( sp == 'air' ) then
	 f:write("db['air'].thermal_conductivity = {\n")
	 f:write("   model = 'Sutherland',\n")
	 f:write(string.format("   k_ref = %.8e,\n", db['air'].sutherlandThermCond.k_ref))
	 f:write(string.format("   T_ref = %.8f,\n", db['air'].sutherlandThermCond.T_ref))
	 f:write(string.format("   S = %.8f,\n", db['air'].sutherlandThermCond.S))
	 f:write("}\n")
      elseif (optsTable and optsTable.gas_giant_trans_props) then
         -- do nothing
      else
	   if optsTable and optsTable.transport_database == "prefer-chemkin" then
	      if ( not db[sp].chemkinThermCond ) then
		 print("WARNING: Chemkin thermal conductivity coefficients have been selected as preferred,")
		 print("WARNING: but they could not be found for species: ", sp)
		 print("WARNING: We default back to using the CEA thermal conductivity coefficients (if available).")
		 print("")
		 writeCeaTransCoeffs(f, sp, db, "ThermCond", "thermal_conductivity")
	      else
		 writeChemkinTransCoeffs(f, sp, db, "ThermCond", "thermal_conductivity")
	      end
	   else
		writeCeaTransCoeffs(f, sp, db, "ThermCond", "thermal_conductivity")
	   end
      end
   end
end

function write2TAir(f, species, db, optsTable)
   f:write("species = {")
   for _,sp in ipairs(species) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n\n")
   f:write("db = {}\n")
   for _,sp in ipairs(species) do
      f:write(string.format("db['%s'] = {}\n", sp))
      f:write(string.format("db['%s'].atomicConstituents = { ", sp))
      for k,v in pairs(db[sp].atomicConstituents) do
         f:write(string.format("%s=%d, ", k, v))
      end
      f:write("}\n")
      f:write(string.format("db['%s'].charge = %d\n", sp, db[sp].charge))
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))
      writeCeaThermoCoeffs(f, sp, db, optsTable)
   end
end

function write2TN2(f, species, db, optsTable)
   -- species for this model are N2 and N in fixed order.
   speciesList = {'N2', 'N'}
   f:write("species = {")
   for _,sp in ipairs(speciesList) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n\n")
   f:write("db = {}\n")
   for _,sp in ipairs(speciesList) do
      f:write(string.format("db['%s'] = {}\n", sp))
      f:write(string.format("db['%s'].atomicConstituents = { ", sp))
      for k,v in pairs(db[sp].atomicConstituents) do
         f:write(string.format("%s=%d, ", k, v))
      end
      f:write("}\n")
      f:write(string.format("db['%s'].charge = %d\n", sp, db[sp].charge))
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))
      writeCeaThermoCoeffs(f, sp, db, optsTable)
   end
end

function write2TGas(f, species, db, optsTable)
    f:write("species = {")
   for _,sp in ipairs(species) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n\n")
   f:write("physical_model = 'two-temperature-gas'\n")
   f:write("db = {}\n")
   for _,sp in ipairs(species) do
      f:write(string.format("db['%s'] = {}\n", sp))
      f:write(string.format("db['%s'].type = '%s'\n", sp, db[sp].type))
      if db[sp].type == "molecule" then
         f:write(string.format("db['%s'].molecule_type = '%s'\n", sp, db[sp].molecule_type))
         f:write(string.format("db['%s'].theta_v = %.3f\n", sp, db[sp].theta_v.value))
      end
      f:write(string.format("db['%s'].atomicConstituents = { ", sp))
      for k,v in pairs(db[sp].atomicConstituents) do
         f:write(string.format("%s=%d, ", k, v))
      end
      f:write("}\n")
      f:write(string.format("db['%s'].charge = %d\n", sp, db[sp].charge))
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))
      if db[sp].Hf then
         f:write(string.format("db['%s'].Hf = %.8e\n", sp, db[sp].Hf.value))
      end
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))
      diffusion_info_missing = false
      if db[sp].sigma then
         sigma = db[sp].sigma.value
      else
         diffusion_info_missing = true
         sigma = db.default.sigma.value
      end
      f:write(string.format("db['%s'].sigma = %.8f\n", sp, sigma))
      if db[sp].epsilon then
         epsilon = db[sp].epsilon.value
      else
          diffusion_info_missing = true
         epsilon = db.default.epsilon.value
      end
      if db[sp].r_eq then
	 f:write(string.format("db['%s'].r_eq = %.8e\n", sp, db[sp].r_eq.value))
      end
      f:write(string.format("db['%s'].epsilon = %.8f\n", sp, epsilon))
      -- Ionised species have a different potentials to LJ, so we don't mind them being missing (NNG)
      if ((db[sp].charge == 0) and diffusion_info_missing) then
          print("------------------------------------------------------------------------------------------")
          print("WARNING: Lennard-Jones potential data could not be found for species: ", sp)
          print("WARNING: As a substitute the values from 'defaults.lua' will be used.")
          print("WARNING: This may affect a multi-species calculation with mass diffusion.")
          print("------------------------------------------------------------------------------------------")
      end

      Lewis = db.default.Lewis.value
      if db[sp].Lewis then
         Lewis = db[sp].Lewis.value
      end
      f:write(string.format("db['%s'].Lewis = %.8f\n", sp, Lewis))
      writeCeaThermoCoeffs(f, sp, db, optsTable)
   end
   writeCollisionIntegrals(f, species, db, optsTable)
end

function write3TGas(f, species, db, optsTable)
    f:write("species = {")
   for _,sp in ipairs(species) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n\n")
   f:write("physical_model = 'three-temperature-gas'\n")
   f:write("db = {}\n")
   for _,sp in ipairs(species) do
      f:write(string.format("db['%s'] = {}\n", sp))
      f:write(string.format("db['%s'].type = '%s'\n", sp, db[sp].type))
      if db[sp].type == "molecule" then
         f:write(string.format("db['%s'].molecule_type = '%s'\n", sp, db[sp].molecule_type))
         f:write(string.format("db['%s'].vib_data = {\n", sp))
         f:write("  model = 'harmonic',\n")
         f:write(string.format("  theta_v = %.3f\n", db[sp].theta_v.value))
         f:write("}\n")
      end
      if db[sp].type ~= "electron" then
         f:write(string.format("db['%s'].electronic_levels = {\n", sp))
         f:write("  model = 'multi-level',\n")
         f:write("  g  = {")
         for _,g in pairs(db[sp].electronic_levels.g.value) do
            f:write(string.format("%d, ", g))
         end
         f:write("},\n")
         f:write("  Te = {")
         for _,Te in pairs(db[sp].electronic_levels.Te.value) do
            f:write(string.format("%.3f, ", Te))
         end
         f:write("},\n")
         f:write("}\n")
      end
      f:write(string.format("db['%s'].atomicConstituents = { ", sp))
      for k,v in pairs(db[sp].atomicConstituents) do
         f:write(string.format("%s=%d, ", k, v))
      end
      f:write("}\n")
      f:write(string.format("db['%s'].charge = %d\n", sp, db[sp].charge))
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))
      if db[sp].Hf then
         f:write(string.format("db['%s'].Hf = %.8e\n", sp, db[sp].Hf.value))
      end
      if db[sp].ionisation_energy then
         f:write(string.format("db['%s'].ionisation_energy = %.8e\n", sp, db[sp].ionisation_energy.value))
      end
      diffusion_info_missing = false
      if db[sp].sigma then
         sigma = db[sp].sigma.value
      else
         diffusion_info_missing = true
         sigma = db.default.sigma.value
      end
      f:write(string.format("db['%s'].sigma = %.8f\n", sp, sigma))
      if db[sp].epsilon then
         epsilon = db[sp].epsilon.value
      else
          diffusion_info_missing = true
         epsilon = db.default.epsilon.value
      end
      f:write(string.format("db['%s'].epsilon = %.8f\n", sp, epsilon))
      if db[sp].r_eq then
	 f:write(string.format("db['%s'].r_eq = %.8e\n", sp, db[sp].r_eq.value))
      end
      -- Ionised species have a different potentials to LJ, so we don't mind them being missing (NNG)
      if ((db[sp].charge == 0) and diffusion_info_missing) then
          print("------------------------------------------------------------------------------------------")
          print("WARNING: Lennard-Jones potential data could not be found for species: ", sp)
          print("WARNING: As a substitute the values from 'defaults.lua' will be used.")
          print("WARNING: This may affect a multi-species calculation with mass diffusion.")
          print("------------------------------------------------------------------------------------------")
      end

      Lewis = db.default.Lewis.value
      if db[sp].Lewis then
         Lewis = db[sp].Lewis.value
      end
      f:write(string.format("db['%s'].Lewis = %.8f\n", sp, Lewis))
      writeCeaThermoCoeffs(f, sp, db, optsTable)
   end
   writeCollisionIntegrals(f, species, db, optsTable)
end

function writeMultiTGas(f, species, db, optsTable)
   f:write("physical_model = 'multi-temperature-gas'\n")
   f:write("species = {")
   for _,sp in ipairs(species) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n")
   f:write("energy_modes = {")
   for _,mode in ipairs(energy_modes) do
      f:write(string.format("'%s', ", mode))
   end
   f:write("}\n\n")

   f:write("db = {}\n")
   f:write("db.modes = {\n")
   for _,mode in ipairs(energy_modes) do
      f:write(string.format("  %s = { ", mode))
      for _,cmpnt in ipairs(mode_components[mode]) do
	 f:write(string.format("'%s', ", cmpnt))
      end
      f:write("  },\n")
   end
   f:write("}\n")


   for _,sp in ipairs(species) do
      f:write(string.format("db['%s'] = {}\n", sp))
      f:write(string.format("db['%s'].type = '%s'\n", sp, db[sp].type))
      if db[sp].type == "molecule" then
         f:write(string.format("db['%s'].molecule_type = '%s'\n", sp, db[sp].molecule_type))
         f:write(string.format("db['%s'].vib_data = {\n", sp))
	 f:write("  model = 'from-cea-thermo-curve',\n")
	 f:write(string.format("  theta_v = %.3f,\n", db[sp].theta_v.value))
	 f:write(string.format("  theta_D = %.3f,\n", db[sp].theta_D.value))
	 f:write("}\n")
      end
      if db[sp].type ~= "electron" then
         f:write(string.format("db['%s'].electronic_levels = {\n", sp))
         f:write("  model = 'multi-level',\n")
         f:write("  g  = {")
         for _,g in pairs(db[sp].electronic_levels.g.value) do
            f:write(string.format("%d, ", g))
         end
         f:write("},\n")
         f:write("  Te = {")
         for _,Te in pairs(db[sp].electronic_levels.Te.value) do
            f:write(string.format("%.3f, ", Te))
         end
         f:write("},\n")
         f:write("}\n")
      end
      f:write(string.format("db['%s'].atomicConstituents = { ", sp))
      for k,v in pairs(db[sp].atomicConstituents) do
         f:write(string.format("%s=%d, ", k, v))
      end
      f:write("}\n")
      f:write(string.format("db['%s'].charge = %d\n", sp, db[sp].charge))
      f:write(string.format("db['%s'].M = %.8e\n", sp, db[sp].M.value))
      if db[sp].Hf then
         f:write(string.format("db['%s'].Hf = %.8e\n", sp, db[sp].Hf.value))
      end
      diffusion_info_missing = false
      if db[sp].sigma then
         sigma = db[sp].sigma.value
      else
         diffusion_info_missing = true
         sigma = db.default.sigma.value
      end
      f:write(string.format("db['%s'].sigma = %.8f\n", sp, sigma))
      if db[sp].epsilon then
         epsilon = db[sp].epsilon.value
      else
          diffusion_info_missing = true
         epsilon = db.default.epsilon.value
      end
      f:write(string.format("db['%s'].epsilon = %.8f\n", sp, epsilon))
      if db[sp].r_eq then
	 f:write(string.format("db['%s'].r_eq = %.8e\n", sp, db[sp].r_eq.value))
      end
      -- Ionised species have a different potentials to LJ, so we don't mind them being missing (NNG)
      if ((db[sp].charge == 0) and diffusion_info_missing) then
          print("------------------------------------------------------------------------------------------")
          print("WARNING: Lennard-Jones potential data could not be found for species: ", sp)
          print("WARNING: As a substitute the values from 'defaults.lua' will be used.")
          print("WARNING: This may affect a multi-species calculation with mass diffusion.")
          print("------------------------------------------------------------------------------------------")
      end

      Lewis = db.default.Lewis.value
      if db[sp].Lewis then
         Lewis = db[sp].Lewis.value
      end
      f:write(string.format("db['%s'].Lewis = %.8f\n", sp, Lewis))
      writeCeaThermoCoeffs(f, sp, db, optsTable)
   end
   writeCollisionIntegrals(f, species, db, optsTable)
end

function writeCollisionIntegrals(f, species, db, optsTable)
   cidb = optsTable.ci_database or "gupta"
   fname = CIDBFileName[cidb]
   if not fname then
      print("Collision integral database is not available: ", cidb)
      print("Bailing out!")
      os.exit(1)
   end
   DGD = os.getenv("DGD")
   dir = DGD.."/data/"
   dbName = dir..fname
   dofile(dbName)
   print("Collision integral database loaded from: ", dbName, "\n")
   
   f:write("db.CIs = {}\n")
   for isp,p in ipairs(species) do
      for jsp=1,isp do
         q = species[jsp]
         key = p .. ":" .. q
         ci = cis[key]
         if not ci then
            -- We'll try a reverse key
            key = q .. ":" .. p
            ci = cis[key]
            if not ci then
               print(string.format("Collision integral data for colliding pair --  %s:%s -- could not be found in database.", p, q))
               if optsTable.CI_default then
                  ci = cis[optsTable.CI_default]
                  print(string.format("Substituting with default: %s.\n", optsTable.CI_default))
               else
                  print("No default selection provided.\n")
                  print("Bailing out!")
                  os.exit(1)
               end
            end
         end
         f:write(string.format("db.CIs['%s'] = {\n", key))
         piList = {"pi_Omega_11", "pi_Omega_22"}
         for _,key in ipairs(piList) do
            f:write(string.format("   %s = {A= % 6.4f, B= % 6.4f, C= % 6.4f, D= % 6.4f},\n",
                                  key, ci[key].A, ci[key].B, ci[key].C, ci[key].D))
         end
         f:write("}\n")
      end
   end
end


gasModels = {}
-- IdealGas and aliases
gasModels["IDEALGAS"] = {writeFn=writeIdealGas, DName="IdealGas"}
gasModels["IDEAL GAS"] = gasModels["IDEALGAS"]
-- Thermally perfect gas and aliases
gasModels["THERMALLYPERFECTGAS"] = {writeFn=writeThermPerfGas, DName="CompositeGas"}
gasModels["THERMALLY PERFECT GAS"] = gasModels["THERMALLYPERFECTGAS"]
gasModels["THERMALLY PERFECT GAS EQUILIBRIUM"] = {writeFn=writeThermPerfGas, DName="ThermallyPerfectGasEquilibrium"}
-- Two-temperature air and aliases
gasModels["TWOTEMPERATUREAIR"] = {writeFn=write2TAir, DName="TwoTemperatureAir"}
gasModels["TWO TEMPERATURE AIR"] = gasModels["TWOTEMPERATUREAIR"]
gasModels["TWO-TEMPERATURE AIR"] = gasModels["TWOTEMPERATUREAIR"]
-- Two-temperature nitogren (N2 and N)
gasModels["TWOTEMPERATUREN2"] = {writeFn=write2TN2, DName="TwoTemperatureDissociatingNitrogen"}
gasModels["TWO TEMPERATURE DISSOCIATING NITROGEN"] = gasModels["TWOTEMPERATUREN2"]
gasModels["TWO-TEMPERATURE DISSOCIATING NITROGEN"] = gasModels["TWOTEMPERATUREN2"]
-- Two-tempeature gas (generalised)
gasModels["TWOTEMPERATUREGAS"] = {writeFn=write2TGas, DName="CompositeGas"}
gasModels["TWO TEMPERATURE GAS"] = gasModels["TWOTEMPERATUREGAS"]
gasModels["TWO-TEMPERATURE GAS"] = gasModels["TWOTEMPERATUREGAS"]
-- Multi temperature gas
gasModels["MULTITEMPERATUREGAS"] = {writeFn=writeMultiTGas, DName="CompositeGas"}
gasModels["MULTI TEMPERATURE GAS"] = gasModels["MULTITEMPERATUREGAS"]
gasModels["MULTI-TEMPERATURE GAS"] = gasModels["MULTITEMPERATUREGAS"]
gasModels["MULTI-T GAS"] = gasModels["MULTITEMPERATUREGAS"]
-- Three-temperature gas
gasModels["THREETEMPERATUREGAS"] = {writeFn=write3TGas, DName="CompositeGas"}
gasModels["THREE TEMPERATURE GAS"] = gasModels["THREETEMPERATUREGAS"]
gasModels["THREE-TEMPERATURE GAS"] = gasModels["THREETEMPERATUREGAS"]

function printHelp()
   print("prep-gas --- Prepares a gas model input file for Eilmer4.")
   print("Usage:")
   print(" > prep-gas input output")
   print("")
   print("   input    : a gas input file with user selections")
   print("   output   : detailed gas file in format ready for Eilmer4.")
   print("")
   print(" > prep-gas --list-available-species")
   print("")
   os.exit(1)
end

function parseSpeciesList(fname)
   f = assert(io.open(fname, 'r'))
   spMap = {}
   while true do
      line = f:read('*line')
      if not line then
	 break
      end
      tks = split_string(line)
      if tks[1] ~= '#' and tks[3] then
	 spSymbol = tks[3]
	 spMap[spSymbol] = true
      end
   end
   f:close()
   return spMap
end


function prepareGasFile(inpFname, outFname)
   dofile(inpFname)
   local inEquilibrium = false
   if equilibrium then
      inEquilibrium = true
   end
   -- Convert species names to database entries
   spNames = {}
   for _,sp in ipairs(species) do
      spNames[#spNames+1] = sp
   end
   -- Locate species list and parse
   listName = dir.."species-list.txt"
   spMap = parseSpeciesList(listName)
   for i,sp in ipairs(spNames) do
      if ( not spMap[sp] ) then
	 print(string.format("The requested species '%s' could not be located in the species list.", sp))
	 print("Exiting without doing anything.")
	 os.exit(1)
      end
      species[i] = sp
   end
   -- Check we have all the species
   for _,sp in ipairs(species) do
      if ( not db[sp] ) then
	 print(string.format("The requested species '%s' could not be located in the species database.", sp))
	 print("Exiting without doing anything.")
	 os.exit(1)
      end
   end
   -- Check we know what to do with the model
   if ( not gasModels[string.upper(model)] ) then
      print(string.format("The requested model '%s' is unknown.", model))
      print("Available models are: ")
      local tmp = {}
      for k,v in pairs(gasModels) do table.insert(tmp, k) end
      table.sort(tmp)
      for _,k in ipairs(tmp) do
	 print(k)
      end
      print("Exiting without doing anything.")
   end
   -- Given we know what to do, write model name then delegate to
   -- specialised function.
   local f = assert(io.open(outFname, 'w'))
   f:write(string.format("-- Auto-generated by prep-gas on: %s\n\n",
			 os.date("%d-%b-%Y %X")))
   modelStr = gasModels[string.upper(model)].DName
   if inEquilibrium then
      modelStr = modelStr .. "Equilibrium"
   end
   f:write(string.format("model = '%s'\n", modelStr))
   options = options or {}
   gasModels[string.upper(model)].writeFn(f, species, db, options)
   f:close()
end

function main()
   local listSpecies = false
   if ( arg[1] == "--help" or arg[1] == "-h" ) then
      printHelp()
   elseif ( arg[1] == "--list-available-species" ) then
      listSpecies = true
   elseif ( #arg < 2 ) then
      print("Not enough arguments or unknown option.")
      print("Exiting program without doing anything.")
      printHelp()
   elseif ( #arg > 3 ) then
      print("Too many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   local inpFname, outFname
   inpFname = arg[1]
   outFname = arg[2]
   
   -- Locate species database and load
   DGD = os.getenv("DGD")
   dir = DGD.."/data/"
   dbName = dir.."species-database.lua"
   dofile(dbName)
   print("Species database loaded from: ", dbName)

   if listSpecies then
      spList = {}
      for sp,_ in pairs(db) do
	 if sp ~= 'default' then
	    spList[#spList+1] = sp
	 end
      end
      table.sort(spList)
      print("")
      for _,sp in ipairs(spList) do
	 print(sp)
      end
      print("")
      return 0
   end

   prepareGasFile(inpFname, outFname)
end

main()
