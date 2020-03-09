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
   f:write("energyModes = {'equilibrium'}\n")
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
      sigma = db.default.sigma.value
      if db[sp].sigma then
	 sigma = db[sp].sigma.value
      end
      f:write(string.format("db['%s'].sigma = %.8f\n", sp, sigma))
      epsilon = db.default.epsilon.value
      if db[sp].epsilon then
	 epsilon = db[sp].epsilon.value
      end
      f:write(string.format("db['%s'].epsilon = %.8f\n", sp, epsilon))
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
      writeCeaThermoCoeffs(f, sp, db, optsTable)
   end
end


function writeCO2Gas(f, sp, db)
   -- This is the Bender Model, entropy Ref values, viscosity etc. are not needed
   if ( #sp > 1 ) then
      print("WARNING: More than one species is listed while trying to prepare")
      print("WARNING: an ideal gas model.")
      print("WARNING: Presently, the ideal gas model is limited to a single species")
      print("WARNING: We will build a file with the first species listed and ignore the rest.")
   end
   local s = sp[1]
   f:write("CO2Gas = {\n")
   f:write(string.format("   speciesName = '%s',\n", s))
   f:write(string.format("   mMass = %.8f,\n", db[s].M.value))
   f:write(string.format("   gamma = %.8f,\n", db[s].gamma.value))
   f:write("   entropyRefValues = {\n")
   f:write(string.format("      s1 = %.8e,\n", db[s].entropyRefValues.s1))
   f:write(string.format("      T1 = %.8f,\n", db[s].entropyRefValues.T1))
   f:write(string.format("      p1 = %.8e,\n", db[s].entropyRefValues.p1))
   f:write("   },\n")
   f:write("   sutherlandVisc = {\n")
   f:write(string.format("      mu_ref = %.8e,\n", db[s].sutherlandVisc.mu_ref))
   f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandVisc.T_ref))
   f:write(string.format("      S = %.8f,\n", db[s].sutherlandVisc.S))
   f:write("   },\n")
   f:write("   sutherlandThermCond = {\n")
   f:write(string.format("      k_ref = %.8e,\n", db[s].sutherlandThermCond.k_ref))
   f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandThermCond.T_ref))
   f:write(string.format("      S = %.8f,\n", db[s].sutherlandThermCond.S))
   f:write("   }\n")
   f:write("}\n")
end

function writeCO2GasSW(f, sp, db)
   -- This is the Span Wagner Model
   if ( #sp > 1 ) then
      print("WARNING: More than one species is listed while trying to prepare")
      print("WARNING: an ideal gas model.")
      print("WARNING: Presently, the ideal gas model is limited to a single species")
      print("WARNING: We will build a file with the first species listed and ignore the rest.")
   end
   local s = sp[1]
   f:write("CO2GasSW = {\n")
   f:write(string.format("   speciesName = '%s',\n", s))
   f:write(string.format("   mMass = %.8f,\n", db[s].M.value))
   f:write(string.format("   gamma = %.8f,\n", db[s].gamma.value))
   f:write("   entropyRefValues = {\n")
   f:write(string.format("      s1 = %.8e,\n", db[s].entropyRefValues.s1))
   f:write(string.format("      T1 = %.8f,\n", db[s].entropyRefValues.T1))
   f:write(string.format("      p1 = %.8e,\n", db[s].entropyRefValues.p1))
   f:write("   },\n")
   f:write("   sutherlandVisc = {\n")
   f:write(string.format("      mu_ref = %.8e,\n", db[s].sutherlandVisc.mu_ref))
   f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandVisc.T_ref))
   f:write(string.format("      S = %.8f,\n", db[s].sutherlandVisc.S))
   f:write("   },\n")
   f:write("   sutherlandThermCond = {\n")
   f:write(string.format("      k_ref = %.8e,\n", db[s].sutherlandThermCond.k_ref))
   f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandThermCond.T_ref))
   f:write(string.format("      S = %.8f,\n", db[s].sutherlandThermCond.S))
   f:write("   },\n")
   f:write("   LUTfilenames = {\n")
   f:write(string.format("      p_rhoe_file = '%s',\n", "../LUT/P_rhoe_Tree.dat"))
   f:write(string.format("      a_rhoe_file = '%s',\n", "../LUT/a_rhoe_Tree.dat"))
   f:write(string.format("      T_rhoe_file = '%s',\n", "../LUT/T_rhoe_Tree.dat"))
   f:write(string.format("      e_rho_sat_file = '%s',\n", "../LUT/e_rho_sat_table.dat"))
   f:write(string.format("      rho_sh_file = '%s',\n", "../LUT/rho_sh_Tree.dat"))
   f:write(string.format("      T_sh_file = '%s',\n", "../LUT/T_sh_Tree.dat"))
   f:write(string.format("      lookup_hsFlag = %s,\n", 1))
   f:write(string.format("      lookup_rhoeFlag = %s,\n", 1))
   f:write("   }\n")
   f:write("}\n")
end

function writeSF6Virial(f, sp, db)
   -- This is the Bender Model
   if ( #sp > 1 ) then
      print("WARNING: More than one species is listed while trying to prepare")
      print("WARNING: an ideal gas model.")
      print("WARNING: Presently, the ideal gas model is limited to a single species")
      print("WARNING: We will build a file with the first species listed and ignore the rest.")
   end
   local s = sp[1]
   f:write("SF6Virial = {\n")
   f:write(string.format("   speciesName = '%s',\n", s))
   f:write(string.format("   mMass = %.8f,\n", db[s].M.value))
   f:write(string.format("   gamma = %.8f,\n", db[s].gamma.value))
   f:write("   entropyRefValues = {\n")
   f:write(string.format("      s1 = %.8e,\n", db[s].entropyRefValues.s1))
   f:write(string.format("      T1 = %.8f,\n", db[s].entropyRefValues.T1))
   f:write(string.format("      p1 = %.8e,\n", db[s].entropyRefValues.p1))
   f:write("   },\n")
   f:write("   sutherlandVisc = {\n")
   f:write(string.format("      mu_ref = %.8e,\n", db[s].sutherlandVisc.mu_ref))
   f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandVisc.T_ref))
   f:write(string.format("      S = %.8f,\n", db[s].sutherlandVisc.S))
   f:write("   },\n")
   f:write("   sutherlandThermCond = {\n")
   f:write(string.format("      k_ref = %.8e,\n", db[s].sutherlandThermCond.k_ref))
   f:write(string.format("      T_ref = %.8f,\n", db[s].sutherlandThermCond.T_ref))
   f:write(string.format("      S = %.8f,\n", db[s].sutherlandThermCond.S))
   f:write("   }\n")
   f:write("}\n")
end

gasModels = {}
-- IdealGas and aliases
gasModels["IDEALGAS"] = {writeFn=writeIdealGas, DName="IdealGas"}
gasModels["IDEAL GAS"] = gasModels["IDEALGAS"]
-- Thermally perfect gas and aliases
gasModels["THERMALLYPERFECTGAS"] = {writeFn=writeThermPerfGas, DName="ThermallyPerfectGas"}
gasModels["THERMALLY PERFECT GAS"] = gasModels["THERMALLYPERFECTGAS"]
gasModels["THERMALLY PERFECT GAS EQUILIBRIUM"] = {writeFn=writeThermPerfGas, DName="ThermallyPerfectGasEquilibrium"}
-- Two-temperature air and aliases
gasModels["TWOTEMPERATUREAIR"] = {writeFn=write2TAir, DName="TwoTemperatureAir"}
gasModels["TWO TEMPERATURE AIR"] = gasModels["TWOTEMPERATUREAIR"]
gasModels["TWO-TEMPERATURE AIR"] = gasModels["TWOTEMPERATUREAIR"]
-- CO2
gasModels["CO2GAS"] = {writeFn=writeCO2Gas, DName = "CO2Gas"}
gasModels["CO2GASSW"] = {writeFn = writeCO2GasSW, DName = "CO2GasSW"}
--SF6
gasModels["SF6VIRIAL"] = {writeFn = writeSF6Virial, DName = "SF6Virial"}

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
   elseif ( #arg > 2 ) then
      print("Too many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end

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

   inpFname = arg[1]
   outFname = arg[2]
   prepareGasFile(inpFname, outFname)
end

main()
