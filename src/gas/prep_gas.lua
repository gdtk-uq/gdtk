#!/usr/bin/env lua
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

function writeIdealGas(f, sp, db)
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

function writeCeaThermoCoeffs(f, sp, db)
   if ( not db[sp].ceaThermoCoeffs ) then
      print("ERROR: The table of CEA coefficients to compute thermodynamic properties")
      print("ERROR: could not be found for species: ", sp)
      print("ERROR: Bailing out!")
      os.exit(1)
   end
   t = db[sp].ceaThermoCoeffs
   f:write(string.format("%s.ceaThermoCoeffs = {\n", sp))
   f:write(string.format("  nsegments = %d, \n", t.nsegments))
   for i=0,t.nsegments-1 do
      seg = "segment"..i
      f:write(string.format("  segment%d = {\n", i))
      f:write(string.format("    T_lower = %.1f,\n", t[seg].T_lower))
      f:write(string.format("    T_upper = %.1f,\n", t[seg].T_upper))
      f:write("    coeffs = {\n")
      for _,c in ipairs(t[seg].coeffs) do
	 f:write(string.format("      % -12.9e,\n", c))
      end
      f:write("    }\n")
      f:write("  },\n")
   end
   f:write("}\n")
end

function writeCeaTransCoeffs(f, sp, db, name)
   secName = "cea"..name
   if ( not db[sp][secName] ) then
      print(string.format("ERROR: The table of CEA coefficients to compute diffusion property: %s", secName))
      print("ERROR: could not be found for species: ", sp)
      print("ERROR: Bailing out!")
      os.exit(1)
   end
   t = db[sp][secName]
   f:write(string.format("%s.%s = {\n", sp, secName))
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

function writeThermPerfGas(f, species, db)
   f:write("species = {")
   for _,sp in ipairs(species) do
      f:write(string.format("'%s', ", sp))
   end
   f:write("}\n\n")
   for _,sp in ipairs(species) do
      f:write(string.format("%s = {}\n", sp))
      f:write(string.format("%s.M = %.8f\n", sp, db[sp].M.value))
      writeCeaThermoCoeffs(f, sp, db)
      writeCeaTransCoeffs(f, sp, db, "Viscosity")
      writeCeaTransCoeffs(f, sp, db, "ThermCond")
   end
end
function writeCO2Gas(f, sp, db)
   -- We can only deal with one species for an ideal gas.
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

gasModels = {}
-- IdealGas and aliases
gasModels["IDEALGAS"] = {writeFn=writeIdealGas, DName="IdealGas"}
gasModels["IDEAL GAS"] = gasModels["IDEALGAS"]
-- Thermally perfect gas and aliases
gasModels["THERMALLYPERFECTGAS"] = {writeFn=writeThermPerfGas, DName="ThermallyPerfectGas"}
gasModels["THERMALLY PERFECT GAS"] = gasModels["THERMALLYPERFECTGAS"]
-- CO2
gasModels["CO2GAS"] = {writeFn=writeCO2Gas, DName = "CO2Gas"}

function printHelp()
   print("prep-gas --- Prepares a gas model input file for Eilmer4.")
   print("Usage:")
   print(" > prep-gas input output")
   print("")
   print("   input    : a gas input file with user selections")
   print("   output   : detailed gas file in format ready for Eilmer4.")
   print("")
   os.exit(1)
end

function prepareGasFile(inpFname, outFname)
   dofile(inpFname)
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
   f:write(string.format("model = '%s'\n", gasModels[string.upper(model)].DName))
   gasModels[string.upper(model)].writeFn(f, species, db)
   f:close()
end

function main()
   if ( arg[1] == "--help" or arg[1] == "-h" ) then
      printHelp()
   end
   if ( #arg < 2 ) then
      print("Not enough arguments or unknown option.")
      print("Exiting program without doing anything.")
      printHelp()
   end
   
   if ( #arg > 2 ) then
      print("Two many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end
   
   -- Locate species database and load
   CFDSRC = os.getenv("CFCFD_SRC") or os.getenv("HOME").."/cfcfd3"
   dir = CFDSRC.."/dlang/gas/species-database/"
   dbName = dir.."species-database.lua"
   dofile(dbName)
   print("Species database loaded from: ", dbName)
   
   inpFname = arg[1]
   outFname = arg[2]
   prepareGasFile(inpFname, outFname)
end

main()
