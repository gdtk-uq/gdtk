#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Date: 2016-01-17
--
-- This script can be used to convert species database
-- files from eilmer3 for use in eilmerD. It only 
-- considers the properties needed for a thermally
-- perfect gas.

function writeSpeciesData(f, spData, spName)
   f:write("db."..spName.." = {}\n")
   -- assemble and write atomic constituents
   str = "{"
   for k,v in pairs(spData.atomic_constituents) do
      str = str .. string.format("%s=%d,", k, v)
   end
   str = str .. "}"
   f:write(string.format("db.%s.atomicConstituents = %s\n", spName, str))
   -- write charge data
   f:write(string.format("db.%s.charge = %d\n", spName, spData.charge))
   -- write out molecular mass data
   f:write("db."..spName..".M = {\n")
   f:write(string.format("   value = %.8f,\n", spData.M.value))
   f:write(string.format("   units = '%s',\n", spData.M.units))
   f:write(string.format("   description = '%s',\n", spData.M.description))
   f:write(string.format("   reference = '%s'\n", spData.M.reference))
   f:write("}\n")
   -- write out gamma value
   f:write("db."..spName..".gamma = {\n")
   f:write(string.format("   value = %.8f,\n", spData.gamma.value))
   f:write(string.format("   units = '%s',\n", spData.gamma.units))
   f:write(string.format("   description = '%s',\n", spData.gamma.description))
   f:write(string.format("   reference = '%s'\n", spData.gamma.reference))
   f:write("}\n")
   -- write out thermo coefficients
   f:write("db."..spName..".ceaThermoCoeffs = {\n")
   f:write(string.format("   nsegments = %d,\n", #spData.CEA_coeffs))
   for i,tab in ipairs(spData.CEA_coeffs) do
      f:write(string.format("   segment%d = {\n", i-1))
      f:write(string.format("      T_lower = %.1f,\n", tab.T_low))
      f:write(string.format("      T_upper = %.1f,\n", tab.T_high))
      f:write("      coeffs = {\n")
      for _,val in ipairs(tab.coeffs) do
	 f:write(string.format("         % -12.9e,\n", val))
      end
      f:write("      }\n")
      f:write("   },\n")
   end
   f:write("}\n")
   -- write out viscosity data
   if ( spData.viscosity and spData.viscosity.model == "CEA" ) then
      f:write("db."..spName..".ceaViscosity = {\n")
      f:write(string.format("   nsegments = %d,\n", #spData.viscosity.parameters))
      for i,tab in ipairs(spData.viscosity.parameters) do
	 f:write(string.format("   segment%d = {\n", i-1))
	 f:write(string.format("      T_lower = %.1f,\n", tab.T_low))
	 f:write(string.format("      T_upper = %.1f,\n", tab.T_high))
	 f:write(string.format("      A = % -10.7e,\n", tab.A))
	 f:write(string.format("      B = % -10.7e,\n", tab.B))
	 f:write(string.format("      C = % -10.7e,\n", tab.C))
	 f:write(string.format("      D = % -10.7e\n", tab.D))
	 f:write("   },\n")
      end
      f:write("}\n")
   else
      print("WARNING: No CEA viscosity data for species: ", spName)
   end
   -- write out thermal conductivity
   if ( spData.thermal_conductivity and spData.thermal_conductivity.model == "CEA" ) then
      f:write("db."..spName..".ceaThermCond = {\n")
      f:write(string.format("   nsegments = %d,\n", #spData.thermal_conductivity.parameters))
      for i,tab in ipairs(spData.thermal_conductivity.parameters) do
	 f:write(string.format("   segment%d = {\n", i-1))
	 f:write(string.format("      T_lower = %.1f,\n", tab.T_low))
	 f:write(string.format("      T_upper = %.1f,\n", tab.T_high))
	 f:write(string.format("      A = % -10.7e,\n", tab.A))
	 f:write(string.format("      B = % -10.7e,\n", tab.B))
	 f:write(string.format("      C = % -10.7e,\n", tab.C))
	 f:write(string.format("      D = % -10.7e\n", tab.D))
	 f:write("   },\n")
      end
      f:write("}\n")
   else
      print("WARNING: No CEA thermal conductivity data for species: ", spName)
   end
end

function main()
   if ( #arg ~= 2 ) then
      print("Incorrect number of arguments.")
      print("Two arguments expected:")
      print("  1. old species database file")
      print("  2. the name of a new species database file")
      print("Exiting without doing anything.")
      os.exit(1)
   end

   dofile(arg[1])
   spName,_ = string.match(arg[1], "([^.]+).([^.]+)")

   f = assert(io.open(arg[2], 'w'))
   writeSpeciesData(f, _G[spName], spName)
   f:close()
   print("Created file: ", arg[2])
   print("Done.")
end

main()
