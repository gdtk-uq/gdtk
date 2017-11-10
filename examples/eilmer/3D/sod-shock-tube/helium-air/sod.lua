-- sod.lua
-- Simple shock tube to test timing of waves.
-- PJ, 08-Sep-2006 adapted from Tcl script to Python
--     11-Feb-2009 ported to Eilmer3 to demonstrate the use of
--        user-supplied functions for geometry and flow conditions.
--     27-Feb-2016 ported to Eilmer4
--     25-May-2017 Helium+nitrogen+oxygen gas mixture for Sam Stennett
config.title = "One-dimensional shock tube with helium driving air."
print(config.title)
config.dimensions = 3

function tube_volume(r, s, t)
   -- User-defined function for the parametric volume maps from
   -- parametric space to physical space.
   -- Note that a table of labelled coordinates is returned.
   -- A simple hexahedron, one unit long in the x-direction
   -- and 0.1 in the x and y directions.
   return {x=1.0*r, y=0.1*s, z=0.1*t}
end
pvol=LuaFnVolume:new{luaFnName="tube_volume"}
grid0 = StructuredGrid:new{pvolume=pvol, niv=101, njv=3, nkv=3}

nsp, nmodes = setGasModel('he-n2-o2-gas-model.lua')
print("GasModel nsp= ", nsp, " nmodes= ", nmodes)

function tube_gas(x, y, z)
   -- User-defined function for the initial gas state works in physical space.
   -- Half the domain has high pressure and the other half low pressure.
   local massFractions, p, T
   if x < 0.5 then
      -- Fill the left-half of the volume with high-pressure helium.
      massFractions = {He=1.0, N2=0.0, O2=0.0}
      p = 1.0e5; T = 348.4
   else
      -- and the right-half with low-pressure air.
      massFractions = {He=0.0, N2=0.7778, O2=0.2222}
      p = 1.0e4; T = 278.8
   end
   -- We use the FlowState object to conveniently set all of
   -- the relevant properties. 
   return FlowState:new{p=p, velx=0.0, vely=0.0, T=T, massf=massFractions}
end

-- Define a single block for the tube and fill it with gas conditions.
FluidBlock:new{grid=grid0, initialState=tube_gas}

config.flux_calculator = "ausmdv"
config.max_time = 0.6e-3  -- seconds
config.max_step = 600
config.dt_init = 1.0e-6

