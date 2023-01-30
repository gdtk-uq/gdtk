-- he-air-1T.lua
-- PJ 2023-01-30
config.title = "One-dimensional shock tube with helium driving air."
print(config.title)

quad0 = CoonsPatch:new{p00={x=-2.0,y=0}, p10={x=9.0,y=0},
                       p11={x=9.0,y=0.1}, p01={x=-2.0,y=0.1}}
grid0 = StructuredGrid:new{psurface=quad0, niv=2001, njv=3}

nsp, nmodes = setGasModel('he-air-5sp-1T-gas-model.lua')
print("GasModel nsp= ", nsp, " nmodes= ", nmodes)
config.reacting = true
config.reactions_file = 'he-air-5sp-1T-chem.lua'

function tube_gas(x, y, z)
   -- User-defined function for the initial gas state works in physical space.
   -- Half the domain has high pressure and the other half low pressure.
   local massFractions, p, T
   if x < 0.0 then
      -- Fill the left-half of the volume with high-pressure helium.
      massFractions = {He=1.0, N2=0.0, O2=0.0}
      p = 0.8e6; T = 800.0
   else
      -- and the right-half with low-pressure air.
      massFractions = {He=0.0, N2=0.7778, O2=0.2222}
      p = 290.0; T = 293.0
   end
   -- We use the FlowState object to conveniently set all of
   -- the relevant properties.
   return FlowState:new{p=p, velx=0.0, vely=0.0, T=T, massf=massFractions}
end

-- Define a single block for the tube and fill it with gas conditions.
blks = FBArray:new{grid=grid0, nib=8, njb=1, initialState=tube_gas}

config.flux_calculator = "ausmdv"
config.max_time = 3.0e-3  -- seconds
config.max_step = 60000
config.dt_init = 1.0e-8
config.dt_plot = 20.0e-6
config.dt_hist = 1.0e-6
-- TODO history points

