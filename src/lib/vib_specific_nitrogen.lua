-- vib_specific_nitrogen.lua
--
-- Service module for the vibrationally-specific nitrogen model in:
-- Giordano, et al. (1997)
-- Vibrationally Relaxing Flow of N2 past an Infinite Cylinder
-- Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27 - 35
--
-- R.Gollan and PJ
-- 2021-11-23 factored out of 2D/cylinder-giordano/ input script.
--

local vsn2 = {}

function vsn2.massf(nsp, T)
   -- Given an equilibrium temperature T,
   -- compute populations in number density for first nsp vibrational levels,
   -- then convert to mass fractions and return that table of mass fractions.
   --
   -- Note that the returned table uses species names that match those defined
   -- in the corresponding Dlang module.
   --
   local massf = {}
   -- Anharmonicity constants
   local w_e = 235857
   local we_xe = 1432.4
   local we_ye = -0.226
   -- Other physical parameters
   local h = 6.626*10^(-34) --m^2.kg/s
   local c = 2.998*10^8     --m/s
   local kb = 1.38*10^(-23) --J/K
   --
   Zsum = 0
   --
   for i=0,nsp-1 do
      E_vib = h*c*(w_e*(i+0.5)-we_xe*(i+0.5)^2+we_ye*(i+0.5)^3)
      Zsum = Zsum + math.exp(-E_vib/(kb*T))
   end
   for i=0,nsp-1 do
      E_vib = h*c*(w_e*(i+0.5)-we_xe*(i+0.5)^2+we_ye*(i+0.5)^3)
      massf["N2-vib-"..i] = math.exp(-E_vib/(kb*T))/Zsum
   end
   return massf
end

return vsn2
