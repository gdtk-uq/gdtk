-- Code for blowing species through a solid wall using a user defined lua function.
-- @author: Nick Gibbons (19/12/22)

function convective_flux(args)

   cell_adjacent = sampleFluidCell(blkId, args.i, args.j, args.k)
   flux = {} -- create empty table for the user flux quantities

   rhow = cell_adjacent.rho
   Tw = 2000.0
   yOw = cell_adjacent.massf.O
   --gamma1 = 0.63*math.exp(-1160.0/Tw)
   gamma1 = 0.3
   MO = 16e-3
   MC = 12e-3
   MCO = 28e-3
   Ru = 8.31446261815324
   nu_O = math.sqrt(Ru*Tw/2.0/3.14159/MO)
   r1 = rhow*yOw*nu_O*gamma1*MC/MO

   mdotCO = r1*MCO/MC
   mdotO = -r1*MO/MC

   -- Let's start by assuming that the momentum and energy are unaffected by a small amount
   -- of surface chemistry. This may not work.
   --flux.momentum_x = -- TODO
   --flux.momentum_y = -- TODO
   --flux.momentum_z = -- TODO
   --flux.total_energy = -- TODO

   -- This needs to be the actual flux of the species densities, in kg/m3/s/m2
   -- The code upstairs will sort out what to do with them, since the conserved quantities
   -- are actually different between the transient solver and the JFNK
   outsign = args.outsign
   flux.species = {}
   flux.species.CO = -1*outsign*mdotCO
   flux.species.O = -1*outsign*mdotO

   return flux
end
