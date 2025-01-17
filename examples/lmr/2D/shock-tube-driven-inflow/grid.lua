print("Sod's Classic Shock tube, just the driven tube.")
-- Adapted from shock-fitting example, PJ 2025-01-17
--
-- We start the calculation with a single region for the driven gas.
-- The effect of the driver gas enters via the inflow boundary condition.
patch_r = CoonsPatch:new{p00={x=0.5, y=0.0}, p01={x=0.5, y=0.1},
                         p10={x=1.0, y=0.0}, p11={x=1.0, y=0.1}}
grid0 = registerFluidGridArray{
   grid=StructuredGrid:new{psurface=patch_r, niv=201, njv=3},
   nib=2, njb=1,
   fsTag="initial",
   bcTags={west="inflow"}
}
