print("Periodic shear layer -- set up transient flow simulation.")
-- PJ 2024-03-26: adapt from Eilmer4 example

config.dimensions = 2

H = 0.010 -- y-direction layer thickness in metres
L = 0.100 -- x-direction wavelength in metres
ymin = -20.0*H; ymax = 20.0*H
xmin = -L; xmax = L

setGasModel('ideal-air.gas')

function initial_flow(x, y, z)
   -- User-defined function for the initial flow state works in physical space.
   local p = 100.0e3 -- Pa
   local T = 300.0 -- K
   local velx0 = 200.0 -- m/s subsonic
   -- The lower half of the domain is flowing left and the upper half, right.
   local velx = -velx0
   if y > H then
      velx = velx0
   elseif y > -H then
      velx = y/H * velx0
   end
   -- Add perturbation that is periodic west to east
   -- but gets smaller toward the north and south boundaries.
   vely = 10.0 * math.exp(-math.abs(y)/H) *
      (math.cos(x/L*math.pi) + math.sin(2*x/L*math.pi))
   -- We use the FlowState object to conveniently set all of
   -- the relevant properties.
   return FlowState:new{p=p, velx=velx, vely=vely, T=T}
end

flowDict = {initial=initial_flow}
bcDict = {}
makeFluidBlocks(bcDict, flowDict)

config.solver_mode = "transient"
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "classic_rk3"
config.max_time = 10.0e-3  -- seconds
config.max_step = 150000
config.dt_init = 1.0e-6
config.cfl_value = 0.8
config.dt_plot = 0.1e-3
