-- psl.lua
-- PJ 2021-08-20 Set up as an example.
--    2021-08-23 Better selection of parameters.
--    2021-09-03 Exercise the option gridArray when building the FBArray
--
config.title = "Periodic shear layer."
print(config.title)
config.dimensions = 2

H = 0.010 -- y-direction layer thickness in metres
L = 0.100 -- x-direction wavelength in metres
ymin = -20.0*H; ymax = 20.0*H
xmin = -L; xmax = L
nib=2; njb = 3
dx = (xmax - xmin)/nib; dy = (ymax - ymin)/njb
factor = 1
niv = math.floor(30*factor)+1
njv = math.floor(40*factor)+1
grids = {}
for ib=1, nib do
   grids[ib] = {}
   for jb=1, njb do
      local domain = CoonsPatch:new{
         p00={x=xmin+(ib-1)*dx, y=ymin+(jb-1)*dy},
         p10={x=xmin+ib*dx,     y=ymin+(jb-1)*dy},
         p01={x=xmin+(ib-1)*dx, y=ymin+jb*dy},
         p11={x=xmin+ib*dx,     y=ymin+jb*dy}}
      grids[ib][jb] = StructuredGrid:new{psurface=domain, niv=niv, njv=njv}
   end
end

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')

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

fba = FBArray:new{gridArray=grids, initialState=initial_flow}
-- We want the domain to be periodic in the x-direction.
for j=1, njb do
   connectBlocks(fba.blockArray[1][j], "west", fba.blockArray[nib][j], "east")
end

config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "classic_rk3"
config.max_time = 10.0e-3  -- seconds
config.max_step = 150000
config.dt_init = 1.0e-6
config.cfl_value = 0.8
config.dt_plot = 0.1e-3
