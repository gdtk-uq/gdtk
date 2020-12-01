-- vortex.lua
-- Peter J. and Lachlan Whyborn
-- A simple version of the advecting subsonic vortex.

config.title = "Inviscid subsonic vortex -- advection test."
print(config.title)
config.dimensions = 2

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
dofile('vortex-flow-spec.lua')
function initial_gas(x, y, z)
   return fillTable({}, x, y)
end

-- Geometry is a simple box.
print("L=", L)
patch = CoonsPatch:new{p00={x=0,y=0}, p10={x=L,y=0},
                       p11={x=L,y=L}, p01={x=0,y=L}}
ncells = math.floor(80*1)
print("ncells=", ncells)
grid0 = StructuredGrid:new{psurface=patch, niv=ncells+1, njv=ncells+1}
-- Flow domain consists of a square grid of blocks and wraps around at
-- the edges of the box so that the flow is effectively periodic in x and y.
nb = 2
fba = FBArray:new{grid=grid0, nib=nb, njb=nb, initialState=initial_gas}
print("Make wrap-around connections")
for i=1,nb do
   connectBlocks(fba.blockArray[i][1], south, fba.blockArray[i][nb], north, 0)
   connectBlocks(fba.blockArray[1][i], west, fba.blockArray[nb][i], east, 0)
end

config.max_time = 50 * tau
config.max_step = 60000
config.dt_plot = tau
config.apply_limiter = false
config.extrema_clipping = false
if false then
   -- 6-cell stencil with Lagrange interpolation
   config.interpolation_order = 3
   config.cfl_value = 0.25
end
if true then
   -- Lachlan's alpha-split scheme
   config.flux_calculator = 'asf'
   config.cfl_value = 0.10
end
