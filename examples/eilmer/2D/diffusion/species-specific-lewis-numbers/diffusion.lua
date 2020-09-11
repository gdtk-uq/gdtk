--
-- This example tests the mass diffusion of two species
-- initially separated by a membrane.
-- 
-- Author: Rowan G.
-- Date: 2017-07-01
--
-- History: This file was ported from the Eilmer3 example.

config.title = "Binary diffusion of N2 into O2 (and vice versa)"

setGasModel("thermally-perfect-N2-O2-mix.lua")
config.mass_diffusion_model = "ficks_first_law"
config.diffusion_coefficient_type = "species_specific_lewis_numbers"

-- initial conditions
-- nitrogen at left end
left = FlowState:new{p=100.0e3, T=273.2, massf={N2=1.0}}
-- oxygen at right end
right = FlowState:new{p=100.0e3, T=273.2, massf={O2=1.0}}

function initialFill(x, y, z)
   if x < 0.0 then
      return left
   else
      return right
   end
end

-- geometry
xL = -2.0e-5
xR = 2.0e-5
ymin = 0.0
ymax = 1.0e-5

a = Vector3:new{x=xL, y=ymin}
b = Vector3:new{x=xR, y=ymin}
c = Vector3:new{x=xL, y=ymax}
d = Vector3:new{x=xR, y=ymax}

quad = CoonsPatch:new{p00=a, p10=b, p01=c, p11=d}
nnx = 100; nny = 10
grid = StructuredGrid:new{psurface=quad, niv=nnx+1, njv=nny+1}
blk = FluidBlock:new{grid=grid,
		     initialState=initialFill}
blk.bcList['north'] = WallBC_NoSlip_Adiabatic:new{}
blk.bcList['east'] = WallBC_NoSlip_Adiabatic:new{}
blk.bcList['south'] = WallBC_NoSlip_Adiabatic:new{}
blk.bcList['west'] = WallBC_NoSlip_Adiabatic:new{}

config.viscous = true
config.flux_calculator = "ausmdv"
config.max_time = 1.0e-6  -- seconds
config.max_step = 200000
config.dt_init = 1.0e-10
config.dt_plot = config.max_time / 10.0

