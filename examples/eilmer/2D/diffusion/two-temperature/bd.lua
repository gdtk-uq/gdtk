--
-- This example tests the mass diffusion of two species
-- initially separated by a membrane. With a separate 
-- vibrational temperature mode.
-- 
-- Author: Nick G.
-- Date: 2022-01-25

config.title = "Binary diffusion of N2 into N+ and e-"

config.mass_diffusion_model = "ficks_first_law"
config.diffusion_coefficient_type = "binary_diffusion"

nsp, nmodes, gmodel = setGasModel('n2-3sp-gm.lua')
config.reacting = true
config.reactions_file = 'n2-blank-rr.lua'
config.energy_exchange_file = 'n2-blank-ee.lua'

leftmf = gmodel:molef2massf({N2=0.9, ['N2+']=0.05, ['e-']=0.05})
rightmf= gmodel:molef2massf({N2=1.0, ['N2+']=0.0, ['e-']=0.0})

-- initial conditions
left = FlowState:new{p=101.325e3, T=2000.0, massf=leftmf, T_modes={2000.0}}
right= FlowState:new{p=101.325e3, T=2000.0, massf=rightmf, T_modes={2000.0}}

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
config.max_time = 8.0e-9  -- seconds
config.max_step = 100000
config.dt_init = 1.0e-12
config.dt_plot = config.max_time
config.cfl_value=0.5

