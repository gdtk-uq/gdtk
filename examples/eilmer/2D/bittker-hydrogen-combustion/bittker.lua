-- Author: Rowan J. Gollan
-- Date: 2016-02-29
--
-- History:
--   2016-02-29 -- Adapted from PJ and Fabs' eilmer3 version.
--   2016-03-19 -- Made a bit tidier (PJ)

config.title = "Test Case 3 from Bittker-Scullin code."
config.dimensions = 2
config.axisymmetric = true
config.stringent_cfl = true

-- Gas model setup
nsp, nmodes, gm = setGasModel("h2-o2-n2-9sp.lua")
config.reacting = true
config.reactions_file = "h2-o2-n2-9sp-18r.lua"

-- Initial flow condition
molefInit = {O2=0.1480, N2=0.5562, H2=0.2958}
massfInit = gm:molef2massf(molefInit)
inflow = FlowState:new{p=96.87e3, T=1559.0, velx=4551.73, massf=massfInit}

-- Geometry and grid setup
L = 0.1; H = 0.01 -- duct is long and skinny
duct = CoonsPatch:new{p00=Vector3:new{x=0.0, y=0.0},
		      p10=Vector3:new{x=L, y=0.0},
		      p11=Vector3:new{x=L, y=H},
		      p01=Vector3:new{x=0.0, y=H}}
grid0 = StructuredGrid:new{psurface=duct, niv=2001, njv=3}

-- Block setup as an array of many blocks
nib = 200; njb = 1
blk0 = SBlockArray{grid=grid0, fillCondition=inflow, label='blk',
		   bcList={east=OutFlowBC_Simple:new{},
			   west=InFlowBC_Supersonic:new{flowCondition=inflow}},
		   nib=nib, njb=njb}

-- Simulation setup
config.block_marching = true
config.nib = nib
config.njb = njb
config.propagate_inflow_data = true
config.flux_calc = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 0.5
config.max_time = 50.0e-6
config.max_step = 1000000
config.dt_plot = 10.0e-6
config.dt_init = 1.0e-10



