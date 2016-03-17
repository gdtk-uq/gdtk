-- Author: Rowan J. Gollan
-- Date: 2016-02-29
--
-- History:
--   2016-02-29 -- Adapted from PJ and Fabs' eilmer3 version.
--

config.title = "Test Case 3 from Bittker-Scullin code."
config.dimensions = 2
config.axisymmetric = true
config.stringent_cfl = true

nxcells = 2000
nycells = 2
nib = 200
njb = 1

-- Gas model setup
nsp, nmodes, gm = setGasModel("h2-o2-n2-9sp.lua")
config.reacting = true
config.reactions_file = "h2-o2-n2-9sp-18r.lua"

-- Initial flow condition
molefInit = {O2=0.1480, N2=0.5562, H2=0.2958}
massfInit = gm:molef2massf(molefInit)
inflow = FlowState:new{p=96.87e3, T=1559.0, velx=4551.73, massf=massfInit}

-- Geometry setup
--[[
   ^ y
   | 
   |
   d----------c
   |          |
   a----------b --> x
--]]

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.1, y=0.0}
c = Vector3:new{x=0.1, y=0.01}
d = Vector3:new{x=0.0, y=0.01}

ab = Line:new{p0=a, p1=b}
ad = Line:new{p0=a, p1=d}
bc = Line:new{p0=b, p1=c}
dc = Line:new{p0=d, p1=c}

-- Grid setup
grid0 = StructuredGrid:new{psurface=makePatch{north=dc, east=bc, south=ab, west=ad},
			   niv=nxcells+1, njv=nycells+1}

-- Block setup
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
config.cfl_value = 1.0
config.max_time = 5.0e-5
config.max_step = 1000000
config.dt_plot = 5.0e-6
config.dt_init = 1.0e-10



