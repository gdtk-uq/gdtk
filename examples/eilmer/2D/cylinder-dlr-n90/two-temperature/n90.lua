-- n90.lua
-- RG & PJ  2015-03-09
-- 2015-04-22 build an original grid in this script
-- 2019-05-25 Billig shock-shape correlation used to shape grid.
--
config.title = "Cylinder in dissociating nitrogen flow."
print(config.title)

nsp, nmodes = setGasModel('nitrogen-gas-model.lua')
inflow = FlowState:new{p=500.0, T=700.0, T_modes={700.0}, velx=5000.0, massf={N2=1.0}}
initial = FlowState:new{p=5.0, T=300.0, T_modes={300.0}, massf={N2=1.0}}
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
Minf = inflow.velx / inflow.a
print("Minf=", Minf)

config.reacting = true
config.reactions_file = '2T-dissociating-N2.lua'
config.energy_exchange_file = 'energy-exchange.lua'

local billig_patch = require "billig_patch"
R = 0.045 -- m
bp = billig_patch.make_patch{Minf=Minf, R=R, xc=R, scale=0.95}
cf = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grid = StructuredGrid:new{psurface=bp.patch, niv=61, njv=41,
                          cfList={west=cf, east=cf}}

-- We can leave east and south as slip-walls.
blk0 = FBArray:new{grid=grid, initialState=initial, label="blk",
		       bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
			       north=OutFlowBC_Simple:new{}},
		       nib=1, njb=4}

-- Set a few more config options.
-- To get a reasonable start, we needed to set dt_init.
config.max_time = 100.0e-6
config.max_step = 40000
config.dt_plot = 20.0e-6
config.dt_init = 5.0e-9
