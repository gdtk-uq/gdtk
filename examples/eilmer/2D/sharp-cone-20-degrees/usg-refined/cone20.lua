-- cone20.lua
-- Unstructured Grid Example -- for use with Eilmer4
-- 2015-11-08  PeterJ, RowanG, KyleD

config.title = "Mach 1.5 flow over a 20 degree cone -- Unstructured Grid."
print(config.title)
config.dimensions = 2
config.axisymmetric = true

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}

-- Define the flow domain using an imported grid.
grid = UnstructuredGrid:new{filename="cone20-refined-near-shock.su2", fmt="su2text"}
my_bcDict = {INFLOW=InFlowBC_Supersonic:new{flowState=inflow},
	     OUTFLOW=OutFlowBC_Simple:new{},
	     SLIP_WALL=WallBC_WithSlip:new{},
	     INTERIOR=ExchangeBC_MappedCell:new{list_mapped_cells=true},
	     CONE_SURFACE=WallBC_WithSlip:new{}
	    }
blk = FluidBlock:new{grid=grid, initialState=inflow,
                     bcDict=my_bcDict}

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3500
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5

setHistoryPoint{x=1.0, y=0.2} -- nose of cone
setHistoryPoint{x=0.201, y=0.001} -- base of cone
