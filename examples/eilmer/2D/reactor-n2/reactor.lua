-- reactor.lua
-- A no-flow chemical reactor.
-- RG & PJ
--   2015-03-10: ported from Eilmer3
--   2018-04-14: refreshed
--   2021-06-24: integral chemical update with the implicit gas-dynamic update
--
config.title = "Nitrogen reactor -- recombination driven."
print(config.title)
config.dimensions = 2
--
setGasModel('nitrogen-2sp.lua')
config.reacting = true
config.reactions_file = 'chem.lua'
config.gasdynamic_update_scheme = 'backward_euler'
config.chemistry_update = 'integral' -- other option is 'split'
initial = FlowState:new{p=1.0e5, T=4000.0, massf={N2=0.8, N=0.2}}
--
-- Geometry is a square
pnts = {}
pnts.a = Vector3:new{x=0.0, y=0.0}
pnts.b = Vector3:new{x=0.0, y=0.01}
pnts.c = Vector3:new{x=0.01, y=0.0}
pnts.d = Vector3:new{x=0.01, y=0.01}
patch0 = CoonsPatch:new{p00=pnts.a, p01=pnts.b, p10=pnts.c, p11=pnts.d}
grid0 = StructuredGrid:new{psurface=patch0, niv=3, njv=3}
blk0 = FluidBlock:new{grid=grid0, initialState=initial}
setHistoryPoint{ib=0, i=0, j=0}
--
-- Finish off config
config.max_time = 2.0e-4
config.max_step = 100000
-- A time step of 1.0e-6 is good for the operator-split chemical update.
-- The backward-euler updates for the chemistry need to take finer steps
-- to get comparable accuracy.
config.dt_init = 0.01e-6
config.dt_history = 5.0e-6
config.dt_plot = 50.0e-6
config.fixed_time_step = true



