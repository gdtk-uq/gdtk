-- reactor.lua
-- A no-flow chemical reactor.
-- RG & PJ  2015-03-10
--

job_title = "Nitrogen reactor -- recombination driven."
print(job_title)

config.dimensions = 2
config.job_title = job_title

setGasModel('nitrogen-2sp.lua')
config.reacting = true
config.reactions_file = 'chem.lua'

initial = FlowState:new{p=1.0e5, T=4000.0, massf={0.8, 0.2}}

-- Geometry is a square
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, 0.01}
c = Vector3:new{0.01, 0.0}
d = Vector3:new{0.01, 0.01}

grid0 = StructuredGrid:new{psurface=makePatch{Line:new{b, d}, Line:new{c, d}, Line:new{a, c}, Line:new{a, b}},
			   niv=3, njv=3}
blk0 = SBlock:new{grid=grid0, fillCondition=initial, label="blk0", hcellList={{0,0}}}

-- Finish off config
config.max_time = 2.0e-4
config.max_step = 100000
config.dt_init = 1.0e-6
config.dt_history = 1.0e-6
config.dt_plot = 1.0e-5
config.fixed_time_step = true



