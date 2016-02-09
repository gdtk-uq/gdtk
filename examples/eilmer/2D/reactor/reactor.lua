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

initial = FlowState:new{p=1.0e5, T=4000.0, massf={N2=0.8, N=0.2}}

-- Geometry is a square
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.0, y=0.01}
c = Vector3:new{x=0.01, y=0.0}
d = Vector3:new{x=0.01, y=0.01}
patch0 = makePatch{north = Line:new{p0=b, p1=d},
		   east = Line:new{p0=c, p1=d},
		   south = Line:new{p0=a, p1=c},
		   west = Line:new{p0=a, p1=b}}

grid0 = StructuredGrid:new{psurface=patch0, niv=3, njv=3}
blk0 = SBlock:new{grid=grid0, fillCondition=initial, label="blk0"}
setHistoryPoint{ib=0, i=0, j=0}
-- Finish off config
config.max_time = 2.0e-4
config.max_step = 100000
config.dt_init = 1.0e-6
config.dt_history = 1.0e-6
config.dt_plot = 1.0e-5
config.fixed_time_step = true



