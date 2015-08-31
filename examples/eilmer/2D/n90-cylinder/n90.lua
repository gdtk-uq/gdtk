-- n90.lua -- Cylinder in dissociating nitrogen flow
-- Rg & PJ  2015-03-09
--          2015-04-22 build an original grid in this script

job_title = "Cylinder in dissociating nitrogen flow."
print(job_title)

config.dimensions = 2
config.title = job_title

nsp, nmodes = setGasModel('nitrogen-2sp.lua')
inflow = FlowState:new{p=500.0, T=700.0, velx=5000.0, massf={1.0, 0.0}}
initial = FlowState:new{p=5.0, T=300.0, massf={1.0, 0.0}}
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)

config.reacting = true
config.reactions_file = 'e4-chem.lua'

print "Building grid."
a = Vector3:new{0.045,    0.0}
b = Vector3:new{0.0,      0.0}
c = Vector3:new{0.013181, 0.031820}
d = Vector3:new{0.045,    0.045}
e = Vector3:new{0.0675,   0.038972}
f = Vector3:new{-0.020,    0.0}
g = Vector3:new{-0.020,    0.050625}
h = Vector3:new{-0.016875, 0.106875}
i = Vector3:new{ 0.045,    0.135}
j = Vector3:new{ 0.07875,  0.095625}
k = Vector3:new{ 0.084375, 0.0675}

bc = Arc:new{b, c, a}
cd = Arc:new{c, d, a}
de = Arc:new{d, e, a}
e_path = Polyline:new{bc, cd, de}
w_path = Bezier:new{points={f, g, h, i}}
s_path = Line:new{f, b}
n_path = Bezier:new{points={i, j, k, e}}

psurf = makePatch{n_path, e_path, s_path, w_path}
grid = StructuredGrid:new{psurface=psurf, niv=61, njv=41}
blk0 = SBlock:new{grid=grid, fillCondition=initial, label="blk-0"}

-- Grid has been built earlier in GridPro
-- gproName = 'n90.gpro'
-- grid = importGridproGrid(gproName)
-- blk0 = SBlock:new{grid=grid[0], fillCondition=initial, label="blk-0"}

-- We can leave east and south as SlipWalls
blk0.bcList[west] = SupInBC:new{flowCondition=inflow}
blk0.bcList[north] = ExtrapolateOutBC:new{}

-- Set a few more config options
config.flux_calc = ADAPTIVE
config.max_time = 100.0e-6
config.max_step = 40000
config.dt_init = 1.0e-9
config.cfl_value = 0.5 
config.dt_plot = 20.0e-6
