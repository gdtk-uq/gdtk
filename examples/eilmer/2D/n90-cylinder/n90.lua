-- n90.lua -- Cylinder in dissociating nitrogen flow
-- Rg & PJ  2015-03-09
--          2015-04-22 build an original grid in this script

job_title = "Cylinder in dissociating nitrogen flow."
print(job_title)

config.dimensions = 2
config.title = job_title

nsp, nmodes = setGasModel('cea-adaptive-lut-air.lua')
inflow = FlowState:new{p=500.0, T=700.0, velx=5000.0, massf={N2=1.0}}
initial = FlowState:new{p=5.0, T=300.0, massf={N2=1.0}}
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)

config.reacting = true
config.reactions_file = 'e4-chem.lua'

print "Building grid."
a = Vector3:new{x=0.045, y=0.0}
b = Vector3:new{x=0.0, y=0.0}
c = Vector3:new{x=0.013181, y=0.031820}
d = Vector3:new{x=0.045, y=0.045}
e = Vector3:new{x=0.0675, y=0.038972}
f = Vector3:new{x=-0.020, y=0.0}
g = Vector3:new{x=-0.020, y=0.050625}
h = Vector3:new{x=-0.016875, y=0.106875}
i = Vector3:new{x=0.045, y=0.135}
j = Vector3:new{x=0.07875, y=0.095625}
k = Vector3:new{x=0.084375, y=0.0675}

bc = Arc:new{p0=b, p1=c, centre=a}
cd = Arc:new{p0=c, p1=d, centre=a}
de = Arc:new{p0=d, p1=e, centre=a}

psurf = makePatch{north=Bezier:new{points={i, j, k, e}},
		  east=Polyline:new{segments={bc, cd, de}},
		  south=Line:new{p0=f, p1=b},
		  west=Bezier:new{points={f, g, h, i}}}
grid = StructuredGrid:new{psurface=psurf, niv=61, njv=41}

-- We can leave east and south as slip-walls
blk0 = SBlockArray{grid=grid, fillCondition=initial, label="blk",
		   bcList={west=InFlowBC_Supersonic:new{flowCondition=inflow},
			   north=OutFlowBC_Simple:new{}}, 
		   nib=1, njb=4}

-- Set a few more config options
config.flux_calc = ADAPTIVE
config.max_time = 100.0e-6
config.max_step = 40000
config.dt_init = 1.0e-9
config.cfl_value = 0.5 
config.dt_plot = 20.0e-6
