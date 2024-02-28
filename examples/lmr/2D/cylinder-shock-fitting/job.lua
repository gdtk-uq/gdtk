-- Cylinder in ideal air flow with shock-fitting boundary.
--
-- Authors: KAD, PAJ and RJG
--
-- History:
--   2015 Jan   -- first version
--   2021       -- clean up
--   2024-02-24 -- reworked for Eilmer 5
--
title = "Cylinder in ideal air flow with shock fitting boundary."
print(title)
config.solver_mode = "transient"
config.dimensions = 2
-- problem parameters
u_inf = 2430.0  -- m/s
radius = 1.0    -- m
--
nsp, nmodes, gm = setGasModel('ideal-air.gas')
initial = FlowState:new{p=100.0e3/3.0, T=200.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=100.0e3, T=300.0, velx=u_inf, vely=0.0}
flowDict = {
   initial=initial,
   inflow=inflow
}
--
print "Building grid."
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-radius, y=0.0}
c = Vector3:new{x=0.0, y=radius}
--
d = Vector3:new{x=-1.5*radius, y=0}
e = Vector3:new{x=-1.5*radius, y=radius}
f = Vector3:new{x=-radius, y=2.0*radius}
g = Vector3:new{x=0.0, y=3.0*radius}
--
psurf = makePatch{
   north=Line:new{p0=g, p1=c},
   east=Arc:new{p0=b, p1=c, centre=a},
   south=Line:new{p0=d, p1=b},
   west=Bezier:new{points={d, e, f, g}}
}
registerGridArray{
   grid=StructuredGrid:new{psurface=psurf, niv=31, njv=41},
   nib=3, njb=2,
   fsTag="initial",
   shock_fitting=true,
   bcTags={west="inflow_sf", north="outflow"}
}
--
bcDict = {
   inflow_sf=InFlowBC_ShockFitting:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- Set a few more config options
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.max_time = (radius*2)/u_inf * 20
config.max_step = 400000
config.cfl_value = 0.5
config.dt_init = 1e-7
config.dt_plot = config.max_time/50
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = (radius*2)/u_inf  -- allow for one flow length
config.max_invalid_cells = 10
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
