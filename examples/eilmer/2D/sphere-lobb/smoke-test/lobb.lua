-- Authors: Rowan G. and Peter J.
-- Date: 2016-01-19
--
-- Ported from eilmer3 examples.
--
--  Lobb, R.K. (1964)
--  Experimental measurement of shock detachment distance on spheres
--  fired in air at hypervelocities.
--  pp 519--527 in
--  The High Temperature Aspects of Hypersonic Flow
--  edited by Nelson, W.C.
--  Pergamon Press, Oxford, 1964
--
-- Description:
--   A nylon sphere, 0.5-inch diameter, fired into a ballistic test range.

config.title = "Sphere fired into air."
print(config.title)

-- This script file acts as a template. It will be configured by a "case.lua" file.
dofile('case.lua')

njb = 4
Db = 0.5 * 0.0254 -- diameter (in m) of ball bearing
Rc = Db/2
setGasModel('air-5sp.lua')
p_inf = 666.0 -- Pa (=5 Torr)
T_inf = 293.0 -- K
u_inf = 4825.0 -- m/s

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, massf={N2=0.78,O2=0.22}}

body_flow_time = Db/u_inf
print("body_flow_time=", body_flow_time)
config.reacting = true
config.reactions_file = 'air-5sp-6r.lua'
config.reaction_time_delay = 0.6*no_flow_times*body_flow_time
print("reaction_time_delay=", config.reaction_time_delay)
config.dimensions = 2
config.axisymmetric = true
config.grid_motion = "shock_fitting"
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.gasdynamic_update_scheme = "backward_euler"
config.max_time = no_flow_times*body_flow_time
print("max_time=", config.max_time)
if useOldSoln then
   config.shock_fitting_delay = 0.0
else
   config.shock_fitting_delay = 2.0*body_flow_time
end
print("shock_ftting_delay=", config.shock_fitting_delay)
config.max_step = 80000
config.dt_init = 1.0e-9
config.cfl_schedule = {{0.0, 0.2}, {3.0*body_flow_time, 0.2}, {5.0*body_flow_time, 0.4}}
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
config.max_invalid_cells = 10
config.dt_plot = config.max_time/10

-- The initial condition may be one of:
-- (a) the inflow state applied everywhere; or
-- (b) a flow field from an old solution
local initial
local fsol
if useOldSoln then
   fsol = FlowSolution:new{jobName=oldSoln_jobName, dir=oldSoln_dir,
			   tindx=oldSoln_tindx, nBlocks=njb}
   function initial(x, y, z)
      cell = fsol:find_nearest_cell_centre{x=x, y=y, z=z}
      cell.fmt = "FlowState"
      return fsol:get_cell_data(cell)
   end
else
   initial = inflow
end


if useOldSoln then
   print("Using MeshPatch for geometry.")
   oldGrid = fsol:get_sgrid{ib=0}
   grid1 = fsol:get_sgrid{ib=1}
   for i=1,njb-1 do
      print("Joining grid= ", i)
      oldGrid:joinGrid(fsol:get_sgrid{ib=i}, "north")
   end
   psurf = MeshPatch:new{sgrid=oldGrid}
else
   print("Building new geometry.")
   a = Vector3:new{x=Rc, y=0.0}
   b = Vector3:new{x=0.0, y=0.0}
   c = Vector3:new{x=Rc, y=Rc}
   d = Vector3:new{x=-0.25*Rc, y=0.0}
   e = Vector3:new{x=-0.25*Rc, y=Rc}
   f = Vector3:new{x=0.5*Rc,  y=2*Rc}
   g = Vector3:new{x=Rc,      y=2*Rc}

   bc = Arc:new{p0=b, p1=c, centre=a}
   dg = Bezier:new{points={d, e, f, g}}
   db = Line:new{p0=d, p1=b}
   gc = Line:new{p0=g, p1=c}

   psurf = makePatch{north=gc, east=bc, south=db, west=dg}
end
grid = StructuredGrid:new{psurface=psurf, niv=ncells+1, njv=ncells+1}
print "Done building grid."

print "Construct block."
blk = FBArray:new{grid=grid, initialState=initial, label='blk',
                  bcList={west=InFlowBC_ShockFitting:new{flowState=inflow},
                          north=OutFlowBC_Simple:new{}},
                  nib=1, njb=njb}
print "Done constructing block."

-- Set some recording points to capture flow data
setHistoryPoint{x=0.0, y=0.0}
setHistoryPoint{x=Rc, y=Rc}
setHistoryPoint{ib=0, i=0, j=0} -- at shock boundary (Which we expect to move)
config.dt_history = 1.0e-7
